//use bio_seq::prelude::*;
use bio_streams::fastq::Fastq;
use flate2::read::MultiGzDecoder;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use crate::Barcodes;

pub fn get_reader(path: &PathBuf) -> Box<dyn BufRead + Send> {
    let filename_str = path.to_str().unwrap();
    let file = match File::open(path) {
        Ok(file) => file,
        Err(error) => panic!("Error opening compressed file: {:?}.", error),
    };
    if filename_str.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

pub fn process_buffer<R: BufRead + Unpin>(
    kmer_limit: Option<u64>,
    barcodes: &Barcodes,
    result_barcodes: &mut HashMap<String, i32>,
    reader: Fastq<R>,
) -> Result<u64, String> {
    let mut kmer_counter: u64 = 0;
    let k = barcodes.k;

    for record in reader {
        // unwrap record (contains name, sequence and quality)
        let record_ready = match record {
            Ok(record) => record,
            Err(err) => {
                return Err(format!("Error in file: {}", err));
            }
        };

        // get sequences and sequence length
        let seq = record_ready.seq;
        //let len_seq = seq.len();

        // only consider sequences long enough to have a kmer
        if seq.len() < k {
            continue;
        }
        // extract kmers (slices from Vect seq)
        for kmer in seq.windows(k) {
            // check if kmer is known -> add to count if yes or create new count if no
            if let Some(id) = barcodes.barcodes.get(kmer) {
                match result_barcodes.get(id) {
                    Some(count) => {
                        result_barcodes.insert(id.to_string(), count + 1);
                    }
                    None => {
                        result_barcodes.insert(id.to_string(), 1);
                    }
                }
            }

            // update kmer counter
            let nb_kmers = (seq.len() - k) as u64;
            kmer_counter += nb_kmers;

            if let Some(max_kmers) = kmer_limit {
                // stop process if number of maximum kmer coverage reached
                if kmer_counter > max_kmers {
                    return Ok(kmer_counter);
                }
            }
        }
    }
    Ok(kmer_counter)
}

pub fn scan_reads(
    mut vect_files: Vec<PathBuf>,
    barcodes: &Barcodes,
    kmer_limit: Option<u64>,
) -> (HashMap<String, i32>, u32, String) {
    // sort vector of paths
    vect_files.sort_by(|a, b| a.file_name().cmp(&b.file_name()));

    let mut result_barcodes: HashMap<String, i32> = HashMap::new();
    let mut kmer_counter: u64 = 0;

    for filename in vect_files {
        // set the reader
        let reader = Fastq::new(get_reader(&filename));
        match process_buffer(kmer_limit, barcodes, &mut result_barcodes, reader) {
            Ok(kmer_count) => {
                kmer_counter += kmer_count;
            }
            Err(err) => {
                return (HashMap::new(), 0, format!("{:?}", err));
            }
        }
    }
    // compute kmer coverage
    let coverage = (kmer_counter as f64 / barcodes.genome_size as f64).round() as u32;

    (result_barcodes, coverage, "".to_string())
}
