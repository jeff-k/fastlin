use flate2::read::MultiGzDecoder;
use seq_io::fasta::{Reader, Record};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::str;

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

pub fn scan_reads(
    mut vect_files: Vec<PathBuf>,
    barcodes: HashMap<String, String>,
    k_size: &u8,
    kmer_limit: Option<u64>,
    genome_size: u64,
) -> (HashMap<String, i32>, u32, String) {
    // initialise kmer size
    let k = *k_size as usize;

    // sort vector of paths
    vect_files.sort_by(|a, b| a.file_name().cmp(&b.file_name()));

    let mut result_barcodes: HashMap<String, i32> = HashMap::new();
    let mut kmer_counter: u64 = 0;
    let mut error_message: String = String::new();

    'label_loop: for filename in vect_files {
        // set the reader
        let buf = get_reader(&filename);
        let mut reader = Reader::new(buf);

        // lookup records 1 by 1
        while let Some(record) = reader.next() {
            // unwrap record (contains name, sequence and quality)
            let record_ready = match record {
                Ok(record) => record,
                Err(err) => {
                    // re-initialise barcode results (-> no result)
                    result_barcodes = HashMap::new();
                    kmer_counter = 0;
                    // save error message
                    error_message = format!("Error in file {:?}: {}", filename, err);
                    // stop reading file(s)
                    break 'label_loop;
                }
            };

            // get sequences and sequence length
            let seq = record_ready.seq();
            //let len_seq = seq.len();

            // only consider sequences long enough to have a kmer
            if seq.len() >= k {
                // extract kmers (slices from Vect seq)
                for n in 0..(seq.len() - k + 1) {
                    // get slice of Vect[u8]
                    let kmer = &seq[n..n + k];

                    // convert Vect[u8] into String
                    let seq_kmer = unsafe { str::from_utf8_unchecked(kmer) };

                    // check if kmer is known -> add to count if yes or create new count if no
                    if let Some(id) = barcodes.get(seq_kmer) {
                        match result_barcodes.get(id) {
                            Some(count) => {
                                result_barcodes.insert(id.to_string(), count + 1);
                            }
                            None => {
                                result_barcodes.insert(id.to_string(), 1);
                            }
                        }
                    }
                }
                // update kmer counter
                let nb_kmers = (seq.len() - k) as u64;
                kmer_counter += nb_kmers;

                if let Some(max_kmers) = kmer_limit {
                    // stop process if number of maximum kmer coverage reached
                    if kmer_counter > max_kmers {
                        break 'label_loop;
                    }
                }
            }
        }
    }
    // compute kmer coverage
    let coverage = (kmer_counter as f64 / genome_size as f64).round() as u32;

    (result_barcodes, coverage, error_message)
}
