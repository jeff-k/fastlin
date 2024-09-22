//use bio_seq::prelude::*;
use bio_streams::fastq::Fastq;
use flate2::read::MultiGzDecoder;

use std::collections::HashMap;
//use std::error::Error;
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

pub struct Lineages(HashMap<String, Vec<u32>>);

impl Lineages {
    fn format_data(&self) -> String {
        // convert hashmap into a string of the following format: key (nb,nb,nb), key2 (nb,nb,nb), ...
        let mut sorted_keys: Vec<String> = self.0.keys().map(|s| s.to_owned()).collect();
        sorted_keys.sort();

        sorted_keys
            .iter()
            .map(|key| {
                let values = self.0.get(key).unwrap();
                let values_string = values
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<String>>()
                    .join(", ");
                format!("{} ({})", key, values_string)
            })
            .collect::<Vec<String>>()
            .join(", ")
    }

    fn filter(self: Self, min_barcodes: usize) -> HashMap<String, u32> {
        // filter lineages with at least min_barcodes barcodes
        let mut filtered_lineages: HashMap<String, u32> = HashMap::new();

        for (lineage_id, nb) in &self.0 {
            if nb.len() >= min_barcodes {
                filtered_lineages.insert(lineage_id.to_string(), median(nb));
            }
        }
        filtered_lineages
    }

    fn non_inclusive(self: Self, min_barcodes: usize) -> Vec<(String, u32)> {
        let filtered: HashMap<String, u32> = self.filter(min_barcodes);
        let all_keys: Vec<String> = filtered.keys().cloned().collect();
        let mut final_vect = vec![];

        for (lin, med_value) in filtered {
            let mut not_included = true;
            for key in all_keys.clone() {
                if key.starts_with(lin.as_str()) && lin != key {
                    not_included = false;
                    break;
                }
            }

            if not_included {
                final_vect.push((lin, med_value));
            }
        }
        final_vect
    }
}

pub struct Analysis {
    pub counts: HashMap<(String, u32), u32>,
    pub coverage: u32,
}

impl Analysis {
    pub fn new() -> Self {
        let counts: HashMap<(String, u32), u32> = HashMap::new();
        Analysis {
            counts,
            coverage: 0,
        }
    }

    pub fn process_buffer<R: BufRead + Unpin>(
        self: &mut Self,
        kmer_limit: Option<u64>,
        barcodes: &Barcodes,
        reader: Fastq<R>,
    ) -> Result<u64, String> {
        let mut kmer_counter: u64 = 0;

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
            if seq.len() < barcodes.k {
                continue;
            }
            // extract kmers (slices from Vect seq)
            for kmer in seq.windows(barcodes.k) {
                // check if kmer is known -> add to count if yes or create new count if no
                if let Some(id) = barcodes.barcodes.get(kmer) {
                    match self.counts.get(id) {
                        Some(count) => {
                            self.counts.insert(id.clone(), count + 1);
                        }
                        None => {
                            self.counts.insert(id.clone(), 1);
                        }
                    }
                }
            }

            // update kmer counter
            let nb_kmers = (seq.len() - barcodes.k) as u64;
            kmer_counter += nb_kmers;

            if let Some(max_kmers) = kmer_limit {
                // stop process if number of maximum kmer coverage reached
                if kmer_counter > max_kmers {
                    return Ok(kmer_counter);
                }
            }
        }
        Ok(kmer_counter)
    }

    pub fn process_barcodes(
        self: &Self,
        min_count: u32,
        min_barcodes: usize,
    ) -> (String, bool, String) {
        // merge barcode IDs to lineages
        let lineages: Lineages = self.merge_barcodes(min_count);

        // save all barcode info into String
        let log_barcodes: String = lineages.format_data();

        // get non-inclusive lineages sorted by nb occurrences
        let lineages: Vec<(String, u32)> = lineages.non_inclusive(min_barcodes);

        // check if mixture of lineages
        let mixture: bool = if lineages.len() > 1 { true } else { false };

        // convert to String
        let formatted_lineages: Vec<String> = lineages
            .iter()
            .map(|(lineage_name, med_value)| format!("{} ({})", lineage_name, med_value))
            .collect();

        let result = formatted_lineages.join(", ");

        (result, mixture, log_barcodes)
    }

    fn merge_barcodes(&self, min_occurences: u32) -> Lineages {
        let mut merged_lineages: HashMap<String, Vec<u32>> = HashMap::new();

        for (barcode_id, nb_occurences) in &self.counts {
            // only consider barcode IDs with abundances >= minimum count
            if *nb_occurences >= min_occurences {
                let lineage = &barcode_id.0;
                match merged_lineages.get(lineage) {
                    Some(_) => {
                        merged_lineages
                            .get_mut(lineage)
                            .unwrap()
                            .push(*nb_occurences);
                    }
                    None => {
                        merged_lineages.insert(lineage.clone(), Vec::new());
                        merged_lineages
                            .get_mut(lineage)
                            .unwrap()
                            .push(*nb_occurences);
                    }
                }
            }
        }
        Lineages(merged_lineages)
    }
}

fn median(values: &[u32]) -> u32 {
    let mut sorted_values = values.to_owned();
    sorted_values.sort();
    let len = sorted_values.len();
    if len % 2 == 0 {
        (sorted_values[len / 2 - 1] + sorted_values[len / 2]) / 2
    } else {
        sorted_values[len / 2]
    }
}

pub fn scan_reads(
    mut files: Vec<PathBuf>,
    barcodes: &Barcodes,
    kmer_limit: Option<u64>,
) -> Result<Analysis, String> {
    // sort vector of paths
    files.sort_by(|a, b| a.file_name().cmp(&b.file_name()));

    let mut analysis = Analysis::new();
    let mut kmer_counter: u64 = 0;

    for filename in files {
        // set the reader
        let reader = Fastq::new(get_reader(&filename));
        match analysis.process_buffer(kmer_limit, &barcodes, reader) {
            Ok(kmer_count) => {
                kmer_counter += kmer_count;
            }
            Err(err) => {
                return Err(err);
            }
        }
    }

    // compute kmer coverage
    analysis.coverage = (kmer_counter as f64 / barcodes.genome_size as f64).round() as u32;

    Ok(analysis)
}
