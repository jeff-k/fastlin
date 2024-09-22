use std::collections::HashMap;
use std::fs::read_to_string;
use std::path::PathBuf;

use bio_seq::prelude::*;

pub struct Barcodes {
    pub barcodes: HashMap<Vec<u8>, (String, u32)>,
    pub genome_size: u64,
    pub k: usize,
}

impl Barcodes {
    pub fn from_file(path: PathBuf, kmer_size: usize) -> Result<Self, String> {
        print!(" . get barcodes and genome size");
        Barcodes::from_string(&read_to_string(path).unwrap(), kmer_size)
    }

    pub fn from_string(barcode_csv: &str, k: usize) -> Result<Self, String> {
        // convert kmer_size to usize and calculate half kmer size
        let half_k_size: usize = (k - 1) / 2;

        // initialise Hashmap and genome size
        let mut barcodes: HashMap<Vec<u8>, (String, u32)> = HashMap::default();
        let mut genome_size: u64 = 0;

        // read barcode file
        let mut counter: u32 = 0;
        for l in barcode_csv.lines() {
            let inserts = l.split('\t');
            let collection = inserts.collect::<Vec<&str>>();

            if collection[0] == "genome_size" || collection[0] == "#genome_size" {
                // convert str to integer
                let parsed_result = collection[1].parse::<u64>();
                // check if the conversion was successful
                genome_size = match parsed_result {
                    Ok(parsed_number) => parsed_number,
                    Err(_) => {
                        return Err("Failed to read the genome size in barcode file".to_string());
                    }
                }
            } else {
                // build id
                let id = (collection[0].to_owned(), counter);

                // parse barcode segments
                let segs: String =
                    format!("{}{}{}", &collection[1], &collection[2], &collection[3]);
                let seg_k = &segs[50 - half_k_size..(50 - half_k_size) + k];

                // build barcode
                let barcode: Seq<Dna> = Seq::try_from(seg_k).map_err(|e| e.to_string())?;

                // save reverse complement
                barcodes.insert(barcode.revcomp().to_string().into(), id.clone());

                // save it in Hashmap
                barcodes.insert(barcode.to_string().into(), id.clone());

                counter += 1;
            }
        }
        // double-check we have the genome size
        if genome_size == 0 {
            return Err("The genome size is missing from the barcode file".to_string());
        }

        println!("	({counter} barcodes)");

        Ok(Barcodes {
            barcodes,
            genome_size,
            k,
        })
    }
}
