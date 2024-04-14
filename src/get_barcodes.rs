use std::collections::HashMap;
use std::fs::read_to_string;
use std::path::PathBuf;

pub fn get_barcodes(file_name: PathBuf, kmer_size: &u8) -> (HashMap<String, String>, u64) {
    print!(" . get barcodes and genome size");
    barcodes(read_to_string(file_name).unwrap(), kmer_size).unwrap()
}

pub fn barcodes(
    barcode_csv: String,
    kmer_size: &u8,
) -> Result<(HashMap<String, String>, u64), String> {
    // convert kmer_size to usize and calculate half kmer size
    let k = *kmer_size as usize;
    let half_k_size: usize = (k - 1) / 2;

    // initialise Hashmap and genome size
    let mut barcodes_id: HashMap<String, String> = HashMap::default();
    let mut genome_size: u64 = 0;

    // read barcode file
    let mut counter = 0;
    for l in barcode_csv.lines() {
        let inserts = l.split('\t');
        let collection = inserts.collect::<Vec<&str>>();

        if collection[0] == "genome_size" {
            // convert str to integer
            let parsed_result = collection[1].parse::<u64>();
            // check if the conversion was successful
            match parsed_result {
                Ok(parsed_number) => genome_size = parsed_number,
                Err(_) => {
                    return Err("Failed to read the genome size in barcode file".to_string());
                }
            }
        } else {
            // build id
            let id = format!("{}__{}", &collection[0], counter);
            // extract both sides
            let left_side = &collection[1][50 - half_k_size..];
            let right_side = &collection[3][..half_k_size];
            // build barcode
            let barcode = left_side.to_owned() + collection[2] + right_side;
            // save it in Hashmap
            barcodes_id.insert(barcode.to_owned(), id.to_owned());
            // build reverse complement and save it
            let rev_comp = revcomp(barcode.as_str());
            barcodes_id.insert(rev_comp.to_owned(), id.to_owned());

            counter += 1;
        }
    }
    // double-check we have the genome size
    if genome_size == 0 {
        panic!("The genome size is missing from the barcode file")
    }

    //println!("	({} barcodes and genome size {})", counter, genome_size);
    println!("	({} barcodes)", counter);

    Ok((barcodes_id, genome_size))
}

fn revcomp(seq: &str) -> String {
    // reverse complement sequence
    let mut rev_compl: String = String::with_capacity(seq.len());

    // iterate through the input sequence
    for c in seq.chars().rev() {
        rev_compl.push(switch_base(c))
    }
    rev_compl
}

fn switch_base(c: char) -> char {
    match c {
        'a' => 'T',
        'c' => 'G',
        't' => 'A',
        'g' => 'C',
        'A' => 'T',
        'C' => 'G',
        'T' => 'A',
        'G' => 'C',
        _ => 'N',
    }
}
