use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use std::fmt;
use std::fs::File;
use std::io::Write;
use std::str;
use std::{path::PathBuf, process};

mod get_barcodes;
use get_barcodes::get_barcodes;

mod input_files;
use input_files::get_input_files;

mod analyse_sample;
use analyse_sample::scan_reads;

mod process_barcodes;
use process_barcodes::process_barcodes;

#[derive(Parser, Debug)]
#[command(author = None, version, about = None, long_about = None)]
struct Args {
    /// directory containing the data files
    #[arg(short, long)]
    dir: String,

    /// file containing the reference barcodes
    #[arg(short = 'b', long)]
    barcodes: String,

    /// output file [out_fastlin.txt]
    #[arg(short = 'o', long, default_value_t = String::from("out_fastlin.txt"))]
    output: String,

    /// kmer size
    #[arg(short, long, default_value_t = 25)]
    kmer_size: u8,

    /// minimum number of kmer occurences
    #[arg(short = 'c', long, default_value_t = 4)]
    min_count: i32,

    /// minimum number of barcodes
    #[arg(short = 'n', long, default_value_t = 3)]
    n_barcodes: usize,

    /// maximum kmer coverage
    #[arg(short = 'x', long)]
    max_cov: Option<u64>,
}

#[derive(PartialEq)]
enum InputType {
    Assembly,
    Single,
    Paired,
}

impl fmt::Display for InputType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            InputType::Assembly => write!(f, "assembly"),
            InputType::Single => write!(f, "single"),
            InputType::Paired => write!(f, "paired"),
        }
    }
}

fn get_data_type(name_sample: String, vec_files: Vec<PathBuf>) -> InputType {
    // depending on the number of files, returns 'single', 'paired' or exit with error message

    let mut count_fasta = 0;
    let mut count_fastq = 0;

    for file_path in vec_files {
        if let Some(file_str) = file_path.to_str() {
            if file_str.ends_with(".fna.gz") || file_str.ends_with(".fas.gz") {
                count_fasta += 1;
            } else if file_str.ends_with(".fq.gz") || file_str.ends_with(".fastq.gz") {
                count_fastq += 1;
            }
        }
    }

    if count_fasta == 1 && count_fastq == 0 {
        InputType::Assembly
    } else if count_fasta == 0 && count_fastq == 1 {
        InputType::Single
    } else if count_fasta == 0 && count_fastq == 2 {
        InputType::Paired
    } else {
        eprintln!(
            "error: the sample {} has {} fasta and {} fastq files",
            name_sample, count_fasta, count_fastq
        );
        process::abort();
    }
}

fn main() {
    println!("\n      fastlin     \n");

    // get command line arguments
    let args = Args::parse();

    // check chosen kmer size
    if args.kmer_size < 11 || args.kmer_size > 99 || args.kmer_size % 2 == 0 {
        // warning message
        eprintln!(" Error: the kmer size should be an odd number between 11 and 99.\n");
        // exit fastlin
        std::process::exit(0);
    }

    // get reference barcodes
    let (barcodes, genome_size) = get_barcodes((&args.barcodes).into(), &args.kmer_size);

    // calculate maximum number of kmers to extract
    let kmer_limit = args.max_cov.map(|limit| limit * genome_size);

    // get samples and input files
    let all_samples = get_input_files(&args.dir);

    // sort samples
    let mut sorted_samples: Vec<_> = all_samples.iter().collect();
    sorted_samples.sort_by_key(|k| k.0);

    // create output file
    let mut output_file =
        File::create(args.output).expect("\n   Warning: couldn't not create output file.\n");
    output_file
        .write_all("#sample	data_type	k_cov	mixture	lineages	log_barcodes	log_errors\n".as_bytes())
        .expect("write failed!");

    // initialise progress bar
    let pb = ProgressBar::new(sorted_samples.len().try_into().unwrap());
    let sty = ProgressStyle::with_template("   {bar:60.cyan/blue} {pos:>7}/{len:7} {msg}")
        .unwrap()
        .progress_chars("##-");
    pb.set_style(sty);

    // process samples 1 by 1
    println!(" . analyse all samples");
    for (sample, list_files) in &sorted_samples {
        // progress bar
        pb.inc(1);

        // get sequencing type ('single' or 'paired' reads)
        let data_type = get_data_type(sample.to_string(), list_files.to_vec());

        let (kmer_limit, min_count) = match &data_type {
            InputType::Assembly => (None, 1),
            InputType::Single | InputType::Paired => (kmer_limit, args.min_count),
        };

        // scan input files
        let (barcode_found, coverage, error_message) = scan_reads(
            list_files.to_vec(),
            barcodes.to_owned(),
            &args.kmer_size,
            kmer_limit,
            genome_size,
        );

        // Note: coverage used to be fixed to 1 for assemblies

        // process barcodes
        let (lineages, mixture, string_occurences) =
            process_barcodes(barcode_found, min_count, args.n_barcodes);

        // write sample info into output file
        writeln!(
            output_file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            sample, data_type, coverage, mixture, lineages, string_occurences, error_message
        )
        .expect("Failed to write to file");
    }

    println!("   done.");
}
