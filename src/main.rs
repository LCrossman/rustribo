use bio::{
    alphabets::dna::revcomp,
    data_structures::{annot_map::AnnotMap, interval_tree::IntervalTree},
    io::fasta,
};
use bio_types::{annot::contig::Contig, strand::ReqStrand};
use clap::{Parser, ValueEnum};
use std::{
    collections::HashMap,
    error::Error,
    fs::{self, File},
    io::{prelude::*, BufReader},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    str,
    vec::Vec,
};

#[derive(Parser, Debug)]
#[clap(
    author = "LCrossman",
    version,
    about = "Ribosomal protein or DNA sequence extraction from DNA sequence file"
)]
struct Arguments {
    #[clap(short, long, default_value_t = 1)]
    threads: u16,
    #[clap(short, long, value_enum)]
    seq_type: SequenceType,
    #[clap(short, long)]
    input_file: PathBuf,
    #[clap(short, long, parse(from_os_str))]
    pathhmm: PathBuf,
}

#[derive(ValueEnum, Clone, Debug)]
pub enum SequenceType {
    Dna,
    Protein,
}

fn validate_package_name(name: &str) -> Result<(), String> {
    if name.trim().len() != name.len() {
        Err(String::from(
            "package name cannot have leading and trailing space",
        ))
    } else {
        Ok(())
    }
}

fn fetch_sequence(cds_char: Vec<u8>, strand: &str) -> String {
    let seq = if strand == "-" {
        revcomp(cds_char)
    } else {
        cds_char.to_vec()
    };
    let sequence = String::from_utf8(seq).unwrap();
    sequence
}
fn fetch_proteinsequence(cds_char: &[u8], strand: &str) -> String {
    let sequence = if strand == "-" {
        protein_translate::translate(&revcomp(cds_char))
    } else {
        protein_translate::translate(&cds_char)
    };
    sequence
}

fn concatenated_vector(
    concatenated_hash: HashMap<String, String>,
    filer: String,
) -> Result<String, ()> {
    let mut concat_vec: Vec<String> = vec![];
    let allkeys: Vec<_> = concatenated_hash.keys().cloned().collect();
    if allkeys.len() != 16 {
        let mut missingnamesfile = fs::OpenOptions::new()
            .write(true)
            .append(true)
            .open("missing_strains.txt")
            .unwrap();
        writeln!(missingnamesfile, "{}", filer).unwrap();
    } else {
        //println!("OK so the length of allkeys is {}", allkeys.len());
        let ribo_list: Vec<String> = vec![
            "rplN".to_string(),
            "rplP".to_string(),
            "rplR".to_string(),
            "rplB".to_string(),
            "rplV".to_string(),
            "rplX".to_string(),
            "rplC".to_string(),
            "rplD".to_string(),
            "rplE".to_string(),
            "rplF".to_string(),
            "rpsJ".to_string(),
            "rpsS".to_string(),
            "rpsC".to_string(),
            "rpsH".to_string(),
            "rpsQ".to_string(),
            "rp10".to_string(),
        ];
        for rib in ribo_list {
            concat_vec.push(concatenated_hash[&rib].to_string());
        }
        let full_join = concat_vec.join("");
        let final_name = format!(">{}\n{}", filer, full_join);
        return Ok(final_name);
    }
    Ok("".to_string())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Arguments::parse();
    let mut seqsannot: AnnotMap<String, String> = AnnotMap::new();
    let mut predictions = IntervalTree::new();
    //let gff = GFFRecord::default();
    let pathhmm = &args.pathhmm;
    let pathhmm = File::open(&pathhmm)?;
    //println!("hmmapth handed");
    let mut nb_bases = 0;
    let mut sequences = HashMap::<String, (String, usize, usize)>::new();
    let reader = fasta::Reader::from_file(args.input_file)?;
    //println!("STUCK!!!?? You need to use the DNA sequence fasta and check if you want protein or DNA sequence (prot) and (dna) and You need to change the file prefix in in the line with strip_prefix");
    for record in reader.records() {
        let record = record?;
        //build the contig bit here
        let lensy = str::from_utf8(record.seq()).unwrap().len();
        //println!("lens, {:?}, nb_bases {:?}, added {:?}", lensy, nb_bases, nb_bases+lensy);
        //let recid = format!("{}_{}",record.id().to_string(),i);
        let recid = format!("{}", record.id().to_string());
        //println!("so the recid will be {:?}", recid);
        let tma24 = Contig::new(
            recid[..].to_string(),
            nb_bases as isize,
            lensy,
            ReqStrand::Forward,
        );
        seqsannot.insert_at(recid.to_string(), &tma24);
        let seqstr = record.seq().to_ascii_uppercase();
        sequences.insert(
            recid,
            ((str::from_utf8(&seqstr)?).to_string(), nb_bases, lensy),
        );
        nb_bases += record.seq().len();
    }

    //println!("file read into records and all into hashmap");

    let input_file_name = args.input_file.file_name().unwrap().to_str().unwrap();
    let outfilename = format!("{}{}", input_file_name, ".faa");
    //println!("running gene predictions");
    let output1 = {
        Command::new("prodigal")
            .arg("-f")
            .arg("gff")
            .arg("-a")
            .arg(&outfilename)
            .arg("-i")
            .arg(args.input_file)
            .arg("-p")
            .arg("meta")
            .output()
            .expect("failed to run gene prediction")
    };
    let mut contigname: Vec<&str> = vec![];
    let flott = str::from_utf8(&output1.stdout).unwrap();
    let floot: Vec<&str> = flott.split('\n').collect();
    let mut count: u32 = 1;
    //println!("going into prodigal loop");
    for flo in floot {
        if !flo.starts_with('#') {
            let prodigal: Vec<&str> = flo.split("\t").collect();
            if flo.contains("Prodigal") {
                let id = prodigal[0];
                count = if !contigname.contains(&id) { 1 } else { count };
                if !contigname.contains(&id) {
                    contigname.push(id);
                };
                let start = prodigal[3].trim().parse::<u32>().unwrap();
                let end = prodigal[4].trim().parse::<u32>().unwrap();
                let strand = prodigal[6];
                let val = (end as isize) - (start as isize);
                //println!("this is id sequences get id {:?} {:?}", id, sequences.get(id));
                for (_x, y, _z) in sequences.get(id) {
                    let newstart = y;
                    let addedstart = (newstart + (start as usize)) as isize;
                    //        println!("new start and end {:?} {:?} {:?}", id, addedstart, (addedstart + val));
                    let creat = format!("{}_{}", id, count);
                    //println!("inserting {:?} {:?} {:?}", addedstart, creat, strand);
                    predictions.insert(addedstart..(addedstart + val), (creat, strand));
                }
            }
        }
        count += 1;
    }
    //println!("finished prodigal loop");
    let pfam_names: HashMap<String, String> = HashMap::from([
        ("1874945at2".to_string(), "rp10".to_string()),
        ("1666043at2".to_string(), "rplV".to_string()),
        ("2040741at2".to_string(), "rplX".to_string()),
        ("1692188at2".to_string(), "rpsH".to_string()),
        ("1963491at2".to_string(), "rplR".to_string()),
        ("1398618at2".to_string(), "rplF".to_string()),
        ("1456375at2".to_string(), "rplE".to_string()),
        ("1799923at2".to_string(), "rplN".to_string()),
        ("2035880at2".to_string(), "rpsQ".to_string()),
        ("1786618at2".to_string(), "rplP".to_string()),
        ("1132353at2".to_string(), "rpsC".to_string()),
        ("2012682at2".to_string(), "rpsJ".to_string()),
        ("1893906at2".to_string(), "rpsS".to_string()),
        ("961486at2".to_string(), "rplB".to_string()),
        ("1572673at2".to_string(), "rplD".to_string()),
        ("1270636at2".to_string(), "rplC".to_string()),
    ]);
    let outfilehmm = format!("{}{}", &input_file_name, "_hmm.out");
    //println!("running hmmer");
    let _output4 = {
        let mut child4 = Command::new("hmmsearch")
            .arg("--domtblout")
            .arg(&outfilehmm)
            .arg("--cpu")
            .arg(&args.threads.into())
            .arg("--notextw")
            .arg("--noali")
            .arg(&args.pathhmm)
            .arg(&outfilename)
            .stdout(Stdio::null())
            .spawn()
            .expect("failed searches");
        //stdout(Stdio::null())
        let _result4 = child4.wait().unwrap();
        //println!("concluded hmm run");
    };
    let hmmfile = File::open(&outfilehmm)?;
    let mut hmm_hash: HashMap<String, (String, String)> = HashMap::new();
    let hmm_reader = BufReader::new(hmmfile);
    let mut seen_before = vec![];
    //println!("collecting hmmer stuff");
    for hmmline in hmm_reader.lines() {
        //  println!("{:?}", hmmline);
        let hmmline = hmmline.unwrap();
        if !hmmline.starts_with("#") {
            let hmmoutput: Vec<&str> = hmmline.split_whitespace().collect();
            let id = hmmoutput[0];
            let pfam = hmmoutput[3];
            let acc = hmmoutput[4];
            //  println!("HMM {:?} {:?} {:?}", id.to_string(), pfam.to_string(), acc.to_string());
            if seen_before.contains(&pfam.to_string()) {
                ()
            } else {
                hmm_hash.insert(id.to_string(), (pfam.to_string(), acc.to_string()));
                seen_before.push(pfam.to_string());
            }
        }
    }
    //println!("concluded hmmer stuff");
    //println!("hmm_hash {:?}", hmm_hash);
    let num_bases = nb_bases as isize;
    let p = Path::new(&args.input_file);
    let mut filename = p.to_str().unwrap().split('/').collect::<Vec<&str>>();
    let filer: String = filename.pop().unwrap().to_string();
    let mut concatenated_hash: HashMap<String, String> = HashMap::new();
    //println!("num_bases is {:?}", num_bases);
    for inty in predictions.find(0..num_bases) {
        let mut sourc: Vec<&str> = inty.data().0.split("_").collect();
        sourc.pop();
        //println!("this is sourc length {:?}", sourc.len());
        //println!("so this means that sourc is {:?}", sourc);
        if inty.data().0.contains("rRNA") {
            sourc.pop()
        } else {
            None
        };
        let src = format!("{}", sourc.join("_"));
        //println!("so now src is {:?}", src);
        let start = inty.interval().start as usize;
        let end = inty.interval().end as usize;
        //println!("this is sourc and src {:?} {:?}", sourc, src);
        let recsequence = sequences[&src].0.to_string();
        let orig_start = sequences[&src].1;
        let sliced_sequence = &recsequence[(start - orig_start) - 1..(end - orig_start)];
        let cds_char = sliced_sequence.as_bytes();
        let dnaprot_seq = match args.seq_type {
            SequenceType::Protein => fetch_sequence(cds_char.to_vec(), inty.data().1),
            SequenceType::Dna => fetch_proteinsequence(&cds_char, inty.data().1),
        };
        if hmm_hash.contains_key(&inty.data().0) {
            let id = inty.data().0.to_string();
            let (pfam, _peef) = hmm_hash
                .entry(id)
                .or_insert(("".to_string(), "".to_string()));
            //println!("identifying pfam {:?}", pfam);
            let geney = &pfam_names[&pfam.to_string()];
            //println!("hmmhash now {:?}", &mut hmm_hash);
            //println!("so p itself is {:?}", p);
            //println!("so filname is now {:?}", filname);
            //let filname: path::PathBuf = p.iter()
            //   .skip_while(|s| *s != &p)
            //   .collect();
            //println!("so filname is {:?}", &filer);
            let fullname = format!("{}|{}", &filer, geney);
            concatenated_hash.insert(geney.to_string(), dnaprot_seq.clone());
            //println!(">{}\n{}",fullname,dnaprot_seq);
        }
    }
    let result = concatenated_vector(concatenated_hash, filer);
    println!("{}", result.unwrap_or("".to_string()));
    Ok(())
}
