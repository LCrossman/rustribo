use clap::Parser;
use anyhow::{Context,Result};
use std::process::{Command, Stdio};
use std::vec::Vec;
use std::fs::File;
use std::fs::OpenOptions;
use std::str;
use std::io::Write;
use std::path::Path;
use std::path;
use std::collections::HashMap;
use protein_translate::translate;
use std::io::{self, prelude::*, BufReader};
use bio::io::fasta::Reader;
use bio::alphabets::dna::revcomp;
use bio::data_structures::interval_tree::IntervalTree;
use bio_types::annot::contig::Contig;
use bio_types::strand::ReqStrand;
use bio::data_structures::annot_map::AnnotMap;

#[derive(Parser,Default,Debug)]
#[clap(author="LCrossman",version,about="Ribosomal protein or DNA sequence extraction from DNA sequence file")]
struct Arguments {
     #[clap(short,long,default_value="1")]
     threads: String,
     #[clap(short,long)]
     seqtype: String,
     #[clap(short,long)]
     filename: String,
     #[clap(short,long,parse(from_os_str))]
     pathhmm: std::path::PathBuf,
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
     let seq = if strand == "-" { revcomp(cds_char) } else { cds_char.to_vec() };
     let sequence = String::from_utf8(seq).unwrap();
     sequence
     }
fn fetch_proteinsequence(cds_char: &[u8], strand: &str) -> String {
     let sequence = if strand == "-" { translate(&revcomp(cds_char)) } else { translate(&cds_char) };
     sequence
     }

fn concatenated_vector(concatenated_hash: HashMap<String,String>, filer: String) -> Result<String, ()> {
      let mut concat_vec: Vec<String> = vec![];
      let allkeys: Vec<_> = concatenated_hash.keys().cloned().collect();
      if allkeys.len() != 16 {
          let mut missingnamesfile = OpenOptions::new()
	     .write(true)
	     .append(true)
	     .open("missing_strains.txt")
	     .unwrap();
	  writeln!(missingnamesfile, "{}", filer).unwrap();
	  }
      else {
          //println!("OK so the length of allkeys is {}", allkeys.len());
          let ribo_list: Vec<String> = vec!["rplN".to_string(),"rplP".to_string(),"rplR".to_string(), "rplB".to_string(),"rplV".to_string(), "rplX".to_string(), "rplC".to_string(), "rplD".to_string(), "rplE".to_string(), "rplF".to_string(), "rpsJ".to_string(), "rpsS".to_string(), "rpsC".to_string(), "rpsH".to_string(),"rpsQ".to_string(),"rp10".to_string()];
          for rib in ribo_list {
              concat_vec.push(concatenated_hash[&rib].to_string());
              }
          let full_join = concat_vec.join("");
          let final_name = format!(">{}\n{}", filer, full_join);
          return Ok(final_name);
	  }
      Ok("".to_string())
      }

fn main() -> io::Result<()> {
   let args = Arguments::parse();
   //println!("args are parsed");
   let dna_protein: &str = &args.seqtype;
   match dna_protein {
      "dna" | "protein" => {
          //println!("type is {}", &dna_protein.to_string());
	  ();
	  }
      _ => {
          panic!("need to specify dna or protein return type");
	  }
   }
   //println!("error handled for the dna protein");
   let threads: &str = &args.threads;
   //println!("threads collected");
   let th = threads.trim_end().parse::<i32>().context("threads must be an integer value, running with default value 1");
   if let Err(err) = th {
      eprintln!("ERROR {}", err);
      err.chain().skip(1).for_each(|cause| eprintln!("because: {}", cause));
      std::process::exit(1);
      }
   //println!("error handled for threads");
   let filename: String = args.filename;
   let mut seqsannot: AnnotMap<String, String> = AnnotMap::new();
   let mut predictions = IntervalTree::new();
   //let gff = GFFRecord::default();
   let file = match File::open(&filename) {
       Ok(file) => file,
       Err(e) =>  panic!("file was found but problem opening the file {:?}", e),
   };
   //println!("file handled");
   let  pathhmm = &args. pathhmm;
   let  pathhmm = match File::open(& pathhmm) {
       Ok( pathhmm) =>  pathhmm,
       Err(e) => panic!("Error opening the hmm file {:?}", e),
   };
   //println!("hmmapth handed");
   let mut nb_bases = 0;
   let mut sequences = HashMap::<String, (String, usize, usize)>::new();
   //println!("initiated hashmap");
   let reader = Reader::new(file);
   //println!("STUCK!!!?? You need to use the DNA sequence fasta and check if you want protein or DNA sequence (prot) and (dna) and You need to change the file prefix in in the line with strip_prefix");
   for record in reader.records() {
       let record = record.unwrap();
       //build the contig bit here
       let lensy = str::from_utf8(record.seq()).unwrap().len();
       //println!("lens, {:?}, nb_bases {:?}, added {:?}", lensy, nb_bases, nb_bases+lensy);
       //let recid = format!("{}_{}",record.id().to_string(),i);
       let recid = format!("{}",record.id().to_string());
       //println!("so the recid will be {:?}", recid);
       let tma24 = Contig::new(recid[..].to_string(), nb_bases as isize, lensy, ReqStrand::Forward);
       seqsannot.insert_at(recid.to_string(), &tma24);
       let seqstr = record.seq().to_ascii_uppercase();
       sequences.insert(recid, ((str::from_utf8(&seqstr).unwrap()).to_string(), nb_bases, lensy));
       nb_bases += record.seq().len();
       }
   //println!("file read into records and all into hashmap");
   let outfilename = format!("{}{}",filename,".faa");
   //println!("running gene predictions");
   let output1 = {
        Command::new("prodigal").arg("-f").arg("gff").arg("-a").arg(&outfilename).arg("-i").arg(&filename).arg("-p").arg("meta").output().expect("failed to run gene prediction")
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
              if !contigname.contains(&id) { contigname.push(id); };
              let start = prodigal[3].trim().parse::<u32>().unwrap();
              let end = prodigal[4].trim().parse::<u32>().unwrap();
              let strand = prodigal[6];
              let val = (end as isize) - (start as isize);
              //println!("this is id sequences get id {:?} {:?}", id, sequences.get(id));
              for (_x,y,_z) in sequences.get(id) {
                  let newstart = y;
                  let addedstart = (newstart + (start as usize)) as isize;
        //        println!("new start and end {:?} {:?} {:?}", id, addedstart, (addedstart + val));
                  let creat = format!("{}_{}", id, count);
                  //println!("inserting {:?} {:?} {:?}", addedstart, creat, strand);
                  predictions.insert(addedstart..(addedstart + val), (creat, strand));
                  }
              }
              }
        count+=1;
        }
   //println!("finished prodigal loop");
   let pfam_names: HashMap<String, String> = HashMap::from([
       ("1874945at2".to_string(),"rp10".to_string()),
       ("1666043at2".to_string(),"rplV".to_string()),
       ("2040741at2".to_string(),"rplX".to_string()),
       ("1692188at2".to_string(),"rpsH".to_string()),
       ("1963491at2".to_string(),"rplR".to_string()),
       ("1398618at2".to_string(),"rplF".to_string()),
       ("1456375at2".to_string(),"rplE".to_string()),
       ("1799923at2".to_string(),"rplN".to_string()),
       ("2035880at2".to_string(),"rpsQ".to_string()),
       ("1786618at2".to_string(),"rplP".to_string()),
       ("1132353at2".to_string(),"rpsC".to_string()),
       ("2012682at2".to_string(),"rpsJ".to_string()),
       ("1893906at2".to_string(),"rpsS".to_string()),
       ("961486at2".to_string(),"rplB".to_string()),
       ("1572673at2".to_string(),"rplD".to_string()),
       ("1270636at2".to_string(),"rplC".to_string()),
   ]);
   let outfilehmm = format!("{}{}",&filename,"_hmm.out");
   //println!("running hmmer");
   let _output4 = {
       let mut child4 = Command::new("hmmsearch").arg("--domtblout").arg(&outfilehmm).arg("--cpu").arg(&threads).arg("--notextw").arg("--noali").arg(&args. pathhmm).arg(&outfilename).stdout(Stdio::null()).spawn().expect("failed searches");
       //stdout(Stdio::null())
       let _result4 = child4.wait().unwrap();
       //println!("concluded hmm run");
       };
   let hmmfile = File::open(&outfilehmm)?;
   let mut hmm_hash: HashMap<String, (String,String)> = HashMap::new();
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
	   if seen_before.contains(&pfam.to_string()) { () } else {
               hmm_hash.insert(id.to_string(), (pfam.to_string(), acc.to_string()));
	       seen_before.push(pfam.to_string());
	       }
        }
	}
       //println!("concluded hmmer stuff");
       //println!("hmm_hash {:?}", hmm_hash);
       let num_bases = nb_bases as isize;
       let p = Path::new(&filename);
       let mut filname: Vec<&str> = p.to_str().unwrap().split('/').collect::<Vec<&str>>();
       let filer: String = filname.pop().unwrap().to_string();
       let mut concatenated_hash: HashMap<String,String> = HashMap::new();
       //println!("num_bases is {:?}", num_bases);
       for inty in predictions.find(0..num_bases)  {
              let mut sourc: Vec<&str> = inty.data().0.split("_").collect();
              sourc.pop();
              //println!("this is sourc length {:?}", sourc.len());
              //println!("so this means that sourc is {:?}", sourc);
              if inty.data().0.contains("rRNA") { sourc.pop() } else { None };
              let src = format!("{}",sourc.join("_"));
              //println!("so now src is {:?}", src);
              let start = inty.interval().start as usize;
              let end = inty.interval().end as usize;
              //println!("this is sourc and src {:?} {:?}", sourc, src);
              let recsequence = sequences[&src].0.to_string();
              let orig_start = sequences[&src].1;
              let sliced_sequence = &recsequence[(start-orig_start)-1..(end-orig_start)];
              let cds_char = sliced_sequence.as_bytes();
              let dnaprot_seq = if args.seqtype == "protein" { fetch_proteinsequence(&cds_char, inty.data().1)} else { fetch_sequence(cds_char.to_vec(), inty.data().1) };
	      if hmm_hash.contains_key(&inty.data().0) {
	          let id = inty.data().0.to_string();
                  let (pfam,_peef) = hmm_hash.entry(id).or_insert(("".to_string(),"".to_string()));
		  //println!("identifying pfam {:?}", pfam);
		  let geney = &pfam_names[&pfam.to_string()];
		  //println!("hmmhash now {:?}", &mut hmm_hash);
		  //println!("so p itself is {:?}", p);
		  //println!("so filname is now {:?}", filname);
		  //let filname: path::PathBuf = p.iter()
		  //   .skip_while(|s| *s != &p)
		  //   .collect();
		  //println!("so filname is {:?}", &filer);
		  let fullname = format!("{}|{}",&filer, geney);
		  concatenated_hash.insert(geney.to_string(), dnaprot_seq.clone());
		  //println!(">{}\n{}",fullname,dnaprot_seq);
                  }
	      }
         let result = concatenated_vector(concatenated_hash, filer);
	 println!("{}", result.unwrap_or("".to_string()));
   Ok(())
}