# rustribo
Rust version of ribosomal_snakemake
Uses HMMER3 dependency which:
1. is faster
2. means it is better for Bacteria than other domains at present
3. can only be run on x86_64 at this time 

Install Rust, details: https://www.rust-lang.org/tools/install

Example installing rustup on MacOSX or Unix OS:

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

check installation with:
rustc --version

You will also need the hmm files from the BUSCO bacteria set (bacteria_odb10.2020-03-06.tar.gz) 

at: https://busco-data.ezlab.org/v4/data/lineages/

*The BUSCO datasets are licensed under the Creative Commons Attribution-NoDerivatives 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nd/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA*

The Rust version takes a DNA sequence fasta file and provides concatenated ribosomal sequences for protein or DNA for each genome sequence.


USAGE:
   rustribo [OPTIONS] --seqtype <SEQTYPE> --filename <FILENAME> --pathhmm <PATHHMM>
   
seqtype can be either dna or protein
pathhmm is the full path for the Downloaded hmm files

