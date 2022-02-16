# filter_combine_reads_v2
This program was designed specifically for Pacbio post-processed reads for use in building a phylogenetic tree. 

Updates to v2:

Fully automated.
File names have changed with updated pacbio software.
Deals with "fasta" formated files from NCBI's Entrez edirect utilities. 
*Please note, this format does not contain "GI" IDs, but rather "NR" IDs. As far as I am aware, and after a few hours of research, there is no way to change this format. Please see "strep_sequence.fasta" file retrieved from ncbi using esearch for specific formatting.

Required Input:

a fasta file of an OTU from Pacbio post-ccs/demultiplexing,
a fasta file containing a genus of interest downloaded from NCBI database.

Output:

"all_filtered_seqs_w_id.fasta" = a combined file of the Pacbio and NCBI sequences filtered to only contain reads where more than one identical read was found (Pacbio only).
"all_reads_info.tsv = contains all of the reads information. This is formatted to be used as annotations for a tree build with the ggtree package in Rstudio.

Please see "filter_combine_reads" for additional information.
