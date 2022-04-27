# README for T2T-chm13v2 HSat123 analyses
#for Altemose N, 2022, Seminars in Cell and Developmental Biology


## Input Tables and Sequences
### In Input_Files you will find 
### 1) CHM13v2.0_HSat123_Strand_Subfam_DistinctArrays.bed
### This contains information on all HSat1-3 arrays in the chm13v2.0 assembly.
### Columns are chr, start, end, family/subfamily, 0, strand, start, end, rgb color, array index (chr_family_number; "small" if <10 kb), number of inversion breakpoints
### 2) chm13v2.0_hsat*_DistinctArrays_FWD.fasta.gz
### This was produced by obtaining sequences from all distinct arrays >10kb, flipping them to the FWD orientation, and concatenating discontiguous subarrays with 50k Ns sandwiched between
### these are the output of Merge_HSat123_Fastas.pl


## Output Files
### 1) Shared_24mer_matrices contains source data used to generate Fig. 2
### Each entry in the matrix represents the proportion of the smaller array 
### "contained in" the larger array, when broken into 24-mers
### Entries along the diagonal represent the 24-mer fold compression of the array, 
### i.e. the ratio of the array length to the total number of unique 24-mers
### these are the output of Count_kmers.pl
### 2) NTRprism_output
### The tophits files produced by NTRprism_ProcessFasta_v0.3e.pl with default parameters and a max span of 20000
### Columns are array index, top repeat unit lengths, column sum values for each length (same order), the best 6-mer for yielding fragments of each length (same order), the proportion of the array covered by fragments of each length when digested by its best kmer (same order)
### 3) Consensus_sequences
### [Pending]

