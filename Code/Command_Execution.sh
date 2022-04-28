#before beginning, install samtools, perl, kalign, hmmer, and optionally stringdecomposer

#download chm13v2.0 reference sequence & unzip
cd Input_Files/
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gunzip chm13v2.0.fa.gz

#extract all HSat1-3 sequences using samtools
awk -v 'OFS=''' '{print $1,":",$2+1,"-",$3}' CHM13v2.0_HSat123_Strand_Subfam_DistinctArrays.bed >CHM13v2.0_HSat123_Strand_Subfam_DistinctArrays.regions
samtools faidx -o chm13v2.0_hsat123.fasta -r CHM13v2.0_HSat123_Strand_Subfam_DistinctArrays.regions chm13v2.0.fa

#merge discontiguous sub-arrays by concatenating with 50000 intervening Ns
perl Merge_HSat123_Fastas.pl

#generate shared 24mer tables for circos plots
perl Count_kmers.pl chm13v2.0_hsat1A_DistinctArrays_FWD.fasta 24
perl Count_kmers.pl chm13v2.0_hsat1B_DistinctArrays_FWD.fasta 24
perl Count_kmers.pl chm13v2.0_hsat2_DistinctArrays_FWD.fasta 24
perl Count_kmers.pl chm13v2.0_hsat3_DistinctArrays_FWD.fasta 24

#run modified version of NTRprism
perl NTRprism_ProcessFasta_v0.3e.pl chm13v2.0_hsat1A_DistinctArrays_FWD.fasta HSat1A_ALL 20000
perl NTRprism_ProcessFasta_v0.3e.pl chm13v2.0_hsat1B_DistinctArrays_FWD.fasta HSat1B_ALL 20000
perl NTRprism_ProcessFasta_v0.3e.pl chm13v2.0_hsat2_DistinctArrays_FWD.fasta HSat2_ALL 20000
perl NTRprism_ProcessFasta_v0.3e.pl chm13v2.0_hsat3_DistinctArrays_FWD.fasta HSat3_ALL 20000

#for each fasta output from NTRprism, align and determine a consensus sequence using kalign and hmmer
for f in HSat*ALL.region*fasta; do
	j=$(perl -lae 'if(m/region_(.+)\.digest.+fraglens\_(.+)\.fasta/){print "$1.$2"}' <(echo $f))
	kalign -f fasta -i $f -o $j.afa; echo $j; grep -c ">" $f
	hmmbuild $j.hmm $j.afa 
	hmmemit -c $j.hmm >$j.consensus.fa
done

cat *consensus.fa >HSat123_consensus_sequences.fa

####applied to a single example: the HSat1B array on chrY
#kalign -f fasta -i HSat1B_ALL.region_Y_1B_1.peak_1.digest_GGGAGG.fraglens_2420bp.fasta -o HSat1B_ALL.region_Y_1B_1.peak_1.digest_GGGAGG.fraglens_2420bp.afa 
#hmmbuild Y_1B_1_1.hmm HSat1B_ALL.region_Y_1B_1.peak_1.digest_GGGAGG.fraglens_2420bp.afa 
#hmmemit -c Y_1B_1_1.hmm >Y_1B_1.consensus1.fa
#
##optional: can also use stringdecomposer to break the array sequence into individual repeat units
#samtools faidx chm13v2.0_hsat1B_DistinctArrays_FWD.fasta Y_1B_1 >Y_1B_1.fa
#stringdecomposer Y_1B_1.fa Y_1B_1.consensus1.fa -o Y_1B_1_SD_output -t 8
#
##can alternatively/additionally use nhmmer to refine the consensus, 
#but in many cases it does not change the consensus appreciably
#nhmmer --notextw -A Y_1B_1_nhmmerout.sto -o Y_1B_1_nhmmerout.txt Y_1B_1_1.hmm Y_1B_1.fa
#hmmbuild Y_1B_1_2.hmm Y_1B_1_nhmmerout.sto
#hmmemit -c Y_1B_1_2.hmm >Y_1B_1.consensus2.fa