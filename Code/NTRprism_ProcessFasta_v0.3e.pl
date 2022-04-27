#Nicolas Altemose, Feb 2022
#input: a fasta file with DNA sequences
#output 1: text file with top repeat period peaks and column sums
#output 2: text file with matrix of k-mer interval information (optional) 

use warnings;
use strict;
srand(0);
my $tic = time;
print "\n\n";

my $usage = "usage: perl NTRprism_ProcessFasta_v0.3e.pl region_fasta_file_path[REQUIRED] output_prefix[NTRspectrum] total_span[30000] kmer_min_count[10] kmer_length[6] bin_size[1] suppress_matrix_output[0] suppress_fasta_output[0] array_coverage_threshold[0.3] column_sum_threshold[0.001] column_sum_window_size[20] fragment_length_tolerance[0.05]\n";

#initialize parameters and hardcode defaults
my $fastafile = "";#fasta file containing sequences to analyze
my $outbase = "NTRprism"; #prefix for output files
my $binsize = 1; #size of interval length bins (highest resolution =1)
my $span = 30000; #max interval size to interrogate, in bp. anything larger will be collapsed into span+1. minimum 1000.
my $countthresh = 10; #minimum number of times a k-mer is seen to count towards heatmap/spectrum. minimum 2
my $klen = 6; #kmer length
my $suppress = 0; #suppress matrix output [set to any value except 0 to suppress]
my $suppressF = 0; #suppress fasta output of tophit repeat instances [set to any value except 0 to suppress]
my $fragpropthresh = 0.3; #threshold for the proportion of the array covered by fragments at each periodicity
my $colsumthresh = 0.001; #threshold for the final column sum to declare a peak significant
my $window = 20; #window size to be used for combining column sums when declaring a peak significant
my $lengthtol = 0.05; #set the threshold for the tolerance of fragment lengths to output (top periodicity +- lengthtolerance * top periodicity)


#read in and check command line arguments
if(defined $ARGV[0]){
	$fastafile = $ARGV[0];
	chomp($fastafile);
}else{
	die "$usage\n";
}
if(defined $ARGV[1]){
	$outbase = $ARGV[1];
	chomp($outbase);
}
if(defined $ARGV[2]){
	$span = $ARGV[2];
	chomp($span);
	if($span<1000){
		$span=1000;
	}
}
if(defined $ARGV[3]){
	$countthresh = $ARGV[3];
	chomp($countthresh);
	if($countthresh<2){
		$countthresh=2;
	}
}
if(defined $ARGV[4]){
	$klen = $ARGV[4];
	chomp($klen);
	if($klen<1){
		$klen=1;
	}
}
if(defined $ARGV[5]){
	$binsize = $ARGV[5];
	chomp($binsize);
	if($binsize<1){
		$binsize=1;
	}
}
if(defined $ARGV[6]){
	$suppress = $ARGV[6];
	chomp($suppress);
	if($suppress != 0){
		$suppress = 1;
	}
}
if(defined $ARGV[7]){
	$suppressF = $ARGV[7];
	chomp($suppressF);
	if($suppressF != 0){
		$suppressF = 1;
	}
}
if(defined $ARGV[8]){
	$fragpropthresh = $ARGV[8];
	chomp($fragpropthresh);
	unless($fragpropthresh>=0 && $fragpropthresh<=1){
		print STDERR "ERROR: array coverage threshold must be between 0 and 1. Resorting to default of 0.1\n";
		$fragpropthresh = 0.3;
	}
}
if(defined $ARGV[9]){
	$colsumthresh = $ARGV[9];
	chomp($colsumthresh);
	unless($colsumthresh>=0 && $colsumthresh<=1){
		print STDERR "ERROR: column sum threshold must be between 0 and 1. Resorting to default of 0.005\n";
		$colsumthresh = 0.005;
	}
}
if(defined $ARGV[10]){
	$window = $ARGV[10];
	chomp($window);
	unless($window>=0 && $window<=1000){
		print STDERR "ERROR: column sum window size must be greater than 1 and less than 1000. Resorting to default of 20 bp\n";
		$window = 20;
	}
}


if($binsize > 1){
	print "NOTE: when binsize is >1, top hits and fragment sequences will not be printed. Instances with binsize >1 are intended only to produce matrix files for heatmap plotting.\n";
	if($suppressF==0){
		die "ERROR: when binsize is set to >1, fragment sequences will not be printed. Re-run with binsize=1 or with fasta output suppressed.\n";
	}
	if($suppress==1){
		die "ERROR: when binsize is set to >1, top hits and fragment sequences will not be printed, i.e. if you suppress matrix output there will be no output at all. Re-run with binsize=1 or without suppressing matrix output.\n";
	}
}

my $tophitsfile = "$outbase.bin"."$binsize.NTRprism_TopHits.txt";

my $maxbin = $span/$binsize+1;
open(TOP,'>'.$tophitsfile);
#open fasta and loop through each sequence one at a time
if(1){
	local $/ = ">"; #change the record separator (locally)-i.e. split files by '>' instead of newlines
	open(FASTA, $fastafile) or die "ERROR: unable to open $fastafile\n"; #open input file for reading

	my $ct2=0; #initialise a counter to count the number of sequences read in
	my $junk = <FASTA>; #remove first > from beginning of file
	while(my $frecord = <FASTA>){ #read in input file one fasta record at a time
		chomp($frecord); #remove any newlines from the end of the record
		$ct2++;
		if($frecord=~/^(\S+?)\n(.*)$/s){ #check whether each record matches the expected format

			my $regionname = $1;
			my $fullseq= uc($2);
			
			$regionname=~s/[\s\:\-]/\./g;
			
			#initialize variables
			my $totcount=0;
			my %kmers; #hash to store the frequency of each kmer
			my %kmerslastpos; #hash to store the last seen position of each k-mer
			my %kmersdist; #hash to store the interval lengths for each k-mer
		
			my $outfile = $outbase.".region_".$regionname.".span".$span.".k".$klen.".mincount".$countthresh.".bin".$binsize.".txt";
			$fullseq=~s/[\n\s]//g;
			my $seqlen = length($fullseq);
			my $seqlenNoN = $fullseq =~ tr/ACGT/ACGT/;

			#count kmers in sequence and record interval lengths
			for(my $i = 0; $i<$seqlen-$klen;$i++){
				my $sub = substr($fullseq,$i,$klen);
				unless($sub=~/N/){
					$kmers{$sub}++;
					if(exists $kmerslastpos{$sub}){
						my $interval = $i-$kmerslastpos{$sub};
						my $bin1 = int($interval/$binsize);
						if($interval>$span){
							$bin1 = $maxbin; #set all interval lengths above the max span to the max span +1
						}
						$kmersdist{$sub}{$bin1}++;
					}
					$kmerslastpos{$sub}=$i;
				}
			}#closes for
			my @orderedkmers = sort {$kmers{$b}<=>$kmers{$a}} keys %kmers;#sort kmers by descending frequency
			
			my @colsums = (0) x ($maxbin+1);
			my @colsums0 = @colsums;
			my @colsums1 = @colsums;
			my @nulldist = @colsums;
			my @windowedcolsums = @colsums;
			my $totalk = 0;
			#print interval distributions for each k-mer to text outfile
			open(OUT2,'>'.$outfile) unless($suppress==1);
	
			foreach my $mer (@orderedkmers){
					my $kcount = $kmers{$mer};
					last if($kcount<$countthresh);
					my $prop = $kcount/$seqlenNoN;
					$totalk++;
					my $freq  = int(1/$prop);
					print OUT2 "$mer\t$kcount\t$prop\t$freq" unless($suppress==1);
					my @bins = (0) x ($maxbin+1);
					foreach my $i (keys %{ $kmersdist{$mer} }){
						$colsums[$i]+=$kmersdist{$mer}{$i}/($kcount-1); #normalize to sum to 1 in each row
						$bins[$i]=$kmersdist{$mer}{$i};
					}
					my $printstr = join("\t",@bins); #convert array to tsv string
					print OUT2 "\t$printstr\n" unless($suppress==1);
			}
			close OUT2 unless($suppress==1);
			
			if($binsize==1){
				if($totalk>0){
					my @nums = 0..($maxbin);
					$colsums[$maxbin]=0; #exclude bins above span from maximum column sum computations
					@colsums0=@colsums;
					@colsums1 = map {$_ / $totalk} @colsums0; #normalize to total number of k-mers above count threshold
					for(my $n=0;$n<scalar(@colsums);$n++){
						if($n<scalar(@colsums)-$window){
							$windowedcolsums[$n]=eval(join('+',@colsums1[$n..$n+$window]));
						}
						@nulldist[$n]=$window*(1/(4**$klen))*exp(-1/(4**$klen)*$n);	
					}
					@colsums=@windowedcolsums;
					my @maxlist = sort {$b<=>$a} @colsums;
					my @whichmaxlist = sort {$colsums[$b]<=>$colsums[$a]} @nums;
				
					#find actual max in each window
					my @whichmaxlist2=@whichmaxlist;
					for(my $n=0;$n<scalar(@whichmaxlist);$n++){
						my $hit = $whichmaxlist[$n];
						my $max = $colsums1[$hit];
						#print "db: $n $hit $max ";
						for(my $m = $hit+1;$m<=$hit+$window;$m++){
							if($m<scalar(@colsums1)){
								if($colsums1[$m]>$max){
									$hit = $m;
									$max = $colsums1[$m];
								}
							}
						}
						#print "$hit\n";
						$whichmaxlist2[$n]=$hit;
						#$colsums[$hit]=$maxlist[$n];
					}
					
					#my $st1 = join(" ",@whichmaxlist[0..20]);
					#my $st2 = join(" ",@whichmaxlist2[0..20]);
					#print "$st1 | $st2\n";
					@whichmaxlist=@whichmaxlist2;
				
					#thin out the list of top periodicities by removing redundancies (remove periods within 10% of each other)
					my @tophits = ();
				
					for(my $i = 0; $i<scalar(@whichmaxlist); $i++){
						my $period = $whichmaxlist[$i];
					
						my $foundflag = 0;
						if($i>0){
							for(my $j=0;$j<scalar(@tophits);$j++){
								my $minsize = int(sprintf("%.0f",$tophits[$j] - $lengthtol*$tophits[$j]));
								my $maxsize = int(sprintf("%.0f",$tophits[$j] + $lengthtol*$tophits[$j]));
								if($period>=$minsize && $period <=$maxsize){
									$foundflag=1;
								}
							}
						}
					
						if($foundflag==0){
							push(@tophits,$period);
						}

					}#closes for
					
					#my $st1 = scalar(@tophits);
					#my $st2 = join(" ",@tophits[0..5]);
					#print "$st1 | $st2\n";

					#for the thinned list of periodicities, identify which ones achieve the best array coverage
					my @finaltophits = ();
					my @finaltopcolsums = ();
					my @bestkmers = ();
					my @fragprops = ();
					
					
					for(my $i = 0;$i<scalar(@tophits);$i++){
						my $period = $tophits[$i];
						my $colsum = $colsums[$period];
						if($colsum>($colsumthresh + 3*$nulldist[$period])){
							my $minsize = int(sprintf("%.0f",$period - $lengthtol*$period));
							my $maxsize = int(sprintf("%.0f",$period + $lengthtol*$period));
							if($maxsize>$span){
								$maxsize=$span;
							}
							my $bestkmer = "NA";
							#my $bestkmercount = 0;
							my $bestkmerfragcov = 0;
							#print "$period\t$colsum\t$minsize\t$maxsize\n";
						
							foreach my $mer (@orderedkmers){
								my $kcount = $kmers{$mer};
								last if($kcount<$countthresh);
								#my $topfragcount = 0;
								my $topfragcoverage=0;
								for(my $j = $minsize;$j<=$maxsize;$j++){
									if(exists $kmersdist{$mer}{$j}){
										#$topfragcount += $kmersdist{$mer}{$j};
										$topfragcoverage += $kmersdist{$mer}{$j}*$j;
									}
								}
								if($topfragcoverage>$bestkmerfragcov){
									$bestkmer = $mer;
									#$bestkmercount = $topfragcount;
									$bestkmerfragcov = $topfragcoverage;
								}
							}
					
							my $fragpropfinal = $bestkmerfragcov/$seqlenNoN;
							if($fragpropfinal>=$fragpropthresh || $i==0){
								push(@finaltophits,$period);
								push(@bestkmers,$bestkmer);
								push(@fragprops, sprintf("%.6f",$fragpropfinal));
								push(@finaltopcolsums,sprintf("%.6f",$colsum));
							}
						}
					}#closes for
	
					#print "$tophits[0]\t$finaltophits[0]\t$finaltopcolsums[0]\t$bestkmers[0]\t$fragprops[0]\n";
					if(scalar(@bestkmers)>0){
						my $testkmer=$bestkmers[0];
						my $lensum=0;
						foreach my $k (keys %{ $kmersdist{$testkmer} }){
							if($k<$maxbin){
								$lensum+=$kmersdist{$testkmer}{$k}*$k;
							}
						}
					}
					#print "$lensum\t$seqlenNoN\n";
					

					
					my $bestkmerstring = join(",",@bestkmers);
					my $fragpropstring = join(",",@fragprops);
					my $colsumstring = join(",",@finaltopcolsums);
					my $tophitstring = join(",",@finaltophits);
					#my $tophitstring2 = join(",",@finaltophits2);
				
					print TOP "$regionname\t$tophitstring\t$colsumstring\t$bestkmerstring\t$fragpropstring\n";
					print "$ct2: printed to $outfile; $seqlen $seqlenNoN; $tophitstring\t$colsumstring\t$bestkmerstring\t$fragpropstring\n";

					####now use the top k-mers to digest the DNA and output the fragments satisfying size range to a fasta
					if($suppressF==0){
				
						my $outcount=1;
						for(my $n=0;$n<scalar(@finaltophits);$n++){
						
								my $motif = $bestkmers[$n];
								my $peak1 = $finaltophits[$n];
						
								my $minsize = int(sprintf("%.0f",$peak1 - $lengthtol*$peak1));
								my $maxsize = int(sprintf("%.0f",$peak1 + $lengthtol*$peak1));
								if($maxsize>$span){
									$maxsize=$span;
								}
				
								my $outfileF = $outbase.".region_".$regionname.".peak_".$outcount.".digest_".$motif.".fraglens_".$peak1."bp.fasta";
								open(OUTF,'>'.$outfileF);
								my $prevpos=0;
								my $start = 0;
								my $end = 0;
								my $len0 =0;
								while($fullseq =~ /$motif/gs){ #loop through all instances of motif, print out fragment sequences
									my $idx = $-[0]; #get index of motif instance
									$start = $prevpos;
									$end = $idx;
									$len0 = $end-$start+1;
									#print "here $len0," if $prevpos<1000;
									if($len0>=$minsize && $len0<=$maxsize){
										my $fragseq = substr($fullseq,$start,$len0);
										my $startp = $start+1;
										my $endp = $end + 1;
										print OUTF ">$startp-$endp\n$fragseq\n";
									}
									$prevpos = $idx;
								}#closes while
									#finish up the last fragment at the end
									$start = $prevpos;
									$end = $seqlen;
									$len0 = $end-$start+1;
									if($len0>=$minsize && $len0<=$maxsize){
										my $fragseq = substr($fullseq,$start,$len0);
										my $startp = $start+1;
										my $endp = $end + 1;
										print OUTF ">$regionname:$startp-$endp\n$fragseq\n";
									}
								close OUTF;
							
								$outcount++;
						
						}#closes for
						
					} #closes if suppressF
				}else{
					print TOP "$regionname\tNA\tNA\tNA\tNA\n";
					print "$ct2: no k-mers satisfy kmer_min_count ($countthresh)\n";
				}# closes if totalk, else
			} #closes if binsize
		}else{
			die "Error reading in header number $ct2\n";
		}#closes if	
	}#closes while
	close FASTA;
}

close TOP;

#calculate and display runtime
my $toc = time;
my $elapsed = $toc-$tic;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));
