use warnings;
use strict;
my $tic = time;
print "\n\n";

#my $infile = "chm13v2.0_DistinctArrays_INFO.tsv";
my $fastafile = "chm13v2.0_hsat1B_DistinctArrays_FWD.fasta";
my $word = 24;

if(defined $ARGV[0]){
	$fastafile = $ARGV[0];
	chomp($fastafile);
}
if(defined $ARGV[1]){
	$word = $ARGV[1];
	chomp($word);
}

my $temp = $fastafile;
$temp=~s/\..?$//;
my $outfile = $temp."_pairwise_".$word."mer_proportions.tsv";


#get all kmers in multidimensional hash
my %kmers;
my @names;
my %tots;
my %unique;

if(1){
	local $/ = ">"; #change the record separator (locally)-i.e. split files by '>' instead of newlines
	open(FASTA, $fastafile);
	my $junk = <FASTA>;
	my $ct = 0;
	while(my $frecord = <FASTA>){ #read in input file one fasta record at a time
		chomp($frecord); #remove any newlines from the end of the record
		if($frecord=~/^(\S+?)\n(.*)$/s){
			my $name = $1;
			$names[$ct]=$name;
			my $seq = uc($2);
			$seq=~s/[\n\s]//g;
			$seq=~s/N+/N/g;
			my $len = length($seq);
			my $kct = 0;
			for(my $n = 0; $n<length($seq)-$word;$n++){
				my $kmer = substr($seq,$n,$word);
				unless($kmer =~ /N/){
					if(exists $kmers{$name}{$kmer}){
						$kmers{$name}{$kmer}++;
					}else{
						$kmers{$name}{$kmer}=1;
						$unique{$name}++;
					}
					$kct++;
				}
			}
			$tots{$name}=$kct;
		}#close if
		$ct++;
	}#close while
	close FASTA;
}#close if

print "Read in all kmers! Now doing pairwise comparisons...\n";

open(OUT,'>'.$outfile);
print OUT "name";
foreach my $name(@names){
	print OUT "\t$name";
}
print OUT "\n";

for(my $i=0;$i<scalar(@names);$i++){
	my $name1 = $names[$i];
	print OUT "$name1";
	for(my $j=0;$j<scalar(@names);$j++){
		if($i>$j){
			print OUT "\t0";
		}elsif($i==$j){
			my $frac = $unique{$name1}/$tots{$name1};
			print OUT "\t$frac";
		}else{
			my $name2 = $names[$j];
			my $larger = $name1;
			my $smaller = $name2;
			if($tots{$name2}>$tots{$name1}){
				$larger = $name2;
				$smaller = $name1;
			}
			my $tot = $tots{$smaller};
			my $found=0;
			foreach my $kmer (keys %{ $kmers{$smaller} }){
				my $smallcount = $kmers{$smaller}{$kmer};
				if(exists $kmers{$larger}{$kmer}){
					my $largecount = $kmers{$larger}{$kmer};
					if($largecount>=$smallcount){
						$found += $smallcount;
					}else{
						$found += $largecount;
					}
				}
			}#close foreach
			my $frac2 = $found/$tot;
			print OUT "\t$frac2";
		}#close else
	}#close for
	print OUT "\n";
}#close for
close OUT;

print "printed distance matrix to $outfile\n";



#calculate and display runtime
my $toc = time;
my $elapsed = $toc-$tic;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));
