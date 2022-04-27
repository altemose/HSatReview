use warnings;
use strict;
my $tic = time;
print "\n\n";


my $infile = "CHM13v2.0_HSat123_Strand_Subfam_DistinctArrays.bed";
my $fastafile = "chm13v2.0_hsat123.fasta";
my $outfile1A = "chm13v2.0_hsat1A_DistinctArrays_FWD.fasta";
my $outfile1B = "chm13v2.0_hsat1B_DistinctArrays_FWD.fasta";
my $outfile2 = "chm13v2.0_hsat2_DistinctArrays_FWD.fasta";
my $outfile3 = "chm13v2.0_hsat3_DistinctArrays_FWD.fasta";
my $outfilenametable = "chm13v2.0_DistinctArrays_INFO.tsv";

open(IN, $infile);
my %coordtoname;
my %coordtostrand;
my %nametosubfam;
my %nametocoords;
my %nametolen;
my %nametoinvs;
while(my $line = <IN>){
	chomp($line);
	my @linearray = split("\t",$line);
	my ($chr, $start, $stop, $subfam, $trash, $strand) = @linearray[0..5];
	my $len = $stop-$start;
	$start=$start+1;
	my $name = $linearray[9];
	my $invs = $linearray[10];
	my $coord = $chr.':'.$start.'-'.$stop;
	$coordtoname{$coord}=$name;
	$coordtostrand{$coord}=$strand;
	$nametosubfam{$name}=$subfam;
	$nametolen{$name}+=$len;
	if(exists $nametocoords{$name}){
		$nametocoords{$name}.=','.$coord.'('.$strand.')';
	}else{
		$nametocoords{$name}=$coord.'('.$strand.')';
		$nametoinvs{$name}=$invs;
	}
}
close IN;

my $nstring = join('',('N'x50000));
#print "$nstring\n";

my @orderednames;
my %nametoseq;
if(1){
	local $/ = ">"; #change the record separator (locally)-i.e. split files by '>' instead of newlines
	open(FASTA, $fastafile);
	my $junk = <FASTA>;
	my $ct = 0;
	while(my $frecord = <FASTA>){ #read in input file one fasta record at a time
		chomp($frecord); #remove any newlines from the end of the record
		if($frecord=~/^(\S+?)\n(.*)$/s){
			my $coord = $1;
			die "cannot find info for $coord\n" unless(exists $coordtoname{$coord});
			my $seq = uc($2);
			$seq=~s/[\n\s]//g;
			my $name = $coordtoname{$coord};
			next if($name eq "small");
			my $strand = $coordtostrand{$coord};
			#my $subfam = $nametosubfam{$name};
			
			if($strand eq '-'){
				$seq = reverse($seq);
				$seq =~ tr/ACGT/TGCA/;
			}
			
			if(exists $nametoseq{$name}){
				$nametoseq{$name}.=$nstring.$seq;
			}else{
				$nametoseq{$name}=$seq;
				$orderednames[$ct]=$name;
				$ct++;
			}
		}#close if
	}#close while
	close FASTA;
}#close if



open(OUT1A, '>'.$outfile1A);
open(OUT1B, '>'.$outfile1B);
open(OUT2, '>'.$outfile2);
open(OUT3, '>'.$outfile3);
open(OUTN, '>'.$outfilenametable);

foreach my $name (@orderednames){
	my @namearray = split('_',$name);
	my $hsat = $namearray[1];
	#print "$hsat\t";
	print OUT1A ">$name\n$nametoseq{$name}\n" if($hsat eq '1A');
	print OUT1B ">$name\n$nametoseq{$name}\n" if($hsat eq '1B');
	print OUT2 ">$name\n$nametoseq{$name}\n" if($hsat eq '2');
	print OUT3 ">$name\n$nametoseq{$name}\n" if($hsat eq '3');
	print OUTN "$name\t$nametolen{$name}\t$nametosubfam{$name}\t$nametoinvs{$name}\t$nametocoords{$name}\n";
}
close OUT1A;
close OUT1B;
close OUT2;
close OUT3;
close OUTN;




#calculate and display runtime
my $toc = time;
my $elapsed = $toc-$tic;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));
