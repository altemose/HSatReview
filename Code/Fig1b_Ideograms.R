library("data.table")
library("glue")
library("karyoploteR")
library("GenomeInfoDb")

setwd("Input_Files/")
FAI   = fread(glue("chm13v2.0.fasta.fai"), col.names = c("chr","chrlen","x","y","z"))
GEN.df = data.table(FAI$chr, 0, FAI$chrlen, gieStain="geng")
GENOME=toGRanges(GEN.df)

convert

# first CYTO option is lifted over from GRCh38, second (uncommented) uses centromere and gaps only
CYTO  = toGRanges(fread(glue("chm13v2.0_cytobands_allchrs.bed"), col.names=c("chr","start","end", "name", "gieStain")))

UN=fread(glue("CHM13v2.0_HSat123_Strand_Subfam_DistinctArrays_plotting.bed"), col.names=c("chr","start", "end","name","score","strand","start2","end2","cols","distinct","inversions","type"))
UN = UN[is.na(UN$type)==F,]
UN$color=UN$type
UN$color[UN$type=="1A"]="#FF9200"
UN$color[UN$type=="1B"]="#00994C"
UN$color[UN$type==2]="purple"
UN$color[UN$type==3]="#0000FA"
UN$name="name"
#UN$type="HSat1-3"
UN$gieStain="gpos100"
UN$status="Unresolved"
UN$version="None"
UN$border="black"

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight=100

kp <- plotKaryotype(genome = GENOME, ideogram.plotter=NULL, cytobands = CYTO, chromosomes=rev(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX","chrY")), plot.type=1,  cex=1,plot.params = pp)
kpAddCytobands(kp, lwd=0.3) # no border so they don't overwhelm small issue colors
#kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 3, tick.col="black", cex=0.5, minor.tick.dist = 1000000, minor.tick.len = 2, minor.tick.col = "gray")

kpPlotRegions(kp, data=UN, r0=0, r1=0.45,lwd=0,col=UN$color)


UNplus = UN[UN$strand=="+",]
UNminus = UN[UN$strand=="-",]
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight=100

kp <- plotKaryotype(genome = GENOME, ideogram.plotter=NULL, cytobands = CYTO, chromosomes=rev(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX","chrY")), plot.type=1,  cex=1,plot.params = pp)
kpAddCytobands(kp, lwd=0.3) # no border so they don't overwhelm small issue colors

autotrack(current.track = 1, total.tracks = 2)
kpPlotRegions(kp, data=UNminus, r0=0.1, r1=0.45,lwd=0,col=UNminus$color,avoid.overlapping=F)

autotrack(current.track = 2, total.tracks = 2)
kpPlotRegions(kp, data=UNplus, r0=0.45, r1=0.8,lwd=0,col=UNplus$color,avoid.overlapping=F)
