library(circlize)
library(reshape2)
setwd("Output_Files/Shared_24mer_matrices/")

heats = c("#FCFCFC","#F9F9AC","#F5E989","#EFCF67","#EAAB47","#DC8031","#B55522","#903317","#6A1B0D","#440905")
heatscale=function(vec){
	heats = c("#FCFCFC","#F9F9AC","#F5E989","#EFCF67","#EAAB47","#DC8031","#B55522","#903317","#6A1B0D","#440905")
	output = rep(heats[0],length(vec))
	ct=1
	for(i in seq(0.1,1,0.1)){
		output[vec>=i-0.1 & vec<i] = heats[ct]
		ct=ct+1
	}
	output[vec==1] = heats[10]
	return(output)
}


data=read.table("../../Input_Files/chm13v2.0_DistinctArrays_INFO_coords.bed")
data$newname = sapply(data$V4,function(x) paste0("A",x))
names(data)=c("name","start","end","subregion","type")
data$newname = sapply(data$subregion,function(x) paste0("A",x))
data$midpoint = round((data$end+data$start)/2)
#colors1 = c("#F09737","#439754","#932CE7","#0000F0")
data$sf ="black"
data$length = data$end - data$start

data$sf[grepl("HSat1B",data$type,perl=T)]="#F09737"
data$sf[grepl("HSat1A",data$type,perl=T)]="#439754"
data$sf[grepl("HSat2A1",data$type,perl=T)]="#011993"
data$sf[grepl("HSat2A2",data$type,perl=T)]="#0088E9"
data$sf[grepl("HSat2B",data$type,perl=T)]="#5CD1DB"
data$sf[grepl("HSat3A1",data$type,perl=T)]="#0A0A91"
data$sf[grepl("HSat3A2",data$type,perl=T)]="#6320EE"
data$sf[grepl("HSat3A3",data$type,perl=T)]="#BA6338"
data$sf[grepl("HSat3A4",data$type,perl=T)]="#0D79D3"
data$sf[grepl("HSat3A5",data$type,perl=T)]="#802368"
data$sf[grepl("HSat3A6",data$type,perl=T)]="#000000"
data$sf[grepl("HSat3B1",data$type,perl=T)]="#FF1B4B"
data$sf[grepl("HSat3B2",data$type,perl=T)]="#C74227"
data$sf[grepl("HSat3B3",data$type,perl=T)]="#EF7F66"
data$sf[grepl("HSat3B4",data$type,perl=T)]="#F7C62F"
data$sf[grepl("HSat3B5",data$type,perl=T)]="#EDC79B"

data$sf[grepl("HSat3A1",data$type,perl=T)]="#3F0171"
data$sf[grepl("HSat3A2",data$type,perl=T)]="#772193"
data$sf[grepl("HSat3A3",data$type,perl=T)]="#774B93"
data$sf[grepl("HSat3A4",data$type,perl=T)]="#9453FF"
data$sf[grepl("HSat3A5",data$type,perl=T)]="#ADA7FF"
data$sf[grepl("HSat3A6",data$type,perl=T)]="#D7C8FF"
data$sf[grepl("HSat3B1",data$type,perl=T)]="#941751"
data$sf[grepl("HSat3B2",data$type,perl=T)]="#FF2F92"
data$sf[grepl("HSat3B3",data$type,perl=T)]="#FF40FF"
data$sf[grepl("HSat3B4",data$type,perl=T)]="#FFA0FF"
data$sf[grepl("HSat3B5",data$type,perl=T)]="#FFC5D8"

###change this to change which family gets plotted!
currenttype = "HSat1A"

data1B = subset(data,type=="HSat1A" & length >= 10000)
propmat = as.matrix(read.table("chm13v2.0_hsat1A_DistinctArrays_FWD.fasta_pairwise_24mer_proportions.tsv",header=T,row.names=1))
thresh=0.1

if(currenttype=="HSat1B"){
	data1B = subset(data,type=="HSat1B" & length >= 10000)
	propmat = as.matrix(read.table("chm13v2.0_hsat1B_DistinctArrays_FWD.fasta_pairwise_24mer_proportions.tsv",header=T,row.names=1))
	thresh=0.1
}

if(currenttype=="HSat2"){
	data1B = data[grepl("HSat2",data$type) & data$length >= 10000,]
	propmat = as.matrix(read.table("chm13v2.0_hsat2_DistinctArrays_FWD.fasta_pairwise_24mer_proportions.tsv",header=T,row.names=1))
	thresh=0.25
}

if(currenttype=="HSat3"){
	data1B = data[grepl("HSat3",data$type) & data$length >= 100000,]
	propmat = as.matrix(read.table("chm13v2.0_hsat3_DistinctArrays_FWD.fasta_pairwise_24mer_proportions.tsv",header=T,row.names=1))
	thresh=0.25
}


circos.clear()
circos.par("track.height" = 0.1,gap.degree=5)
circos.genomicInitialize(data1B,tickLabelsStartFromZero=F)

circos.genomicTrack(data1B[,c(1:3,9)], ylim = c(0.5, 1.5), bg.border = NA, panel.fun = function(region, value, ...) {
     circos.genomicRect(region=region, value=1, col = sapply(value,function(x) as.character(x)), border = NA)
})

colnames(propmat)=sapply(rownames(propmat),function(x) paste0("A",x))
rownames(propmat)=colnames(propmat)
for(i in 1:nrow(propmat)) { propmat[i,i]=0}

propmatmelt = melt(propmat)
propmatmeltsub0=subset(propmatmelt,value>=thresh)
propmatmeltsub=propmatmeltsub0[order(propmatmeltsub0$value),]
df1 = data.frame()
df2 = data.frame()
for(i in 1:dim(propmatmeltsub)[1]){
	name1 = propmatmeltsub[i,1]
	name2 = propmatmeltsub[i,2]
	if(sum(data1B$newname==name1)>0 & sum(data1B$newname==name2)>0){
	midpoint1 = data1B[data1B$newname==name1,"midpoint"]
	midpoint2 = data1B[data1B$newname==name2,"midpoint"]
	chr1 = data1B[data1B$newname==name1,"name"]
	chr2 = data1B[data1B$newname==name2,"name"]
	df1[i,1]=chr1
	df2[i,1]=chr2
	df1[i,2]=midpoint1
	df2[i,2]=midpoint2
	df1[i,3]=5*propmatmeltsub[i,3]
	df2[i,3]=5*propmatmeltsub[i,3]
	df1[i,4]=propmatmeltsub[i,3]
	df2[i,4]=propmatmeltsub[i,3]
	df1[i,5]=rgb(t(col2rgb(heatscale(propmatmeltsub[i,3])))/255,alpha=0.7)
	df2[i,5]=rgb(t(col2rgb(heatscale(propmatmeltsub[i,3])))/255,alpha=0.7)
	circos.link(df1[i,1],df1[i,2],df2[i,1],df2[i,2], col = df1[i,5], border = NA,lwd=df1[i,3])
	}
}
