library(onemap)
library(magrittr)

####read in ALL onemap data:
onemapdata<-read.table("onemap.out3.csv",skip=1,stringsAsFactors=FALSE)

#create a vector of scaffold names corresponding to each locus
scafnames<-onemapdata[,1]
gsub("\\*","",scafnames) %>% gsub("_[0-9]*", "",.) -> scafnames

#create an object containing the names of the 60 scaffolds with the largest numbers of markers
first60<-names(head(sort(table(scafnames),decreasing=TRUE),n=60))

#make a new onemap file containing only markers from these 60 scaffolds
fname<-"onemapsubset.txt"
markersubset<-onemapdata[scafnames%in%first60,]
markersubset<-paste(markersubset[,1],markersubset[,2],markersubset[,3],sep=" ")
cat("90 ", length(markersubset), "\n", file=fname)
cat(markersubset, sep="\n", file=fname, append=TRUE)

nmarkers<-length(markersubset)

#read in new onemap file
onemapsubdata <- read.table(fname,skip=1,stringsAsFactors=FALSE)
crosstypes<-onemapsubdata[,2]
onemapsub <- read.outcross(".",fname)

#calculate 2 point recombination frequencies. 
onemapsub.rec <- rf.2pts(onemapsub)

	#####the output of rf.2pts is a large object. 
	#####$analysis contains a 3 dimensional array
	#####the dimensions are [nmarkers,4,c("Theta","LODs")]
	#####theta is the ML recombination fraction, LODs are the LOD scores 
	
#The next step is to pull out all the LOD scores for the recombination frequencies

#table containing numbers for all marker pairs
pcom3<-(t(combn(1:nmarkers,2)))
lodout3.1<-((apply(test3.rec$analysis[,,"LODs"],MAR=1,FUN=max)))
lodout3.2<-lodout3.1[unlist(lapply(pcom3[,2]-2,FUN=acum))+pcom3[,1]]
lodout3.3<-lodout3.2
lodout3.3[lodout3.2>10]<-10
lodout3.3[lodout3.2<3]<-0
