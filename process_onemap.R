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
pcom<-(t(combn(1:nmarkers,2)))

##internal onemap function used to retrieve data for marker pairs
acum<-function (w){
    if (w < 0) 
        stop("'w' should be equal to or higher than zero")
    w * (w + 1)/2
	}	


#each pair of markers has 4 LOD scores for the four possible phases. get the highest of these four for each marker pair
lodout.1<-((apply(onemapsub.rec$analysis[,,"LODs"],MAR=1,FUN=max)))

#reorder and truncate LOD scores
lodout.2<-lodout.1[unlist(lapply(pcom[,2]-2,FUN=acum))+pcom[,1]]
lodout.3<-lodout.2
lodout.3[lodout.2>10]<-10
lodout.3[lodout.2<3]<-0

#plot LOD scores
plot(NULL,xlim=c(1,nmarkers),ylim=c(1,nmarkers))
points(pcom[,1],pcom[,2],pch=15,col=rgb(1,0,0,((lodout.3)/10)),cex=.2)
abline(v=(1:nmarkers)[!duplicated(gsub("(?<=[0-9]_).*","",onemapsub.rec$marnames,perl=TRUE))])

######
pcomscaf<-cbind(onemapsub.rec$marnames[pcom[,1]],onemapsub.rec$marnames[pcom[,2]])
pcomscaf<-gsub("_[0-9]*$","",pcomscaf)

gsub("^","_",first60)->first60.2
pairscafs<-t(combn(first60.2,2))

meanlod<-c()
for(i in 1:length(pairscafs[,1])){
	
	lvec1<-pcomscaf[,1]==pairscafs[i,1]&pcomscaf[,2]==pairscafs[i,2]
	lvec2<-pcomscaf[,1]==pairscafs[i,2]&pcomscaf[,2]==pairscafs[i,1]
	meanlod<-c(meanlod, mean(lodout.3[lvec1|lvec2]))
	if((i%%10)==0){print(i)}
	}

