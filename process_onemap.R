###In this script we are testing whether or not we can leverage linkage information from markers already clustered into scaffolds
###to generate linkage groups. our radseq data are too sparse and error prone to effectively use joinmap or onemap naively

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

#table containing numbers for all marker pairs. pcom = "pairwise comparisons of markers"
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

#create a square matrix containing LOD scores. dunno if this one works. didn't run it yet. 
LODmat<-matrix(nrow=nmarkers,ncol=nmarkers)
for(i in 1:length(lodout.3)){
	LODmat[pcom[i,1],pcom[i,2]]<-lodout.3[i]
	LODmat[pcom[i,2],pcom[i,1]]<-lodout.3[i]
	}
dimnames(LODmat)<-list(onemapsub.rec$marnames,onemapsub.rec$marnames)

colfunc<- colorRampPalette(c("white", "red"))
image(LODmat,zlim=c(3,10),col=colfunc(7))

#####now we calculate the mean LOD score between pairs of scaffolds

######pcomscaf enumerates pairwise comparisons of markers, but contains only scaffold names, not marker positions
pcomscaf<-cbind(onemapsub.rec$marnames[pcom[,1]],onemapsub.rec$marnames[pcom[,2]])
pcomscaf<-gsub("_[0-9]*$","",pcomscaf)

###pairscafs contains all pairwise comparisons of scaffolds in the data
gsub("^","_",first60)->first60.2
pairscafs<-t(combn(first60.2,2))

####iterate over all pairwise comparisons of scaffolds, calculating mean lod scores. 
### the number of comparisons should be (nscafs*nscafs-1)/2

meanlod<-c()
for(i in 1:length(pairscafs[,1])){
	
	lvec1<-pcomscaf[,1]==pairscafs[i,1]&pcomscaf[,2]==pairscafs[i,2]
	lvec2<-pcomscaf[,1]==pairscafs[i,2]&pcomscaf[,2]==pairscafs[i,1]
	meanlod<-c(meanlod, mean(lodout.3[lvec1|lvec2]))
	if((i%%10)==0){print(i)}
	}

#### to visualize scaffold linkages... need to play with the rgb settings and/or scale meanlod to get the best results. 
sn<-(t(combn(1:60,2)))
plot(sn[,1],sn[,2],pch=15,col=rgb(1,0,0,meanlod/9))

#####Now create clusters of scaffolds based on the mean LOD scores connecting them. 

###given a scaffold, or set of scaffolds, mean lod scores, and a Xx2 matrix of scaffold pairs, 
	###this recursive function generates a cluster of all scaffolds connected to at greater than a threshold

growclust<-function(clust, pairwise, lods,thresh=.1){
	
	startl<-length(clust)
	scafs<-rbind(pairwise[pairwise[,1]%in%clust&lods>thresh,], pairwise[pairwise[,2]%in%clust&lods>thresh,])
	scafs<-unique(as.vector(scafs))
	if(length(scafs)==0){return(clust)}	
	clust<-scafs
	endl<-length(clust)
	
	if(startl<endl){growclust(clust,pairwise,lods,thresh)}
	
	else{return(clust)}
	
	}

lgs<-cbind(gsub("^","_",first60),NA)
tout<-c()
lgroup<-1

while(sum(is.na(lgs[,2]))>0){
	
	tout<-lgs[is.na(lgs[,2]),1][1]
	tout<-growclust(tout,pairscafs,meanlod,thresh=.2)
	lgs[lgs[,1]%in%tout,2]<-lgroup
	lgroup<-lgroup+1
	
	}

lgs[order(lgs[,2]),]->lgs

#####Now we want to reorder the square matrix of LOD scores to see if our scaffold clustering is working well. 

#####this function reorders the square LODmat based on the linkage groups assigned in scaf_links (lgs above)
reorder.LODmat<-function(marker_LODs, scaf_links){
	
	loci<-rownames(marker_LODs)
	lgs<-unique(scaf_links[,2])
	ngroups<-length(lgs)
	
	neworder<-c()
	
	for(i in lgs){
		
		scafs<-scaf_links[scaf_links[,2]==i,1]
		
		for(j in scafs){
			
			neworder<-c(neworder,grep(paste(j,"_",sep=""),loci))
			
			}
		
		}
	marker_LODs[neworder,neworder]
	}

LODmat2<-reorder.LODmat(LODmat,lgs)
colfunc<- colorRampPalette(c("white", "red"))
image(LODmat2,zlim=c(1,10),col=colfunc(10),x=1:nmarkers,y=1:nmarkers)
abline(v=(1:nmarkers)[!duplicated(gsub("(?<=[0-9]_).*","",rownames(LODmat2),perl=TRUE))])
