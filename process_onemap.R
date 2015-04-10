library(onemap)

####read in ALL onemap data:
onemapdata<-read.table("onemap.out3.csv",skip=1,stringsAsFactors=FALSE)

scafnames<-onemapdata[,1]
gsub("\\*","",locnames) %>% gsub("_.*", "",.) -> scafnames

