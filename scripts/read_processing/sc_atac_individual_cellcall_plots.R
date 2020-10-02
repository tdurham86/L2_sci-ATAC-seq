#This script was written by Darren Cusanovich (c) 2018

library(mclust)
args = commandArgs(TRUE)

#Prefix of current experiment
currexpt = args[1]
#cellfloor = No. reads to be considered a cell
#Either specify a number or "mclust" to calculate automatically
cellfloor = args[2]

report2 = read.table(paste0(currexpt,".report.txt"),header=T)
bkgd.ind = grep("bkgd",report2$Tag)
if(length(bkgd.ind) > 0){
	nobkgdmat = report2[-bkgd.ind,]
} else {
	nobkgdmat = report2
}
cellcall = Mclust(data.frame(log10(nobkgdmat$Total)),G=2)
if(cellfloor == "mclust"){
  cellfloor = min(nobkgdmat[which(cellcall$classification == 2 & cellcall$uncertainty < 0.05),2])
}else{
  cellfloor = as.numeric(cellfloor)
}
#cellfloor = 1000

subsetmat = nobkgdmat[which(nobkgdmat$Total >= cellfloor),]

write.table(cbind(as.character(rownames(subsetmat)),as.character(subsetmat[,1])),paste0(currexpt,".readdepth.cells.indextable.txt"),row.names=F,col.names=F,sep="\t",quote=F)

pdf(paste0(currexpt,".results.pdf"),height=8,width=8)
par(mfrow=c(2,2))
subsamples = levels(nobkgdmat$Tag)

for(i in 1:length(subsamples)){
  if(subsamples[i] == "bkgd"){next}
  currind = grep(subsamples[i],	nobkgdmat$Tag)
  currsubind = grep(subsamples[i], subsetmat$Tag)
  currsub = subset(nobkgdmat,Tag == subsamples[i])
  currsubcells = which(currsub$Total >= cellfloor)
  currcellcall = Mclust(data.frame(log10(currsub$Total)),G=2)
  currcellfloor = min(currsub[which(currcellcall$classification == 2 & currcellcall$uncertainty < 0.05),2])
  hist(log10(currsub$Total)[currsub$Total > 4],breaks=60,col="mediumseagreen",main=subsamples[i],
       xlab="Number of Reads (log10)",las=1,xlim=c(0,7))
  abline(v=log10(cellfloor),lty="dashed",  lwd=2,col="blue")
  abline(v=log10(currcellfloor),lty="dashed",  lwd=2,col="red")
  legend("topright",legend = c(paste0("global ",cellfloor),paste0("current ",currcellfloor)),fill=c("blue","red"))
  boxplot(subset(subsetmat,Tag == subsamples[i])[,3]/subset(subsetmat,Tag == subsamples[i])[,2],
        ylim=c(0,1),ylab="Fraction of Reads Mapping to DHS",col="dodgerblue2",lwd=2,pch=20,las=1)
  abline(h=0.5,lwd=2,lty="dashed")
  boxplot(subset(subsetmat,Tag == subsamples[i])[,3]/subset(subsetmat,Tag == subsamples[i])[,2],
        col="dodgerblue2",lwd=2,pch=20,add=T,las=1)
}
dev.off()

pdf(paste0(currexpt,".results.hists.pdf"),height=12,width=12)
par(mfrow=c(2,2))
for(i in 1:length(subsamples)){
  if(subsamples[i] == "bkgd"){next}
  currind = grep(subsamples[i], nobkgdmat$Tag)
  currsubind = grep(subsamples[i], subsetmat$Tag)
  currsub = subset(nobkgdmat,Tag == subsamples[i])
  currsubcells = which(currsub$Total >= cellfloor)
  hist(log10(currsub$Total)[currsub$Total > 4],breaks=60,col="mediumseagreen",main=subsamples[i],
       xlab="Number of Reads (log10)",las=1,xlim=c(0,7))
  abline(v=log10(cellfloor),lty="dashed",  lwd=2)
  legend("topright",c(paste0("Total Reads: ",sum(currsub$Total)),
                    paste0("\n Total Reads (cells only): ",sum(currsub[currsubcells,2])),
                    paste0("\n Total Barcodes: ",length(currsub$Total)),
                    paste0("\n Number of Cells: ",length(subsetmat$Total[currsubcells])),
                    paste0("\n Fraction HS (cells only): ",round(sum(currsub[currsubcells,3])/sum(currsub[currsubcells,2]),3)),
                    paste0("\n Median Reads/Cell: ",median(currsub[currsubcells,2])),
                    paste0("\n Range of Reads/Cell: ",min(currsub[currsubcells,2])," - ",max(currsub[currsubcells,2]))),bty="n")
}

hist(log10(nobkgdmat$Total)[nobkgdmat$Total > 4],breaks=60,col="mediumseagreen",main="Overall",
	xlab="Number of Reads (log10)",las=1,xlim=c(0,7))
abline(v=log10(cellfloor),lty="dashed",  lwd=2)
legend("topright",c(paste0("Total Reads: ",sum(nobkgdmat$Total)),
	paste0("Total Reads (cells only): ",sum(subsetmat$Total)),
	paste0("\n Total Barcodes: ",length(nobkgdmat$Total)),
	paste0("\n Number of Cells: ",length(subsetmat$Total)),
	paste0("\n Fraction HS (cells only): ",round(sum(subsetmat[,3])/sum(subsetmat[,2]),3)),
	paste0("\n Median Reads/Cell: ",median(subsetmat[,2])),
	paste0("\n Range of Reads/Cell: ",min(subsetmat[,2])," - ",max(subsetmat[,2]))),bty="n")
dev.off()
