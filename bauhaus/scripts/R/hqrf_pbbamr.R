#
# R analysis script.
#
# Used for generating HQRF metrics to guage the performance. The output is a PDF file.
#
# Before using, the dependent libraries need to be installed. See the "install.r" script.

options(error = traceback)


# library(hexbin)
library(ggplot2)
require(gridExtra)
library(plyr)
library(pbbamr)
library(dplyr)

baseAND_ <- function(amin,amax,bmin,bmax) {
  x=c(amin,amax,bmin,bmax)
  x <- sort(x)
  x[3] - x[2] # take two middle numbers
}

baseAND <- function(amin,amax,bmin,bmax) {
  ifelse ((!is.na(amin)) & (!is.na(bmin)) & (amin < bmax) & (amax > bmin), 
          baseAND_(amin,amax,bmin,bmax) ,0)
}



if (F) {
  setwd("/home/UNIXHOME/mlakata/hqrm/m54003_161110_020238__2_B01__Sequel_3.0.17_rfSfy4")
  origbam     <- "m54003_161110_020238.subreads.bam"
  alignedbam  <- "m54003_161110_020238_nohq.aligned.bam"
  pdffile     <- "analysis.pdf"
  csvfile     <- "analysis.csv"
}

if (F) {
  setwd("/home/UNIXHOME/mlakata/hqrm/m54003_161110_020238__2_B01__Sequel_3.0.17_vMsVME")
  origbam     <- "m54003_161110_020238_poly.zmws.bam"
  alignedbam  <- "m54003_161110_020238_nohq.aligned.bam"
  pdffile     <- "m54003_161110_020238.pdf"
  csvfile     <- "m54003_161110_020238_metrics.csv"
} 

if (F)
{
  setwd("/home/UNIXHOME/mlakata/hqrm/emilia_3150469_4.0_32K")
  origbam     <- "/home/UNIXHOME/mlakata/hqrm/emilia_3150469_4.0_32K/m54003_161119_020419.subreads.bam" 
  scrapsbam   <- "/home/UNIXHOME/mlakata/hqrm/emilia_3150469_4.0_32K/m54003_161119_020419.scraps.bam" 
  polyzmqbam  <- ""# fixme
  polyscrapsbam   <- "" # fixme
  alignedbam <- "m54003_161119_020419_nohq.aligned.bam"
  pdffile <- "m54003_161119_020419.pdf" 
  csvfile <- "m54003_161119_020419_metrics.csv"

} else {
  args <- commandArgs(trailingOnly = T)
  origbam     <- args[1]
  scrapsbam   <- args[2] # ignored
  alignedbam  <- args[3]
  pdffile     <- args[4]
  csvfile     <- args[5]
}

# load the HQ regions first
hqr1 <- loadRegionsTable(origbam)
hqr <- filter(hqr1, RegionType=="HQREGION")
hqr$RegionStart[hqr$RegionStart==hqr$RegionEnd] = NA


ind1=loadPBI(alignedbam,loadSNR=T)
ind1$identity <- ind1$matches/(ind1$aend - ind1$astart)
ind1$asize <- ind1$aend - ind1$astart
ind1$qsize <- ind1$qend - ind1$qstart

ind2 <- merge(x=ind1,y=hqr,by.x="hole",by.y="HoleNumber",all=TRUE)

throwAwayNonoverlappingMappedReads = F
if (throwAwayNonoverlappingMappedReads) {
  ind3 <- ind2[(ind2$astart < ind2$RegionEnd & ind2$aend > ind2$RegionStart) | is.na(ind2$RegionStart) | is.na(ind2$astart),]
} else {
  
  ind2$c0bases <- mapply(baseAND,ind2$astart,ind2$aend,ind2$RegionStart, ind2$RegionEnd)
  
  ind3 <- ind2
}

# ind <- aggregate(subreads ~ hole, FUN = "sum", data = ind1)
ind <- ddply(ind3, .(hole), summarize, 
             subreads = length(hole), 
             astart = min(astart), aend = max(aend), 
             aSize  = sum(asize),
             c0bases = sum(c0bases),
             snrA=mean(snrA),
             snrC=mean(snrC),snrG=mean(snrG),snrT=mean(snrT),
             RegionStart = mean(RegionStart),
             RegionEnd = mean(RegionEnd),
             identity=mean(identity) )

qq1 <- loadPBI(origbam)
qq1$qsize <- qq1$qend - qq1$qstart
qq2 <- ddply(qq1, .(hole), summarize, qSize  = sum(qsize))
qq3 <- loadPBI(scrapsbam)
qq3$qsize <- qq3$qend - qq3$qstart
qq4 <- ddply(qq3, .(hole), summarize, qSize  = sum(qsize))
qq <- merge(x=qq2,y=qq4, by="hole",all=TRUE)
qq$qSize.x[is.na(qq$qSize.x)]=0
qq$qSize.y[is.na(qq$qSize.y)]=0
qq$qSize <- qq$qSize.x + qq$qSize.y
             



#hqr  <- mutate(hqr2, 
#  RegionStart = ifelse(RegionStart==RegionEnd, NA , RegionStart),
#  RegionEnd   = ifelse(RegionStart==RegionEnd, NA , RegionEnd)
#)

d_orig  <- merge(x=ind, y=qq, by="hole", all=TRUE)

d_orig$minSnr <- pmin( d_orig$snrA,d_orig$snrC,d_orig$snrG,d_orig$snrT )
d_orig$maxSnr <- pmax( d_orig$snrA,d_orig$snrC,d_orig$snrG,d_orig$snrT )
d_orig$meanSnr <- ( d_orig$snrA + d_orig$snrC + d_orig$snrG + d_orig$snrT ) / 4
# d_orig$medianSnr <- pmedian( d_orig$snrA,d_orig$snrC,d_orig$snrG,d_orig$snrT )
#d_orig$qSize  <- d_orig$qend -d_orig$qstart
#d_orig$aSize  <- d_orig$aend -d_orig$astart
d_orig$hqSize <- d_orig$RegionEnd-d_orig$RegionStart
head(d_orig)
minSnr <- 0
d <- d_orig[(d_orig$minSnr >= minSnr) | is.na(d_orig$snrA),]

d$snrBin = floor(d$minSnr)
d$class  = ifelse(is.na(d$astart),2,0) + ifelse(is.na(d$RegionStart),1,0)

d$c1bases <- ifelse(is.na(d$astart)     ,0,d$aSize  - d$c0bases)
d$c2bases <- ifelse(is.na(d$RegionStart),0,d$hqSize - d$c0bases)
d$c3bases <- d$qSize - d$c0bases - d$c1bases - d$c2bases  


# Calculate the sizes of the various "reads".
#  q=query read
#  a=aligned read
#  hq=HQ region


#  cut the data on the number of subreads
d0 <- d[d$subreads==0 & d$hqSize > 0,]
# d0 <- d[d$subreads==0,]
d1 <- d[d$subreads==1,]
d012 <- d[d$subreads<=2,]


pdf(pdffile)


#---------------------
# page 1. Title page, with statistics.

plot.new()
field <- c("bamname","host","dir","rerundir","num ZMWs","reference","comment")
value <- c(alignedbam , Sys.getenv("HOSTNAME"), getwd(), Sys.getenv("rerundir"), nrow(d), Sys.getenv("reference"),Sys.getenv("COMMENT"))
summary <- data.frame(field,value)
grid.table(summary)

mtext(getwd(), side=1, line=3, outer=F, adj=0, cex=0.7) 


#---------------------
# page 2. Pie chart showing distribution of ZMW

map_hq     <- nrow(d[!is.na(d$astart) & !is.na(d$RegionStart),])
map_nohq   <- nrow(d[!is.na(d$astart) &  is.na(d$RegionStart),])
nomap_hq   <- nrow(d[ is.na(d$astart) & !is.na(d$RegionStart),])
nomap_nohq <- nrow(d[ is.na(d$astart) &  is.na(d$RegionStart),])

zmws <- c(
  map_hq,
  map_nohq,
  nomap_hq,
  nomap_nohq
)

percents <- round(zmws * 100 / sum(zmws) , 1)

labels <- c(
  "mapped, HQ",
  "mapped, no HQ",
  "not mapped, HQ",
  "not mapped, no HQ"
)

mydf <- data.frame(labels,zmws,percents)

# gp <- pie(slices,labels=labels, main = "ZMW classification", clockwise=T)
# gp<- ggplot(data=mydf,aes(x=factor(1),y=slices,fill=factor(labels))) +  geom_bar(width=1, stat = "identity")
gp<- ggplot(data=mydf,aes(x=factor(1),y=zmws,fill=factor(labels))) +  geom_bar(width=1, stat = "identity") + coord_polar(theta = "y") + 
  ylab('') + labs(fill='labels')


gt <- tableGrob(mydf)

grid.arrange(gp,gt)


#---------------------
# page 3 Pie chart showing distribution of Bases


if (F) {
  map_hq     <- sum(as.double(d[d$class == 0,]$hqSize))
  map_nohq   <- sum(as.double(d[d$class == 1,]$aSize))
  nomap_hq   <- sum(as.double(d[d$class == 2,]$hqSize))
  allBases   <- sum(as.double(d$qSize))
  if (T) {
    nomap_nohq <- allBases - map_hq - map_nohq - nomap_hq
  } else {
    nomap_nohq <- 0 #allBases - map_hq - map_nohq - nomap_hq
  }
} else {
  map_hq     <- sum(d$c0bases)
  map_nohq   <- sum(d$c1bases)
  nomap_hq   <- sum(d$c2bases)
  nomap_nohq <- sum(d$c3bases)
}
bases <- c(
  map_hq,
  map_nohq,
  nomap_hq,
  nomap_nohq
)

percents <- round(bases * 100 / sum(bases) , 1)

labels <- c(
  "mapped, HQ",
  "mapped, no HQ",
  "not mapped, HQ",
  "not mapped, no HQ"
)

mydf <- data.frame(labels,bases,percents)

gp<- ggplot(data=mydf,aes(x=factor(1),y=bases,fill=factor(labels))) +  geom_bar(width=1, stat = "identity") + coord_polar(theta = "y") + 
  ylab('') + labs(fill='labels')


gt <- tableGrob(mydf)

grid.arrange(gp,gt)

#---------
# page 4
# looking at classification verses minSNR


dhq <- d[d$class != 3,]
p1 <- ggplot(data=d,aes(x=factor(snrBin),fill=factor(class,levels = c(0,1,2,3),
                                                     labels = labels))) + geom_bar( )

p2 <- ggplot(data=dhq,aes(x=factor(snrBin),fill=factor(class,levels = c(0,1,2,3),
                                                       labels = labels))) + geom_bar( )

p3 <- ggplot(d,aes(identity,..density..,fill=factor(class,levels = c(0,1,2,3),
                                                    labels = labels))) + 
  geom_histogram(bins=50,pos="identity",alpha=0.5) 

p4 <- ggplot(d,aes(aSize,..density..,fill=factor(class,levels = c(0,1,2,3),
                                                 labels = labels))) + 
  geom_histogram(bins=50,pos="identity",alpha=0.5)  + scale_x_log10()


grid.arrange(p1,p2,p3, p4, ncol=1)

#--------------
# page 5

unmapped <- d_orig[is.na(d_orig$aSize) & !is.na(d_orig$hqSize) ,]
mapped <- d_orig[!is.na(d_orig$aSize),]
p1<-qplot(unmapped$minSnr)+ xlim(minSnr,18)
p2<-qplot(mapped$minSnr)+ xlim(minSnr,18)
p3<-qplot(unmapped$maxSnr)+ xlim(minSnr,18)
p4<-qplot(mapped$maxSnr)+ xlim(minSnr,18)
p5<-qplot(unmapped$maxSnr,unmapped$minSnr,size=I(.1))+ xlim(minSnr,18) + ylim(minSnr,18)
p6<-qplot(mapped$maxSnr  ,  mapped$minSnr,size=I(.1))+ xlim(minSnr,18) + ylim(minSnr,18)
grid.arrange(p1,p3,p2,p4,p5,p6)


# page 6
p1 <- qplot(unmapped$maxSnr-unmapped$minSnr) + xlim(0,10)
p2 <- qplot(mapped$maxSnr-mapped$minSnr)+ xlim(0,10)

p3<- qplot(mapped$maxSnr-mapped$minSnr, mapped$aSize,size=I(.1)) + scale_y_log10()
p4<- qplot(mapped$maxSnr-mapped$minSnr, mapped$hqSize,size=I(.1)) + scale_y_log10()

# p5<- qplot(unmapped$maxSnr-unmapped$minSnr, unmapped$hqSize,size=I(.1)) + scale_y_log10()

p5 <- ggplot(unmapped, aes(maxSnr-minSnr, hqSize)) + geom_point() + scale_y_log10() + xlim(0,10)

grid.arrange(p1,p2,p3,p4,p5)


# ---
# page 7


p3<- qplot(mapped$minSnr, mapped$aSize,size=I(.1)) + scale_y_log10()
p4<- qplot(mapped$minSnr, mapped$hqSize,size=I(.1)) + scale_y_log10()

p5<- qplot(unmapped$minSnr, unmapped$hqSize,size=I(.1)) + scale_y_log10() + xlim(0,10)

grid.arrange(p3,p4,p5)

#----------------------------------------------------------------
# page 8
# scatter plot of HQR size versus mapped read size, for facets
# of 1,2,3 or 4 subreads per ZMW

npoints <- function(df){
  dim(df)[1]
}


dd = d[d$subreads>=1 & d$subreads<=4,]
dd5 = d[d$subreads>=0 & d$subreads<=5,]

eq <- ddply(dd,.(subreads), npoints)
meanBasesBadHQ <- function(df){
  c(dim(df)[1] , mean(df$hqSize) )
}
eq5 <- ddply(dd5,.(subreads), meanBasesBadHQ)

#maxReadLength=55000
ggplot(dd ,aes(aSize,hqSize)) + geom_point(alpha=1/20, size=0.01) + scale_x_log10() + scale_y_log10() + facet_wrap( ~ subreads, labeller = label_both)
#p2 <- p2 + geom_point(alpha=1/20, size=0.01, aes(colour= factor(pmin(subreads,5)))) + scale_x_log10() + scale_y_log10()
#p2 + geom_point(size=0.01)+ xlim(0,maxReadLength) + ylim(0,maxReadLength)
#p3 <- p3 + geom_point()

# ------------------------------
# page 9
#  difference in HQR size and mapped read size


ggplot(dd , aes(x=hqSize - aSize)) + geom_histogram(binwidth=100) +
  facet_wrap( ~ subreads, scale="free_y", labeller = label_both) + scale_y_log10() +
  geom_text(data=eq,aes(label=sprintf("%d reads",V1)),x=Inf,y=Inf,vjust=1,hjust=1)


# ------------------------------
# page 10
# size of HQR, binned in facets of num mapped subreads (0,1,2,3,4,5)

# meanBasesBadHQ = mean(d[d$subreads == 0,"hqSize"])

ggplot(dd5, aes(hqSize)) + geom_histogram(binwidth=100) + 
  #  annotate("text",label=sprintf("%d reads\n%f mean bases",dim(d[d$subreads == 0,])[1],meanBasesBadHQ),x=Inf,y=Inf,vjust=1,hjust=1) +
  geom_text(data=eq5,aes(label=sprintf("%d reads\n%.1f mean HQ bases",V1,V2)),x=Inf,y=Inf,vjust=1,hjust=1) +
  facet_wrap( ~ subreads, scale="free_y", labeller = label_both) + scale_y_log10()



# ----------------------------------------------
# page 11
# num mapped subreads per ZMW

ggplot(d,aes(subreads)) + geom_histogram(color="black",fill="white") + 
  scale_y_log10() +
  ggtitle("Mapped subreads per ZMW")

# ----------------------------------------------
# page 12
# number and size of mapped reads, with and without HQ

p1 <- ggplot(d[ is.na(d$RegionStart),],aes(subreads)) + geom_histogram(color="black",fill="white") + 
  scale_y_log10() +
  ggtitle("Mapped subreads per ZMW\nwith no HQ")
p2 <- ggplot(d[ is.na(d$RegionStart),],aes(aSize)) + geom_histogram(color="black",fill="white",binwidth=100) + 
  ggtitle("subread length\nwith no HQ")

d_nohq <- d[ is.na(d$RegionStart),]
bases_mapped_nohq <- sum(d_nohq$aSize)
mean_bases_mapped_nohq <- mean(d_nohq$aSize)

# geom_text(data=eq5,aes(label=sprintf("%d reads\n%.1f mean HQ bases",V1,V2)),x=Inf,y=Inf,vjust=1,hjust=1) +


p3 <- ggplot(d[ !is.na(d$RegionStart),],aes(subreads)) + geom_histogram(color="black",fill="white") + 
  scale_y_log10() +
  ggtitle("Mapped subreads per ZMW\nwith HQ")
p4 <- ggplot(d[ !is.na(d$RegionStart),],aes(aSize)) + geom_histogram(color="black",fill="white",binwidth=100) + 
  ggtitle("subread length\nwith HQ")

d_whq <- d[! is.na(d$RegionStart) & !is.na(d$aSize),]
bases_mapped_whq <- sum(d_whq$aSize)
mean_bases_mapped_whq <- mean(d_whq$aSize)

nrow(d_nohq)/ (nrow(d_nohq) + nrow(d_whq) )

print(bases_mapped_nohq / (bases_mapped_nohq + bases_mapped_whq ))


grid.arrange(p1,p2, p3,p4,ncol=2)

# ----------------------------------------------
# page 13
# comparison of HQR size distribution, for 0,1 or 2 mapped subreads

ggplot(d012,aes(hqSize, fill=factor(subreads))) + geom_histogram(alpha=0.5,binwidth=100,pos="identity") + 
  xlim(0,25000)  +
  scale_fill_discrete(name="Mapped Subreads\nper ZMW") +
  theme(legend.position=c(.5, .5))
# see http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/

nrow(d0)/nrow(d1)

#-----------------------------------------
# page 14
# distributions of advanced and retarded regions of HQRs

# d1$overlap <- ! ( d1$aend < d1$RegionStart || d1$RegionEnd < d1$astart)

if (T) {
  d$hqadv <- d$RegionStart - d$astart
  d$hqret <- d$RegionEnd   - d$aend
  d1234 <- d[!is.na(d$subreads) & d$subreads <=4,]
  p1 <- ggplot(d1234, aes(hqadv)) + geom_histogram(binwidth=100) + xlim(-10000,10000) + scale_y_log10() + facet_wrap( ~ subreads, labeller = label_both)
  p2 <- ggplot(d1234, aes(hqret)) + geom_histogram(binwidth=100) + xlim(-10000,10000) + scale_y_log10() + facet_wrap( ~ subreads, labeller = label_both)
  qadv <- quantile(d1234$hqadv,probs=c(0,0.1,0.5,0.9,1.0),na.rm=T)
  qret <- quantile(d1234$hqret,probs=c(0,0.1,0.5,0.9,1.0),na.rm=T)
} else {  
  d1$hqadv <- d1$RegionStart - d1$astart
  d1$hqret <- d1$RegionEnd - d1$aend
  p1 <- ggplot(d1, aes(hqadv)) + geom_histogram(binwidth=100) + xlim(-10000,10000) + scale_y_log10()
  p2 <- ggplot(d1, aes(hqret)) + geom_histogram(binwidth=100) + xlim(-10000,10000) + scale_y_log10()
  qadv <- quantile(d1$hqadv,probs=c(0,0.1,0.5,0.9,1.0),na.rm=T)
  qret <- quantile(d1$hqret,probs=c(0,0.1,0.5,0.9,1.0),na.rm=T)
}

#pg1 <- ggplot_build(p1)
#g1 <- pg1$data[[1]]
#fit <- nls(ndensity ~ (C1 * exp(-(x-mean1)**2/(2 * sigma1**2)) +
#                         C2 * exp(-(x-mean2)**2/(2 * sigma2**2))
#), data=g1,
#start=list(C1=6345151, mean1=-30, sigma1=53,
#           C2= 150000, mean2=-1000, sigma2=1000
#))
#
#
#dffit <- data.frame(x=seq(-10000,10000,100))
#dffit$ndensity <- predict(fit, newdata=dffit)
#p1 + geom_smooth(data=dffit, aes(ndensity, x),stat="identity", color="red", size=1.5)


qnames <- c("0%","10%","50%","90%","100%")

qqdf <- data.frame(qnames,qadv,qret)

#names(qadv) <- c("0%","10%","50%","90%","100%")
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(qqdf, rows=NULL, theme=tt)
grid.arrange(p1,p2,tbl, as.table=T, ncol=2)

  

#----------------------------------------------------
if(FALSE) {
  C1 <- 100000
  C2 <-  30
  mean1<- 10
  mean2<- 100
  sigma1=20
  sigma2=4000
  #mycurve <- function(x) C1 * exp(-(x-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(x-mean2)**2/(2 * sigma2**2))
  mycurve <- function(x) {
    y1 <-C2 * exp(-(x-mean2)**2/(2 * sigma2**2))
    #  y0 <- ifelse(abs(x-mean1)<sigma1,C1*(1-abs(x-mean1)/sigma1),0)
    y0 <- C1*(1/(sigma1* (((x-mean1)/sigma1)**2+1) ))
    return( y1+y0)
  }
  p1 <- ggplot(d1, aes(hqret)) + geom_histogram(binwidth=10) + xlim(-2000,2000) + scale_y_log10() + stat_function(fun=mycurve,color="red")
  p1
  pg1 <- ggplot_build(p1)
  g1 <- pg1$data[[1]]
  fit <- nls(count ~ mycurve(x), data=g1,
             start=list(C1=100000, mean1=10, sigma1=20,
                        C2= 30, mean2= 100, sigma2=4000))
}

#-----------------------------------------------------
# page 15
# HQR ZMWs with no mapped subreads

d0_big <- d0[!is.na(d0$hqSize) & d0$hqSize > 10000,]
print("d0_big summary")
summary(d0_big)
print("head d0_big")
head(d0_big$hole)
xx <- d0_big$hole %/% 65536
yy <- d0_big$hole %% 65536
print("xx")
head(xx)
locations <- data.frame(d0_big$hole,xx,yy,round(d0_big$minSnr,1),d0_big$RegionStart)
plot.new()
grid.table(head(locations,n=20)) + title("HQR ZMWs with no mapped subreads\n(first 20)")
write.csv(locations,"hq_nomapped_zmws.csv")

#-----------------------------------------------------
# page 16
# no HQR with mapped subreads

d1_big    <- d1[is.na(d1$hqSize),]
zmwNumber <- d1_big$hole
xx        <- d1_big$hole %/% 65536
yy        <- d1_big$hole %% 65536
aSize     <- d1_big$aSize
identity  <- d1_big$identity

locations <- data.frame(zmwNumber,xx,yy,aSize,identity)
plot.new()
grid.table(head(locations,n=20)) + title("NO HQR ZMWs with mapped subreads\n(first 20)")
write.csv(locations,"nohq_mapped_zmws.csv")

#-----------------------------------------------------
# page 17
# HQR with mapped subreads (good ones)

good <- d1[d1$class==0,]
zmwNumber <- good$hole
xx        <- good$hole %/% 65536
yy        <- good$hole %% 65536
aSize     <- good$aSize
identity  <- good$identity

locations <- data.frame(zmwNumber,xx,yy,aSize,identity)
plot.new()
grid.table(head(locations,n=20)) + title("HQR ZMWs with mapped subreads\n(first 20)")
write.csv(locations,"hq_mapped_zmws.csv")
#-------------------------------------------------------
# page 18
# summary page

plot.new()

metric <- c(
  "NUM_ZMWs",
  "NUM_Mapped_ZMWs",
  "NUM_HQ_ZMWS",
  "HQR_FP_BW",
  "HQR_FP_ZW",
  "HQR_FN_BW",
  "HQR_FN_ZW",
  "HQR_ADV",
  "HQR_RET"
)
description <- c(
  "Interesting ZMWs\nwith HQ or mapped reads",
  "Count of ZMWs\nwith good subreads",
  "Count of ZMWS\nwith HQ",
  "false positive\nbase weighted",
  "false positive\nZMW weighted",
  "false negative\nbase weighted",
  "false negative\nZMW weighted",
  "10,50,90 percentiles",
  "10,50,90 percentiles"
)
good         <- d[d$class==0,]
fn           <- d[d$class==1,]
fp           <- d[d$class==2,]
hqr          <- d[!is.na(d$hqSize),]
mapped       <- d[!is.na(d$aSize),]
d12345 <- d[d$subreads>0 & d$hqSize > 0,]

value = c(
  nrow(d),
  nrow(mapped),
  nrow(hqr),
  round(sum(fp$hqSize)/sum(hqr$hqSize),3),  # false pos
  round(nrow(fp)/nrow(hqr),3),  # false pos
  round(sum(fn$aSize)/sum(mapped$aSize),3),  # false neg
  round(nrow(fn)/nrow(mapped),3),  # false neg
  paste(qadv[2:4],collapse=","),
  paste(qret[2:4],collapse=",")
)

ideal <- c(
  "1M",
  "1M",
  "1M",
  " -> 0.0",
  " -> 0.0",
  " -> 0.0",
  " -> 0.0",
  "-eps,0,+eps",
  "-eps,0,+eps"
)

summary <- data.frame(metric,description,value,ideal)

write.csv(summary,csvfile)

grid.table(summary)

#-----------------------------------------------------

dev.off()

