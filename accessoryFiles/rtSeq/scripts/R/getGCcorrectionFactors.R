#Sys.setenv(RTSEQ_RSCRIPTDIR='/data/RDCO/kevbrick/GL_Sorting_Paper/finalPipe/output/')
#Sys.setenv(KBPIPEWORKDIR='/data/RDCO/kevbrick/GL_Sorting_Paper/')

## Load custom functions & H3 data
accessorydir <- paste0(Sys.getenv('NXF_PIPEDIR'),'/accessoryFiles/rtSeq/')

imgdir  <- './'

source(paste0(accessorydir,'/scripts/R/genericFunctions.R'))

library(data.table)
library(ggplot2)
library(ggpubr)

theme7point()

z <- fread('wholeGenome.DS.tab',header=FALSE)
names(z) <- c('cs','from','to','real','GC')
z$csType <- 'sex-CS'
z$csType[z$cs %in% paste0('chr',1:100)] <- 'autosomes'

## Get TRUE GC percent (original in 101bp windows)
z$GC <- round(z$GC/101*100);

## Get median cover per GC
zAutosomes <- z[z$csType == 'autosomes',]
df <-as.data.frame(aggregate(zAutosomes$real,by=list(GC=zAutosomes$GC),FUN=median))
names(df) <- c('GC','median')

## Express median coverage as a multiple of the average coverage
## Allow simple normalization on autosomes and on X/Y
df$median <- df$median / mean(df$median)

## Add median to original data for plot
zDF <- plyr:::join(z,df,by='GC')

zDF$real[zDF$csType == 'sex-CS'] <- zDF$real[zDF$csType == 'sex-CS']*1.5
zDF$norm<-zDF$real/zDF$median

############################################################ Plot before v after figure
#######################################################################################
loB4 <- quantile(zDF$real[zDF$csType == 'autosomes' & zDF$norm > 0],0.05)
hiB4 <- quantile(zDF$real[zDF$csType == 'autosomes' & zDF$norm > 0],0.95)

gBefore <- ggplot(zDF,aes(x=GC,y=real/101,group=GC)) +
  geom_boxplot(fill=alpha('forestgreen',.3),
               color=alpha('black',.7),
               outlier.size=.1,
               outlier.alpha=.3,
               lwd=.1) +
  coord_cartesian(ylim=c(loB4/101,hiB4/101)) +
  xlab('GC content (%)') + ylab('Read coverage (per bp)')+
  ggtitle('Before GC-normalization') + facet_wrap(~csType, ncol=2)

loAft <- quantile(zDF$norm[zDF$csType == 'autosomes' & zDF$norm > 0],0.05)
hiAft <- quantile(zDF$norm[zDF$csType == 'autosomes' & zDF$norm > 0],0.95)

gAfter <- ggplot(zDF,aes(x=GC,y=norm,group=GC)) +
  geom_boxplot(fill=alpha('firebrick',.3),
               color=alpha('black',.7),
               outlier.size=.1,
               outlier.alpha=.3,
               lwd=.1) +
  coord_cartesian(ylim=c(loAft,hiAft)) +
  xlab('GC content (%)') + ylab('GC-normalized coverage (per bp)') +
  ggtitle('After GC-normalization') + facet_wrap(~csType, ncol=2)

gOut <- ggarrange(gBefore,gAfter,
                  ncol=1,nrow=2,
                  labels=c('a','b'),
                  hjust=0,vjust=1,
                  font.label = list(size=8,font='bold'))

ggsave(getIMGname('rtSeq_GCCorrection_Before_v_After',saveLocation = imgdir, type = 'PNG'),
       plot = gOut,
       width = 4,
       height = 5)

ggsave(getIMGname('rtSeq_GCCorrection_Before_v_After',saveLocation = imgdir, type = 'PDF'),
       plot = gOut,
       width = 4,
       height = 5)

######################################################## Finalize and export normalization table
################################################################################################
df$loQ <- loAft
df$hiQ <- hiAft

fwrite(df, file='GCcorrectiontable.tab',append=FALSE,quote=FALSE,sep="\t")

############################################################### Get chromosome & genome averages
################################################################################################
zDF$OK <- FALSE
zDF$OK[zDF$norm > loAft & zDF$norm < hiAft] <- TRUE

## ONLY APPLY LIMITS TO AUTOSOMES
## OTHERWISE THIS BIASES SEX CS RT DOWN IN MALES
zDF$OK[!(zDF$cs %in% paste0("chr",1:100))] <- TRUE

zOK <- zDF[zDF$OK,]

csStat <- aggregate(zOK$norm,
                    by=list(cs=zOK$cs),
                    FUN=mean)

names(csStat) <- c('cs','mean')

csMedian <- aggregate(zOK$norm,
                    by=list(cs=zOK$cs),
                    FUN=median)

csStat$median <- csMedian$x

zOKAuto <- zOK[zOK$cs %in% paste0("chr",1:100),]
dfMean <- as.data.frame(rbind(csStat,
                              data.frame(cs='all',
                                         mean=mean(zOKAuto$norm),
                                         median=median(zOKAuto$norm))))

names(dfMean) <- c('cs','mean','median')
fwrite(dfMean, file='chromCorrectiontable.tab',append=FALSE,quote=FALSE,sep="\t")
