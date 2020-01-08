## Load custom functions & H3 data
source(paste0(Sys.getenv('NXF_PIPEDIR'),'/accessoryFiles/rtSeq/scripts/R/genericFunctions.R'))

library(png)
library(ggpubr)
library(ggplot2)
library(data.table)

theme7point()

covData <- fread('GCdata.tab',header=FALSE)
names(covData) <- c('GC','raw','s1','coverage')

gB4    <- ggplot(covData) + 
  geom_boxplot(fill=alpha('firebrick',.3),
               aes(x=GC,y=raw,group=GC),
               outlier.size  = 0.1,
               outlier.alpha  = 0.1,
               lwd=.2) + 
  coord_cartesian(xlim=c(0,100),
                  ylim=c(0,quantile(covData$raw,.999)),
                  expand = FALSE) + 
  xlab('GC content (%; 101bp wins)') + 
  ylab('Raw Coverage')

gAfter <- ggplot(covData) + 
  geom_boxplot(fill=alpha('forestgreen',.3),
               aes(x=GC,y=coverage,group=GC),
               outlier.size  = 0.1,
               outlier.alpha  = 0.1,
               lwd=.2) + 
  coord_cartesian(xlim=c(0,100),
                  ylim=c(0,quantile(covData$coverage,.999)),
                  expand = FALSE) + 
  xlab('GC content (%; 101bp wins)') + 
  ylab('Adjusted Coverage')

gN <- ggplot(covData) + 
  geom_histogram(aes(x=GC, y = (..count..)/sum(..count..)*100),
                 binwidth = 1,
                 color='darkorange',fill='grey68',
                 lwd=.2) +
  coord_cartesian(xlim=c(0,100),
                  expand = FALSE) + 
  xlab('GC content (%; 101bp wins)') + 
  ylab('Bins (%)')

ggarrange(gB4,gAfter,gN,align='v',ncol=1,nrow=3)

ggsave('./gcCorrection.png',width=5,height=5)
ggsave('./gcCorrection.pdf',width=5,height=5)