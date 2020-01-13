## Inferring DNA replication timing from WGS (RT-Seq)

### Sample execution:

nextflow run rtSeq_V1.4.groovy \
    --bam bwaAligned.bam \
    --gcCorrectionFile accessoryFiles/rtSeq/gcCalibration/Admera_HiSeqX.tab \
    --genome mm10 \
    --name sampleRTseq \
    --outdir sampleRTseq 

### Command line args:
arg                | required | detail
------------------ | -------- | --------------------------------
--bam              | *        | BWA aligned BAM file
--genome           | *        | Genome name (mm10 / hg38)
--gcCorrectionFile | *        | GC calibration file (choose one from accessoryFiles/rtSeq/gcCalibration)
..                 |          | If omitted, pipeline will generate a calibration file from the BAM
..                 |          | provided. This is a time-consuming step. 
--name             | *        | Output file name stem
--outdir           |          | Output folder (default = ".")


### Bash environment variables: 
$NXF_PIPEDIR \
Path to a the folder containing the pipeline .groovy file and the accessoryFiles folder. 

### Requirements: 
R (3.5.2) : NOT compatible with later versions of R \
BEDtools (2.27.1) \
Picard tools (2.9.2) \
SAMtools (1.9) \

### R packages: 
data.table \
dplyr \
extrafont \
factoextra \
ggplot2 \
ggpmisc \
ggpubr \
grid \
gridExtra \
lsr \
numform \
pROC \
plyr \
png \
preprocessCore \
reshape2 \
tictoc 






