nextflow run rtSeq_V1.4.groovy \
    --bam bwaAligned.bam \
    --gcCorrectionFile accessoryFiles/rtSeq/gcCalibration/Admera_HiSeqX.tab \
    --genome mm10 \
    --name sampleRTseq \
    --outdir sampleRTseq 
