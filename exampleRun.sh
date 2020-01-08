nextflow run -work-dir rtseq_work rtSeq_V1.3.groovy --bam bwaAligned.bam --skipAlignment true --gcCorrectionFile GCcorrectiontable.tab --name testRTseq  --project testRTseq  --threads  16  --sample_name testSample1 --library testLib --rundate  20200101 --outdir rtseq  --mem 32G --genome mm10 -with-dag rtseq.nextflowDAG.dot -with-report rtseq.nextflowPipelineReport.html -with-trace rtseq.nextflowPipelineTrace.html -with-timeline rtseq.nextflowPipelineTimeline.html 