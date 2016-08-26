#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-m", "--mode", type="integer", default=1, help="run mode: 1 for all analysis; 2 for reads preprocessing only; 3 for alignment-post analysis  [default %(default)s]")
parser$add_argument("-t", "--thread", type="integer", default=8, help="number of threads [default %(default)s]")
parser$add_argument("-B", "--BatchMode", type="integer", default=4, help="gsnap's \"-B\" [default %(default)s]")
parser$add_argument("-n", "--nCore", type="integer", default=8, help="number of cores the machine have [default %(default)s]")
parser$add_argument("-a", "--trim", type="integer", default=0, help="if specified do reads trimming to INT bp")
parser$add_argument("-p", "--prefq", action="store_true", default=FALSE, help="if specified preserve processed fq")
parser$add_argument("-d", "--DB", default="/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap", help="database directory [default %(default)s]")
parser$add_argument("-s", "--splicesites", default="hg19_knownGene.splicesites", help="splicesite index [default %(default)s]")
parser$add_argument("-r", "--rRNA", default="human_rRNA", help="rRNA [default %(default)s]")
parser$add_argument("-g", "--genome", default="human_hg19", help="genome [default %(default)s]")
parser$add_argument("-e", "--geneModel", default="hg19.knownGene.RData", help="gene model [default %(default)s]")
parser$add_argument("outDir", help="output directory")
parser$add_argument("sampleID", help="sample's id")
parser$add_argument("fq1", help="fq1 file")
parser$add_argument("fq2", help="fq2 file")
args <- parser$parse_args()

dir.create(args$outDir, showWarnings = FALSE)

library("HTSeqGenie")
dbDir <- args$DB
#dbDir <- "/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap"
#dbDir <- "/lustre1/zeminz_pkuhpc/00.database/broad/bundle/2.8/hg19/gmap-gsnap"
rp <- list(
  ## aligner
  path.gsnap_genomes=dbDir,
  alignReads.genome=args$genome,
  alignReads.static_parameters=paste("-B ",args$BatchMode," --novelsplicing 1 -n 10 -i 1 -M 2 ",sep=""),
  alignReads.nbthreads_perchunk=args$thread, 
  num_cores=args$nCore,
  alignReads.splice_index=args$splicesites,
  ## trim
  trimReads.do = as.logical(args$trim),
  trimReads.length = args$trim,
  ## rRNA contamination
  detectRRNA.do = T,
  detectRRNA.rrna_genome = args$rRNA,
  ## gene model
  path.genomic_features=dbDir,
  countGenomicFeatures.gfeatures=args$geneModel,
  #for debug and test
  #max_nbchunks=2,
  alignReads.sam_id=args$sampleID,
  #markDuplicates.do = TRUE,
  #path.picard_tools = "/WPS/BP/zhenglt/01.bin/picard/current",
  alignReads.use_gmapR_gsnap = FALSE,
  ## input
  input_file=args$fq1,
  input_file2=args$fq2,
  paired_ends=TRUE,
  #quality_encoding="illumina1.8",
  analyzeVariants.do=FALSE,
  ## output
  save_dir=args$outDir,
  prepend_str=args$sampleID,
  overwrite_save_dir="erase",
  remove_processedfastq=(!args$prefq)
)

print(rp)

if(args$mode==1)
{
    resdir <- do.call(runPipeline, c(rp))
}else if(args$mode==2)
{
    initPipelineFromConfig(config_update=c(rp,remove_processedfastq=FALSE))
    preprocessReads()

}else if(args$mode==3)
{
    initPipelineFromSaveDir(save_dir=args$outDir,config_update=c(overwrite_save_dir="overwrite"))
    if (getConfig.logical("countGenomicFeatures.do")) {
        countGenomicFeatures()
     }
     if (getConfig.logical("coverage.do")) calculateCoverage()

     if (getConfig.logical("analyzeVariants.do")) analyzeVariants()

     #removeChunkDir()
     loginfo("Pipeline run successful.")
     invisible(getConfig("save_dir"))
}
