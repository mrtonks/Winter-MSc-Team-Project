#!/home/yuan03/R-3.2.0/bin/Rscript
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
icgc <- as.character(args[1])
sampleName <- as.character(args[2])
prefix <- as.character(args[3])
resultsFolder <- as.character(args[4])
burnIn <- as.integer(args[5])
iteration <- as.integer(args[6])
cellularity <- as.double(args[7])
inputFile <- paste0(prefix,".tsv")

dataFolder <- paste0(icgc, "/inputTmp/",sampleName)

#if (is.na(jobList)) {
#} else {
#    sampleList <- read.table(jobList, sep="\t", header = T,  stringsAsFactors=FALSE)$x
#}
ssm <- read.csv(paste0(dataFolder, "/", inputFile), header=T, sep='\t',na.strings=c("NA"))

pycloneFolder <- paste0(dataFolder, "/pyclone")
unlink(pycloneFolder, recursive = T, force = T)
dir.create(pycloneFolder, recursive = T)

maxSnv <- 5e3
if (nrow(ssm) > maxSnv) {
    ssm <- sample_n(ssm, maxSnv)
}
pycloneData <- data.frame(mutation_id = ssm$mutation_id, ref_counts = ssm$ref_counts,
var_counts = ssm$var_counts,
normal_cn = ssm$normal_cn, minor_cn =ssm$minor_cn, major_cn = ssm$major_cn )
write.table(pycloneData, file=paste0(pycloneFolder, "/pyclone_data.tsv"), quote=FALSE, sep='\t', row.names = F)
rm(pycloneData)

save.image("inputTmp/myfile")
