rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)

args <- commandArgs(T)
path <- as.character(args[1])
name <- as.character(args[2])
num1 <- as.integer(args[3])
num2 <- as.integer(args[4])
num3 <- as.integer(args[5])
num4 <- as.integer(args[6])
num5 <- as.integer(args[7])
registerDoParallel(cores=num2)

set.seed(num1)
numOfClusterPool <- 1:num4
numOfRepeat <- num3
resultFolder1 =paste0(path , "/results/ccube/multiplicity")
dir.create(resultFolder1, recursive = T)
resultFolder2 =paste0(path, "/results/ccube/mutation_assignments")
dir.create(resultFolder2, recursive = T)
resultFolder3 =paste0(path, "/results/ccube/subclonal_structure")
dir.create(resultFolder3,recursive = T)
read_path = paste0(path,"/inputTmp/",name,"/",name,".tsv")
mydata <- read.delim(read_path, stringsAsFactors = F)
#f <- function(start_time) {
# start_time <- as.POSIXct(start_time)
# dt <- difftime(Sys.time(), start_time, units="secs")
# format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
#}
time1<-Sys.time()
results <- RunCcubePipeline(sampleName = name,
                            writeOutput=F,
                            ssm = mydata, numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                            runAnalysis = T, runQC = F, multiCore = T,maxiter=num5,
                            basicFormats = F, allFormats = F, returnAll = T)

#change

res=results$res
ssm=results$ssm
sampleName = name
uniqLabels <- unique(res$label)
id <- do.call(rbind, strsplit(as.character(ssm$mutation_id), ":", fixed = T))

{
  mult <- data.frame(chr = id[,3], pos = id[,4])
  mult$tumour_copynumber <- ssm$major_cn+ssm$minor_cn
  mult$multiplicity <- ssm$ccube_mult
  fn <- paste0(resultFolder1, "/",
               sampleName, "_multiplicity.txt")
  write.table(mult, file = fn, sep = "\t", row.names = F, quote = F)
  rm(mult)
}

{ cellularity <- unique(ssm$purity)
    clusterCertainty <- as.data.frame(table(res$label), stringsAsFactors = F)
    clusterCertainty <- dplyr::rename(clusterCertainty, cluster = Var1, n_ssms = Freq)
    clusterCertainty$proportion <- res$full.model$ccfMean[as.integer(clusterCertainty$cluster)] * cellularity
    clusterCertainty$cluster <- seq_along(uniqLabels)
    mutAssign <- data.frame(chr = id[,3], pos = id[,4])
    if (length(uniqLabels) == 1) {
        mutAssign$cluster = 1
    } else {
        mutAssign$cluster <- res$label
    }
    mutAssign <- mutAssign[order(mutAssign[,3]),]
    k =1
    m1=list()
    m4=list()
    for(i in clusterCertainty$n_ssms){
        m2 <- matrix(clusterCertainty$proportion[k],ncol=i)
        m1 <- append(m1,m2)
        m3 <- matrix(res$full.model$ccfMean[as.integer(clusterCertainty$cluster)][k],ncol=i)
        m4 <- append(m4,m3)
        k=k+1
    }
    mutAssign$proportion <- m1
    mutAssign$ccfmean <- m4
    mutAssign <- as.matrix(mutAssign)
    fn <- paste0(resultFolder2, "/",
    sampleName, "_mutation_assignments.txt")
    write.table(mutAssign, file = fn, sep = "\t", row.names = F, quote = F)
    rm(mutAssign)
    rm(m1)
}
{
  cellularity <- unique(ssm$purity)
  clusterCertainty <- as.data.frame(table(res$label), stringsAsFactors = F)
  clusterCertainty <- dplyr::rename(clusterCertainty, cluster = Var1, n_ssms = Freq)
  clusterCertainty$proportion <- res$full.model$ccfMean[as.integer(clusterCertainty$cluster)] * cellularity
  clusterCertainty$cluster <- seq_along(uniqLabels)
  fn <- paste0(resultFolder3, "/",
               sampleName, "_subclonal_structure.txt")
  write.table(clusterCertainty, file = fn, sep = "\t", row.names = F, quote = F)
  
}
sign <-  paste0(sampleName," finished")
print(sign)

writepath <- paste0(path,"/results/ccube/runtimes.tsv")
time2 <- Sys.time()-time1
cat(c(name,time2*60,"\n"), file= writepath,append= TRUE)
print(time2)
rm(time1)
rm(time2)
