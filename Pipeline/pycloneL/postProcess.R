library(dplyr)
library(mcclust)
compute_mpear_label <- function(label_traces){
    ltmat <- as.matrix(label_traces)
    ltmat <- ltmat + 1
    psm <- comp.psm(ltmat)
    mpear <- maxpear(psm)
    mpear_label <- mpear$cl
    return(mpear_label)
}

gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}

#' bsxfun {pracma} with single expansion (Matlab style)
#' @param func the function used by bsxfun
#' @param x a matrix
#' @param y a vector need to be expanded
#' @param  expandByRow applies only when x is a square matrix
#' @return value of func
bsxfun.se <- function(func, x, y, expandByRow=TRUE) {
    
    if(length(y) == 1) return(pracma::arrayfun(func, x, y)) else
    stopifnot(nrow(x) == length(y) || ncol(x) == length(y))
    
    expandCol <- nrow(x) == length(y)
    expandRow <- ncol(x) == length(y)
    
    if(expandCol & expandRow & expandByRow) expandCol <- FALSE
    if(expandCol & expandRow & !expandByRow) expandRow <- FALSE
    
    # repeat row (if dim2expand = 1, then length(y) = ncol(x))
    if(expandRow) y.repmat <- matrix(rep(as.numeric(y), each=nrow(x)), nrow=nrow(x))
    
    # repeat col (if dim2expand = 2, then length(y) = nrow(x))
    if(expandCol) y.repmat <- matrix(rep(as.numeric(y), ncol(x)), ncol=ncol(x))
    
    pracma::bsxfun(func, x, y.repmat)
}

# matlab style helper
repmat = function(X,m,n){
    ##R equivalent of repmat (matlab)
    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

#' Compute log(sum(exp(x),dim)) while avoiding numerical underflow
#' @param x a matrix
#' @param margin used for apply
#' @return log(sum(exp(x),dim)) a matrix of the sample size as x
logsumexp <- function(x, margin=1) {
    
    if ( ! is.matrix(x) ) {
        x <- as.matrix(x)
    }
    
    # subtract the largest in each column
    y <- apply(x, margin, max)
    
    if (nrow(x) == ncol(x)) {
        x <- bsxfun.se("-", x, y, expandByRow = F)
    } else {
        x <- bsxfun.se("-", x, y)
    }
    
    s <- y + log(apply(exp(x), margin, sum))
    
    i <- which(!is.finite(s))
    
    if(length(i) > 0) s[i] <- y[i]
    
    s
}

Assign <- function(x, centers, s) {
    n <- length(x)
    k <- length(centers)
    logRho <- array(0, dim= c(n ,k))
    
    for (ii in 1:k) {
        logRho[,ii] = bsxfun.se("-", -(x-centers[ii])^2/(2*s[ii]), log(s[ii]))
    }
    
    if (n==k) {
        logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)  # 10.49
    } else {
        logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1)) # 10.49
    }
    
    R <- exp(logR)
    return(list(R=R, logR=logR))
}

GetMultFromCcf <- function(bn, dn, ccf, major_cn, minor_cn, purity, epi = 1e-3) {
    total_cn = major_cn + minor_cn
    z = (1-purity)*2 + purity*total_cn
    k = max(major_cn)
    multPool <- seq(0, k, by = 1)
    aa = repmat(as.matrix(purity * ccf * (1 - epi) / z), 1, k+1)
    bb = repmat(as.matrix(epi * ( z - purity * (1-ccf) * total_cn)), 1, k+1)
    pp = bsxfun.se("*", aa, multPool) + bb
    ll = bsxfun.se("*", log(pp), bn) + bsxfun.se("*", log(1-pp), dn-bn)
    return(multPool[apply(ll, 1, which.max)])
}

GetCcfFromLabel <- function(ccfTraces, idx, label) {
    
    allData <- data.frame(mutation_id = idx, cluster_id = label,
    ccf = apply(ccfTraces, 1, median, na.rm = T), stringsAsFactors = F)
    idxDf <- data.frame(mutation_id = idx, stringsAsFactors = F)
    
    tt <- table(allData$cluster_id)
    clusterMean <- vector(mode = "numeric", length = length(tt))
    clusterSd <- clusterMean
    for (ii in seq_along(tt)) {
        
        dataIdx <- which(allData$cluster_id %in% as.integer(names(tt[ii])))
        
        clusterMean[ii] <- if (tt[ii]/sum(tt) > 0.01 &
        length(dataIdx) > 1
        ) {
            median(c( ccfTraces[dataIdx, ]), na.rm = T )
        } else {NA}
        
        clusterSd[ii] <- if (tt[ii]/sum(tt) > 0.01 &
        length(dataIdx) > 1
        ) {
            if (sd(c( ccfTraces[dataIdx, ]),  na.rm = T) ==0) {
                1e-20
            } else { sd( c( ccfTraces[dataIdx, ]), na.rm = T)}
        } else {NA}
    }
    
    clusterDf <- data.frame(cluster_id=as.integer(names(tt)),
    average_ccf = clusterMean,
    lower_95_ci = clusterMean - 2*clusterSd,
    upper_95_ci = clusterMean + 2*clusterSd)
    
    allData <- left_join(allData, clusterDf, by="cluster_id" )
    clusterDf <- filter(clusterDf, !is.na(average_ccf) & !is.na(lower_95_ci))
    ## Reassign low support data
    lowSupportData <- dplyr::filter(allData, is.na(average_ccf) | is.na(lower_95_ci))
    allData <- dplyr::filter(allData, !is.na(average_ccf))
    lowSupportDataFlag <- F
    if (nrow(lowSupportData) > 0) {
        lowSupportDataFlag <- T
        lowSupportR <- Assign(lowSupportData$ccf, clusterDf$average_ccf,
        ((clusterDf$upper_95_ci - clusterDf$average_ccf)/2)^2)$R
        lowSupportLabel <- apply(lowSupportR, 1, which.max)
        lowSupportData$cluster_id <- clusterDf$cluster_id[lowSupportLabel]
        lowSupportData$average_ccf <- clusterDf$average_ccf[lowSupportLabel]
        lowSupportData$lower_95_ci <- clusterDf$lower_95_ci[lowSupportLabel]
        lowSupportData$upper_95_ci <- clusterDf$upper_95_ci[lowSupportLabel]
        allData <- rbind(allData, lowSupportData)
        rm(lowSupportR)
    }
    
    allData$cluster_id <- match(allData$average_ccf, sort(unique(allData$average_ccf)))
    return(left_join(idxDf, allData, by = "mutation_id"))
}

load("inputTmp/myfile")
id <- Reduce(rbind, strsplit(as.character(ssm$mutation_id), ":", fixed = T), c())

## write results
if (dir.exists(resultsFolder)) {
    "Folder exists and continue."
} else {
    dir.create(resultsFolder, recursive = T)
}


# load trace and mpear label
traceFile <- dir(paste0(pycloneFolder, "/trace"),
                 pattern = "cellular_prevalence", full.names = T)
paramsTrace <- read.delim(traceFile, stringsAsFactors = F, header = F)
idx <- as.character(paramsTrace[1,])
#burnIn <- 20
#paramsTrace <- as.matrix( paramsTrace[-1:-burnIn, ]) #!!!!!!!!need to be change
paramsTrace <- as.matrix(paramsTrace[-1:-(burnIn+1), ])
class(paramsTrace) <- "numeric"

mpearFile <- paste0(pycloneFolder, "/pyclone_mpear.tsv")

if (file.exists(mpearFile)) {
  mpear <- read.delim(paste0(pycloneFolder, "/pyclone_mpear.tsv"), stringsAsFactors = F)
} else {

  labelTrace <- read.delim(paste0(pycloneFolder, "/trace/labels.tsv.bz2"), stringsAsFactors = F)
  #burnIn <- 20
  #mpearLabels <- compute_mpear_label(labelTrace[-1:-burnIn,])
  mpearLabels <- compute_mpear_label(labelTrace[-1:-burnIn,])
  mpear <- data.frame(mutation_id = colnames(labelTrace), cluster_id = mpearLabels)
  write.table(mpear, file = paste0(pycloneFolder, '/pyclone_mpear.tsv'),
              row.names = F, sep = "\t", quote = F)
  rm(labelTrace)
}

allData <- GetCcfFromLabel(t(paramsTrace), idx, mpear$cluster_id)
rm(paramsTrace)

# load tsv
sampleTsvFile <- paste0(pycloneFolder, "/pyclone_data.tsv")
sampleTsv <- read.delim(sampleTsvFile, stringsAsFactors = F)
sampleTsv <- mutate(sampleTsv, vaf = var_counts/(ref_counts+var_counts))
allData <- left_join(allData, sampleTsv, by = "mutation_id")
rm(sampleTsv)

mutAssign <- data.frame(mutation_id = ssm$mutation_id, chr = id[,3], pos = id[,4])

allData <- left_join(allData, mutAssign, by = "mutation_id")
allData <- dplyr::filter(allData, !is.na(average_ccf))

rm(mutAssign)

# co-assignments matrix file
#labelFile <- dir(paste0(pycloneFolder, "/trace/"), pattern = "labels", full.names = T)
#labelTrace <- read.delim(labelFile, stringsAsFactors = F)
#labelTrace <- as.matrix(labelTrace[-1:-(burnIn+1), ]) + 1
#labelTrace <- as.matrix(labelTrace[-1:-(burnIn+1), ]) + 1
#psm <- comp.psm(labelTrace)
#fn <- paste0(resultsFolder, "/",
#             sampleName, "_coassignment_probabilities.txt")
#write.table(psm, file = fn, sep = "\t", row.names = F, quote = F)
#shellCommand <- paste0("gzip -f ", fn)
#system(shellCommand, intern = TRUE)
#rm(psm)

# index file
#index <- allData[, c("chr", "pos")]
#index$col <- seq_along(allData$mutation_id)
#fn <- paste0(resultsFolder, "/",
#             sampleName, "_index.txt")
#write.table(index, file = fn, sep = "\t", row.names = F, quote = F)
#shellCommand <- paste0("gzip -f ", fn)
#system(shellCommand, intern = TRUE)

# cluster certainty file
clusterCertainty <- subset(allData,
                           select = c("chr", "pos", "cluster_id",
                                      "average_ccf", "lower_95_ci", "upper_95_ci"))
clusterCertainty <- rename(clusterCertainty, most_likely_assignment = cluster_id)
clusterCertainty$most_likely_assignment <-
  match(clusterCertainty$most_likely_assignment,
        sort(unique(clusterCertainty$most_likely_assignment)))

tmp11 <- clusterCertainty[, c("chr", "pos", "most_likely_assignment","average_ccf")]
tmp11 <- rename(tmp11, cluster = most_likely_assignment)
tmp11 <- mutate(tmp11,average_ccf1=as.numeric(average_ccf) * cellularity)
tmp11 <- rename(tmp11, proportion = average_ccf1)
dir.create(paste0(resultsFolder, "pyclone/mutation_assignments"), recursive = T)
fn <- paste0(resultsFolder, "pyclone/mutation_assignments/",
             prefix, "_mutation_assignments.txt")
write.table(tmp11, file = fn, sep = "\t", row.names = F, quote = F)
#shellCommand <- paste0("gzip -f ", fn)
#system(shellCommand, intern = TRUE)

# multiplicity
allData <- mutate(allData, multiplicity = GetMultFromCcf(bn = ssm$var_counts,
                                             dn = var_counts + ref_counts,
                                             ccf = ccf,
                                             major_cn = major_cn,
                                             minor_cn = minor_cn,
                                             purity = cellularity))

mult <- allData[, c("chr", "pos", "multiplicity")]
mult$tumour_copynumber <- allData$major_cn+allData$minor_cn
dir.create(paste0(resultsFolder, "pyclone/multiplicity"), recursive = T)
fn <- paste0(resultsFolder, "pyclone/multiplicity/",
             prefix, "_multiplicity.txt")
write.table(mult, file = fn, sep = "\t", row.names = F, quote = F)
#shellCommand <- paste0("gzip -f ", fn)
#system(shellCommand, intern = TRUE)

# subclonal_structure file
tmp1 <- as.data.frame(table(clusterCertainty$most_likely_assignment), stringsAsFactors = F)
tmp2 <- as.data.frame(table(clusterCertainty$average_ccf), stringsAsFactors = F)
tmp <- left_join(tmp1, tmp2, by ="Freq")
delx <- duplicated(tmp$Var1.x)
dely <- duplicated(tmp$Var1.y)
dels <- xor(delx,dely)
delI <- which(dels)
if (!length(delI)==0){
    tmp <- tmp[-delI,]
}
tmp <- rename(tmp, cluster = Var1.x, n_ssms = Freq, proportion = Var1.y)
tmp <- mutate(tmp, proportion = as.numeric(proportion) * cellularity)
dir.create(paste0(resultsFolder, "pyclone/subclonal_structure"), recursive = T)
fn <- paste0(resultsFolder, "pyclone/subclonal_structure/",
             prefix, "_subclonal_structure.txt")
write.table(tmp, file = fn, sep = "\t", row.names = F, quote = F)
#shellCommand <- paste0("gzip -f ", fn)
#system(shellCommand, intern = TRUE)

# save sample summary
#allData$purity = ssm$purity
#fn <- paste0(pycloneFolder, "/", sampleName, "_pyclone_results_table.csv")
#write.csv(allData, file = fn , row.names = F)

# graph summary
#fn = paste0(resultsFolder, "/",
#           sampleName, "_pyclone_results_summary.pdf")
#pdf(fn, width=8, height=8)
#myColors <- gg_color_hue(n_distinct(allData$cluster_id))
#par(mfrow=c(2,2))
#plot(allData$ccf, allData$vaf, col = myColors[allData$cluster_id],
#     xlab = "cancer cell fraction", ylab = "variant allele frequecy",
#     main = "ccf vs vaf (colored by cluster memebership)")
#hist(allData$ccf, density=20, breaks=20, prob=TRUE,
#     main = "ccf histogram",
#     xlab = "cancer cell fraction")
#clusterSize <- table(allData$average_ccf)/nrow(allData)
#names(clusterSize) <- as.character(format(round(as.numeric(names(clusterSize)), 2), nsmall = 2))
#tmp1 <- as.data.frame(table(allData$average_ccf))
#tmp2 <- as.data.frame(table(allData$cluster_id))
#tmp3 <- left_join(tmp1, tmp2, by ="Freq")
#barplot(clusterSize, las = 2, col = myColors[tmp3$Var1.y], xlab = "cluster mean", ylab="mutation proportions",
#        main = "cluster sizes")
#dev.off()

