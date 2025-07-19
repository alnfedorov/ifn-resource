library("tximport")
library("DESeq2")

args = commandArgs(TRUE)

SAMPLES = args[1]
TX2GROUP = args[2]
BASELINE = args[3]
LOG2FC = as.numeric(args[4])
PADJ = as.numeric(args[5])
SAVETO = args[6]
COMPARISONS = unlist(strsplit(args[7], "$", fixed = TRUE))

########################################################################################################################
# Prepare samples table
samples = read.csv(SAMPLES, header = TRUE, row.names = "file")
paths = rownames(samples)
names(paths) = rownames(samples)

# Load mapping: transcript ID -> group ID
TX2GROUP = read.csv(TX2GROUP, sep = "\t", header = FALSE)
TX2GROUP = TX2GROUP[, 1:2]
colnames(TX2GROUP) = c("tx", "gene_id") # For compatibility with tximport

# Import Salmon estimations
txi <- tximport(paths, type = "salmon", tx2gene = TX2GROUP)
stopifnot(rownames(samples) == colnames(txi$counts))
########################################################################################################################
# DE analysis
samples$condition = relevel(factor(samples$condition), ref = BASELINE)
dds = DESeqDataSetFromTximport(txi, colData = samples, design = ~donor + condition)
dds = estimateSizeFactors(dds)

# Filter lowly expressed transcript groups
min_samples = 3
min_counts = 20
cnts = counts(dds, normalize=FALSE)

keep = rep(FALSE, nrow(cnts))
for (condition in unique(samples$condition)) {
    condition = rownames(samples)[samples$condition == condition]
    condition = cnts[, condition, drop = FALSE]
    keep = keep | (rowSums(condition >= min_counts) >= min_samples)
}
print(paste("Total transcript groups with low expression:", sum(!keep), "( out of", nrow(dds), ")"))
dds = dds[keep,]

# Run DE analysis
dds = DESeq(dds)
for (coef in COMPARISONS) {
    name = paste("condition", coef, "vs", BASELINE, sep = "_")
    res = results(dds, name = name, alpha = PADJ, lfcThreshold = 0, altHypothesis = 'greaterAbs')
    res.shrunk = lfcShrink(dds, coef = name, res = res, type = "apeglm")
    res.shrunk$padj[is.na(res.shrunk$padj)] = 1

    name = paste(coef, "vs", BASELINE, sep = "_")
    name = paste(name, "csv.gz", sep = ".")
    name = file.path(SAVETO, name)
    write.csv(res.shrunk, gzfile(name))
}
########################################################################################################################
# Variance stabilizing (vs) transformation and rlog
if (BASELINE == "mock") {
    vsd = varianceStabilizingTransformation(dds, blind = TRUE)
    name = file.path(SAVETO, "vsd.csv.gz")
    write.csv(assay(vsd), gzfile(name))

    rld = rlog(dds, blind = TRUE)
    name = file.path(SAVETO, "rld.csv.gz")
    write.csv(assay(rld), gzfile(name))
}
