## Average transcript length normalization for count based methods

The count of RNA-seq reads which can be assigned to a gene depends on
a number of factors, including the abundance of that gene, the library
size, the length of the gene, sequence biases, the distribution of
fragment lengths and other factors. Comparing counts across samples and only
accounting for library size assumes that the gene-specific factors
other than gene abundance do not change across samples: for example
the length of the gene or the gene-specific effects of sequence
biases. These factors can also be included as model parameters into a
count based method, through the use of a gene-by-sample matrix of
normalization factors, or an offset matrix (on the log scale).

Here we show how to use the estimates of average transcript length for each
gene and each sample from the output of RSEM and Cufflinks as
normalization for count base methods DESeq2 and edgeR.

By gene length, we refer to a gene's average transcript
length, averaging with respect to the abundance of each
transcript. For example, if a gene has two transcripts, one of length
1000 which accounts for 1/4 of the abundance, and one of length 2000
which accounts for 3/4 of the abundance, the average transcript length, would be:


```r
0.25 * 1000 + 0.75 * 2000
```

```
## [1] 1750
```

## Experiment files

The files are RNA-seq runs from the UCSD Human Reference Epigenome
Mapping Project. Note that these are individuals runs from experiments
which might have had multiple runs per experiment.
The experiments are two samples of each: left ventricle and sigmoid colon.


```r
samples <- list.files("rsem_out")
condition <- scan("condition.txt",what="char")
ord <- order(condition)
(samples <- samples[ord])
```

```
## [1] "SRR577587" "SRR578627" "SRR577595" "SRR578635"
```

```r
(condition <- condition[ord])
```

```
## [1] "LV" "LV" "SC" "SC"
```

## RSEM output

The RSEM `.genes.results` files contain a column
`effective_length`. The effective length as defined by the RSEM
authors as "transcript length - mean fragment length + 1", and for the
gene results, the reported effective length is averaged over the
transcripts, weighting by the percent of total abundance, as described above.
All we have to do is collect this column over the samples.


```r
library(data.table)
rsem_list <- list()
for (i in seq_along(samples)) {
  cat(i)
  rsem_raw <- fread(paste0("rsem_out/",samples[i],"/",samples[i],".genes.results"))
  rsem_list[[i]] <- data.frame(eff_length=rsem_raw$effective_length, row.names=rsem_raw$gene_id)
}
```

```
## 1234
```

```r
head(rsem_list[[1]])
```

```
##                 eff_length
## ENSG00000000003    1917.31
## ENSG00000000005     786.96
## ENSG00000000419     891.14
## ENSG00000000457    5934.49
## ENSG00000000460    2617.80
## ENSG00000000938    1964.99
```

We combine these columns into one data frame, checking that the
`gene_id` is consistent.


```r
for (i in seq_along(samples)[-1]) {
  stopifnot(all(rownames(rsem_list[[i]]) == rownames(rsem_list[[1]])))
}
rsem <- do.call(cbind, rsem_list)
rsem <- as.matrix(rsem)
colnames(rsem) <- samples
head(rsem)
```

```
##                 SRR577587 SRR578627 SRR577595 SRR578635
## ENSG00000000003   1917.31   2058.40   1880.45   2073.04
## ENSG00000000005    786.96   1170.59    781.80    374.96
## ENSG00000000419    891.14    907.13    916.23    908.85
## ENSG00000000457   5934.49   3243.15   3251.06   3240.75
## ENSG00000000460   2617.80   1491.74   3568.30   2433.27
## ENSG00000000938   1964.99   2225.01   2290.79   2192.96
```

## Cufflinks output

For Cufflinks output files, we do not have the exact same
quantity, the average effective length, but the `isoforms.fpkm_tracking` file does
report the length of each isoform and the FPKM estimate for each
isoform. We can calculate the average gene length by averaging over
the isoforms, weighting by the ratio of each isoform's FPKM to the total
FPKM for all isoforms.


```r
library(dplyr)
cuff_list <- list()
for (i in seq_along(samples)) {
  cat(i)
  suppressWarnings({cuff_raw <- fread(paste0("cufflinks_out/",samples[i],"/isoforms.fpkm_tracking"))})
  res <- cuff_raw %>% group_by(gene_id) %>% summarize(avg_length=sum(length*FPKM)/sum(FPKM))
  cuff_list[[i]] <- data.frame(avg_length=res$avg_length, row.names=res$gene_id)
}
```

```
## 1234
```

```r
head(cuff_list[[1]])
```

```
##                 avg_length
## ENSG00000000003  2205.0650
## ENSG00000000005        NaN
## ENSG00000000419   996.7745
## ENSG00000000457  6364.0000
## ENSG00000000460  2653.6558
## ENSG00000000938  2316.1061
```

We combine these columns into one data frame, checking that the
`gene_id` is consistent.


```r
for (i in seq_along(samples)[-1]) {
  stopifnot(all(rownames(cuff_list[[i]]) == rownames(cuff_list[[1]])))
}
cuff <- do.call(cbind, cuff_list)
cuff <- as.matrix(cuff)
colnames(cuff) <- samples
head(cuff)
```

```
##                 SRR577587 SRR578627 SRR577595 SRR578635
## ENSG00000000003 2205.0650  2284.611  2222.429 2263.9102
## ENSG00000000005       NaN  1339.000  1339.000  677.3134
## ENSG00000000419  996.7745  1075.000  1075.000 1075.0000
## ENSG00000000457 6364.0000  4307.202  4670.106 4997.9242
## ENSG00000000460 2653.6558  1751.074  3824.666 2690.7123
## ENSG00000000938 2316.1061  2360.064  2507.561 2343.4121
```

## Genes with zero length

We have to remove gene's with zero length to continue.


```r
rsem.nz <- rsem[apply(rsem, 1, function(x) all(x > 0)),]
cuff.nz <- cuff[apply(cuff, 1, function(x) all(x > 0)),]
```

## Dividing out the geometric mean

As the count based methods have an intercept on the log scale, and 
the effect of gene length on counts is multiplicative, we can
divide each row by it's geometric median to produce a matrix which is
centered on zero in the log scale.


```r
norm.mat <- rsem.nz / exp(rowMeans(log(rsem.nz)))
head(norm.mat)
```

```
##                 SRR577587 SRR578627 SRR577595 SRR578635
## ENSG00000000003 0.9680983 1.0393382 0.9494867 1.0467303
## ENSG00000000005 1.0916744 1.6238477 1.0845165 0.5201462
## ENSG00000000419 0.9838251 1.0014782 1.0115246 1.0033771
## ENSG00000000457 1.5726358 0.8594325 0.8615287 0.8587966
## ENSG00000000460 1.0848411 0.6181912 1.4787374 1.0083702
## ENSG00000000938 0.9076811 1.0277912 1.0581767 1.0129864
```

## Input to DESeq2

Here we show the code for including these gene lengths in a
DESeq2 analysis. The `estimateSizeFactors` function corrects for
library size after taking into account the information in `norm.mat`.
Here we demonstrate using the normalization matrix on a matrix of
random counts, where the actual read counts for each gene and sample
would go.


```r
n <- nrow(norm.mat)
m <- length(samples)
library(DESeq2)
cts <- matrix(rpois(n*m,lambda=100),ncol=m)
coldata <- data.frame(condition=factor(condition), row.names=samples)
dds <- DESeqDataSetFromMatrix(cts, coldata, ~ condition)
dds <- estimateSizeFactors(dds, normMatrix=norm.mat)
```

```
## adding normalization factors which account for library size
```

```r
head(normalizationFactors(dds))
```

```
##      SRR577587 SRR578627 SRR577595 SRR578635
## [1,] 0.9386213 1.0654183 0.9386604 1.0653221
## [2,] 1.0584347 1.6645950 1.0721505 0.5293849
## [3,] 0.9538692 1.0266083 0.9999909 1.0211989
## [4,] 1.5247517 0.8809983 0.8517053 0.8740504
## [5,] 1.0518095 0.6337035 1.4618764 1.0262807
## [6,] 0.8800437 1.0535816 1.0461111 1.0309790
```

## Input to edgeR

One can provide the normalization factors as an offset to edgeR, which
needs to be on the natural log scale.  The offset matrix, like the
`normalizationFactors` in DESeq2, needs to account for library size.


```r
library(edgeR)
o <- log(calcNormFactors(cts/norm.mat)) + log(colSums(cts/norm.mat))
y <- DGEList(cts)
y$offset <- t(t(log(norm.mat)) + o)
```

## Exploring RSEM and Cufflinks estimates

As a caveat, we note that the RSEM column we used is the effective
length (where the average fragment length has been subtracted),
whereas the Cufflinks column we used was the transcript length alone.
This might explain some discrepancy across different software.


```r
idx <- intersect(rownames(rsem), rownames(cuff))
length(idx)
```

```
## [1] 60083
```

```r
rsem <- rsem[idx,]
cuff <- cuff[idx,]
stopifnot(all.equal(rownames(rsem), rownames(cuff)))
colnames(rsem) <- condition
colnames(cuff) <- condition
```


```r
library(rafalib)
mypar(2,2)
for (i in 1:4) {
  suppressWarnings(plot(rsem[,i], cuff[,i],log="xy",cex=.3,main=samples[i],
       xlab="rsem", ylab="cufflinks"))
  abline(0,1,col="blue")
}
```

![plot of chunk scatter](figure/scatter-1.png) 


```r
suppressWarnings(pairs(rsem,log="xy",cex=.1))
```

![plot of chunk rsem_pairs](figure/rsem_pairs-1.png) 


```r
suppressWarnings(pairs(cuff,log="xy",cex=.1))
```

![plot of chunk cuff_pairs](figure/cuff_pairs-1.png) 


```r
mypar()
rsem_a <- rowMeans(rsem)
rsem_m <- log2(rowMeans(rsem[,condition == "SC"])/rowMeans(rsem[,condition == "LV"]))
plot(rsem_a, rsem_m, log="x", xlab="A", ylab="M", main="rsem", cex=.3)
```

```
## Warning in xy.coords(x, y, xlabel, ylabel, log): 2444 x values <= 0
## omitted from logarithmic plot
```

```r
abline(h=0,col="blue")
```

![plot of chunk rsem_ma](figure/rsem_ma-1.png) 


```r
mypar()
cuff_a <- rowMeans(cuff)
cuff_m <- log2(rowMeans(cuff[,condition == "SC"])/rowMeans(cuff[,condition == "LV"]))
plot(cuff_a, cuff_m, log="x", xlab="A", ylab="M", main="cufflinks", cex=.3)
abline(h=0,col="blue")
```

![plot of chunk cuff_ma](figure/cuff_ma-1.png) 


```r
mypar()
plot(rsem_m, cuff_m, xlab="rsem", ylab="cufflinks", main="M vs M", cex=.3)
legend("bottomright",legend=paste0("pearson=",round(cor(rsem_m,cuff_m,use="complete"),2)),inset=.05)
abline(0,1,col="blue")
```

![plot of chunk m_vs_m](figure/m_vs_m-1.png) 


```r
mypar()
plot((rsem_a + cuff_a)/2, rsem_m - cuff_m, xlab="(rsem A + cufflinks A)/2", ylab="rsem M - cufflinks M",
     main="M minus M", cex=.3, log="x")
abline(h=0,col="blue")
```

![plot of chunk m_minus_m](figure/m_minus_m-1.png) 
