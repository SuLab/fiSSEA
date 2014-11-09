# Hackathon R script
library(VariantAnnotation)
library(mygene)
#Extract CDS from gene model

txdb <- makeTxDbFromMyGene("GDF11", scopes="symbol", species="human")
chrom <- sapply(as.list(seqnames(cds(txdb))), as.integer)
cds.start <- start(cds(txdb))
cds.end <- end(cds(txdb))
#ranges <- paste("chr", chrom, ":", cds.start, "-", cds.end, sep="")
## cds.start and cds.end input to VCF query need optimization
param <- GRanges(seqnames=chrom, ranges=IRanges(start=cds.start, end=cds.start))
tab <- TabixFile("SWGR_GDF11.vcf.gz")
vcf_rng <- readVcf(tab, "hg19", param=param)
vcf_df <- rowData(vcf_rng)
chrom <- sapply(as.list(seqnames(vcf_df)), as.integer)
vcf.start <- start(vcf_df)[1]
vcf.end <- tail(end(vcf_df), n=1)

ranges <- paste("chr", chrom[1], ":", vcf.start, "-", vcf.end, sep="")
## the [1] above is a quick hack-optimize

# returns data.frame of id, mutation assessor score, and gene name from ranged queries.
getMAS <- function(genomic.position){
  # make query
  res <- queryVariant(genomic.position, return.as="records")$hits$hits
  # retrieve HGVS id
  id <- unlist(lapply(res, function(i) i$`_id`))
  # retrieve mutation assessor score
  ma <- unlist(lapply(res, function(i) i$`_source`$dbnsfp$mutationassessor$score))
  # retrieve European allele frequencies
  #af <- unlist(lapply(res, function(i) i$`_source`$dbnsfp$1000gp1$eur_af))
  # retrieve gene name
  gene <- unlist(lapply(res, function(i) i$`_source`$dbnsfp$genename))
  # create boolean vector of results which have mutation assessor scores
  scores.exist <- unlist(lapply(res, function(i) !is.null(i$`_source`$dbnsfp$mutationassessor$score)))
  # subset vector of genes in which scores exist 
  gene <- gene[scores.exist]
  id <- id[scores.exist]
  # create data.frame
  df <- data.frame(HGVS=id, mutation.assessor.score=ma, gene=gene)
}

merge.df <- function(ranges){
  df <- lapply(ranges, getMAS)
  merged <- do.call(rbind, df)
  merged
}
merged <- merge.df(ranges)

getAllGenes <- function(){
  i <- 0
  genes <- list()
  while (i < 21000){
    q <- query(q="type_of_gene:protein-coding", fields="entrezgene", species="human", skip=i, size=i+1000)
    genes <- append(genes, q$hits$`_id`)
    i <- i + 1000}
  genes <- unlist(unique(genes))
  
}
