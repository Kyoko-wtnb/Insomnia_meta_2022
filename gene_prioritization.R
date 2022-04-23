##### Insomnia manuscript gene prioritization #####
#
# 13 Oct 2021: updated for sharing
# 10 Nov 2021: finalized
###################################################

setwd("/path/to/example")
library(data.table)
library(dplyr)
library(GenomicRanges)

##### create gene mapping summary #####
for(phe in c("ISM")){
  print(paste0("Phenotype: ", phe))
  ### set variable
  gwas_dir <- paste0(phe, "/gwas")
  cred_dir <- paste0(phe, "/cred")
  cred_file <- paste0(phe, "/finemap/cred_snps.txt.gz")
  coloc_file <- paste0(phe, "/coloc/coloc_eqtls.txt")
  out_dir <- phe
  
  min_PIP = 0.1
  coloc_nSNPs = 10
  max_gwasP = 5e-8
  min_pp4 = 0.9
  ensg <- fread("data/ENSG.genes.txt", data.table = F) %>% filter(gene_biotype=="protein_coding")
  ex <- c()
  if(file.exists(paste0(phe, "/exclude_loci.txt"))){
    ex <- fread(paste0(phe, "/exclude_loci.txt"), header=F)[[1]]
  }
  
  ### map loci
  ### This is because credible SNPs might be outside of the loci defined by LD clump
  ### Assign loci based on the credible SNPs to the closest LD clumped loci
  loci = fread(paste0(gwas_dir, "/GenomicRiskLoci.txt"), data.table = F)
  cred_loci = fread(paste0(cred_dir, "/GenomicRiskLoci.txt"), data.table = F)
  loci_gr = with(loci, GRanges(seqnames = chr, IRanges(start=start, end=end), idx=GenomicLocus))
  cred_loci_gr = with(cred_loci, GRanges(seqnames = chr, IRanges(start=start, end=end), idx=GenomicLocus))
  nearest_loci = nearest(cred_loci_gr, loci_gr, ignore.strand=T)
  cred_loci_map = data.frame(origin=cred_loci$GenomicLocus, mapped_loci=nearest_loci)
  rm(loci_gr, cred_loci_gr, nearest_loci)
  
  ### credible SNPs
  credSNPs <- fread(cmd=paste0("gzip -cd ", cred_file), data.table = F) %>% filter(!locus %in% ex)
  credSNPs <- credSNPs %>% filter(prob>min_PIP)
  credSNPs$uniqID <- sub("_", ":", credSNPs$uniqID)
  snps <- fread(paste0(cred_dir, "/snps.txt"), data.table = F) %>% filter(uniqID %in% credSNPs$uniqID)
  snps$GenomicLocus = cred_loci_map$mapped_loci[match(snps$GenomicLocus, cred_loci_map$origin)]
  ## deleterious coding
  annov <- fread(cmd=paste0("gzip -cd ", cred_dir, "/annov.txt.gz"), data.table = F) %>% filter(uniqID %in% credSNPs$uniqID) %>% dplyr::select(c(1,2,3,5))
  annov <- unique(annov %>% filter(annot %in% c("exonic", "splicing"))) %>% filter(annot=="splicing" | exonic_func!="synonymous SNV")
  geneTable <- data.frame(ensg=unique(annov$gene), stringsAsFactors = F) 
  tmp <- as.data.frame(table(annov$gene))
  geneTable$credDC <- tmp$Freq[match(geneTable$ensg, tmp$Var1)]
  tmp_gene_snps <- annov[,c("uniqID", "gene")]
  tmp_gene_snps$CADD <- snps$CADD[match(tmp_gene_snps$uniqID, snps$uniqID)]
  tmp_gene_snps$PIP <- credSNPs$prob[match(tmp_gene_snps$uniqID, credSNPs$uniqID)]
  tmp_gene_snps$gwasP <- snps$gwasP[match(tmp_gene_snps$uniqID, snps$uniqID)]
  tmp_gene_snps$locus <- snps$GenomicLocus[match(tmp_gene_snps$uniqID, snps$uniqID)]
  tmp <- with(tmp_gene_snps, aggregate(CADD, list(gene), max))
  geneTable$credDC.maxCADD <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(tmp_gene_snps, aggregate(PIP, list(gene), max))
  geneTable$credDC.maxPIP <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(tmp_gene_snps, aggregate(gwasP, list(gene), min))
  geneTable$credDC.minP <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(tmp_gene_snps, aggregate(locus, list(gene), function(x){paste(unique(x), collapse = ":")}))
  geneTable$credDC.loci <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  rm(tmp, tmp_gene_snps)
  
  ## coloc eQTL + active chromatin
  coloc <- fread(coloc_file, data.table = F) %>% filter(pp.h4>min_pp4 & nSNPs>=coloc_nSNPs & gene %in% ensg$ensembl_gene_id)
  eqtl <- fread(cmd=paste0("gzip -cd ", cred_dir, "/eqtl.txt.gz"), data.table=F) %>% filter(gene %in% coloc$gene & uniqID %in% credSNPs$uniqID)
  eqtl$tissue <- sub(".sig", "", eqtl$tissue)
  coloc$label <- with(coloc, paste(ds, gene, sep="/"))
  eqtl$label <- with(eqtl, paste(tissue, gene, sep="/"))
  eqtl <- eqtl %>% filter(label %in% coloc$label)
  
  genes <- unique(eqtl$gene)
  genes <- data.frame(ensg=genes, ensg[match(genes,ensg$ensembl_gene_id), 3:6], stringsAsFactors = F)
  colnames(genes) <- c("ensg", "chr", "start", "end", "strand")
  genes <- (genes %>% filter(ensg %in% eqtl$gene & !ensg %in% geneTable$ensg))$ensg
  genes_add <- matrix(nrow=length(genes), ncol=ncol(geneTable))
  colnames(genes_add) <- colnames(geneTable)
  genes_add <- data.frame(genes_add)
  genes_add$ensg <- genes
  geneTable <- rbind(geneTable, genes_add)
  
  eqtl$PIP <- credSNPs$prob[match(eqtl$uniqID, credSNPs$uniqID)]
  eqtl$gwasP <- snps$gwasP[match(eqtl$uniqID, snps$uniqID)]
  eqtl$locus <- snps$GenomicLocus[match(eqtl$uniqID, snps$uniqID)]
  eqtl$tissue <- as.character(apply(eqtl[,2:3], 1, function(x){if(grepl(x[1], x[2])){x[2]}else{paste0(x[1],"/",x[2])}}))
  tmp <- with(eqtl, aggregate(uniqID, list(gene), function(x){length(unique(x))}))
  geneTable$cred_eQTL <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(eqtl, aggregate(tissue, list(gene), function(x){paste(unique(x), collapse=":")}))
  geneTable$cred_eQTL.ts <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(eqtl, aggregate(PIP, list(gene), max))
  geneTable$cred_eQTL.maxPIP <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(eqtl, aggregate(gwasP, list(gene), min))
  geneTable$cred_eQTL.minP <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(eqtl, aggregate(locus, list(gene), function(x){paste(unique(x), collapse = ":")}))
  geneTable$cred_eQTL.loci <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  
  mapped_loci <- as.numeric(unique(unlist(strsplit(c(geneTable$credDC.loci, geneTable$cred_eQTL.loci), ":"))))
  mapped_loci <- mapped_loci[!is.na(mapped_loci)]
  
  ### GWS SNPs
  snps <- fread(paste0(gwas_dir, "/snps.txt"), data.table = F) %>% filter(gwasP<max_gwasP & !GenomicLocus %in% mapped_loci & !GenomicLocus %in% ex)
  ## deleterious coding
  annov <- fread(cmd=paste0("gzip -cd ", gwas_dir, "/annov.txt.gz"), data.table = F) %>% filter(uniqID %in% snps$uniqID) %>% dplyr::select(c(1,2,3,5))
  annov <- unique(annov %>% filter(annot %in% c("exonic", "splicing"))) %>% filter(annot=="splicing" | exonic_func %in% c("nonsynonymous SNV", "stopgain")) %>% filter(gene %in% ensg$ensembl_gene_id)
  genes <- unique(annov$gene[!annov$gene %in% geneTable$ensg])
  genes_add <- matrix(nrow=length(genes), ncol=ncol(geneTable))
  colnames(genes_add) <- colnames(geneTable)
  genes_add <- as.data.frame(genes_add)
  genes_add$ensg <- genes
  geneTable <- rbind(geneTable, genes_add)
  tmp <- as.data.frame(table(annov$gene))
  geneTable$gwsDC <- tmp$Freq[match(geneTable$ensg, tmp$Var1)]
  tmp_gene_snps <- annov[,c("uniqID", "gene")]
  tmp_gene_snps$CADD <- snps$CADD[match(tmp_gene_snps$uniqID, snps$uniqID)]
  tmp_gene_snps$gwasP <- snps$gwasP[match(tmp_gene_snps$uniqID, snps$uniqID)]
  tmp_gene_snps$locus <- snps$GenomicLocus[match(tmp_gene_snps$uniqID, snps$uniqID)]
  tmp <- with(tmp_gene_snps, aggregate(CADD, list(gene), max))
  geneTable$gwsDC.maxCADD <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(tmp_gene_snps, aggregate(gwasP, list(gene), min))
  geneTable$gwsDC.minP <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(tmp_gene_snps, aggregate(locus, list(gene), function(x){paste(unique(x), collapse = ":")}))
  geneTable$gwsDC.loci <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  rm(tmp, tmp_gene_snps)
  
  ## eqtl
  coloc <- fread(coloc_file, data.table = F) %>% filter(pp.h4>min_pp4 & nSNPs>=coloc_nSNPs & gene %in% ensg$ensembl_gene_id)
  eqtl <- fread(cmd=paste0("gzip -cd ", gwas_dir, "/eqtl.txt.gz"), data.table=F) %>% filter(gene %in% coloc$gene & uniqID %in% snps$uniqID)
  eqtl$tissue <- sub(".sig", "", eqtl$tissue)
  coloc$label <- with(coloc, paste(ds, gene, sep="/"))
  eqtl$label <- with(eqtl, paste(tissue, gene, sep="/"))
  eqtl <- eqtl %>% filter(label %in% coloc$label)
  
  genes <- unique(eqtl$gene)
  genes <- data.frame(ensg=genes, ensg[match(genes,ensg$ensembl_gene_id), 3:6], stringsAsFactors = F)
  colnames(genes) <- c("ensg", "chr", "start", "end", "strand")
  genes <- (genes %>% filter(ensg %in% eqtl$gene & !ensg %in% geneTable$ensg))$ensg
  genes_add <- matrix(nrow=length(genes), ncol=ncol(geneTable))
  colnames(genes_add) <- colnames(geneTable)
  genes_add <- data.frame(genes_add)
  genes_add$ensg <- genes
  geneTable <- rbind(geneTable, genes_add)
  
  eqtl$gwasP <- snps$gwasP[match(eqtl$uniqID, snps$uniqID)]
  eqtl$locus <- snps$GenomicLocus[match(eqtl$uniqID, snps$uniqID)]
  eqtl$tissue <- as.character(apply(eqtl[,2:3], 1, function(x){if(grepl(x[1], x[2])){x[2]}else{paste0(x[1],"/",x[2])}}))
  tmp <- with(eqtl, aggregate(uniqID, list(gene), function(x){length(unique(x))}))
  geneTable$gws_eQTL <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(eqtl, aggregate(tissue, list(gene), function(x){paste(unique(x), collapse=":")}))
  geneTable$gws_eQTL.ts <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(eqtl, aggregate(gwasP, list(gene), min))
  geneTable$gws_eQTL.minP <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  tmp <- with(eqtl, aggregate(locus, list(gene), function(x){paste(unique(x), collapse = ":")}))
  geneTable$gws_eQTL.loci <- tmp$x[match(geneTable$ensg, tmp$Group.1)]
  
  tmp <- as.numeric(unique(unlist(strsplit(c(geneTable$gwsDC.loci, geneTable$gws_eQTL.loci), ":"))))
  tmp <- tmp[!is.na(tmp)]
  mapped_loci <- c(mapped_loci, tmp)
  geneTable = geneTable[geneTable$ensg!="UNKNOWN",]
  
  ### label HC-1 or HC-m
  geneTable$HC_cred_loci <- apply(geneTable[,c("credDC.loci", "cred_eQTL.loci")], 1, function(x){l <- unique(unlist(strsplit(gsub(" ", "", x), ":"))); paste(l[!is.na(l)], collapse = ":")})
  geneTable$HC_cred_loci[geneTable$HC_cred_loci==""] <- NA
  tmp <- geneTable[which(geneTable$credDC>0), c("ensg", "credDC.maxPIP", "credDC.loci")]
  colnames(tmp) <- c("ensg", "PIP", "loci")
  if(length(which(grepl(":", tmp$loci)))>0){
    tmp <- do.call(rbind, apply(tmp, 1, function(x){l<-unique(unlist(strsplit(x[3],":"))); cbind(rep(x[1], length(l)), rep(x[2], length(l)),l)}))
    colnames(tmp) <- c("ensg", "PIP", "loci")
    tmp <- data.frame(tmp, stringsAsFactors = F)
    rownames(tmp) <- NULL
  }
  genes <- tmp
  
  tmp <- geneTable[which(geneTable$cred_eQTL>0), c("ensg", "cred_eQTL.maxPIP", "cred_eQTL.loci")]
  colnames(tmp) <- c("ensg", "PIP", "loci")
  if(length(which(grepl(":", tmp$loci)))>0){
    tmp <- do.call(rbind, apply(tmp, 1, function(x){l<-unique(unlist(strsplit(x[3],":"))); cbind(rep(x[1], length(l)), rep(x[2], length(l)),l)}))
    colnames(tmp) <- c("ensg", "PIP", "loci")
    tmp <- data.frame(tmp, stringsAsFactors = F)
    rownames(tmp) <- NULL
  }
  genes <- rbind(genes, tmp)
  
  genes$label <- paste(genes$ensg, genes$loci, sep=":")
  genes$PIP <- as.numeric(genes$PIP)
  genes <- genes %>% arrange(-PIP) %>% filter(!duplicated(label))

  geneTable$type <- NA
  geneTable$type[geneTable$ensg %in% genes$ensg] <- "HC_cred"
  geneTable$score <- genes$score[match(geneTable$ensg, genes$ensg)]
  
  geneTable$HC_gws_loci <- apply(geneTable[,c("gwsDC.loci", "gws_eQTL.loci")], 1, function(x){l <- unique(unlist(strsplit(gsub(" ", "", x), ":"))); paste(l[!is.na(l)], collapse = ":")})
  geneTable$HC_gws_loci[geneTable$HC_gws_loci==""] <- NA
  
  tmp <- geneTable[which(geneTable$gwsDC>0), c("ensg","gwsDC.loci")]
  colnames(tmp) <- c("ensg", "loci")
  if(length(which(grepl(":", tmp$loci)))>0){
    tmp <- do.call(rbind, apply(tmp, 1, function(x){l<-unique(unlist(strsplit(x[2],":"))); cbind(rep(x[1], length(l)),l)}))
    colnames(tmp) <- c("ensg", "loci")
    tmp <- data.frame(tmp, stringsAsFactors = F)
    rownames(tmp) <- NULL
  }
  genes <- tmp
  
  tmp <- geneTable[which(geneTable$gws_eQTL>0), c("ensg","gws_eQTL.loci")]
  colnames(tmp) <- c("ensg", "loci")
  if(length(which(grepl(":", tmp$loci)))>0){
    tmp <- do.call(rbind, apply(tmp, 1, function(x){l<-unique(unlist(strsplit(x[2],":"))); cbind(rep(x[1], length(l)),l)}))
    colnames(tmp) <- c("ensg", "loci")
    tmp <- data.frame(tmp, stringsAsFactors = F)
    rownames(tmp) <- NULL
  }
  genes <- rbind(genes, tmp)
  
  tmp <- as.data.frame(table(genes$loci))
  tmp$Var1 <- as.character(tmp$Var1)
  geneTable$type[geneTable$ensg %in% genes$ensg] <- "HC_gws"
  
  geneTable$HC_loci <- with(geneTable, ifelse(is.na(HC_cred_loci), HC_gws_loci, HC_cred_loci))
  hc_loci <- table(unlist(strsplit(geneTable$HC_loci, ":")))
  hc1_loci <- names(hc_loci)[hc_loci==1]
  geneTable$final_gene_label <- "LC"
  geneTable$final_gene_label[!is.na(geneTable$HC_loci)] <- "HC-m"
  geneTable$final_gene_label[sapply(geneTable$HC_loci, function(x){any(unlist(strsplit(x, ":")) %in% hc1_loci)})] <- "HC-1"
  
  ### prioritize by PPI
  inweb <- fread(cmd="gzip -cd data/InWeb_20160912_ensg.txt.gz", data.table = F)
  genes <- list()
  genes[["HC1"]] <- geneTable$ensg[geneTable$final_gene_label=="HC-1"]
  genes[["HCm"]] <- geneTable$ensg[geneTable$final_gene_label=="HC-m"]
  tmp <- inweb %>% filter(ensg1 %in% genes$HC1 & ensg2 %in% genes$HCm)
  tmp <- rbind(tmp, inweb %>% filter(ensg2 %in% genes$HC1 & ensg1 %in% genes$HCm))
  genes[["HC_inweb"]] <- unique(c(tmp$ensg1, tmp$ensg2))
  geneTable$HC_prioritized = with(geneTable, ifelse(final_gene_label=="HC-1", "Yes", "No"))
  geneTable$HC_prioritized[geneTable$ensg %in% genes$HC_inweb] = "Yes"
  write.table(geneTable, paste0(out_dir, "/mapped_genes.txt"), quote=F, row.names=F, sep="\t")
  save(genes, file=paste0(out_dir, "/genes.list.RData"))
}
