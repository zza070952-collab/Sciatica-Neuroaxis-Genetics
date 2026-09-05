#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
root <- "."
.libPaths(c(file.path(root,"work/software/phase3_coloc_Rlib"), .libPaths()))
suppressPackageStartupMessages({library(coloc); library(data.table)})
out <- file.path(root,"results/phase3_final/06_coloc")
inp <- file.path(root,"work/phase3_final/coloc/inputs")
harmdir <- file.path(root,"work/phase3_final/coloc/harmonized"); dir.create(harmdir,recursive=TRUE,showWarnings=FALSE)
inv <- fread(file.path(out,"coloc_input_inventory.tsv"))
inv <- inv[status=="AVAILABLE"]
res <- list(); harmon <- list()
for(i in seq_len(nrow(inv))){
  x <- inv[i]
  q <- fread(file.path(root,x$qtl_file)); g <- fread(file.path(root,x$gwas_file))
  numq <- c("position","maf","pvalue","beta","se"); numg <- c("position","z","beta","se","p","EAF")
  for(v in intersect(numq,names(q))) q[,(v):=as.numeric(get(v))]
  for(v in intersect(numg,names(g))) g[,(v):=as.numeric(get(v))]
  q <- q[!is.na(rsid)&rsid!=""&rsid!="NA"]; g <- g[!is.na(rsid)&rsid!=""&rsid!="NA"]
  setorder(q,pvalue); q <- unique(q,by="rsid"); setorder(g,p); g <- unique(g,by="rsid")
  d <- merge(g,q,by="rsid",suffixes=c("_gwas","_qtl"))
  same <- d$ref_gwas==d$ref_qtl & d$effect_allele==d$alt
  swap <- d$ref_gwas==d$alt & d$effect_allele==d$ref_qtl
  d <- d[same|swap]
  d[swap[same|swap],beta_qtl:=-beta_qtl]
  d <- d[is.finite(beta_gwas)&is.finite(se_gwas)&se_gwas>0&is.finite(p)&
         is.finite(beta_qtl)&is.finite(se_qtl)&se_qtl>0&is.finite(pvalue)&
         is.finite(EAF)&EAF>0&EAF<1&is.finite(maf)&maf>0&maf<1]
  d[,maf_gwas:=pmin(EAF,1-EAF)]
  fwrite(d,file.path(harmdir,paste0(x$gene_symbol,"__",x$dataset_id,".tsv.gz")),sep="\t",na="NA")
  harmon[[length(harmon)+1]] <- data.table(gene_symbol=x$gene_symbol,context=x$context,dataset_id=x$dataset_id,
    qtl_rows=nrow(q),gwas_rows=nrow(g),shared_rsid_before_alleles=nrow(merge(g[,.(rsid)],q[,.(rsid)],by="rsid")),
    harmonized_snps=nrow(d),same_alleles=sum(same),swapped_alleles=sum(swap))
  if(nrow(d)<50){
    res[[length(res)+1]] <- data.table(gene_symbol=x$gene_symbol,ensembl_gene=x$ensembl_gene,context=x$context,
      dataset_id=x$dataset_id,qtl_sample_size=x$qtl_sample_size,n_snps=nrow(d),status="NOT_ENOUGH_OVERLAP",
      PP0=NA_real_,PP1=NA_real_,PP2=NA_real_,PP3=NA_real_,PP4=NA_real_,lead_gwas=NA_character_,lead_gwas_p=NA_real_,
      lead_qtl=NA_character_,lead_qtl_p=NA_real_,lead_shared_direction=NA_character_)
    next
  }
  ds1 <- list(beta=d$beta_gwas,varbeta=d$se_gwas^2,snp=d$rsid,position=d$position_gwas,
              type="cc",s=28094/(28094+347768),N=28094+347768,MAF=d$maf_gwas)
  ds2 <- list(beta=d$beta_qtl,varbeta=d$se_qtl^2,snp=d$rsid,position=d$position_qtl,
              type="quant",N=as.numeric(x$qtl_sample_size),MAF=d$maf)
  z <- tryCatch(coloc.abf(ds1,ds2,p1=1e-4,p2=1e-4,p12=1e-5),error=function(e)e)
  if(inherits(z,"error")){
    res[[length(res)+1]] <- data.table(gene_symbol=x$gene_symbol,ensembl_gene=x$ensembl_gene,context=x$context,
      dataset_id=x$dataset_id,qtl_sample_size=x$qtl_sample_size,n_snps=nrow(d),status=paste0("FAILED: ",conditionMessage(z)),
      PP0=NA_real_,PP1=NA_real_,PP2=NA_real_,PP3=NA_real_,PP4=NA_real_,lead_gwas=NA_character_,lead_gwas_p=NA_real_,
      lead_qtl=NA_character_,lead_qtl_p=NA_real_,lead_shared_direction=NA_character_)
    next
  }
  lg <- d[which.min(p)]; lq <- d[which.min(pvalue)]
  atlead <- d[rsid==lg$rsid][1]
  direction <- if(nrow(atlead)) ifelse(sign(atlead$beta_gwas)==sign(atlead$beta_qtl),"concordant","discordant") else NA_character_
  s <- z$summary
  res[[length(res)+1]] <- data.table(gene_symbol=x$gene_symbol,ensembl_gene=x$ensembl_gene,context=x$context,
    dataset_id=x$dataset_id,qtl_sample_size=x$qtl_sample_size,n_snps=nrow(d),status="COMPLETED",
    PP0=unname(s[["PP.H0.abf"]]),PP1=unname(s[["PP.H1.abf"]]),PP2=unname(s[["PP.H2.abf"]]),
    PP3=unname(s[["PP.H3.abf"]]),PP4=unname(s[["PP.H4.abf"]]),lead_gwas=lg$rsid,lead_gwas_p=lg$p,
    lead_qtl=lq$rsid,lead_qtl_p=lq$pvalue,lead_shared_direction=direction)
}
r <- rbindlist(res,fill=TRUE); h <- rbindlist(harmon,fill=TRUE)
r[,evidence:=fcase(PP4>=.8,"STRONG",PP4>=.6,"MODERATE",status=="COMPLETED","NO_SUPPORT",default=status)]
setorder(r,-PP4)
fwrite(r,file.path(out,"coloc_all.tsv"),sep="\t",na="NA")
fwrite(r[PP4>=.8],file.path(out,"coloc_strong.tsv"),sep="\t",na="NA")
fwrite(h,file.path(out,"coloc_harmonization.tsv"),sep="\t",na="NA")
cat("completed=",sum(r$status=="COMPLETED")," strong=",sum(r$PP4>=.8,na.rm=TRUE)," moderate=",sum(r$PP4>=.6&r$PP4<.8,na.rm=TRUE),"\n",sep="")
