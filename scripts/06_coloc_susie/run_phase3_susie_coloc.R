#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
root <- "."
.libPaths(c(file.path(root,"work/software/phase3_coloc_Rlib"),.libPaths()))
suppressPackageStartupMessages({library(coloc);library(susieR);library(data.table)})
out <- file.path(root,"results/phase3_final/06_coloc")
m <- fread(file.path(out,"susie_ld_manifest.tsv")); ans <- list()
for(i in seq_len(nrow(m))){
 x<-m[i]; pref<-file.path(root,x$plink_prefix); h<-fread(file.path(root,x$harmonized_file))
 bim<-fread(paste0(pref,".bim"),header=FALSE); setnames(bim,c("chr","rsid","cm","bp","A1","A2"))
 ld<-as.matrix(fread(cmd=paste("gzip -cd",shQuote(paste0(pref,".ld.gz"))),header=FALSE))
 if(nrow(ld)!=nrow(bim)){ans[[length(ans)+1]]<-data.table(gene_symbol=x$gene_symbol,context=x$context,dataset_id=x$dataset_id,status="LD_DIMENSION_MISMATCH");next}
 j<-match(bim$rsid,h$rsid); keep<-which(!is.na(j)); bim<-bim[keep]; h<-h[j[keep]]; ld<-ld[keep,keep,drop=FALSE]
 same<-h$effect_allele==bim$A1; swap<-h$effect_allele==bim$A2; ok<-same|swap
 bim<-bim[ok];h<-h[ok];ld<-ld[ok,ok,drop=FALSE]; flip<-swap[ok]; h[flip,beta_gwas:=-beta_gwas];h[flip,beta_qtl:=-beta_qtl]
 rownames(ld)<-h$rsid; colnames(ld)<-h$rsid
 diagv<-diag(ld); symmetry<-max(abs(ld-t(ld))); mineig<-tryCatch(min(eigen(ld,symmetric=TRUE,only.values=TRUE)$values),error=function(e) NA_real_)
 if(nrow(h)<100){ans[[length(ans)+1]]<-data.table(gene_symbol=x$gene_symbol,context=x$context,dataset_id=x$dataset_id,status="NOT_ENOUGH_REFERENCE_SNPS",n_snps=nrow(h),ld_min_eigen=mineig,ld_max_asymmetry=symmetry);next}
 d1<-list(beta=h$beta_gwas,varbeta=h$se_gwas^2,snp=h$rsid,position=h$position_gwas,type="cc",s=28094/(28094+347768),N=375862,MAF=h$maf_gwas,LD=ld)
 d2<-list(beta=h$beta_qtl,varbeta=h$se_qtl^2,snp=h$rsid,position=h$position_qtl,type="quant",N=as.numeric(x$qtl_sample_size),MAF=h$maf,LD=ld)
 z<-tryCatch({s1<-runsusie(d1,maxit=100);s2<-runsusie(d2,maxit=100);cs<-coloc.susie(s1,s2);list(s1=s1,s2=s2,cs=cs)},error=function(e)e)
 if(inherits(z,"error")){ans[[length(ans)+1]]<-data.table(gene_symbol=x$gene_symbol,context=x$context,dataset_id=x$dataset_id,status=paste0("FAILED: ",conditionMessage(z)),n_snps=nrow(h),ld_min_eigen=mineig,ld_max_asymmetry=symmetry);next}
 sm<-as.data.table(z$cs$summary)
 if(!nrow(sm)){ans[[length(ans)+1]]<-data.table(gene_symbol=x$gene_symbol,context=x$context,dataset_id=x$dataset_id,status="COMPLETED_NO_SIGNAL_PAIR",n_snps=nrow(h),gwas_converged=z$s1$converged,qtl_converged=z$s2$converged,ld_min_eigen=mineig,ld_max_asymmetry=symmetry);next}
 sm[,`:=`(gene_symbol=x$gene_symbol,context=x$context,dataset_id=x$dataset_id,status="COMPLETED",n_snps=nrow(h),gwas_converged=z$s1$converged,qtl_converged=z$s2$converged,ld_min_eigen=mineig,ld_max_asymmetry=symmetry,ABF_PP4=x$PP4_abf)]
 ans[[length(ans)+1]]<-sm
}
r<-rbindlist(ans,fill=TRUE);fwrite(r,file.path(out,"susie_coloc.tsv"),sep="\t",na="NA")
print(r)
