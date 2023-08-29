# Modified from /project2/guiming/xsun/proteome_alzheimer/codes/locus_plot_fusion.R
locus_plot <- function(study_name, tissue, ctwas_res, chrom=22, region_tag2=5,
                       eqtl_or_sqtl, gwide_tiss_sig_thresh,
                       xlim=NULL, return_table=F, focus=NULL,
                       label_panel="TWAS", label_genes=NULL, use_gname=T, label_pos=NULL,
                       plot_eqtl=NULL, rerun_ctwas=F, rerun_load_only=F,
                       legend_side="right", legend_panel="cTWAS",
                       twas_ymax=NULL,draw_gene_track = draw_gene_track,
                       
                       results_dir="data/ctwas_rss/",
                       zscore_dir="data/impute_expr_z/",
                       calced_ld_dir="/scratch/jmcclellan/scratch_pdb_filt_ctwas_gtex_v8_eur_ld/",
                       by_chrom_fix=T # Think my hacky way of doing by chromosome ctwas requires a fix somewhere, regionlist seems to have numbers as characters at second list layer
                       ){
    require(glue)
    require(tidyverse)
    study_prefix <- switch(study_name, "meta_analysis_BCAC_UKB_ovr_brca" = "meta_bcac_ukb_ovr")
    region_tag1 <- chrom

    
    a <- ctwas_res 
    ctwas_gene_res <- a %>%
      filter(type == "gene")
    
    ##added for fusion weights
    a$genename[a$type == "gene"] <- unlist(strsplit(a$id[a$type == "gene"],split = "[.]"))[seq(2,2*nrow(a[a$type == "gene",]), by = 2)]
    a$id_ensembl[a$type == "gene"] <- unlist(strsplit(a$id[a$type == "gene"],split = "[.]"))[seq(1,2*nrow(a[a$type == "gene",]), by = 2)]
    
    
    regionlist <- readRDS(file.path(results_dir, glue("{study_prefix}__{tissue}__{eqtl_or_sqtl}_chr{chrom}.regionlist.RDS")))

    if (isTRUE(by_chrom_fix)){
      region <- regionlist[[as.numeric(region_tag1)]][[as.character(region_tag2)]]
    } else {
      region <- regionlist[[as.numeric(region_tag1)]][[region_tag2]]
    }
  
    R_snp_info <- do.call(rbind, lapply(region$regRDS,
                                        function(x){data.table::fread(paste0(file.path(calced_ld_dir,
                                                                                       tools::file_path_sans_ext(basename(x))),
                                                                             ".Rvar"))}))
  
    if (isTRUE(rerun_ctwas)){
        ld_exprfs <- paste0(results_dir, "/", analysis_id, "_expr_chr", 1:22, ".expr.gz")
        temp_reg <- data.frame("chr" = paste0("chr",region_tag1), "start" = region$start, "stop" = region$stop)
  
        write.table(temp_reg, 
              #file= paste0(results_dir, "/", analysis_id, "_ctwas.temp.reg.txt") , 
              file= "temp_reg.txt",
              row.names=F, col.names=T, sep="\t", quote = F)
  
        load(paste0(results_dir, "/", analysis_id, "_expr_z_snp.Rd"))
  
        z_gene_temp <-  z_gene[z_gene$id %in% a$id[a$type=="gene"],]
        z_snp_temp <-  z_snp[z_snp$id %in% R_snp_info$id,]
    
        if (!rerun_load_only){
            ctwas::ctwas_rss(z_gene_temp, z_snp_temp, ld_exprfs, ld_pgenfs = NULL, 
                      ld_R_dir = dirname(region$regRDS)[1],
                      ld_regions_custom = "temp_reg.txt", thin = 1, 
                      outputdir = ".", outname = "temp", ncore = 1, ncore.rerun = 1, prob_single = 0,
                      group_prior = estimated_group_prior, group_prior_var = estimated_group_prior_var,
                      estimate_group_prior = F, estimate_group_prior_var = F)
        }
    
        a_bkup <- a         
        a <- as.data.frame(data.table::fread("temp.susieIrss.txt", header = T))
    
        rownames(z_snp_temp) <- z_snp_temp$id
        z_snp_temp <- z_snp_temp[a$id[a$type=="SNP"],]
        z_gene_temp <- z_gene_temp[a$id[a$type=="gene"],]
    
        a$genename <- NA
        a$gene_type <- NA

        a[a$type=="gene",c("genename", "gene_type")] <- a_bkup[match(a$id[a$type=="gene"], a_bkup$id),c("genename","gene_type")]
    
        a$z <- NA
        a$z[a$type=="SNP"] <- z_snp_temp$z
        a$z[a$type=="gene"] <- z_gene_temp$z
    }
  
    a_pos_bkup <- a$pos
    #a$pos[a$type=="gene"] <- G_list$tss[match(sapply(a$genename[a$type=="gene"], function(x){unlist(strsplit(x, "[.]"))[1]}) ,G_list$hgnc_symbol)]
    a$pos[is.na(a$pos)] <- a_pos_bkup[is.na(a$pos)]
    a$pos <- a$pos/1000000
  
    if (!is.null(xlim)){
    
        if (is.na(xlim[1])){
            xlim[1] <- min(a$pos)
        }
    
        if (is.na(xlim[2])){
            xlim[2] <- max(a$pos)
        }
    
        a <- a[a$pos>=xlim[1] & a$pos<=xlim[2],,drop=F]
    }
  
    if (is.null(focus)){
        focus <- a$id[which.max(abs(a$z)[a$type=="gene"])]
    }
  
    if (is.null(label_genes)){
        label_genes <- focus
    }
  
    
    if (is.null(label_pos)){
        label_pos <- rep(3, length(label_genes))
    }
  
    if (is.null(plot_eqtl)){
        plot_eqtl <- focus
    }
  
    focus <- a$id[which(a$id==focus)]
    a$focus <- 0
    a$focus <- as.numeric(a$id==focus)
    
    a$PVALUE <- (-log(2) - pnorm(abs(a$z), lower.tail=F, log.p=T))/log(10)
  
    # results_dir="data/ctwas_rss/"
    R_gene <- readRDS(file.path(results_dir, glue("{study_prefix}__{tissue}__{eqtl_or_sqtl}_chr{chrom}_LDR"), basename(region$R_g_file)))
    R_snp_gene <- readRDS(file.path(results_dir, glue("{study_prefix}__{tissue}__{eqtl_or_sqtl}_chr{chrom}_LDR"), basename(region$R_sg_file)))
    R_snp <- as.matrix(Matrix::bdiag(lapply(region$regRDS,
                                            function(x){readRDS(file.path(calced_ld_dir,
                                                                          basename(x)))})))
    rownames(R_gene) <- region$gid
    colnames(R_gene) <- region$gid
    rownames(R_snp_gene) <- R_snp_info$id
    colnames(R_snp_gene) <- region$gid
    rownames(R_snp) <- R_snp_info$id
    colnames(R_snp) <- R_snp_info$id
  
    a$r2max <- NA
    a$r2max[a$type=="gene"] <- R_gene[focus,a$id[a$type=="gene"]]
    a$r2max[a$type=="SNP"] <- R_snp_gene[a$id[a$type=="SNP"],focus]
  
    r2cut <- 0.4
    colorsall <- c("#7fc97f", "#beaed4", "#fdc086")
  
    start <- min(a$pos)
    end <- max(a$pos)
  
    #layout(matrix(1:3, ncol = 1), widths = 1, heights = c(1.5,1.75,0.75), respect = FALSE)
    
    if (draw_gene_track){
      layout(matrix(1:4, ncol = 1), widths = 1, heights = c(1.5,0.25,1.5,0.75), respect = FALSE)
    } else {
      layout(matrix(1:3, ncol = 1), widths = 1, heights = c(1.5,0.25,1.5), respect = FALSE)
    }
    
  
    par(mar = c(0, 4.1, 0, 2.1))
  
    if (is.null(twas_ymax)){
        twas_ymax <- max(a$PVALUE)*1.1
    }
  
    plot(a$pos[a$type=="SNP"], a$PVALUE[a$type == "SNP"], pch = 21, xlab=paste0("Chromosome ", region_tag1, " position (Mb)"), frame.plot=FALSE, bg = colorsall[1], ylab = "-log10(p value)", panel.first = grid(), ylim =c(0, twas_ymax), xaxt = 'n', xlim=c(start, end))
  
    # abline(h=-log10(alpha/nrow(ctwas_gene_res)), col ="red", lty = 2)
    # gwide_tiss_sig_thresh
    abline(h=-log10(gwide_tiss_sig_thresh), col ="red", lty = 2)
    points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$PVALUE[a$type == "SNP"  & a$r2max > r2cut], pch = 21, bg = "purple")
    points(a$pos[a$type=="SNP" & a$focus == 1], a$PVALUE[a$type == "SNP" & a$focus == 1], pch = 21, bg = "salmon")
    points(a$pos[a$type=="gene" & a$group=="Expression"], a$PVALUE[a$type == "gene" & a$group=="Expression"], pch = 22, bg = colorsall[1], cex = 2)
    points(a$pos[a$type=="gene" & a$group=="Splicing"], a$PVALUE[a$type == "gene" & a$group=="Splicing"], pch = 23, bg = colorsall[1], cex = 2)
    points(a$pos[a$type=="gene" & a$group=="Methylation"], a$PVALUE[a$type == "gene" & a$group=="Methylation"], pch = 24, bg = colorsall[1], cex = 2)

    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Expression"], a$PVALUE[a$type == "gene"  & a$r2max > r2cut & a$group=="Expression"], pch = 22, bg = "purple", cex = 2)
    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Splicing"], a$PVALUE[a$type == "gene"  & a$r2max > r2cut & a$group=="Splicing"], pch = 23, bg = "purple", cex = 2)
    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Methylation"], a$PVALUE[a$type == "gene"  & a$r2max > r2cut & a$group=="Methylation"], pch = 24, bg = "purple", cex = 2)

    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Expression"], a$PVALUE[a$type == "gene" & a$focus == 1 & a$group=="Expression"], pch = 22, bg = "salmon", cex = 2)
    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Splicing"], a$PVALUE[a$type == "gene" & a$focus == 1 & a$group=="Splicing"], pch = 23, bg = "salmon", cex = 2)
    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Methylation"], a$PVALUE[a$type == "gene" & a$focus == 1 & a$group=="Methylation"], pch = 24, bg = "salmon", cex = 2)
  
    if (legend_panel=="TWAS"){
        x_pos <- ifelse(legend_side=="right", max(a$pos)-0.2*(max(a$pos)-min(a$pos)), min(a$pos))
        legend(x_pos, y= twas_ymax*0.95, c("Expression","Splicing","Methylation","SNP","Lead TWAS Gene", "R2 > 0.4", "R2 <= 0.4"), pch = c(22,23,24,21,19,19,19), col = c("black","black","black","black", "salmon", "purple", colorsall[1]), cex=0.7, title.adj = 0)
    }
    
    og_label_genes <- label_genes
    label_genes <- a %>%
      filter(id %in% og_label_genes) %>%
      pull(id)
    
    if (label_panel=="TWAS" | label_panel=="both"){
      # if (length(label_genes) == 0){
      #   browser()
      # }
      for (i in 1:length(label_genes)){
        if (isTRUE(use_gname)){
          text(a$pos[a$id==label_genes[i]], a$PVALUE[a$id==label_genes[i]], labels=a$genename[a$id==label_genes[i]], pos=label_pos[i], cex=0.7)
        } else {
          text(a$pos[a$id==label_genes[i]], a$PVALUE[a$id==label_genes[i]], labels=a$id[a$id==label_genes[i]], pos=label_pos[i], cex=0.7)
        }
      }
    }

    par(mar = c(0.25, 4.1, 0.25, 2.1))
    
    plot(NA, xlim = c(start, end), ylim = c(0, length(plot_eqtl)), frame.plot = F, axes = F, xlab = NA, ylab = NA)
    
    for (i in 1:length(plot_eqtl)){
      cgene <- a$id[which(a$id==plot_eqtl[i])]
      
      if (isTRUE(by_chrom_fix)){
        load(file.path(zscore_dir, glue("{study_prefix}__{tissue}__{eqtl_or_sqtl}_chr{chrom}_chr{chrom}.exprqc.Rd"))) 
      } else {
        load(file.path(zscore_dir, glue("{study_prefix}__{tissue}__{eqtl_or_sqtl}_chr{chrom}.exprqc.Rd"))) 
      }
      # load(paste0(results_dir, "/",analysis_id, "_expr_chr", region_tag1, ".exprqc.Rd"))
      eqtls <- rownames(wgtlist[[cgene]])
      eqtl_pos <- a$pos[a$id %in% eqtls]
      
      col="grey"
      
      rect(start, length(plot_eqtl)+1-i-0.8, end, length(plot_eqtl)+1-i-0.2, col = col, border = T, lwd = 1)
      
      if (length(eqtl_pos)>0){
        for (j in 1:length(eqtl_pos)){
          segments(x0=eqtl_pos[j], x1=eqtl_pos[j], y0=length(plot_eqtl)+1-i-0.2, length(plot_eqtl)+1-i-0.8, lwd=1.5)  
        }
      }
    }
    
    text(start, length(plot_eqtl)-(1:length(plot_eqtl))+0.5,  
         labels = plot_eqtl, srt = 0, pos = 2, xpd = TRUE, cex=0.7)
    
  
    par(mar = c(4.1, 4.1, 0, 2.1))
  
    plot(a$pos[a$type=="SNP"], a$susie_pip[a$type == "SNP"], pch = 19, xlab=paste0("Chromosome ", region_tag1, " position (Mb)"),frame.plot=FALSE, col = "white", ylim= c(0,1.1), ylab = "cTWAS PIP", xlim = c(start, end))
  
    grid()
    points(a$pos[a$type=="SNP"], a$susie_pip[a$type == "SNP"], pch = 21, xlab="Genomic position", bg = colorsall[1])
    points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$susie_pip[a$type == "SNP"  & a$r2max >r2cut], pch = 21, bg = "purple")
    points(a$pos[a$type=="SNP" & a$focus == 1], a$susie_pip[a$type == "SNP" & a$focus == 1], pch = 21, bg = "salmon")
    points(a$pos[a$type=="gene" & a$group=="Expression"], a$susie_pip[a$type == "gene" & a$group=="Expression"], pch = 22, bg = colorsall[1], cex = 2)
    points(a$pos[a$type=="gene" & a$group=="Splicing"], a$susie_pip[a$type == "gene" & a$group=="Splicing"], pch = 23, bg = colorsall[1], cex = 2)
    points(a$pos[a$type=="gene" & a$group=="Methylation"], a$susie_pip[a$type == "gene" & a$group=="Methylation"], pch = 24, bg = colorsall[1], cex = 2)
    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Expression"], a$susie_pip[a$type == "gene"  & a$r2max > r2cut & a$group=="Expression"], pch = 22, bg = "purple", cex = 2)
    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Splicing"], a$susie_pip[a$type == "gene"  & a$r2max > r2cut & a$group=="Splicing"], pch = 23, bg = "purple", cex = 2)
    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Methylation"], a$susie_pip[a$type == "gene"  & a$r2max > r2cut & a$group=="Methylation"], pch = 24, bg = "purple", cex = 2)

    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Expression"], a$susie_pip[a$type == "gene" & a$focus == 1 & a$group=="Expression"], pch = 22, bg = "salmon", cex = 2)
    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Splicing"], a$susie_pip[a$type == "gene" & a$focus == 1 & a$group=="Splicing"], pch = 23, bg = "salmon", cex = 2)
    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Methylation"], a$susie_pip[a$type == "gene" & a$focus == 1 & a$group=="Methylation"], pch = 24, bg = "salmon", cex = 2)
  
    if (legend_panel=="cTWAS"){
        x_pos <- ifelse(legend_side=="right", max(a$pos)-0.2*(max(a$pos)-min(a$pos)), min(a$pos))
        legend(x_pos, y= 1 ,c("Expression","Splicing","Methylation","SNP","Lead TWAS Gene", "R2 > 0.4", "R2 <= 0.4"), pch = c(22,23,24,21,19,19,19), col = c("black","black","black", "black", "salmon", "purple", colorsall[1]), cex=0.7, title.adj = 0)
    }


  #####modified for fusion weights
    # if (label_panel=="cTWAS" | label_panel=="both"){
    #     for (i in 1:length(label_genes)){
    #     text(a$pos[a$id==label_genes[i]], a$susie_pip[a$id==label_genes[i]], labels=a$genename[a$id==label_genes[i]], pos=label_pos[i], cex=0.7)
    #     }
    # }
    if (label_panel=="cTWAS" | label_panel=="both"){
      for (i in 1:length(label_genes)){
        if (isTRUE(use_gname)){
          text(a$pos[a$id==label_genes[i]], a$susie_pip[a$id==label_genes[i]], labels=a$genename[a$id==label_genes[i]], pos=label_pos[i], cex=0.7)
        } else {
          text(a$pos[a$id==label_genes[i]], a$susie_pip[a$id==label_genes[i]], labels=a$id[a$id==label_genes[i]], pos=label_pos[i], cex=0.7)
        }
      }
    }
  
    
    if (draw_gene_track){
      source("/project2/guiming/xsun/proteome_alzheimer/codes/trackplot.R")
      
      query_ucsc = TRUE
      build = "hg19"
      col = "gray70"
      txname = NULL
      genename = NULL
      collapse_txs = TRUE
      gene_model = "/project2/guiming/xsun/proteome_alzheimer/data_others/hg19.ensGene.gtf.gz"          #######
      isGTF = T
      
      ##########
      start <- min(a$pos)*1000000
      end <- max(a$pos)*1000000
      chr <- paste0("chr",as.character(unique(a$chrom)))
      
      #collect gene models
      if(is.null(gene_model)){
        if(query_ucsc){
          message("Missing gene model. Trying to query UCSC genome browser..")
          etbl = .extract_geneModel_ucsc(chr, start = start, end = end, refBuild = build, txname = txname, genename = genename)
        } else{
          etbl = NULL
        }
      } else{
        if(isGTF){
          etbl = .parse_gtf(gtf = gene_model, chr = chr, start = start, end = end, txname = txname, genename = genename)  
        } else{
          etbl = .extract_geneModel(ucsc_tbl = gene_model, chr = chr, start = start, end = end, txname = txname, genename = genename)  
        }
      }
      
      #draw gene models
      if(!is.null(etbl)){
        if(collapse_txs){
          etbl = .collapse_tx(etbl)
        }
        
        #subset to protein coding genes in ensembl and lincRNAs included in the analysis
        #etbl <- etbl[names(etbl) %in% G_list$ensembl_gene_id] 
        
        ######modified for fusion weights
        etbl <- etbl[names(etbl) %in% c(G_list$ensembl_gene_id[G_list$gene_biotype=="protein_coding"], sapply(a$id[a$type=="gene"], function(x){unlist(strsplit(x, split="[.]"))[1]}))]
        
        for (i in 1:length(etbl)){
          ensembl_name <- attr(etbl[[i]], "gene")
          gene_name <- G_list$hgnc_symbol[match(ensembl_name, G_list$ensembl_gene_id)]
          if (gene_name==""){
            gene_name <- a$genename[sapply(a$id, function(x){unlist(strsplit(x,split="[.]"))[1]})==ensembl_name]
          }
          if (length(gene_name)>0){
            attr(etbl[[i]], "gene") <- gene_name
          }
          
          attr(etbl[[i]], "imputed") <- ensembl_name %in% sapply(ctwas_gene_res$id, function(x){unlist(strsplit(x, split="[.]"))[1]})
          
        }
        
        par(mar = c(0, 4.1, 0, 2.1))
        
        plot(NA, xlim = c(start, end), ylim = c(0, length(etbl)), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
        
        etbl <- etbl[order(-sapply(etbl, function(x){attr(x, "start")}))]
        
        for(tx_id in 1:length(etbl)){
          txtbl = etbl[[tx_id]]
          
          if (attr(txtbl, "imputed")){
            exon_col = "#192a56"
          } else {
            exon_col = "darkred"
          }
          
          segments(x0 = attr(txtbl, "start"), y0 = tx_id-0.45, x1 = attr(txtbl, "end"), y1 = tx_id-0.45, col = exon_col, lwd = 1)
          
          if(is.na(attr(txtbl, "tx"))){
            text(x = start, y = tx_id-0.45, labels = paste0(attr(txtbl, "gene")), cex = 0.7, adj = 0, srt = 0, pos = 2, xpd = TRUE)
          } else {
            text(x = start, y = tx_id-0.45, labels = paste0(attr(txtbl, "tx"), " [", attr(txtbl, "gene"), "]"), cex = 0.7, adj = 0, srt = 0, pos = 2, xpd = TRUE)
          }
          
          rect(xleft = txtbl[[1]], ybottom = tx_id-0.75, xright = txtbl[[2]], ytop = tx_id-0.25, col = exon_col, border = NA)
          if(attr(txtbl, "strand") == "+"){
            dirat = pretty(x = c(min(txtbl[[1]]), max(txtbl[[2]])))
            dirat[1] = min(txtbl[[1]]) #Avoid drawing arrows outside gene length
            dirat[length(dirat)] = max(txtbl[[2]])
            points(x = dirat, y = rep(tx_id-0.45, length(dirat)), pch = ">", col = exon_col)
          }else{
            dirat = pretty(x = c(min(txtbl[[1]]), max(txtbl[[2]])))
            dirat[1] = min(txtbl[[1]]) #Avoid drawing arrows outside gene length
            dirat[length(dirat)] = max(txtbl[[2]])
            points(x = dirat, y = rep(tx_id-0.45, length(dirat)), pch = "<", col = exon_col)
          }
        }
      }
    }
    
    
    
    if (return_table){
        return(a)
    }
}








