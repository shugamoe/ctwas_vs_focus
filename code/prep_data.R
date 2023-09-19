main <- function(study_name,
                 tiss_list="data/trait_brca/tissue_lists/11tiss.txt",
                 focus_no_pval_filt=T,
                 eqtl_or_sqtl="eqtl",
                 focus_low_cutoff=0.5,
                 ctwas_high_cutoff=0.9,
                 trait="brca",
                 want_genes=NULL
                 ){
  require(data.table)
  require(tidyverse)
  require(glue)

  tiss_fname_only <- basename(tiss_list)
  tiss_name <- str_split_1(tiss_fname_only, "\\.")[1]

  want_tissues <- fread(tiss_list, header=F)[["V1"]]

  gene_key_df <- fread(glue("data/trait_{trait}/reporting_friendly_tables_output/gencode_26_40_conversion_and_distances_to_known_{trait}_index_snps.tsv")) %>%
    select(all_of(c("gene_id.v26", "gene_name.v40", "gene_type.v40", "gene_id_no_decimal")))

  # Compile ctwas results as expected by locus_plot function
  print("Getting ctwas results.")
  zscore_dir <-glue("data/trait_{trait}/impute_expr_z/")
  ctwas_dir <- glue("data/trait_{trait}/ctwas_rss/")
  ctwas_res <- get_ctwas_res(study_name, want_tissues, eqtl_or_sqtl, gene_key_df,
                             zscore_dir=zscore_dir, ctwas_dir=ctwas_dir)
  print("Got ctwas results.")

  focus_type_tag <- switch(focus_no_pval_filt,
                           T = "no_pval_filt",
                           F = "")

  # Recreate what Dezheng did, get data on FOCUS and cTWAS to see which
  # high susie pip from cTWAS correspond with low FOCUS pips
  if (trait == "ad"){
  } else {
    focus_rolled_df <- fread(glue("data/trait_{trait}/reporting_friendly_tables_output/gene_rolled_focus_{eqtl_or_sqtl}_{tiss_name}_{study_name}_{focus_type_tag}_by_tissue.focus.tsv")) %>%
      select(-all_of(c("chrom")))
  }


  twas_df <- tibble()
  for (cur_tiss in want_tissues){
    cur_tiss_z_name <- glue("{tolower(cur_tiss)}_highest_abs_twas_z")
      # rename(
      #        !!glue("{cur_tiss}_num_in_focus") := !!glue("{tolower(cur_tiss)}_num_in_focus"),
      #        !!glue("{cur_tiss}_num_in_cred_set") := !!glue("{tolower(cur_tiss)}_num_in_cred_set"),
      #        !!glue("{cur_tiss}_highest_pip") := !!glue("{tolower(cur_tiss)}_highest_pip"),
      #        !!glue("my_cur_z") := !!glue("{tolower(cur_tiss)}_highest_abs_twas_z")
      #        ) %>%
      # mutate(!!glue("{cur_tiss}_highest_abs_twas_z_pvalue") := 2 * dnorm(-abs(my_cur_z))) %>%
      # rename(!!glue("{cur_tiss}_highest_abs_twas_z") := my_cur_z)
    if (trait == "ad"){
    } else {
      focus_rolled_df <- focus_rolled_df %>%
        rename(
               !!glue("{cur_tiss}_num_in_focus") := !!glue("{tolower(cur_tiss)}_num_in_focus"),
               !!glue("{cur_tiss}_num_in_cred_set") := !!glue("{tolower(cur_tiss)}_num_in_cred_set"),
               !!glue("{cur_tiss}_highest_pip") := !!glue("{tolower(cur_tiss)}_highest_pip"),
               ) %>%
        mutate(!!glue("{cur_tiss}_highest_abs_twas_z_pvalue") := 2 * dnorm(-abs(get(cur_tiss_z_name)))) %>%
        rename(!!glue("{cur_tiss}_highest_abs_twas_z") := !!cur_tiss_z_name)
    }

    # zscore to pvalue for TWAS:
    # 2 * dnorm(-abs(zscore)) 
    cur_tiss_twas_df <- fread(glue("data/trait_{trait}/spredixcan_{eqtl_or_sqtl}_mashr/spredixcan_igwas_gtexmashrv8_{study_name}__PM__{cur_tiss}.csv")) %>%
                                                rename(!!glue("{cur_tiss}_spxcan_twas_z")
                                                       := zscore,
                                                       !!glue("{cur_tiss}_spxcan_twas_pvalue")
                                                       := pvalue)

    cur_tiss_twas_df <- cur_tiss_twas_df %>%
          select(-all_of(c("gene_name", "effect_size", "var_g", "pred_perf_r2",
                           "pred_perf_pval", "n_snps_used", "n_snps_in_cov",
                           "n_snps_in_model", "best_gwas_p", "largest_weight",
                           "pred_perf_qval")))
    if (nrow(twas_df) == 0){
      twas_df <- cur_tiss_twas_df 
    } else {
      twas_df <- full_join(twas_df, cur_tiss_twas_df, by = "gene")
    }
  }
  twas_df <- twas_df %>%
    mutate(ens_gene_id_no_decimal = str_extract(gene, "(ENSG\\d{1,})"))

  ctwas_rolled_df <- fread(glue("data/trait_{trait}/reporting_friendly_tables_output/gene_rolled_ctwas_{study_name}_eqtl.susieIrss.txt")) %>%
    mutate(ens_gene_id_no_decimal = str_extract(id, "(ENSG\\d{1,})"))

  if (trait == "ad"){
    gene_df_master <- full_join(ctwas_rolled_df, twas_df, by = c("ens_gene_id_no_decimal"))
  } else {
    gene_df_master <- full_join(focus_rolled_df, ctwas_rolled_df, by = c("ens_gene_id_no_decimal")) %>%
      full_join(twas_df, by = "ens_gene_id_no_decimal") 
  }

  # Go through each tissue and recreate Dezheng's "321" for breast
  # 321 reached by seeing how many genes were > 0.9 cTWAS susie PIP & < 0.5 FOCUS pip
                 # focus_low_cutoff=0.5,
                 # ctwas_high_cutoff=0.9
  
  mismatch_gene_tissue_tracker <- tibble()
  for (cur_tiss in want_tissues){
    if (trait == "ad"){
      cur_tiss_names <- paste0(cur_tiss, c("_spxcan_twas_pvalue", "_highest_susie_pip"))
    } else {
      cur_tiss_names <- paste0(cur_tiss, c("_spxcan_twas_pvalue", "_highest_susie_pip", "_highest_pip"))
    }
    ct_num_genes <- gene_df_master %>%
      filter(!is.na(get(glue("{cur_tiss}_spxcan_twas_pvalue")))) %>%
      nrow()
    ct_bfr_thresh <- .05 / ct_num_genes

    if (is.null(want_genes)){
      ct_mismatch_df <- gene_df_master %>%
        filter(get(glue("{cur_tiss}_highest_susie_pip")) > ctwas_high_cutoff,
               get(glue("{cur_tiss}_highest_pip")) < focus_low_cutoff
               )
    } else {
      ct_mismatch_df <- gene_df_master %>%
        filter(ens_gene_id_no_decimal %in% want_genes)
    }

    ct_mismatch_df <- ct_mismatch_df %>%
      select(all_of(c("ens_gene_id_no_decimal", "chrom", cur_tiss_names))) %>% as_tibble() %>%
      mutate(`twas_sig` = get(cur_tiss_names[1]) < ct_bfr_thresh,
             tissue=cur_tiss,
             tiss_bfr_thresh=ct_bfr_thresh
             ) %>%
      left_join(gene_key_df, by = c("ens_gene_id_no_decimal" = "gene_id_no_decimal"))

    if (eqtl_or_sqtl == "eqtl"){
      ct_mismatch_df <- ct_mismatch_df %>%
        left_join(ctwas_res[[cur_tiss]] %>% select(all_of(c("id", "region_tag2"))), by = c("gene_id.v26" = "id"))
    } else {
      ct_mismatch_df <- ct_mismatch_df %>%
        left_join(ctwas_res[[cur_tiss]] %>% select(all_of(c("gene_id", "id", "region_tag2"))), by = c("gene_id.v26" = "gene_id"))
    }

    gene_df_master <- gene_df_master %>%
      mutate(!!glue("{cur_tiss}_twas_sig") := get(cur_tiss_names[1]) < ct_bfr_thresh,
             !!glue("{cur_tiss}_pip_mismatch") := ens_gene_id_no_decimal %in% ct_mismatch_df$ens_gene_id_no_decimal)
    
    num_bfr_sig <- sum(ct_mismatch_df %>% pull(get("twas_sig")))

    mismatch_gene_tissue_tracker <- bind_rows(mismatch_gene_tissue_tracker,
                                              ct_mismatch_df %>%
                                                select(any_of(c("id", "gene_id.v26",
                                                                "gene_name.v40",
                                                                "chrom",
                                                                "twas_sig",
                                                                "tiss_bfr_thresh",
                                                                "region_tag2",
                                                                "tissue")))
                                                  )
    print(glue("{cur_tiss} | {nrow(ct_mismatch_df)} ({num_bfr_sig} bfr sig in that tissue) with FOCUS pip < {focus_low_cutoff} and cTWAS pip > {ctwas_high_cutoff}\n"))
  }

  top_mismatch_genes <- mismatch_gene_tissue_tracker %>%
    group_by(`gene_id.v26`, `gene_name.v40`, `chrom`) %>%
    summarize(n_tiss_mismatch = n()) %>%
    arrange(desc(n_tiss_mismatch))
  print(glue("Got information on mismatches. (FOCUS pip < {focus_low_cutoff} and cTWAS pip > {ctwas_high_cutoff})"))
  
  rlist <- list(top_mismatch_genes = top_mismatch_genes,
                mismatch_gene_tissue_tracker = mismatch_gene_tissue_tracker,
                want_tissues=want_tissues,
                ctwas_res=ctwas_res,
                focus_low_cutoff=focus_low_cutoff,
                ctwas_high_cutoff=ctwas_high_cutoff
                )
  return(rlist)
}

get_ctwas_res <- function(study_name, want_tissues, eqtl_or_sqtl, gene_key_df,
                          zscore_dir="data/impute_expr_z/",
                          ctwas_dir="data/ctwas_rss/"){
  if (eqtl_or_sqtl == "sqtl"){
    gene_intron_df <- fread("data/combine_phenotype_groups.txt.gz") %>%
      select(all_of(c("gene_id", "intron_id"))) %>% as_tibble()
    gene_intron_df[c("tissue", "intron")] <- stringr::str_split_fixed(gene_intron_df$intron_id, pattern="\\.", 2)
  }
  gene_key_df <- gene_key_df %>%
    select(all_of(c("gene_id.v26", "gene_name.v40", "gene_type.v40"))) %>%
    rename(genename=`gene_name.v40`,
           gene_type=`gene_type.v40`)

  # Taking from /project2/guiming/xsun/proteome_alzheimer/8.brain_jansen_fusionwgts_ukbb_lasso/process_results.R
  study_prefix <- switch(study_name, "meta_analysis_BCAC_UKB_ovr_brca" = "meta_bcac_ukb_ovr",
                                     "AD_Bellenguez_GWAS_NG_2022" = "bellen",
                                     "AD_Wightmen_NG_2021" = "wight")
  sample_size <- switch(study_name, "meta_analysis_BCAC_UKB_ovr_brca" = 195650 + 228951,
                                     "AD_Bellenguez_GWAS_NG_2022" = 111326 + 677663,
                                     "AD_Wightmen_NG_2021" = 90338 + 1036225)

  ctwas_res <- list()
  num_want_tiss <- length(want_tissues)
  num_tiss <- 0
  for (cur_tiss in want_tissues){
    num_tiss <- num_tiss + 1
    print(glue("CTWAS: Tiss {num_tiss}/{num_want_tiss}"))
    cur_tiss_df <- tibble()
    for (cur_chr in 1:22){
      # Get ctwas results first
      cur_ctwas_res <- fread(file.path(ctwas_dir, glue("{study_prefix}__{cur_tiss}__{eqtl_or_sqtl}_chr{cur_chr}.susieIrss.txt"))) %>%
        mutate(PVE = susie_pip * mu2/sample_size,
               id_ensembl = str_extract(id, "(ENSG\\d{1,})"))
      if (eqtl_or_sqtl == "sqtl"){
        cur_tiss_key <- gene_intron_df %>%
          filter(tissue == !!cur_tiss)
        cur_ctwas_res <- cur_ctwas_res %>%
          left_join(cur_tiss_key %>% select(all_of(c("gene_id", "intron"))), by = c("id" = "intron"))
      } 
      cur_ctwas_res <- cur_ctwas_res %>%
        left_join(gene_key_df, by = c("id" = "gene_id.v26"))
      # Add zscore
      cur_zscore_res <- readRDS(file.path(zscore_dir, glue("{study_prefix}__{cur_tiss}__{eqtl_or_sqtl}_chr{cur_chr}.rds")))

      cur_ctwas_res <- left_join(cur_ctwas_res, cur_zscore_res$z_snp %>% select(all_of(c("id", "z"))), by = "id")
      cur_ctwas_res <- left_join(cur_ctwas_res, cur_zscore_res$z_gene %>% select(all_of(c("id", "z"))), by = "id") %>%
        mutate(z = pmax(z.x, z.y, na.rm=T))

      cur_ctwas_res <- cur_ctwas_res %>%
        select(any_of(c("chrom", "id", "gene_id", "pos", "type", "region_tag2",
                        "cs_index", "susie_pip", "mu2", "PVE", "genename", "id_ensembl", "gene_type", "z")))

      cur_tiss_df <- bind_rows(cur_tiss_df, cur_ctwas_res)
    }
      if (eqtl_or_sqtl == "eqtl"){
        cur_tiss_df <- cur_tiss_df %>%
          mutate(group = ifelse(type == "gene", "Expression", NA))
      } else if (eqtl_or_sqtl == "sqtl"){
        cur_tiss_df <- cur_tiss_df %>%
          mutate(group = ifelse(type == "gene", "Splicing", NA))
      }
    ctwas_res[[cur_tiss]] <- cur_tiss_df
  }


  # > names(tt)
  # [1] "chrom"       "id"          "pos"         "type"        "region_tag1"
  # [6] "region_tag2" "cs_index"    "susie_pip"   "mu2"         "region_tag" 
  # [11] "PVE"         "genename"    "id_ensembl"  "gene_type"   "z"     
  return(ctwas_res)
}

make_diagnostic_graphs <- function(study_prefix,
                                   tiss_list="data/trait_brca/tissue_lists/11tiss.txt",
                                   eqtl_or_sqtl="eqtl",
                                   chroms=1:22,
                                   plot_output_dir="docs/assets/",
                                   cache_output_dir="output/cache/",
                                   trait="brca",
                                   sample_size=NULL
                                   ){
  require(tidyverse)
  require(data.table)

  param_convergence <- function(study_prefix, tissue, chr,  ctwas_results_dir="data/ctwas_rss/",
                                thin=0.1, sample_size=NULL, plot_output_dir=NULL, cache_output_dir=NULL){
  library(ggplot2)
  library(cowplot)
  library(glue)
  library(data.table)
  library(tidyverse)
  study_long <- switch(study_prefix, 
                       "meta_bcac_ukb_ovr" = "meta_analysis_BCAC_UKB_ovr_brca",
                       "bellen" = "AD_Bellenguez_GWAS_NG_2022",
                       "wight" = "AD_Wightmen_NG_2021",
                       study_prefix)
  cur_p_grid <- runonce::save_run({    
    if (is.null(sample_size)){
      sample_size <- 195650 + 228951
    }
    
    load(file.path(ctwas_results_dir, glue("{study_prefix}__{tissue}__{eqtl_or_sqtl}_chr{chr}.s2.susieIrssres.Rd")))
    
    ctwas_res_s1 <- fread(file.path(ctwas_results_dir, glue("{study_prefix}__{tissue}__{eqtl_or_sqtl}_chr{chr}.s1.susieIrss.txt")))
    n_snps <- sum(ctwas_res_s1$type=="SNP")/thin
    
    ctwas_gene_res <- fread(file.path(ctwas_results_dir, glue("{study_prefix}__{tissue}__{eqtl_or_sqtl}_chr{chr}.susieIrss.txt"))) %>%
      filter(type == "gene")
    
    #estimated group prior (all iterations)
    if (!("SNP" %in% rownames(group_prior_rec))){
      rownames(group_prior_rec) <- c("gene", "SNP")
    }
    estimated_group_prior_all <- group_prior_rec
    estimated_group_prior_all["SNP",] <- estimated_group_prior_all["SNP",]*thin #adjust parameter to account for thin argument
    
    #estimated group prior variance (all iterations)
    if (!("SNP" %in% rownames(group_prior_var_rec))){
      rownames(group_prior_var_rec) <- c("gene", "SNP")
    }
    estimated_group_prior_var_all <- group_prior_var_rec
    
    #set group size
    group_size <- c(table(ctwas_gene_res$type), structure(n_snps, names="SNP"))
    group_size <- group_size[rownames(estimated_group_prior_all)]
    
    #estimated group PVE (all iterations)
    estimated_group_pve_all <- estimated_group_prior_var_all*estimated_group_prior_all*group_size/sample_size #check PVE calculation
    
    #estimated enrichment of genes (all iterations)
    estimated_enrichment_all <- t(sapply(rownames(estimated_group_prior_all)[rownames(estimated_group_prior_all)!="SNP"], function(x){estimated_group_prior_all[rownames(estimated_group_prior_all)==x,]/estimated_group_prior_all[rownames(estimated_group_prior_all)=="SNP"]}))
    
    title_size <- 11
    
    df <- data.frame(niter = rep(1:ncol(estimated_group_prior_all), nrow(estimated_group_prior_all)),
                     value = unlist(lapply(1:nrow(estimated_group_prior_all), function(x){estimated_group_prior_all[x,]})),
                     group = rep(rownames(estimated_group_prior_all), each=ncol(estimated_group_prior_all)))
    
    # df$group[df$group=="Liver"] <- "Liver_Expression"
    df$group <- as.factor(df$group)
    
    p_pi <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
      xlab("Iteration") + ylab(bquote(pi)) +
      ggtitle(glue("chr{chr}: Proportion Causal")) +
      theme_cowplot()
    
    p_pi <- p_pi + theme(plot.title=element_text(size=title_size)) + 
      expand_limits(y=0) + 
      guides(color = guide_legend(title = "Group")) + theme (legend.title = element_text(size=12, face="bold"))
    
    df <- data.frame(niter = rep(1:ncol(estimated_group_prior_var_all), nrow(estimated_group_prior_var_all)),
                     value = unlist(lapply(1:nrow(estimated_group_prior_var_all), function(x){estimated_group_prior_var_all[x,]})),
                     group = rep(rownames(estimated_group_prior_var_all), each=ncol(estimated_group_prior_var_all)))
    
    # df$group[df$group=="Liver"] <- "Liver_Expression"
    df$group <- as.factor(df$group)
    
    p_sigma2 <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
      xlab("Iteration") + ylab(bquote(sigma^2)) +
      ggtitle("Effect Size") +
      theme_cowplot()
    
    p_sigma2 <- p_sigma2 + theme(plot.title=element_text(size=title_size)) + 
      expand_limits(y=0) + 
      guides(color = guide_legend(title = "Group")) + theme (legend.title = element_text(size=12, face="bold"))
    
    df <- data.frame(niter = rep(1:ncol(estimated_group_pve_all), nrow(estimated_group_pve_all)),
                     value = unlist(lapply(1:nrow(estimated_group_pve_all), function(x){estimated_group_pve_all[x,]})),
                     group = rep(rownames(estimated_group_pve_all), each=ncol(estimated_group_pve_all)))
    
    # df$group[df$group=="Liver"] <- "Liver_Expression"
    df$group <- as.factor(df$group)
    
    p_pve <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
      xlab("Iteration") + ylab(bquote(h[G]^2)) +
      ggtitle("PVE") +
      theme_cowplot()
    
    p_pve <- p_pve + theme(plot.title=element_text(size=title_size)) + 
      expand_limits(y=0) + 
      guides(color = guide_legend(title = "Group")) + theme (legend.title = element_text(size=12, face="bold"))
    
    df <- data.frame(niter = rep(1:ncol(estimated_enrichment_all), nrow(estimated_enrichment_all)),
                     value = unlist(lapply(1:nrow(estimated_enrichment_all), function(x){estimated_enrichment_all[x,]})),
                     group = rep(rownames(estimated_enrichment_all), each=ncol(estimated_enrichment_all)))
    
    # df$group[df$group=="Liver"] <- "Liver_Expression"
    df$group <- as.factor(df$group)
    
    p_enrich <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
      xlab("Iteration") + ylab(bquote(pi[V]/pi[G])) +
      ggtitle("Enrichment") +
      theme_cowplot()
    
    p_enrich <- p_enrich + theme(plot.title=element_text(size=title_size)) + 
      expand_limits(y=0) + 
      guides(color = guide_legend(title = "Group")) + theme (legend.title = element_text(size=12, face="bold"))
 
    cur_p_grid <- plot_grid(p_pi, p_sigma2, p_enrich, p_pve, labels="AUTO")
    cur_p_grid},
    file.path(cache_output_dir, glue("ctwas_diag_{study_long}__{tissue}__{eqtl_or_sqtl}_chr{chr}.rds"))
  )
  save_plot(file.path(cache_output_dir, glue("ctwas_diag_{study_long}__{tissue}__{eqtl_or_sqtl}_chr{chr}.png")),
            cur_p_grid, ncol=2)
  # return(cur_p_grid)
  }
  want_tissues <- fread(tiss_list, header=F)[["V1"]]
  for (cur_tiss in want_tissues){
    for (chr in chroms){
      param_convergence(study_prefix, cur_tiss, chr, plot_output_dir=plot_output_dir,
                        cache_output_dir=cache_output_dir, ctwas_results_dir=glue("data/trait_{trait}/ctwas_rss/"),
                        sample_size=sample_size)
    }
  }
}

make_ctwas_vs_twas_graph <- function(study,
                                     tiss_list="data/tissue_lists/11tiss.txt"
                                     ){
  want_tissues <- fread(tiss_list, header=F)[["V1"]]
}


ad_bell_eqtl_want <- c("ENSG00000189298",
"ENSG00000164938",
"ENSG00000174226",
"ENSG00000159314",
"ENSG00000225190")

ad_bell_sqtl_want <- c("ENSG00000197062",
"ENSG00000229391",
"ENSG00000156170", 
"ENSG00000149932", 
"ENSG00000214425", 
"ENSG00000204650", 
"ENSG00000266903")

ad_wight_sqtl_want <- c("ENSG00000196821")


library(glue)
library(tidyverse)
WANT_STUDIES <- c("meta_analysis_BCAC_UKB_ovr_brca")
if (interactive()){
  # make_diagnostic_graphs("meta_bcac_ukb_ovr")
  # prepped_data <- list()
  # for (study in WANT_STUDIES){
  #   # make_ctwas_vs_twas_graph(study)
  #   prepped_data[[study]] <- runonce::save_run({main(study)}, glue("output/cache/{study}__eqtl.rds"))
  # }

  # AD
  # study <- "AD_Bellenguez_GWAS_NG_2022" 
  # runonce::save_run({main(study,
  #     tiss_list="data/trait_ad/tissue_lists/12_of_13_brain_tissues.txt",
  #     trait="ad", want_genes=ad_bell_eqtl_want)}, glue("output/cache/{study}__eqtl.rds"))
  # make_diagnostic_graphs("bellen", tiss_list="data/trait_ad/tissue_lists/12_of_13_brain_tissues.txt",
  #                        trait="ad", sample_size=111326 + 677663)
  # study <- "AD_Bellenguez_GWAS_NG_2022" 
  # runonce::save_run({main(study,
  #     tiss_list="data/trait_ad/tissue_lists/12_of_13_brain_tissues.txt",
  #     eqtl_or_sqtl="sqtl", 
  #     trait="ad", want_genes=ad_bell_sqtl_want)}, glue("output/cache/{study}__sqtl.rds"))
  # make_diagnostic_graphs("bellen", eqtl_or_sqtl="sqtl", tiss_list="data/trait_ad/tissue_lists/12_of_13_brain_tissues.txt",
  #                        trait="ad", sample_size=111326 + 677663)

  study <- "AD_Wightmen_NG_2021" 
  runonce::save_run({main(study,
      tiss_list="data/trait_ad/tissue_lists/12_of_13_brain_tissues.txt",
      eqtl_or_sqtl="sqtl", 
      trait="ad", want_genes=ad_wight_sqtl_want)}, glue("output/cache/{study}__sqtl.rds"))
  make_diagnostic_graphs("wight", eqtl_or_sqtl="sqtl", tiss_list="data/trait_ad/tissue_lists/12_of_13_brain_tissues.txt",
                         trait="ad", sample_size=90338 + 1036225)
} else {
}
