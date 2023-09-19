do_tissue <- function(wt_idx, study, eqtl_or_sqtl="eqtl"){
  library(tidyverse)
  library(glue)
  library(cowplot)
  library(purrr)
  # First a table showing all mismatching genes
  do_chrom_region <- function(chrom, region_tag2){
    # First the locus plot
    cur_tiss_chrom_region_mismatch_tracker <- cur_tiss_chrom_mismatch_tracker %>%
      filter(region_tag2 == !!region_tag2)

    ctwas_res_filt <- ctwas_res[[want_tissues[wt_idx]]] %>%
      filter(chrom == !!chrom, region_tag2 == !!region_tag2)
    region_tag_pos_min <- ctwas_res_filt %>% pull(pos) %>% min(., na.rm=T)
    region_tag_pos_max <- ctwas_res_filt %>% pull(pos) %>% max(., na.rm=T)
    
    pre_fig_table <- DT::datatable(cur_tiss_chrom_region_mismatch_tracker %>% 
                  select(all_of(c("gene_id.v26", "gene_name.v40", "twas_sig", "tiss_bfr_thresh"))) %>% unique(),
              caption=glue("Chr{chrom} | ~ {region_tag_pos_min}-{region_tag_pos_max}"))

    cat("\n\n")
    cat(knitr::knit_print(pre_fig_table))
    cat("\n\n")
  
    if (eqtl_or_sqtl == "sqtl"){
      label_genes <- cur_tiss_chrom_region_mismatch_tracker %>%
        pull(`id`) %>% unique()
    } else {
      label_genes <- cur_tiss_chrom_region_mismatch_tracker %>%
        pull(`gene_id.v26`) %>% unique()
    }
  
    if (eqtl_or_sqtl == "sqtl"){
      plot_eqtl <- c()
    } else {
      plot_eqtl <- NULL
    }
      a <- locus_plot(study,
                      want_tissues[wt_idx],
                      ctwas_res_filt,
                      chrom=chrom,
                      region_tag2=region_tag2,
                      eqtl_or_sqtl=eqtl_or_sqtl,
                      gwide_tiss_sig_thresh=gwide_tiss_sig_thresh,
                      return_table=F,
                      focus=NULL,
                      label_genes=label_genes,
                      plot_eqtl = plot_eqtl,
                      rerun_ctwas=F,
                      rerun_load_only=F,
                      label_panel="both",
                      legend_side="left",
                      legend_panel="",
                      draw_gene_track = F)
    cat("\n\n")
  }
  cur_tiss_mismatch_tracker <- mismatch_gene_tissue_tracker %>%
    filter(tissue == want_tissues[wt_idx])
  
  gwide_tiss_sig_thresh <- cur_tiss_mismatch_tracker %>%
    filter(tissue == want_tissues[wt_idx]) %>%
    pull(tiss_bfr_thresh) %>%
    max(., na.rm=T)
  
  cat("\n\n")
  cat(knitr::knit_print(DT::datatable(cur_tiss_mismatch_tracker %>%
                  arrange(chrom, region_tag2) %>%
                  select(all_of(c("gene_id.v26", "gene_name.v40", "chrom", "twas_sig"))) %>% unique(),
                            caption=glue("{want_tissues[wt_idx]}: FOCUS PIP < {focus_low_cutoff} & cTWAS PIP > {ctwas_high_cutoff}")
              )))
  cat("\n\n")

  cur_tiss_chrom_region_df <- cur_tiss_mismatch_tracker %>%
    group_by(chrom, region_tag2) %>%
    summarize()

  for (cur_chr in cur_tiss_chrom_region_df %>% pull(chrom) %>% unique() %>% sort()){
    cur_diag_graph <- readRDS(glue("output/cache/ctwas_diag_{study}__{want_tissues[wt_idx]}__{eqtl_or_sqtl}_chr{cur_chr}.rds"))
    # print(glue("Chromosome {cur_chr}. {want_tissues[wt_idx]}"))
    cat("\n")
    cat(glue("### Chromosome {cur_chr}"))
    cat("\n")
    print(cur_diag_graph)
  
    cur_tiss_chrom_mismatch_tracker <- cur_tiss_mismatch_tracker %>%
      filter(chrom == !!cur_chr)
  
    cur_tiss_chrom_region_df %>%
      filter(chrom == !!cur_chr, !is.na(region_tag2)) %>%
      pwalk(do_chrom_region)
  }
}
