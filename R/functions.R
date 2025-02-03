#' Helping function to get started 
#'
#' @returns
#' @export
#'
#' @examples
hello_I_want_to_analyse_my_PELSA_data <- function() {
  
  message("\nHi, no worries! We'll get this done in no time.")
  
  Sys.sleep(3)
  
  message("\nLet's start with selecting the DIA-NN report output file.")
  
  Sys.sleep(3.5)
  
  message("\nWe'll start in")
  
  Sys.sleep(1.5)
  
  message("3")
  
  Sys.sleep(1)
  
  message("2")
  
  Sys.sleep(0.75)
  
  message("1")
  
  Sys.sleep(0.5)
  
  file <- choose.files(multi = F)
  
  message("\nGood job! Now you can copy the code below and continue with importing your data.\n")
  
  Sys.sleep(2)
  to_cat <- paste0('data_raw <- import_diann_report("', 
                   gsub("\\\\", "\\\\\\\\", file), 
                   '")')
  
  for (i in unlist(strsplit(to_cat, split = ""))) {
    cat(i)
    Sys.sleep(0.001)
  }
  # cat(paste0('data_raw <- import_diann_report("', 
  #            gsub("\\\\", "\\\\\\\\", file), 
  #            '")'))
  
}


#' Import DIA-NN report 
#'
#' @param file 
#' 
#' @importFrom magrittr %>%
#' @returns
#' @export
#'
#' @examples
import_diann_report <- function(file) {
  
  # Check file 
  if (!hasArg(file)) 
    file <- choose.files()
  
  # Load .parquet or .tsv
  if (tools::file_ext(file) == "parquet") 
    data_import <- arrow::read_parquet(file)
  else if (tools::file_ext(file) == "tsv") 
    data_import <- vroom::vroom(file)
  
  
  # Identify file names 
  sample_names <- data_import$Run %>% 
    unique() %>% 
    sort()
  
  # Remove longest common prefix 
  sample_names <- sample_names %>% 
    setNames(., stringr::str_remove(., .longest_common_prefix(.)))
  
  sample_groups <- sample_names %>% 
    setNames(., rep("Control", length(sample_names)))
  
  # Suggest process function
  samples_n <- paste0(names(sample_names), '" = "', sample_names, '"')
  samples_v <- paste0('c("', paste(samples_n, collapse = ',\n\t"'), ')')
  
  groups_n <- paste0(names(sample_groups), '" = "', sample_groups, '"')
  groups_v <- paste0('c("', paste(groups_n, collapse = ',\n\t"'), ')')
  
  message("Your data was imported successfully.")
  
  Sys.sleep("3")
  
  message("\nPlease copy the code below to prepare the precursor data for limma analysis.\n")
  
  Sys.sleep("2")
  
  to_cat <- paste0('data_processed <- process_diann_report(data_raw,\n', 
                   "\t# You can change the name for each sample here below\n", 
                   '\tsample_names = ', samples_v, ',\n', 
                   "\t# You need to change the groups for each sample here below\n", 
                   '\tsample_groups = ', groups_v, ',\n',
                   '\tquant_column = "Precursor.Normalised",\n', 
                   '\tmin_val_per_group = 1,\n', 
                   '\tnorm.method = "none",\n', 
                   '\tLib.Q.Value = 0.01,\n', 
                   '\tLib.PG.Q.Value = 0.01)')
  
  for (i in unlist(strsplit(to_cat, split = ""))) {
    cat(i)
    Sys.sleep(0.001)
  }
  
  # cat(paste0('data_processed <- process_diann_report(data_raw,\n', 
  #            '\tsample_names = ', samples_v, ',\n', 
  #            '\tsample_groups = ', groups_v, ',\n',
  #            '\tquant_column = "Precursor.Normalised",\n', 
  #            '\tmin_val_per_group = 1,\n', 
  #            '\tnorm.method = "none",\n', 
  #            '\tLib.Q.Value = 0.01,\n', 
  #            '\tLib.PG.Q.Value = 0.01)'))
  
  message("\n\nYou can change the names of the proposed <sample_names> vector. The names of the <group_names> vector HAVE TO be changed to match the experimental design.")
  
  return(data_import)
  
}




#' Processs the imported DIA-NN precursor data 
#'
#' @param data_raw 
#' @param sample_names 
#' @param sample_groups 
#' @param peptide.id 
#' @param quant_column 
#' @param min_val_per_group 
#' @param norm.method 
#' @param proteotypic.only 
#' @param Q.Value 
#' @param PG.Q.Value 
#' @param Lib.Q.Value 
#' @param Lib.PG.Q.Value 
#' @param protein.q 
#' @param gg.q 
#'
#' @importFrom magrittr %>%
#' @returns
#' @export
#'
#' @examples
process_diann_report <- function(data_raw, 
                                 sample_names, 
                                 sample_groups, 
                                 peptide.id = "Stripped.Sequence", 
                                 quant_column = "Precursor.Normalised", 
                                 min_val_per_group = 1, 
                                 norm.method = "none", 
                                 proteotypic.only = F, 
                                 Q.Value = 1, 
                                 PG.Q.Value = 1, 
                                 Lib.Q.Value = 1, 
                                 Lib.PG.Q.Value = 1, 
                                 protein.q = 1, 
                                 gg.q = 1) {
  
  # Check input
  if (!hasArg(data_raw)) 
    stop("Please provide an imported DIA-NN report as <data_raw>.")
  
  # Check input
  if (all(c(Q.Value, 
            PG.Q.Value, 
            Lib.Q.Value, 
            Lib.PG.Q.Value, 
            protein.q, 
            gg.q) == 0)) 
    stop("All your q-value cutoffs are 1, please change them according to your study.")
  
  # Check groups 
  if (length(unique(names(sample_groups))) < 2) 
    stop("Please define the sample groups by changing the names of the <sample_groups> argument vector. There must by minimum two groups.")
  
  
  # Filter precursors by proteotypicity
  if (proteotypic.only) 
    data_raw_filtered <- data_raw %>% 
      dplyr::filter(Proteotypic != 0)
  else 
    data_raw_filtered <- data_raw
  
  # Filter precursors by Q-values
  data_raw_filtered <- data_raw_filtered %>% 
    dplyr::filter(Q.Value <= Q.Value, 
                  PG.Q.Value <= PG.Q.Value, 
                  Lib.Q.Value <= Lib.Q.Value, 
                  Lib.PG.Q.Value <= Lib.PG.Q.Value)
  
  # Extract quantitative data
  
  # Rename runs 
  data_raw_filtered <- data_raw_filtered %>% 
    dplyr::mutate(Run = setNames(names(sample_names), unname(sample_names))[Run]) %>% 
    dplyr::arrange(Run)
  
  # Summarise precursors to modified peptides 
  data_quant <- data_raw_filtered %>%  
    dplyr::summarise(!!quant_column := sum(!!rlang::sym(quant_column)), 
                     .by = c("Run", peptide.id)) 
  
  # Add group information 
  data_quant <- data_quant %>% 
    dplyr::mutate(Group = setNames(sample_groups, names(sample_names))[Run], 
                  .after = "Run") %>% 
    dplyr::mutate(n = sum(!!rlang::sym(quant_column) > 0), .by = c(peptide.id, "Group")) %>% 
    dplyr::mutate(p = n / max(n), .by = "Group") %>% 
    dplyr::filter(all(p >= min_val_per_group), .by = peptide.id) %>% ##### Check this 
    dplyr::mutate(n_peptide = length(Run), .by = peptide.id) %>% 
    dplyr::filter(n_peptide == max(n_peptide)) %>% 
    tidyr::pivot_wider(id_cols = "Run", 
                       names_from = peptide.id, 
                       values_from = quant_column)
  
  
  # Check sample groups 
  sample_groups_m <- names(sample_groups)[match(sample_names[data_quant$Run], sample_groups)]
  
  
  list_output <- list(data_raw_filtered = data_raw_filtered, 
                      data_processed = data_quant, 
                      sample_groups = sample_groups_m)
  
  
  
  # Message 
  message("/nYour data was successfully processed and can be used for the limma analysis. Please copy the code below./n/n")
  
  to_cat <- paste0('data_limma <- limma_diann_report(data_processed,\n', 
                   'conditions = c("', 
                   unique(sample_groups_m)[1], 
                   '", "', 
                   unique(sample_groups_m)[2], 
                   '"),\n', 
                   '\tp.adjust.method = "BH",\n', 
                   '\tp.threshold = 0.01,\n', 
                   '\tfc.threshold = log2(1.2))')
  
  for (i in unlist(strsplit(to_cat, split = ""))) {
    cat(i)
    Sys.sleep(0.001)
  }
  
  # Return output list 
  return(list_output)
}



#' Limma function 
#'
#' @param data_processed 
#' @param conditions 
#' @param p.adjust.method 
#' @param p.threshold 
#' @param fc.threshold 
#'
#' @importFrom magrittr %>%
#' @returns
#' @export
#'
#' @examples
limma_diann_report <- function(data_processed, 
                               conditions = c("Control", "Treatment"), 
                               p.adjust.method = "BH", 
                               p.threshold = 0.01, 
                               fc.threshold = log2(1.2)) {
  
  if (hasArg(data_processed)) {
    data_quant <- data_processed[["data_processed"]]
    sample_groups <- data_processed[["sample_groups"]]
  } else {
    stop("Please use the output of the process_diann_report() function as input for <data_processed>.")
  }
  
  
  # Add limma input 
  results_list <- list()
  
  ## Add base data frame 
  results_list[["data"]] <- data_quant %>% 
    tidyr::pivot_longer(-1) %>% 
    dplyr::mutate(value = log2(value)) %>% 
    tidyr::pivot_wider()
  
  ## Add eset 
  results_list[["eset"]] <- results_list[["data"]] %>%
    dplyr::rename(rowname = 1) %>% 
    tibble::column_to_rownames() %>%
    as.matrix() %>% 
    t() %>% 
    Biobase::ExpressionSet()
  
  ## Describe experimental groups
  sample_groups_mod <- setNames(c("Control", "Treatment"), conditions)[sample_groups]
  
  results_list[["design"]] <- 
    model.matrix(
      ~0+factor(sample_groups_mod, 
                unique(sample_groups_mod)))
  
  ## Correct column names
  colnames(results_list[["design"]]) <- unique(sample_groups_mod)
  
  
  # Do first linear model fit
  results_list[["fit"]] <- limma::lmFit(results_list[["eset"]], results_list[["design"]])
  
  # Describe comparisons (contrasts)
  results_list[["contrast.matrix"]] <- limma::makeContrasts(Treatment-Control, 
                                                            levels=results_list[["design"]])
  
  # Compute contrasts between groups
  results_list[["fit2"]] <- limma::contrasts.fit(fit = results_list[["fit"]], 
                                                 contrasts = results_list[["contrast.matrix"]])
  
  # Apply empirical ebayes to estimate t-statistic and calculate p-values
  results_list[["fit2_eBayes"]] <- limma::eBayes(results_list[["fit2"]])
  
  
  # Define thresholds 
  results_list[["par"]] <- list(p.threshold = p.threshold, 
                                p.adjust.method = p.adjust.method, 
                                fc.threshold = fc.threshold)
  
  # Extract results 
  results_list[["results"]] <- biobroom::tidy.MArrayLM(results_list[["fit2_eBayes"]]) %>% 
    # rename contrasts 
    dplyr::mutate(term = paste(rev(conditions), collapse = " - ")) %>% 
    dplyr::rename(id = gene) %>% 
    dplyr::arrange(p.value) %>% 
    # Adjust p-values
    dplyr::mutate(p.adjust = 
                    p.adjust(p.value, 
                             method = results_list[["par"]]$p.adjust.method), 
                  .after = "p.value") %>% 
    # Apply thresholds
    dplyr::mutate(regulation = dplyr::case_when(
      p.adjust < results_list[["par"]]$p.threshold & 
        estimate >= results_list[["par"]]$fc.threshold ~ "up", 
      p.adjust < results_list[["par"]]$p.threshold & 
        estimate <= - results_list[["par"]]$fc.threshold ~ "down", 
      .default = "none"
    ))
  
  require(ggplot2)
  
  results_list[["p"]] <- results_list[["results"]] %>%
    dplyr::filter(!is.na(regulation)) %>%
    ggplot(aes(x = estimate,
               y = -log10(p.value),
               color = regulation)) +
    geom_point(shape = 16, alpha = 0.5) +
    scale_color_manual(values = c(none = "grey",
                                  up = "red",
                                  down = "blue")) +
    theme_classic() +
    coord_cartesian(expand = F)
  
  print(results_list[["p"]])
  
  return(results_list)
  
}




#' Title
#'
#' @param strs 
#'
#' @returns
#' @export
#'
#' @examples
.longest_common_prefix <- function(strs) {
  if (length(strs) == 0) return("")
  
  min_len <- min(nchar(strs))
  if (min_len == 0) return("")
  
  prefix <- character(0)
  
  for (i in 1:min_len) {
    current_char <- substr(strs[1], i, i)
    if (all(substr(strs, i, i) == current_char)) {
      prefix <- c(prefix, current_char)
    } else {
      break
    }
  }
  
  return(paste(prefix, collapse = ""))
}
