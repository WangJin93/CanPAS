#' @title Print Method for cpas_COX_by_genes Class
#' @description Custom print method for the result object returned by COX_by_genes().
#' @param x An object of class cpas_COX_by_genes.
#' @param ... Additional arguments passed to print.
#' @return No return value, prints a summary of the univariable Cox regression results.
#' @export
#' @examples
#' \dontrun{
#'    ## A real per-gene screen on a real cohort.
#'    d <- cohort_merged("GSE14814", c("GAPDH", "ACTB", "TP53"), type = "OS")
#'    r <- COX_by_genes(d, type = "OS", genes = c("GAPDH", "ACTB", "TP53"))
#'    print(r)
#' }
print.cpas_COX_by_genes <- function(x, ...) {
  cat("\nUnivariate Cox Regression Analysis Results")
  cat("\n=======================================")
  
  # Print analysis metadata
  cat("\n\nAnalysis Information:")
  cat(sprintf("\n- Survival Type: %s", x$metadata$survival_type))
  cat(sprintf("\n- Number of Genes Analyzed: %d", x$metadata$gene_count))
  cat(sprintf("\n- Sample Size: %d", x$metadata$sample_size))
  cat(sprintf("\n- Number of Events: %d", x$metadata$events))
  cat(sprintf("\n- Analysis Time: %s", format(x$input_params$analysis_time)))
  
  # Print input parameters
  cat("\n\nInput Parameters:")
  cat(sprintf("\n- Data Frame: %s", x$input_params$df_name))
  cat(sprintf("\n- Genes: %s", paste(x$input_params$markers, collapse = ", ")))
  
  # Print main results
  cat("\n\nMain Results (Top 10 Genes by P-value):")
  cat("\n=====================================")
  
  # Sort results by p-value
  sorted_results <- x$results_table[order(x$results_table$Pvalue), ]
  
  # Print only the first 10 results
  print(utils::head(sorted_results, 10), row.names = FALSE, ...)
  
  # If there are more than 10 results, indicate that
  if (nrow(sorted_results) > 10) {
    cat(sprintf("\n[... %d more genes not shown]", nrow(sorted_results) - 10))
  }
  
  # Print additional information about the object structure
  cat("\n\nObject Structure:")
  cat("\n- $input_params: Input parameters used for analysis")
  cat("\n- $processed_data: Processed data frame used for analysis")
  cat("\n- $individual_models: List of individual Cox regression models")
  cat("\n- $individual_summaries: List of model summaries")
  cat("\n- $results_table: Complete results table")
  cat("\n- $metadata: Additional analysis metadata")
  
  cat("\n\nUse str(object) to view full object structure.")
  cat("\n=======================================")
  
  invisible(x)
}

#' @title Print Method for cpas_COX Class
#' @description Custom print method for the result object returned by COX_analysis function.
#' @param x An object of class cpas_COX.
#' @param ... Additional arguments passed to print.
#' @return No return value, prints a summary of the COX regression results.
#' @export
#' @examples
#' \dontrun{
#'    ## A real Cox fit: expression + survival + clinical covariates.
#'    d <- cohort_merged("GSE13507", "GAPDH", type = "OS", clin = TRUE)
#'    r <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
#'                       cate_Variates = c("grade", "N"), method = "uni")
#'    print(r)
#' }
print.cpas_COX <- function(x, ...) {
  cat("\nCOX Regression Analysis Results")
  cat("\n================================")
  
  # Print analysis metadata
  cat("\n\nAnalysis Information:")
  cat(sprintf("\n- Analysis Type: %s", x$metadata$analysis_type))
  cat(sprintf("\n- Survival Type: %s", x$metadata$survival_type))
  cat(sprintf("\n- Sample Size: %d", x$metadata$sample_size))
  cat(sprintf("\n- Number of Events: %d", x$metadata$events))
  cat(sprintf("\n- Analysis Time: %s", format(x$input_params$analysis_time)))
  
  # Print input parameters
  cat("\n\nInput Parameters:")
  cat(sprintf("\n- Data Frame: %s", x$input_params$df_name))
  cat(sprintf("\n- Method: %s", x$input_params$method))
  
  if (length(x$input_params$cont_Variates) > 0) {
    cat(sprintf("\n- Continuous Variables: %s", paste(x$input_params$cont_Variates, collapse = ", ")))
  }
  if (length(x$input_params$cate_Variates) > 0) {
    cat(sprintf("\n- Categorical Variables: %s", paste(x$input_params$cate_Variates, collapse = ", ")))
  }
  
  # Print main results
  cat("\n\nMain Results:")
  cat("\n=============")
  print(x$results_table, row.names = FALSE, ...)
  
  # Print additional information about the object structure
  cat("\n\nObject Structure:")
  cat("\n- $input_params: Input parameters used for analysis")
  cat("\n- $processed_data: Processed data frame used for analysis")
  cat("\n- $models: List of COX regression models")
  cat("\n- $summaries: List of model summaries")
  cat("\n- $results_table: Complete results table")
  cat("\n- $metadata: Additional analysis metadata")
  
  cat("\n\nUse str(object) to view full object structure.")
  cat("\n================================")
  
  invisible(x)
}

#' @title Print Method for cpas_get Class
#' @description Custom print method for the result object returned by get_data function.
#' @param x An object of class cpas_get.
#' @param ... Additional arguments passed to print.
#' @return No return value, prints a summary of the API response.
#' @export
#' @examples
#' \dontrun{
#'    ## The raw API accessor, printed.
#'    r <- get_data("GSE14814", "surv_data")
#'    print(r)
#' }
print.cpas_get <- function(x, ...) {
  cat("\nAPI Data Retrieval Results")
  cat("\n========================")
  
  # Print request metadata
  cat("\n\nRequest Information:")
  cat(sprintf("\n- Table: %s", x$input_params$table))
  cat(sprintf("\n- Action: %s", x$input_params$action))
  cat(sprintf("\n- Request Time: %s", format(x$input_params$request_time)))
  cat(sprintf("\n- Response Time: %s", format(x$metadata$response_time)))
  
  # Print API URL
  cat("\n\nAPI URL:")
  cat(sprintf("\n%s", x$url))
  
  # Print response information
  cat("\n\nResponse Information:")
  cat(sprintf("\n- API Version: %s", x$metadata$api_version))
  cat(sprintf("\n- Response Format: %s", x$metadata$response_format))
  cat(sprintf("\n- Result Count: %d", x$metadata$result_count))
  
  # Print response preview
  cat("\n\nResponse Preview:")
  cat("\n================")
  print(utils::head(x$response), ...)
  
  # Print additional information about the object structure
  cat("\n\nObject Structure:")
  cat("\n- $input_params: Input parameters used for API request")
  cat("\n- $response: Raw response from the API")
  cat("\n- $url: Constructed API URL")
  cat("\n- $metadata: Request metadata")
  
  cat("\n\nUse str(object) to view full object structure.")
  cat("\n========================")
  
  invisible(x)
}

#' @title Print Method for cpas_get_expr Class
#' @description Custom print method for the result object returned by get_expr_data function.
#' @param x An object of class cpas_get_expr.
#' @param ... Additional arguments passed to print.
#' @return No return value, prints a summary of the gene expression data retrieval.
#' @export
#' @examples
#' \dontrun{
#'    e <- get_expr_data("GSE14814", c("GAPDH", "ACTB"), process_duplicates = "max")
#'    print(e)
#' }
print.cpas_get_expr <- function(x, ...) {
  cat("\nGene Expression Data Retrieval Results")
  cat("\n=====================================")
  
  # Print analysis metadata
  cat("\n\nDataset Information:")
  cat(sprintf("\n- Dataset Accession: %s", x$metadata$dataset_accession))
  cat(sprintf("\n- Platform: %s", x$metadata$platform))
  cat(sprintf("\n- Analysis Time: %s", format(x$input_params$analysis_time)))
  
  # Print gene information
  cat("\n\nGene Information:")
  cat(sprintf("\n- Total Genes Requested: %d", x$metadata$total_genes_requested))
  cat(sprintf("\n- Genes Found: %d", x$metadata$genes_found))
  cat(sprintf("\n- Genes Matched to Platform: %d", x$metadata$genes_matched_platform))
  cat(sprintf("\n- Duplicate Handling: %s", x$input_params$process_duplicates))
  
  # Print sample information
  cat("\n\nSample Information:")
  cat(sprintf("\n- Number of Samples: %d", x$metadata$sample_count))
  
  # Print expression data preview
  cat("\n\nExpression Data Preview:")
  cat("\n======================")
  print(utils::head(x$expr_data), ...)
  
  # Print reference IDs preview
  cat("\n\nReference IDs Preview:")
  cat("\n====================")
  print(utils::head(x$ref_ids), ...)
  
  # Print additional information about the object structure
  cat("\n\nObject Structure:")
  cat("\n- $input_params: Input parameters used for analysis")
  cat("\n- $raw_ids: Raw reference IDs retrieved from ID_map")
  cat("\n- $platform_info: Platform information for the dataset")
  cat("\n- $ref_ids: Processed reference IDs matching genes and platform")
  cat("\n- $expr_data: Processed expression data")
  cat("\n- $metadata: Analysis metadata")
  
  cat("\n\nUse str(object) to view full object structure.")
  cat("\n=====================================")
  
  invisible(x)
}

#' @title Print Method for cpas_signature Class
#' @description Custom print method for the result object returned by get_signature_value function.
#' @param x An object of class cpas_signature.
#' @param ... Additional arguments passed to print.
#' @return No return value, prints a summary of the gene signature calculation.
#' @export
#' @examples
#' \dontrun{
#'    s <- get_signature_value("0.5*GAPDH + 0.5*ACTB", "GSE14814")
#'    print(s)
#' }
print.cpas_signature <- function(x, ...) {
  cat("\nGene Signature Calculation Results")
  cat("\n================================")
  
  # Print calculation metadata
  cat("\n\nCalculation Information:")
  cat(sprintf("\n- Dataset: %s", x$metadata$dataset))
  cat(sprintf("\n- Analysis Time: %s", format(x$input_params$analysis_time)))
  
  # Print signature information
  cat("\n\nSignature Information:")
  cat(sprintf("\n- Original Signature: %s", x$signature_info$original_signature))
  cat(sprintf("\n- Number of Genes: %d", x$metadata$signature_length))
  cat(sprintf("\n- Duplicate Handling: %s", x$input_params$process_duplicates))
  
  # Print weight-gene pairs
  cat("\n\nWeight-Gene Pairs:")
  cat("\n=================")
  print(x$signature_info$weight_gene_pairs, row.names = FALSE, ...)
  
  # Print signature values preview
  cat("\n\nSignature Values Preview:")
  cat("\n=======================")
  print(utils::head(x$results_table), ...)
  
  # Print additional information about the object structure
  cat("\n\nObject Structure:")
  cat("\n- $input_params: Input parameters used for analysis")
  cat("\n- $signature_info: Detailed information about the gene signature")
  cat("\n- $expr_data: Processed expression data used for calculation")
  cat("\n- $results_table: Sample IDs and their corresponding signature values")
  cat("\n- $metadata: Analysis metadata")
  
  cat("\n\nUse str(object) to view full object structure.")
  cat("\n================================")
  
  invisible(x)
}

#' @title Print Method for cpas_merge Class
#' @description Custom print method for the result object returned by merge_surv_expr function.
#' @param x An object of class cpas_merge.
#' @param ... Additional arguments passed to print.
#' @return No return value, prints a summary of the merged data.
#' @export
#' @examples
#' \dontrun{
#'    e <- get_expr_data("GSE14814", c("GAPDH", "ACTB"), process_duplicates = "max")
#'    m <- merge_surv_expr("GSE14814", e)
#'    print(m)
#' }
print.cpas_merge <- function(x, ...) {
  cat("\nSurvival and Expression Data Merge Results")
  cat("\n=========================================")
  
  # Print merge metadata
  cat("\n\nMerge Information:")
  cat(sprintf("\n- Dataset Accession: %s", x$metadata$dataset_accession))
  cat(sprintf("\n- Merge Time: %s", format(x$input_params$merge_time)))
  
  # Print sample matching information
  cat("\n\nSample Matching:")
  cat(sprintf("\n- Samples in Survival Data: %d", x$metadata$samples_in_survival))
  cat(sprintf("\n- Samples in Expression Data: %d", x$metadata$samples_in_expression))
  cat(sprintf("\n- Samples After Merge: %d", x$metadata$samples_merged))
  
  # Print column information
  cat("\n\nColumn Information:")
  cat(sprintf("\n- Survival Columns: %s", paste(x$metadata$survival_columns, collapse = ", ")))
  cat(sprintf("\n- Genes in Expression: %d", x$metadata$genes_in_expression))
  cat(sprintf("\n- Total Columns After Merge: %d", length(x$metadata$merged_columns)))
  
  # Print merged data preview
  cat("\n\nMerged Data Preview:")
  cat("\n===================")
  print(utils::head(x$merged_data), ...)
  
  # Print additional information about the object structure
  cat("\n\nObject Structure:")
  cat("\n- $input_params: Input parameters used for merging")
  cat("\n- $raw_surv_data: Raw survival data retrieved from the database")
  cat("\n- $raw_expr_data: Raw expression data provided as input")
  cat("\n- $merged_data: Merged survival and gene expression data")
  cat("\n- $metadata: Merge metadata")
  
  cat("\n\nUse str(object) to view full object structure.")
  cat("\n=========================================")
  
  invisible(x)
}

#' @title Print Method for cpas_COX_by_datasets Class
#' @description Custom print method for the result object returned by COX_by_datasets().
#' @param x An object of class cpas_COX_by_datasets.
#' @param ... Additional arguments passed to print.
#' @return No return value, prints a summary of the multi-dataset COX analysis.
#' @export
#' @examples
#' \dontrun{
#'    r <- COX_by_datasets(c("GSE14814", "GSE31210"), gene = "GAPDH", type = "OS")
#'    print(r)
#' }
print.cpas_COX_by_datasets <- function(x, ...) {
  cat("\nMulti-Dataset COX Regression Analysis Results")
  cat("\n==========================================")
  
  # Print analysis metadata
  cat("\n\nAnalysis Information:")
  cat(sprintf("\n- Number of Datasets: %d", x$metadata$total_datasets))
  cat(sprintf("\n- Successful Datasets: %d", x$metadata$successful_datasets))
  cat(sprintf("\n- Failed Datasets: %d", x$metadata$failed_datasets))
  cat(sprintf("\n- Gene Analyzed: %s", x$metadata$gene_analyzed))
  cat(sprintf("\n- Survival Type: %s", x$metadata$survival_type))
  cat(sprintf("\n- Analysis Time: %s", format(x$input_params$analysis_time)))
  
  # Print input parameters
  cat("\n\nInput Parameters:")
  cat(sprintf("\n- Datasets: %s", paste(x$input_params$tables, collapse = ", ")))
  cat(sprintf("\n- Precision: %d decimal places", x$metadata$precision))
  
  # Print failed datasets if any
  if (x$metadata$failed_datasets > 0) {
    cat("\n\nFailed Datasets:")
    cat("\n==============")
    for (ds in names(x$errors)) {
      cat(sprintf("\n- %s: %s", ds, x$errors[[ds]]))
    }
  }
  
  # Print combined results
  cat("\n\nCombined Results:")
  cat("\n================")
  if (nrow(x$combined_results) > 0) {
    print(x$combined_results, row.names = FALSE, ...)
  } else {
    cat("\nNo results available.")
  }
  
  # Print additional information about the object structure
  cat("\n\nObject Structure:")
  cat("\n- $input_params: Input parameters used for analysis")
  cat("\n- $individual_results: Results for each successful dataset")
  cat("\n- $combined_results: Merged results from all datasets")
  cat("\n- $metadata: Analysis metadata")
  cat("\n- $errors: Errors encountered during analysis")
  
  cat("\n\nUse str(object) to view full object structure.")
  cat("\n==========================================")
  
  invisible(x)
}