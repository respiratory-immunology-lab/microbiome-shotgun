phyloseq_limma <- function(phyloseq_object, metadata_vars = NULL, metadata_condition = NULL, model_matrix = NULL,
                           model_formula_as_string = NULL, use_contrast_matrix = TRUE, coefficients = NULL,
                           factor_reorder_list = NULL, continuous_modifier_list = NULL,
                           contrast_matrix = NULL, adjust_method = 'BH', rownames = NULL, 
                           tax_id_col = 'taxa', adj_pval_threshold = 0.05, logFC_threshold = 1, 
                           legend_metadata_string = NULL, volc_plot_title = NULL, volc_plot_subtitle = NULL,
                           volc_plot_xlab = NULL, volc_plot_ylab = NULL) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'phyloseq', 'limma', 'ggplot2', 'stringr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Set 'use_contrast_matrix' to FALSE if coefficients are provided
  if (!is.null(coefficients)) {
    use_contrast_matrix <- FALSE
  }
  
  # Filter phyloseq object by the conditional statement if present
  if(!is.null(metadata_condition)) {
    phyloseq_object <- prune_samples(metadata_condition, phyloseq_object)
  }
  
  # Prepare limma data.frame
  ps_limma_data <- data.frame(otu_table(phyloseq_object))
  
  # Set rownames to tax_id_col unless provided another vector
  if (is.null(rownames)) {
    rownames(ps_limma_data) <- make.unique(tax_table(phyloseq_object)[,tax_id_col])
  } else {
    rownames(ps_limma_data) <- rownames
  }
  
  # Ensure that NA values are actually NA, and not character type 'NA'
  if (!is.null(metadata_vars)) {
    if (!is.null(metadata_condition)) {
      test_df <- data.frame(sample_data(phyloseq_object)[, metadata_vars])[metadata_condition,]
    } else {
      test_df <- data.frame(sample_data(phyloseq_object)[, metadata_vars])
    }
  } else {
    if (!is.null(metadata_condition)) {
      test_df <- data.frame(sample_data(phyloseq_object))[metadata_condition,]
    } else {
      test_df <- data.frame(sample_data(phyloseq_object))
    }
  }
  
  # Define a function to ensure that character 'NA' values are switched with real 'NA' values
  ensure_NA <- function(vector) {
    replace(vector, vector == 'NA', NA)
  }
  
  # Define a function to determine whether any 'NA' values exist on a given row
  any_NA <- function(data) {
    keep_vector <- c()
    for (i in 1:nrow(data)) {
      if (sum(is.na(data[i, ])) > 0) {
        keep_vector <- c(keep_vector, TRUE)
      } else {
        keep_vector <- c(keep_vector, FALSE)
      }
    }
    return(keep_vector)
  }
  
  # Ensure that character 'NA' values are switched with real 'NA' values
  for (i in 1:ncol(test_df)) {
    if (class(test_df[, i]) == 'character' | class(test_df[, i]) == 'factor') {
      test_df[, i] <- ensure_NA(test_df[, i])
    }
  }
  
  # Attempt to determine which columns are being referenced by the formula
  if (!is.null(model_formula_as_string)) {
    split_formula <- unlist(str_split(model_formula_as_string, pattern = '[+\\d~ ]* '))
    split_formula <- split_formula[!split_formula == ''] # Remove any empty character vector remnants
    
    if (length(split_formula) > 1) {
      test_df <- test_df[, split_formula] # Subset the identified columns from the test data.frame
    } else if (length(split_formula == 1)) {
      test_df <- data.frame(test_df[, split_formula])
      colnames(test_df) <- split_formula
    } else {
      stop('Incorrect dimensions for the data.frame used to generate the model matrix.')
    }
  }
  
  # Reorder factors where requested
  if (!is.null(factor_reorder_list)) {
    for (i in 1:length(factor_reorder_list)) {
      test_df[, names(factor_reorder_list)[i]] <- factor(test_df[, names(factor_reorder_list)[i]],
                                                         levels = factor_reorder_list[[i]])
    }
  }
  
  # Modify continuous variables where requested
  if (!is.null(continuous_modifier_list)) {
    for (i in 1:length(continuous_modifier_list)) {
      test_df[, names(continuous_modifier_list)[i]] <- sapply(test_df[, names(continuous_modifier_list)[i]],
                                                              FUN = continuous_modifier_list[[i]])
    }
  }
  
  # Filter limma data to remove NA values
  ps_limma_data <- ps_limma_data[,!any_NA(test_df)]
  
  # Create design matrix if none provided
  if (is.null(model_matrix)) {
    if (is.null(model_formula_as_string)) {
      if (!is.null(metadata_vars)) {
        if (length(metadata_vars == 1)) {
          ps_limma_design <- model.matrix(~ 0 + test_df[,1])
          colnames(ps_limma_design) <- levels(test_df[,1])
        } else {
          stop('Please provide a single character value that matches one of the variables in the phyloseq sample data object.')
        }
      } else {
        stop('Please provide either a model matrix or a formula (as a string - using a combination of the metadata variables you selected) for use with limma.
             Alternatively, please provide a single character value to the metadata_vars paramter that matches one of the variables in the phyloseq sample data object.')
      }
    } else {
      ps_limma_design <- model.matrix(as.formula(model_formula_as_string), data = test_df)
    }
  } else {
    ps_limma_design <- model_matrix
  }
  
  # Fit the expression matrix to a linear model
  fit <- lmFit(ps_limma_data, ps_limma_design)
  
  # Define a function to handle creation of the contrasts matrix
  make_contrasts_vector <- function(design_matrix) {
    # Define levels and make new lists
    levels <- colnames(design_matrix)
    num_levels <- length(levels)
    ci <- list()
    cj <- list()
    
    # For loops for generate contrast statements
    for (i in 2:num_levels) {
      ci[length(ci) + 1] <- paste0(levels[1], '-', levels[i])
      for (j in (i + 1):num_levels) {
        cj[length(cj) + 1] <- paste0(levels[i], '-', levels[j])
      }
    }
    
    # Unlist elements and remove NA values, then remove the last element (contrasts itself)
    cx <- c(unlist(ci), unlist(cj)); cx <- cx[!str_detect(cx, 'NA')]
    cx <- cx[1:length(cx) - 1]
    
    # Generate contrasts matrix and return it
    cont_matrix <- makeContrasts(contrasts = cx, levels = levels)
    cont_matrix
  }
  
  # Create the contrasts matrix if none provided
  if (is.null(contrast_matrix)) {
    cont_matrix <- make_contrasts_vector(ps_limma_design)
  } else {
    cont_matrix <- contrast_matrix
  }
  
  # Bayes statistics of differential expression
  if (use_contrast_matrix) {
    fit2 <- contrasts.fit(fit, contrasts = cont_matrix)
  } else {
    if (is.null(coefficients)) {
      stop('As you have selected to use coefficients instead of a contrast matrix, please assign indices to the 
           coefficients parameter to select coefficients you want to retain for Bayes statistics.')
    } else {
      fit2 <- contrasts.fit(fit, coefficients = coefficients)
    }
  }
  fit2 <- eBayes(fit2, robust = TRUE, trend = FALSE)
  
  # Generate a list of the differentially abundant metabolites
  get_all_topTables <- function(fit2, top = TRUE) {
    # Prepare an empty list object
    all_topTables <- list()
    
    # Run through each of the coefficients and add topTable to the list
    for (coef in 1:dim(fit2$contrasts)[2]) {
      coef_name <- colnames(fit2$contrasts)[coef]
      metab_topTable <- topTable(fit2, number = dim(fit2)[1], adjust.method = adjust_method, coef = coef)
      if (top == TRUE) {
        metab_topTable <- metab_topTable %>%
          filter(adj.P.Val < adj_pval_threshold) %>%
          filter(abs(logFC) >= logFC_threshold)
      }
      all_topTables[[coef_name]] <- metab_topTable
    }
    
    # Return list object
    all_topTables
  }
  
  # Get both significant results and all results lists
  ps_limma_signif <- get_all_topTables(fit2, top = TRUE)
  ps_limma_all <- get_all_topTables(fit2, top = FALSE)
  
  # Generate a venn diagram of the results
  results <- decideTests(fit2, adjust.method = adjust_method, p.value = adj_pval_threshold, lfc = logFC_threshold)
  if (dim(results)[2] <= 4) {
    venn_diagram <- vennDiagram(results)
  } else {
    venn_diagram <- 'No Venn diagram. Greater than 4 comparisons.'
  }
  
  # Add a direction column to the 'metab_limma_all' topTables
  limma_add_categorical_direction <- function(ps_limma_all) {
    for (i in 1:length(ps_limma_all)) {
      ps_limma_all[[i]] <- data.frame(ps_limma_all[[i]]) %>%
        mutate(direction = case_when(
          adj.P.Val < adj_pval_threshold & logFC <= -logFC_threshold ~ 'Decreased',
          adj.P.Val < adj_pval_threshold & logFC >= logFC_threshold ~ 'Increased',
          adj.P.Val >= adj_pval_threshold | abs(logFC) < logFC_threshold ~ 'NS'
        ))
    }
    ps_limma_all
  }
  
  ps_limma_all <- limma_add_categorical_direction(ps_limma_all)
  
  # Assign test_string if provided
  if (!is.null(legend_metadata_string)) {
    test_string <- legend_metadata_string
  } else {
    test_string <- 'Test Variable'
  }
  
  # Create a vector of strings for color legend
  plot_color_names <- function(test_string, ps_limma_all) {
    vector = ''
    if (use_contrast_matrix) {
      for (i in 1:length(ps_limma_all)) {
        contrast_name <- names(ps_limma_all[i])
        contrast_name <- gsub('-', ' vs. ', contrast_name)
        contrast_name <- paste0(test_string, ':\n', contrast_name)
        vector <- c(vector, contrast_name)
      }
    } else {
      for (i in 1:length(ps_limma_all)) {
        contrast_name <- names(ps_limma_all[i])
        contrast_name <- paste0(test_string, ':\n', contrast_name)
        vector <- c(vector, contrast_name)
      }
    }
    vector <- vector[2:length(vector)]
    vector
  }
  
  plot_color_labels <- plot_color_names(test_string, ps_limma_all)
  
  # Define function to plot multiple volcano plots
  limma_volcano_plots <- function(topTable, plot_title, plot_subtitle, plot_xlab, plot_ylab, plot_color_labels) {
    # Create a blank list to hold plots
    plots <- list()
    
    # Make volcano plot for each topTable
    for (i in 1:length(ps_limma_all)) {
      plot <- ggplot(ps_limma_all[[i]], aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(aes(color = direction)) +
        geom_hline(yintercept = -log10(adj_pval_threshold), linetype = 2) +
        geom_vline(xintercept = -1, linetype = 2) +
        geom_vline(xintercept = 1, linetype = 2) +
        scale_color_manual(values = c('Decreased' = 'blue', 'Increased' = 'red', 'NS' = 'grey70'), name = plot_color_labels[i]) +
        geom_text_repel(data = ps_limma_all[[i]][ps_limma_all[[i]]$adj.P.Val < adj_pval_threshold, ],
                        label = rownames(ps_limma_all[[i]][ps_limma_all[[i]]$adj.P.Val < adj_pval_threshold, ]),
                        size = 2) +
        labs(title = plot_title,
             subtitle = plot_subtitle,
             x = plot_xlab,
             y = plot_ylab)
      
      # Add plot to list
      var_name <- names(ps_limma_all[i])
      plots[[var_name]] <- plot
    }
    
    # Return plots
    plots
  }
  
  # Prepare labelling
  if (is.null(volc_plot_title)) {
    plot_title <- 'Differentially Abundant Taxa'
  } else {
    plot_title <- volc_plot_title
  }
  
  if (is.null(volc_plot_subtitle)) {
    plot_subtitle <- NULL
  } else {
    plot_subtitle <- volc_plot_subtitle
  }
  
  if (is.null(volc_plot_xlab)) {
    plot_xlab <- 'log2FC'
  } else {
    plot_xlab <- volc_plot_xlab
  }
  
  if (is.null(volc_plot_ylab)) {
    plot_ylab <- 'Significance\n-log10(adjusted p-value)'
  } else {
    plot_ylab <- volc_plot_ylab
  }
  
  # Plot all volcano plots
  volcano_plots <- limma_volcano_plots(metab_limma_all, plot_title, plot_subtitle, plot_xlab, plot_ylab, plot_color_labels)
  
  # Make a list of all components above to return from function
  if (use_contrast_matrix) {
    return_list <- list(input_data = ps_limma_data,
                        test_variables = test_df,
                        model_matrix = ps_limma_design,
                        contrast_matrix = cont_matrix,
                        limma_significant = ps_limma_signif,
                        limma_all = ps_limma_all,
                        volcano_plots = volcano_plots,
                        venn_diagram = venn_diagram)
  } else {
    return_list <- list(input_data = ps_limma_data,
                        test_variables = test_df,
                        model_matrix = ps_limma_design,
                        coefficients = coefficients,
                        limma_significant = ps_limma_signif,
                        limma_all = ps_limma_all,
                        volcano_plots = volcano_plots,
                        venn_diagram = venn_diagram)
  }
  
  
  # Return the list
  return(return_list)
}