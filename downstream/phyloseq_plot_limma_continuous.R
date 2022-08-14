# Plot continuous metadata variable vs taxa
phyloseq_limma_plot_indiv_continuous <- function(ps_limma_cont_object, comparison = NULL, feature_num, plot_subtitle = NULL, 
                                                 plot_x_label = NULL, plot_y_label = NULL, text_size = 8, geom_point_fill = 'black', 
                                                 geom_point_alpha = 0.7, geom_point_size = 2) {
  if (is.null(comparison) & length(ps_limma_cont_object) == 6) {
    tax_id <- rownames(ps_limma_cont_object$limma_significant)[feature_num]
    
    logFC <- ps_limma_cont_object$limma_significant$logFC[feature_num]
    direction <- ifelse(logFC > 0, 'Increased', 'Decreased')
    
    ps_df <- data.frame(t(ps_limma_cont_object$input_data[tax_id,]),
                        ps_limma_cont_object$test_variable[!is.na(ps_limma_cont_object$test_variable)],
                        direction)
    
    colnames(ps_df) <- c('taxa', 'var', 'direction')
  } else {
    tax_id <- rownames(ps_limma_cont_object$limma_significant[[comparison]])[feature_num]
    
    logFC <- ps_limma_cont_object$limma_significant[[comparison]]$logFC[feature_num]
    direction <- ifelse(logFC > 0, 'Increased', 'Decreased')
    
    ps_df <- data.frame(t(ps_limma_cont_object$input_data[tax_id,]),
                        ps_limma_cont_object$test_variables[[comparison]][!is.na(ps_limma_cont_object$test_variables[[comparison]])],
                        direction)
    
    colnames(ps_df) <- c('taxa', 'var', 'direction')
  }
  
  if (is.null(plot_x_label)) {
    plot_x_label <- metadata_var
  }
  
  if (is.null(plot_y_label)) {
    plot_y_label <- 'Abundance'
  }
  
  plot <- ggplot(ps_df, aes(x = var, y = taxa)) +
    geom_point(shape = 21, fill = 'black', alpha = 0.7, size = geom_point_size) +
    geom_smooth(aes(color = direction), method = 'gam', formula = y ~ s(x, bs = 'cs')) +
    labs(title = tax_id,
         subtitle = plot_subtitle,
         x = plot_x_label,
         y = plot_y_label) +
    theme(legend.position = 'NONE',
          text = element_text(size = text_size)) +
    scale_color_manual(values = c('Increased' = 'red', 'Decreased' = 'blue'))
  
  return(plot)
}

# Plot continuous metadata variable vs all significant taxa
phyloseq_limma_plot_all_continuous <- function(ps_limma_cont_object, comparison = NULL, plot_subtitle = NULL, 
                                               plot_x_label = NULL, plot_y_label = NULL, text_size = 8, 
                                               geom_point_fill = NULL, geom_point_alpha = NULL, geom_point_size = 2) {
  # Create an empty list to hold the plots
  plot_list <- list()
  
  # Throw error if the name of comparison/test_variable isn't provided
  if (is.null(comparison)) {
    if (length(ps_limma_cont_object) != 6) { # If the length is 6, then the object was created by the old limma_continuous function
      stop('Please provide the name of the comparison or test variable you want to plot as a string value.')
    } else {
      # Create plots for each significant taxa
      for (taxon in 1:dim(ps_limma_cont_object$limma_significant)[1]) {
        tax_id <- rownames(ps_limma_cont_object$limma_significant)[taxon]
        
        p <- phyloseq_limma_plot_indiv_continuous(ps_limma_cont_object = ps_limma_cont_object,
                                                  feature_num = taxon,
                                                  comparison = comparison,
                                                  plot_subtitle = plot_subtitle,
                                                  plot_x_label = plot_x_label,
                                                  plot_y_label = plot_y_label,
                                                  text_size = text_size,
                                                  geom_point_fill = geom_point_fill,
                                                  geom_point_alpha = geom_point_alpha,
                                                  geom_point_size = geom_point_size)
        
        plot_list[[tax_id]] <- p
      }
    }
  } else {
    # Create plots for each significant taxa
    for (taxon in 1:dim(ps_limma_cont_object$limma_significant[[comparison]])[1]) {
      tax_id <- rownames(ps_limma_cont_object$limma_significant[[comparison]])[taxon]
      
      p <- phyloseq_limma_plot_indiv_continuous(ps_limma_cont_object = ps_limma_cont_object,
                                                feature_num = taxon,
                                                comparison = comparison,
                                                plot_subtitle = plot_subtitle,
                                                plot_x_label = plot_x_label,
                                                plot_y_label = plot_y_label,
                                                text_size = text_size,
                                                geom_point_fill = geom_point_fill,
                                                geom_point_alpha = geom_point_alpha,
                                                geom_point_size = geom_point_size)
      
      plot_list[[tax_id]] <- p
    }
  }
  
  return(plot_list)
}









