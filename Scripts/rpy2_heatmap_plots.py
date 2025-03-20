import streamlit as st
from rpy2.rinterface_lib import openrlib
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects import conversion, default_converter
import os
os.environ['R_HOME'] = '/usr/lib/R'
os.environ['R_USER'] = '/usr/lib/R'


@st.cache_resource(show_spinner="R Libraries are Loading...")
def r_libraries():
    """
    Load necessary R libraries for plotting and analysis.
    """
    with openrlib.rlock:

        with conversion.localconverter(default_converter):
            rscript = """
            library(ComplexHeatmap)
            library(survival)
            library(survminer)
            library(circlize)
            library(devtools)
            library(grDevices)
            library(base)
            library(utils)
            library(ggplot2)
            library(ggbeeswarm)
            library(ggpubr)
            library(surviplot)
            library(survivalAnalysis)
            library(data.table)
            library(dplyr)
            library(clusterProfiler)
            library(org.Hs.eg.db)
            library(enrichplot)
            library(DOSE)
            library(compiler)
            library(extrafont)
            loadfonts(device = "pdf")
            library(Cairo)
            library(biomaRt)
            """

            # Execute the R script
            robjects.r(rscript)


def surv_plot(surv_dataframe_path, surv_png, surv_pdf, user_id, dataset_option, surv_color_path):
    """
    Generate a survival plot using R libraries.

    Arguments:
    surv_dataframe_path -- Path to the survival data CSV file.
    surv_png -- Path to save the PNG output of the survival plot.
    surv_pdf -- Path to save the PDF output of the survival plot.
    user_id -- User ID for file naming and tracking.
    dataset_option -- Dataset option to determine the time unit (e.g., 'TCGA').
    surv_color_path -- Path to the CSV file containing color definitions for clusters.
    """
    with openrlib.rlock:

        with conversion.localconverter(default_converter):
            pandas2ri.activate()

            if dataset_option == 'TCGA':
                surv_time_type = "Months"
            else:
                surv_time_type = "Months"

            robjects.r(f'surv_color <- "{surv_color_path}"')
            robjects.r(f'srv_dataframe <- "{surv_dataframe_path}"')
            robjects.r(f'srv_png <- "{surv_png}"')
            robjects.r(f'srv_pdf <- "{surv_pdf}"')
            robjects.r(f'srv_time_type <- "{surv_time_type}"')
            robjects.r(f'user_id <- "{user_id}"')

            # Define the R script
            rscript = """
            ######################
            # Survival functions #
            ######################
             read_cluster_colors = function(infile = NULL) {
                  print(paste("Reading cluster colors from:", infile))
                  if(endsWith(infile, ".tsv")) {
                    indf = read.table(infile, header = T, sep = "\t", comment.char = "")
                  } else if(endsWith(infile, ".csv")) {
                    indf = read.table(infile, header = T, sep = ",", comment.char = "")
                  } else {
                    print("Trying as tsv")
                    indf = read.table(infile, header = T, sep = "\t", comment.char = "")
                  }
                  
                  colnames(indf)[1:2] = c("cluster", "color")
                  indf = indf[complete.cases(indf[,1:2]),]
                  indf$cluster = as.character(indf$cluster)
                  indf = indf[order(indf$cluster, decreasing = F),]
                  return(indf[,1:2])
                }
                
                read_survival_table <- function(infile = NULL) {
                  print(paste("Reading survival table from:", infile))
                  indf <- read.table(infile, header = T, 
                                     sep = ifelse(grepl(".csv$", infile), ",", "\t"),
                                     comment.char = "")
                  
                  colnames(indf)[1:4] <- c("patient", "time", "status", "cluster")
                  indf <- indf[complete.cases(indf[,1:4]),]
                  indf$cluster <- as.factor(indf$cluster)  # Ensure cluster is a factor
                  print("Survival table read successfully.")
                  return(indf)
                }
                
                PlotSurvival = function(survival_df_file = NULL,
                                        cluster_color_df_file = NULL,
                                        sur_x_time = "Months",
                                        output_png = NULL,
                                        output_pdf = NULL) {
                  print("Starting PlotSurvival function")
                  surv_df = read_survival_table(infile = survival_df_file)
                  cluster_color_df = data.frame()
                  if(!is.null(cluster_color_df_file)) {
                    cluster_color_df = read_cluster_colors(infile = cluster_color_df_file)
                  }
                  
                  print("Preparing survival fit and plot")
                  if(nrow(cluster_color_df) > 0) {
                    surv_df$cluster = factor(surv_df$cluster, levels = cluster_color_df$cluster)
                    fit <- survfit(Surv(time, status) ~ cluster, data=surv_df)
                
                    ggsurv_plot = ggsurvplot(fit, data = surv_df, risk.table = TRUE, conf.int = FALSE,
                                             legend.title = "Cluster",
                                             legend.labs = levels(surv_df$cluster), # Use cluster levels as legend labels
                                             ylim = c(0, 1),
                                             palette = cluster_color_df$color,  # Change to palette
                                             xlab = paste0("Time in ", sur_x_time),
                                             ylab = "Survival Probability",
                                             ggtheme = theme_bw(base_family = "Arial", base_size = 15),
                                             main = "Kaplan-Meier Plot by Cluster", dpi = 500)
                  } else {
                    fit <- survfit(Surv(time, status) ~ cluster, data=surv_df)
                
                    ggsurv_plot = ggsurvplot(fit, data = surv_df, risk.table = TRUE, conf.int = FALSE, 
                                             legend.title = "Cluster",
                                             legend.labs = levels(surv_df$cluster), # Use cluster levels as legend labels
                                             ylim = c(0, 1),
                                             xlab = paste0("Time in ", sur_x_time),
                                             ylab = "Survival Probability",
                                             ggtheme = theme_bw(base_family = "Arial", base_size = 15),
                                             main = "Kaplan-Meier Plot by Cluster", dpi = 500)
                  }
                
                  plot_text = ""
                
                  print("Calculating hazard ratios and p-values")
                  if(length(unique(as.character(surv_df$cluster))) == 2) {
                    # Define the method based on the number of rows in surv_df
                    method_type <- if (nrow(surv_df) < 25) "exact" else "efron"
                    coxres = summary(coxph(Surv(time, status) ~ cluster, data=surv_df, method = method_type))
                    surv_p_value = coxres$coefficients[[5]]
                
                    if(surv_p_value < 0.01) {
                      surv_p_value = formatC(surv_p_value, format = "e", digits = 2)
                    } else {
                      surv_p_value = format(round(surv_p_value, 2), nsmall = 2)
                    }
                
                    hr = coxres$conf.int[[1]]
                    hr_lower = coxres$conf.int[[3]]
                    hr_upper = coxres$conf.int[[4]]
                    
                    # Ensure HR is always >= 1
                    if (hr < 1) {
                        hr = 1 / hr
                        temp = hr_lower
                        hr_lower = 1 / hr_upper
                        hr_upper = 1 / temp
                    }
                    
                    hr = format(hr, scientific = FALSE, digits = 2)
                    hr_lower = format(round(hr_lower, 2), nsmall = 2)
                    hr_upper = format(round(hr_upper, 2), nsmall = 2)
                    
                    plot_text = paste0("HR = ", hr, " (", hr_lower, "-", hr_upper, ")\n",
                                       "logrank p = ", surv_p_value)
                  } else if (length(unique(as.character(surv_df$cluster))) > 2) {
                    surv_p_value = surv_pvalue(surv_fit(Surv(time, status) ~ cluster, data=surv_df))$pval
                
                    if(surv_p_value < 0.01) {
                      surv_p_value = formatC(surv_p_value, format = "e", digits = 2)
                    } else {
                      surv_p_value = format(round(surv_p_value, 2), nsmall = 2)
                    }
                
                    plot_text = paste0("logrank p = ", surv_p_value)
                  }
                  
                  print("Customizing the plot")
                  # Ensure plot_text is only added if plot_text is not empty
                  if(plot_text != "") {
                    # Further customization for bigger axis texts, axis numbers, legend, and p-value annotation
                    ggsurv_plot$plot <- ggsurv_plot$plot + 
                      theme(
                        axis.title = element_text(size = 20, face = "bold", color = "black"),  # Axis titles in black
                        axis.text = element_text(size = 18, color = "black"),  # Axis text in black
                        legend.title = element_text(size = 20, face = "bold"),  # Legend title
                        legend.text = element_text(size = 18),  # Legend text
                        plot.title = element_text(size = 24, face = "bold"),  # Plot title
                        text = element_text(family = "Arial"),  # Arial everywhere
                        panel.grid.major = element_blank(),  # Remove major grid lines
                        panel.grid.minor = element_blank(),  # Remove minor grid lines
                        panel.background = element_blank(),  # Remove panel background
                        panel.border = element_blank(),  # Remove panel border
                        plot.background = element_blank(),  # Remove plot background
                        axis.line = element_line(color = "black", linewidth = 1),  # Add black border to axes
                        plot.margin = margin(t = 10, r = 10, b = 10, l = 10)  # Add margins around the plot
                      ) 
                    
                    ggsurv_plot$table <- ggsurv_plot$table + 
                      theme(
                        axis.title = element_text(size = 20, face = "bold", color = "black"),  # Axis titles for risk table in black
                        axis.text = element_text(size = 15, color = "black"),  # Axis text for risk table in black
                        text = element_text(family = "Arial"),  # Arial everywhere
                        panel.grid.major = element_blank(),  # Remove major grid lines
                        panel.grid.minor = element_blank(),  # Remove minor grid lines
                        panel.background = element_blank(),  # Remove panel background
                        panel.border = element_blank(),  # Remove panel border
                        plot.background = element_blank(),  # Remove plot background
                        axis.line = element_line(color = "black", linewidth = 1),  # Add black border to axes
                        plot.margin = margin(t = 10, r = 10, b = 10, l = 10)  # Add margins around the plot
                      ) 
                    
                    ggsurv_plot$pval <- ggsurv_plot$pval + 
                      theme(
                        text = element_text(size = 18, face = "bold", family = "Arial")  # P-value annotation
                      )
                    
                    # Adding the annotated text with Arial bold and larger size
                    ggsurv_plot$plot <- ggsurv_plot$plot +
                      ggplot2::annotate(
                        "text",
                        x = Inf, y = Inf,
                        vjust = 1, hjust = 1,
                        label = plot_text,
                        size = 6,  # Size for annotate (use 6 for large text in ggplot2)
                        fontface = "bold",  # Bold text
                        family = "Arial"  # Arial font
                      )
                  }
                
                  print("Saving the plot")
                  if(!is.null(output_png)) {
                    png(output_png, width = 1800, height = 1800,
                        units = "px", pointsize = 6,
                        res = 150)
                    print(ggsurv_plot)
                    dev.off()  
                  }
                
                  if (!is.null(output_pdf)) {
                    CairoPDF(output_pdf, width = 14, height = 14, family = "Arial")
                    print(ggsurv_plot)
                    dev.off()
                  }
                  print("Plot saved successfully")
                }
                
                #################
                # Analysis part #
                #################
                
                print("Starting analysis")
                
                if(!is.null(srv_dataframe)) {
                  if(is.null(srv_png) & is.null(srv_pdf)) {
                    stop("[ERROR] no output survival pdf or png specified")
                  }
                  
                  PlotSurvival(survival_df_file = srv_dataframe, 
                               cluster_color_df_file = surv_color, 
                               sur_x_time = srv_time_type, 
                               output_png = srv_png, output_pdf = srv_pdf)
                }
                print("Analysis completed")
                

            """

            # Execute the R script
            robjects.r(rscript)


def multi_var(user_id, value_list, multi_info_name, forest_pdf, forest_png, multi_db_name):
    """
    Perform multivariate analysis and generate a forest plot using R libraries.

    Arguments:
    user_id -- User ID for file naming and tracking.
    value_list -- DataFrame containing values for multivariate analysis.
    multi_info_name -- Path to save the multivariate analysis information.
    forest_pdf -- Path to save the PDF output of the forest plot.
    forest_png -- Path to save the PNG output of the forest plot.
    multi_db_name -- Path to the CSV file containing clinical data for analysis.
    """
    value_list = value_list.reset_index()

    with openrlib.rlock:

        with conversion.localconverter(default_converter + pandas2ri.converter):
            pandas2ri.activate()

            # r_dataframes contains references of annotations
            r_dataframe = robjects.conversion.get_conversion().py2rpy(value_list)
            robjects.r.assign("r_dataframe", r_dataframe)

            # Get user_id
            robjects.r(f'user_id <- "{user_id}"')
            robjects.r(f'multi_db_name <- "{multi_db_name}"')
            robjects.r(f'multi_info_name <- "{multi_info_name}"')
            robjects.r(f'forest_pdf <- "{forest_pdf}"')
            robjects.r(f'forest_png <- "{forest_png}"')

            # define the R script you want to run
            rscript = """
            clinical_file <- multi_db_name
            reference_table_file <- r_dataframe
            pdf_out <- forest_pdf
            png_out <- forest_png
            tsv_out <- multi_info_name
            
            clinical.df <- read.csv(clinical_file, header = TRUE, check.names = FALSE)
            rownames(clinical.df) <- clinical.df$`Patient ID`
            clinical.df$`Patient ID` <- NULL
            clinical.df$Status <- as.numeric(clinical.df$Status)
            clinical.df$Time <- as.numeric(clinical.df$Time)
            
            clinical.df <- clinical.df[!is.na(clinical.df$Time) & !is.na(clinical.df$Status), ]
            
            reference_table <- r_dataframe
            colnames(reference_table) <- c("Column", "Selected_reference")
            reference_table <- reference_table[!reference_table$Column %in% c("Time", "Status"), ]
            rownames(reference_table) <- reference_table$Column
            reference_table <- reference_table[rownames(reference_table) %in% colnames(clinical.df), ]
            covariate_columns <- reference_table$Column
            
            for (i in 1:nrow(reference_table)) {
              reference_table$Column[i]
              if (reference_table$Column[i] == "age" | reference_table$Column[i] == "Age") {
                clinical.df[, reference_table$Column[i]] <- as.numeric(clinical.df[, reference_table$Column[i]])
              } else {
                clinical.df[, reference_table$Column[i]] <- ifelse(clinical.df[, reference_table$Column[i]] == "nan", 
                NA, clinical.df[, reference_table$Column[i]])
                clinical.df[, reference_table$Column[i]] <- ifelse(is.nan(clinical.df[, reference_table$Column[i]]), 
                NA, clinical.df[, reference_table$Column[i]])
                clinical.df[, reference_table$Column[i]] <- factor(clinical.df[, reference_table$Column[i]])
                clinical.df[, reference_table$Column[i]] <- relevel(clinical.df[, reference_table$Column[i]], 
                                                                    ref = reference_table$Selected_reference[i])
              }
            }
            
            tmp_form <- paste0("Surv(", "Time", ",", "Status", ") ~ ")
            multivar_formula <- as.formula(paste(tmp_form, paste(covariate_columns, collapse = "+")))
            clinical.df <- as.data.table(clinical.df)
            model_multivar <- coxph(multivar_formula, data = clinical.df, na.action = na.exclude)
            
            # Customize the ggforest plot
            forest_plot <- ggforest(model_multivar, data = clinical.df, fontsize = 1) +
              theme(
                text = element_text(size = 14, family = "Arial"),
                plot.title = element_text(size = 18, family = "Arial", face = "bold"),
                axis.title = element_text(size = 16, family = "Arial", face = "bold"),
                axis.text = element_text(size = 14, family = "Arial")
              )
            
            clinical.df %>%
              analyse_multivariate(vars(Time, Status), covariates = covariate_columns,
                                   covariate_name_dict = covariate_columns,
                                   covariate_label_dict = covariate_columns) -> result
            
            # Save the plot as PNG
            forest_plot_path <- png_out
            png(forest_plot_path, width = 720, height = 720)
            print(forest_plot)
            dev.off()
            
            # Save the plot as PDF
            forest_plot_path_pdf <- pdf_out
            pdf(forest_plot_path_pdf, width = 12, height = 12)
            print(forest_plot)
            dev.off()
            dev.off()
            
            forest.df <- result$summaryAsFrame
            forest.df$factor.name <- factor(forest.df$factor.name, levels = covariate_columns)
            forest.df <- forest.df[order(forest.df$factor.name, forest.df$HR), ]
            write.table(forest.df, tsv_out, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
            """

            # Run the script
            robjects.r(rscript)


def py_hm(ht_max_perc, ht_mid_perc, ht_min_perc, ht_color_max, ht_color_mid, ht_color_min, ht_name,
          cluster_color_path, ht_top_annotation_color, ht_expression_path, column_cluster_path,
          ht_top_annotation_path, surv_dataframe_path, ht_png, ht_pdf, surv_png, surv_pdf, user_id,
          ht_row_cluster_path, dataset_option, dendrogram_display, dendrogram_display_gene):
    """
    Generate a heatmap and survival plot using R libraries.

    Arguments:
    ht_max_perc -- Maximum percentile for heatmap color scaling.
    ht_mid_perc -- Mid percentile for heatmap color scaling.
    ht_min_perc -- Minimum percentile for heatmap color scaling.
    ht_color_max -- Color for the maximum value in the heatmap.
    ht_color_mid -- Color for the mid value in the heatmap.
    ht_color_min -- Color for the minimum value in the heatmap.
    ht_name -- Name for the heatmap matrix.
    cluster_color_path -- Path to the CSV file containing cluster color definitions.
    ht_top_annotation_color -- Path to the CSV file for top annotation colors.
    ht_expression_path -- Path to the CSV file containing expression data.
    column_cluster_path -- Path to the CSV file for column clustering information.
    ht_top_annotation_path -- Path to the CSV file for top annotation data.
    surv_dataframe_path -- Path to the CSV file containing survival data.
    ht_png -- Path to save the PNG output of the heatmap.
    ht_pdf -- Path to save the PDF output of the heatmap.
    surv_png -- Path to save the PNG output of the survival plot.
    surv_pdf -- Path to save the PDF output of the survival plot.
    user_id -- User ID for file naming and tracking.
    ht_row_cluster_path -- Path to the CSV file for row clustering information.
    dataset_option -- Dataset option to determine the time unit (e.g., 'TCGA').
    """
    with openrlib.rlock:

        with conversion.localconverter(default_converter):
            pandas2ri.activate()

            ht_max_row_name = '100'
            ht_max_column_name = '75'

            if dataset_option == 'TCGA':
                surv_time_type = "Months"
            else:
                surv_time_type = "Months"

            robjects.r(f'cluster_colors <- "{cluster_color_path}"')
            robjects.r(f'dendrogram_display <- "{dendrogram_display}"')
            robjects.r(f'dendrogram_display_gene <- "{dendrogram_display_gene}"')
            robjects.r(f'ht_expression <- "{ht_expression_path}"')
            robjects.r(f'ht_column_clusters <- "{column_cluster_path}"')
            robjects.r(f'ht_png <- "{ht_png}"')
            robjects.r(f'ht_matrix_name <- "{ht_name}"')
            robjects.r(f'ht_top_annotation <- "{ht_top_annotation_path}"')
            robjects.r(f'ht_row_clusters <- "{ht_row_cluster_path}"')
            robjects.r(f'ht_top_annotation_color <- "{ht_top_annotation_color}"')
            robjects.r(f'ht_color_max_percentile <- "{ht_max_perc}"')
            robjects.r(f'ht_color_mid_percentile <- "{ht_mid_perc}"')
            robjects.r(f'ht_color_min_percentile <- "{ht_min_perc}"')
            robjects.r(f'ht_color_max <- "{ht_color_max}"')
            robjects.r(f'ht_color_mid <- "{ht_color_mid}"')
            robjects.r(f'ht_color_min <- "{ht_color_min}"')
            robjects.r(f'ht_max_row_name <- "{ht_max_row_name}"')
            robjects.r(f'ht_max_column_name <- "{ht_max_column_name}"')
            robjects.r(f'ht_pdf <- "{ht_pdf}"')
            robjects.r(f'srv_dataframe <- "{surv_dataframe_path}"')
            robjects.r(f'srv_png <- "{surv_png}"')
            robjects.r(f'srv_pdf <- "{surv_pdf}"')
            robjects.r(f'srv_time_type <- "{surv_time_type}"')

            # Define the R script
            rscript = """
    
            #####################
            # Heatmap functions #
            #####################
            
            compiler::enableJIT(3)
    
            read_expression_table = function(infile = NULL) {
              if(endsWith(infile, ".tsv")) {
                indf = read.table(infile, header = T, sep = "\t", comment.char = "",
                                  row.names = 1, check.names = F)
              } else if(endsWith(infile, ".csv")) {
                indf = read.table(infile, header = T, sep = ",", comment.char = "",
                                  row.names = 1, check.names = F)
              } else {
                print("Trying as tsv")
                indf = read.table(infile, header = T, sep = "\t", comment.char = "",
                                  row.names = 1, check.names = F)
              }
    
              return(indf)
            }
    
            read_clustering_table = function(infile = NULL) {
              if(endsWith(infile, ".tsv")) {
                indf = read.table(infile, header = T, comment.char = "",
                                  sep = "\t", check.names = F)
              } else if(endsWith(infile, ".csv")) {
                indf = read.table(infile, header = T,  comment.char = "",
                                  sep = ",", check.names = F)
              } else {
                print("Trying as tsv")
                indf = read.table(infile, header = T,  comment.char = "",
                                  sep = "\t", check.names = F)
              }
    
              rownames(indf) = indf[,1]
              return(indf)
            }
    
            read_annotation_color_table = function(infile = NULL) {
              if(endsWith(infile, ".tsv")) {
                indf = read.table(infile, header = T,  comment.char = "",
                                  sep = "\t", check.names = F)
              } else if(endsWith(infile, ".csv")) {
                indf = read.table(infile, header = T,  comment.char = "",
                                  sep = ",", check.names = F)
              } else {
                print("Trying as tsv")
                indf = read.table(infile, header = T, comment.char = "",
                                  sep = "\t", check.names = F)
              }
    
              return(indf)
            }
    
    
            PlotHeatmap = function(expression_file = NULL,
                                   out_png = NULL,
                                   out_pdf = NULL,
                                   annotation_file = NULL,
                                   ht_color_file = NULL,
                                   column_cluster_file = NULL,
                                   row_cluster_file = NULL,
                                   max_row_name = 75,
                                   max_column_name = 0,
                                   column_ann_border_max_sample = 100,
                                   ht_annotation_border_color = "#444444",
                                   matrix_name = "Z-scored",
                                   missing_ann_color = "lightgrey",
                                   ht_max_percentile = 99,
                                   ht_mid_percentile = 50,
                                   ht_min_percentile = 1,
                                   ht_max_color = "red",
                                   ht_mid_color = "white",
                                   ht_min_color = "blue") {
              expression_df = read_expression_table(infile = expression_file)
              ht_colsplit = rep("Patients", ncol(expression_df))
              ht_rowsplit = rep("Genes", nrow(expression_df))
              annotation_color_df = data.frame()
              ht_top_annotation = NULL
    
              if(!is.null(column_cluster_file)) {
                column_cluster_df = read_clustering_table(infile = column_cluster_file)
                expression_df = expression_df[,colnames(expression_df) %in% rownames(column_cluster_df)]
                expression_df = expression_df[,rownames(column_cluster_df)]
    
                if(nrow(column_cluster_df) > 0) {
                  ht_colsplit = as.character(column_cluster_df[,2])
                }
              }
    
              if(!is.null(row_cluster_file)) {
                row_cluster_df = read_clustering_table(infile = row_cluster_file)
                rownames(row_cluster_df) = row_cluster_df[,1]
                expression_df = expression_df[rownames(expression_df) %in% rownames(row_cluster_df),]
                expression_df = expression_df[rownames(row_cluster_df),]
                if(nrow(row_cluster_df) > 0) {
                  ht_rowsplit = as.character(row_cluster_df[,2])
                }
              }
    
              ann_border_width = 0.5 + (100 - ncol(expression_df))/100
    
              if(ncol(expression_df) > column_ann_border_max_sample) {
                ht_annotation_border_color = NULL
              }
    
              if(!is.null(annotation_file)) {
                annotation_df = read_expression_table(infile = annotation_file)
                annotation_df = annotation_df[colnames(expression_df), , drop = FALSE]
                annotation_df = annotation_df
                if("Cluster" %in% colnames(annotation_df)) {
                  annotation_df$Cluster = as.character(annotation_df$Cluster)
                }
    
                anno_legend_param_list = list()
                break_size = 4
                for(i in 1:ncol(annotation_df)) {
                  n_breaks = ceiling(length(unique(annotation_df[,i]))/break_size)
                  anno_legend_param_list[[i]] = list(ncol = 2,
                                                     title = colnames(annotation_df)[i])#,
                                                     #labels_gp = gpar(fontsize = 18),
                                                     #title_gp = gpar(fontsize = 22, fontface = "bold"))
                }
    
                names(anno_legend_param_list) = colnames(annotation_df)
                annotation_size_par = unit(8, "mm")
                annot_name_font = 12
    
                if(!is.null(ht_color_file)) {
                  annotation_color_df = read_annotation_color_table(infile = ht_color_file)
                  color_tmp = annotation_color_df$color
                  names(color_tmp) = annotation_color_df$type
                  ht_color_list = split(color_tmp, annotation_color_df$column)
                  annotation_df = annotation_df[,colnames(annotation_df) %in% names(ht_color_list)]
    
    
                  if(!is.null(ht_annotation_border_color)) {
                    ht_top_annotation = HeatmapAnnotation(df = annotation_df,
                                                          col = ht_color_list, border = TRUE,
                                                          gp = gpar(col = ht_annotation_border_color, 
                                                          lwd = 1*ann_border_width),
                                                          annotation_name_gp = gpar(fontsize = annot_name_font, 
                                                          fontface = "bold"),
                                                          na_col = missing_ann_color,
                                                          annotation_name_side = "left",
                                                          simple_anno_size = annotation_size_par,
                                                          annotation_legend_param = anno_legend_param_list)
                  } else {
                    ht_top_annotation = HeatmapAnnotation(df = annotation_df,
                                                          col = ht_color_list, border = TRUE,
                                                          annotation_name_gp = gpar(fontsize = annot_name_font, 
                                                          fontface = "bold"),
                                                          na_col = missing_ann_color,
                                                          annotation_name_side = "left",
                                                          simple_anno_size = annotation_size_par,
                                                          annotation_legend_param = anno_legend_param_list)
                  }
                } else {
                  if(!is.null(ht_annotation_border_color)) {
                    ht_top_annotation = HeatmapAnnotation(df = annotation_df,
                                                          gp = gpar(col = ht_annotation_border_color,
                                                                    lwd = 1*ann_border_width),
                                                          annotation_name_gp = gpar(fontsize = annot_name_font, 
                                                          fontface = "bold", fontfamily = "Arial"),
                                                          na_col = missing_ann_color,
                                                          annotation_name_side = "left",
                                                          simple_anno_size = annotation_size_par,
                                                          annotation_legend_param = anno_legend_param_list)
                  } else {
                    ht_top_annotation = HeatmapAnnotation(df = annotation_df, border = TRUE,
                                                          annotation_name_gp = gpar(fontsize = annot_name_font, 
                                                          fontface = "bold", fontfamily = "Arial"),
                                                          na_col = missing_ann_color,
                                                          annotation_name_side = "left",
                                                          simple_anno_size = annotation_size_par,
                                                          annotation_legend_param = anno_legend_param_list)
                  }
                }
              }
    
    
              row_name_font_size = 20 * (100 - nrow(expression_df))/100
              
              ht_color_limits = c(quantile(as.matrix(expression_df),
                                                         as.numeric(ht_min_percentile)/100,
                                                         na.rm = T)[1],
                                  (quantile(as.matrix(expression_df),
                                                         as.numeric(ht_min_percentile)/100,
                                                         na.rm = T)[1] + quantile(as.matrix(expression_df),
                                                         as.numeric(ht_max_percentile)/100,
                                                         na.rm = T)[1])/2,
                                  quantile(as.matrix(expression_df),
                                                         as.numeric(ht_max_percentile)/100,
                                                         na.rm = T)[1]
                                  )
                                              
              if(min(expression_df, na.rm=T) < 0) {
                ht_color_limits[2] = 0
              }
              
              heatmap_color_func = colorRamp2(ht_color_limits,
                                              c(ht_min_color, ht_mid_color, ht_max_color))
              
                # Customize the heatmap
              logical_value_dendro <- ifelse(dendrogram_display == "True", TRUE, FALSE)
              logical_value_dendro_gene <- ifelse(dendrogram_display_gene == "True", TRUE, FALSE)
              heatmap_result = Heatmap(as.matrix(expression_df),
                                       top_annotation = ht_top_annotation,
                                       show_row_names = ifelse(nrow(expression_df) > as.numeric(max_row_name), F, T),
                                       #show_column_names = ifelse(ncol(expression_df) > as.numeric(max_column_name), F, T),
                                       show_column_names =F,
                                       name = matrix_name,
                                       column_split = ht_colsplit,
                                       row_title_rot = 0,
                                       use_raster = T,
                                       show_column_dend = logical_value_dendro,
                                       show_row_dend = logical_value_dendro_gene,
                                       row_split = ht_rowsplit,
                                       row_names_gp = gpar(fontsize = row_name_font_size, fontfamily = "Arial"),
                                       column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
                                       row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
                                       col = heatmap_color_func,
                                       border_gp=grid::gpar(col = "black", lwd = 1),
                                       border = TRUE,
                                       heatmap_legend_param = list(
                                           title_gp = gpar(fontsize = 12, fontfamily = "Arial", fontface = "bold"),
                                           labels_gp = gpar(fontsize = 10, fontfamily = "Arial"),
                                           legend_direction = "vertical"
                                       ))
                
              if(!is.null(out_png)) {
                #png(out_png, width = 2750, height = 2750,
                #    units = "px", pointsize = 6,
                #    res = 175)
                png(out_png, width = 2000, height = 2000, units = "px", pointsize = 32, res = 215)
                print(heatmap_result)
                dev.off()
              }
              if(!is.null(out_pdf)) {
                CairoPDF(out_pdf, width = 8, height = 8, family = "Arial")
                print(heatmap_result)
                dev.off()
              }
            }
            ######################
            # Survival functions #
            ######################
    
    
            read_cluster_colors = function(infile = NULL) {
              if(endsWith(infile, ".tsv")) {
                indf = read.table(infile, header = T, sep = "\t", comment.char = "")
              } else if(endsWith(infile, ".csv")) {
                indf = read.table(infile, header = T, sep = ",", comment.char = "")
              } else {
                print("Trying as tsv")
                indf = read.table(infile, header = T, sep = "\t", comment.char = "")
              }
    
              colnames(indf)[1:2] = c("cluster", "color")
              indf = indf[complete.cases(indf[,1:2]),]
              indf$cluster = as.character(indf$cluster)
              indf = indf[order(indf$cluster, decreasing = F),]
              return(indf[,1:2])
            }
    
            read_survival_table = function(infile = NULL) {
              if(endsWith(infile, ".tsv")) {
                indf = read.table(infile, header = T, sep = "\t", comment.char = "")
              } else if(endsWith(infile, ".csv")) {
                indf = read.table(infile, header = T, sep = ",", comment.char = "")
              } else {
                print("Trying as tsv")
                indf = read.table(infile, header = T, sep = "\t", comment.char = "")
              }
    
              colnames(indf)[1:4] = c("patient", "time", "status", "cluster")
              indf = indf[complete.cases(indf[,1:4]),]
              indf$cluster = as.character(indf$cluster)
              return(indf[,1:4])
            }
    
    
            PlotSurvival = function(survival_df_file = NULL,
                                    cluster_color_df_file = NULL,
                                    sur_x_time = "Months",
                                    output_png = NULL,
                                    output_pdf = NULL) {
              surv_df = read_survival_table(infile = survival_df_file)
              cluster_color_df = data.frame()
              if(!is.null(cluster_color_df_file)) {
                cluster_color_df = read_cluster_colors(infile = cluster_color_df_file)
              }
    
              if(nrow(cluster_color_df) > 0) {
                surv_df$cluster = factor(surv_df$cluster, levels = cluster_color_df$cluster)
                fit <- survfit(Surv(time, status) ~ cluster, data=surv_df)
                ggsurv_plot = ggsurvplot(fit, data = surv_df, risk.table = TRUE, conf.int = FALSE,
                                         legend.title = "Cluster",
                                         ylim = c(0, 1),
                                         color = "cluster",
                                         palette = cluster_color_df$color,
                                         xlab = paste0("Time in ", sur_x_time),
                                         ylab = "Survival Probability",
                                         ggtheme = theme_bw(base_family = "Arial", base_size = 15),
                                         main = "Kaplan-Meier Plot by Cluster", dpi = 500)
                                         
                                            
              } else {
                fit <- survfit(Surv(time, status) ~ cluster, data=surv_df)
    
                ggsurv_plot = ggsurvplot(fit, data = surv_df, risk.table = TRUE, conf.int = FALSE, legend.title = "Cluster",
                                         ylim = c(0, 1),
                                         xlab = paste0("Time in ", sur_x_time),
                                         ylab = "Survival Probability",
                                         ggtheme = theme_bw(base_family = "Arial", base_size = 15),
                                         main = "Kaplan-Meier Plot by Cluster", dpi = 500)
              }
    
              plot_text = ""
    
              if(length(unique(as.character(surv_df$cluster))) == 2) {
                # Define the method based on the number of rows in surv_df
                method_type <- if (nrow(surv_df) < 25) "exact" else "efron"
                coxres = summary(coxph(Surv(time, status) ~ cluster, data=surv_df, method = method_type))
                surv_p_value = coxres$coefficients[[5]]
    
                if(surv_p_value < 0.01) {
                  surv_p_value = formatC(surv_p_value, format = "e", digits = 2)
                } else {
                  surv_p_value = format(round(surv_p_value, 2), nsmall = 2)
                }
                
                                    hr = coxres$conf.int[[1]]
                    hr_lower = coxres$conf.int[[3]]
                    hr_upper = coxres$conf.int[[4]]
                    
                    # Ensure HR is always >= 1
                    if (hr < 1) {
                        hr = 1 / hr
                        temp = hr_lower
                        hr_lower = 1 / hr_upper
                        hr_upper = 1 / temp
                    }
                    
                    hr = format(hr, scientific = FALSE, digits = 2)
                    hr_lower = format(round(hr_lower, 2), nsmall = 2)
                    hr_upper = format(round(hr_upper, 2), nsmall = 2)
                    
                    plot_text = paste0("HR = ", hr, " (", hr_lower, "-", hr_upper, ")\n",
                                       "logrank p = ", surv_p_value)
    
              } else if (length(unique(as.character(surv_df$cluster))) > 2) {
                surv_p_value = surv_pvalue(surv_fit(Surv(time, status) ~ cluster, data=surv_df))$pval
    
                if(surv_p_value < 0.01) {
                  surv_p_value = formatC(surv_p_value, format = "e", digits = 2)
                } else {
                  surv_p_value = format(round(surv_p_value, 2), nsmall = 2)
                }
    
                plot_text = paste0("logrank p = ", surv_p_value)
                
              }
                # Further customization for bigger axis texts, axis numbers, legend, and p-value annotation
                ggsurv_plot$plot <- ggsurv_plot$plot + 
                  theme(
                    axis.title = element_text(size = 20, face = "bold", color = "black"),  # Axis titles in black
                    axis.text = element_text(size = 18, color = "black"),  # Axis text in black
                    legend.title = element_text(size = 20, face = "bold"),  # Legend title
                    legend.text = element_text(size = 18),  # Legend text
                    plot.title = element_text(size = 24, face = "bold"),  # Plot title
                    text = element_text(family = "Arial"),  # Arial everywhere
                    panel.grid.major = element_blank(),  # Remove major grid lines
                    panel.grid.minor = element_blank(),  # Remove minor grid lines
                    panel.background = element_blank(),  # Remove panel background
                    panel.border = element_blank(),  # Remove panel border
                    plot.background = element_blank(),  # Remove plot background
                    axis.line = element_line(color = "black", size = 1),  # Add black border to axes
                    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)  # Add margins around the plot
                  ) 
                
                ggsurv_plot$table <- ggsurv_plot$table + 
                  theme(
                    axis.title = element_text(size = 20, face = "bold", color = "black"),  # Axis titles for risk table in black
                    axis.text = element_text(size = 15, color = "black"),  # Axis text for risk table in black
                    text = element_text(family = "Arial"),  # Arial everywhere
                    panel.grid.major = element_blank(),  # Remove major grid lines
                    panel.grid.minor = element_blank(),  # Remove minor grid lines
                    panel.background = element_blank(),  # Remove panel background
                    panel.border = element_blank(),  # Remove panel border
                    plot.background = element_blank(),  # Remove plot background
                    axis.line = element_line(color = "black", size = 1),  # Add black border to axes
                    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)  # Add margins around the plot
                  ) 
                
                ggsurv_plot$pval <- ggsurv_plot$pval + 
                  theme(
                    text = element_text(size = 18, face = "bold", family = "Arial")  # P-value annotation
                  )
                
                # Adding the annotated text with Arial bold and larger size
                ggsurv_plot$plot <- ggsurv_plot$plot +
                  ggplot2::annotate(
                    "text",
                    x = Inf, y = Inf,
                    vjust = 1, hjust = 1,
                    label = plot_text,
                    size = 6,  # Size for annotate (use 6 for large text in ggplot2)
                    fontface = "bold",  # Bold text
                    family = "Arial"  # Arial font
                  )
    
              if(!is.null(output_png)) {
                png(output_png, width = 1800, height = 1800,
                    units = "px", pointsize = 6,
                    res = 150)
                print("ok_saving1")
                print(ggsurv_plot)
                print("ok_saving2")
                dev.off()  
              }
    
              if(!is.null(output_pdf)) {
                CairoPDF(output_pdf, width = 14, height = 14, family = "Arial")
                print(ggsurv_plot)
                dev.off()
              }
            }
    
            #################
            # Analysis part #
            #################
            
            if(!is.null(ht_expression)) {
              if(is.null(ht_png) & is.null(ht_pdf)) {
                stop("[ERROR] no output heatmap pdf or png specified")
              }
              PlotHeatmap(expression_file = ht_expression,
                          out_png = ht_png,
                          out_pdf = ht_pdf,
                          annotation_file = ht_top_annotation,
                          ht_color_file = ht_top_annotation_color,
                          column_cluster_file = ht_column_clusters, 
                          row_cluster_file = ht_row_clusters, 
                          max_row_name = ht_max_row_name, 
                          max_column_name = ht_max_column_name, 
                          column_ann_border_max_sample = 100, 
                          ht_annotation_border_color = "black",
                          matrix_name = ht_matrix_name,
                          missing_ann_color = "lightgrey", 
                          ht_max_percentile = ht_color_max_percentile, 
                          ht_mid_percentile = ht_color_mid_percentile, 
                          ht_min_percentile = ht_color_min_percentile, 
                          ht_max_color = ht_color_max,
                          ht_mid_color = ht_color_mid, 
                          ht_min_color = ht_color_min)
            }
            
            tryCatch({
                      if (!is.null(srv_dataframe)) {
                        if (is.null(srv_png) & is.null(srv_pdf)) {
                          stop("[ERROR] no output survival pdf or png specified")
                        }
                        
                        PlotSurvival(
                          survival_df_file = srv_dataframe, 
                          cluster_color_df_file = cluster_colors, 
                          sur_x_time = srv_time_type, 
                          output_png = srv_png, 
                          output_pdf = srv_pdf
                        )
                      }
                    }, 
                    error = function(e) {
                      cat("[ERROR] An error occurred:\n", e$message, "\n")
                    })

    
            """

            # Execute the R script
            robjects.r(rscript)


def corr_heatmap(user_id, selected_genes_dm_name, cluster_name, title, anno_colors_path, e, cluster_toggle):
    """
    Generate a correlation heatmap using R libraries.

    Arguments:
    corr_df_path -- Path to the CSV file containing correlation data.
    ht_corr -- Name for the correlation matrix.
    ht_corr_png -- Path to save the PNG output of the correlation heatmap.
    ht_corr_pdf -- Path to save the PDF output of the correlation heatmap.
    user_id -- User ID for file naming and tracking.
    """
    with openrlib.rlock:

        with conversion.localconverter(default_converter):
            pandas2ri.activate()
            robjects.r(f'user_id_r <- "{user_id}"')
            robjects.r(f'selected_genes_dm_name <- "{selected_genes_dm_name}"')
            robjects.r(f'cluster_name <- "{cluster_name}"')
            robjects.r(f'title <- "{title}"')
            robjects.r(f'anno_colors_path <- "{anno_colors_path}"')

            r_script = '''
            annotation_color_df = read.csv(anno_colors_path, header = TRUE, check.names = F, comment.char = "")
            color_tmp = annotation_color_df$color
            names(color_tmp) = annotation_color_df$type
            ht_color_list = split(color_tmp, annotation_color_df$column)
            
            mat_corr <- read.csv(selected_genes_dm_name, header = TRUE, check.names = F, comment.char = "")
            numCols <- ncol(mat_corr)
            row.names(mat_corr) <- mat_corr[, numCols]
            mat_corr <- mat_corr[, -numCols]
    
            mat_cluster <- read.csv(cluster_name, header = TRUE, check.names = F, comment.char = "", row.names = 1)
            
            cl_names = rownames(mat_cluster)
            cl_names = cl_names[cl_names %in% colnames(mat_corr)]
            
            mat_corr = mat_corr[,cl_names]
            mat_cluster = mat_cluster[cl_names,,drop=F]
            
            cl_names = cl_names[cl_names %in% rownames(mat_corr)]
            mat_corr = mat_corr[cl_names,]
            
            ht_top_annotation = HeatmapAnnotation(df = mat_cluster, col = ht_color_list, border = TRUE, 
            annotation_name_side = "left", gp = grid::gpar(col = "black", lwd = 0.5, fontsize = 20, fontface = "bold"), 
            show_annotation_name = TRUE)
            
            # Create the heatmap
            heatmap <- Heatmap(mat_corr, 
                               name = "Correlation",
                               show_row_names = F,
                               show_column_names = F,
                               column_title_side = "top",
                               row_title = title,
                               row_title_gp = gpar(fontsize = 20, fontface = "bold"),
                               border = T,
                               column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                               column_split = mat_cluster$Cluster,
                               row_split = mat_cluster$Cluster,
                               column_names_gp = gpar(fontsize = 20, fontface = "bold"),
                               border_gp = gpar(col = "black", lwd = 1),
                               column_dend_gp = gpar(fontsize = 20, fontface = "bold"),
                               column_gap = unit(3, "mm"),
                               row_gap = unit(3, "mm"),
                               show_column_dend =T,
                               show_row_dend = F,
                               top_annotation = ht_top_annotation
            )
            
            if (title == "Genes") {
              corr_heatmap_path <- paste0("result_data/result_data_" ,user_id_r, "/corr_hm_genes", user_id_r, ".png")
            } else {
              corr_heatmap_path <- paste0("result_data/result_data_" ,user_id_r, "/corr_hm_samples", user_id_r, ".png")
            }   
            
            png(corr_heatmap_path, height = 800, width = 800, pointsize = 12)
            print(heatmap)
            dev.off()
    
            if (title == "Genes") {
              corr_heatmap_path_pdf <- paste0("result_data/result_data_" ,user_id_r, "/corr_hm_genes", user_id_r, ".pdf")
            } else {
              corr_heatmap_path_pdf <- paste0("result_data/result_data_" ,user_id_r, "/corr_hm_samples", user_id_r, ".pdf")
            }   
            
            pdf(corr_heatmap_path_pdf, width=7, height=7, title="Correlation Heatmap")
            print(heatmap)
            dev.off()
    
            '''

            if not cluster_toggle:
                try:
                    robjects.r(r_script)
                    img_cols = st.columns(3)
                    img_cols[1].image("result_data/result_data_" + user_id + "/corr_hm_samples" + user_id + ".png")

                except:
                    img_cols[1].image("style_items/error.svg", use_container_width=True)
                    img_cols[1].warning(
                        ":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                        ":blue[or ] :red[refresh ] :blue[page]")
            else:
                    if title == "Genes":
                        try:
                            robjects.r(r_script)
                            e[1].image("result_data/result_data_" + user_id + "/corr_hm_genes" + user_id + ".png")
                        except:
                            e[1].image("style_items/error.svg", use_container_width=True)
                            e[1].warning(
                                ":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                                ":blue[or ] :red[refresh ] :blue[page]")
                    else:
                        try:
                            robjects.r(r_script)
                            e[0].image("result_data/result_data_" + user_id + "/corr_hm_samples" + user_id + ".png")
                        except:
                            e[1].image("style_items/error.svg", use_container_width=True)
                            e[1].warning(
                                ":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                                ":blue[or ] :red[refresh ] :blue[page]")

    # Cleanup R environment
    robjects.r('rm(list = ls())')  # Remove all R variables
    robjects.r('graphics.off()')  # Close all open graphics devices


def py_hm_ssgsea_surv(ht_name, ht_top_annotation_color, ht_expression_path, column_cluster_path,
                      ht_top_annotation_path, ht_png, ht_pdf, user_id, ht_row_cluster_path, dataset_option,
                      dendrogram_display, dendrogram_display_gene):
    """
    Generate an SS-GSEA heatmap with survival analysis using R libraries.

    Arguments:
    user_id_r -- User ID for file naming and tracking.
    """
    with openrlib.rlock:

        with conversion.localconverter(default_converter):
            pandas2ri.activate()

            ht_max_row_name = '100'
            ht_max_column_name = '75'

            if dataset_option == 'TCGA':
                surv_time_type = "Months"
            else:
                surv_time_type = "Months"

            robjects.r(f'ht_matrix_name <- "{ht_name}"')
            robjects.r(f'ht_expression <- "{ht_expression_path}"')
            robjects.r(f'dendrogram_display <- "{dendrogram_display}"')
            robjects.r(f'dendrogram_display_gene <- "{dendrogram_display_gene}"')
            robjects.r(f'ht_column_clusters <- "{column_cluster_path}"')
            robjects.r(f'ht_png <- "{ht_png}"')
            robjects.r(f'ht_matrix_name <- "{ht_name}"')
            robjects.r(f'ht_top_annotation <- "{ht_top_annotation_path}"')
            robjects.r(f'ht_row_clusters <- "{ht_row_cluster_path}"')
            robjects.r(f'ht_top_annotation_color <- "{ht_top_annotation_color}"')
            robjects.r(f'ht_pdf <- "{ht_pdf}"')
            robjects.r(f'ht_max_row_name <- "{ht_max_row_name}"')
            robjects.r(f'ht_max_column_name <- "{ht_max_column_name}"')
            robjects.r(f'surv_time_type <- "{surv_time_type}"')

            # Define the R script
            rscript = """
    
            #####################
            # Heatmap functions #
            #####################
    
            read_expression_table = function(infile = NULL) {
              if(endsWith(infile, ".tsv")) {
                indf = read.table(infile, header = T, sep = "\t", comment.char = "",
                                  row.names = 1, check.names = F)
              } else if(endsWith(infile, ".csv")) {
                indf = read.table(infile, header = T, sep = ",", comment.char = "",
                                  row.names = 1, check.names = F)
              } else {
                print("Trying as tsv")
                indf = read.table(infile, header = T, sep = "\t", comment.char = "",
                                  row.names = 1, check.names = F)
              }
    
              return(indf)
            }
    
            read_clustering_table = function(infile = NULL) {
              if(endsWith(infile, ".tsv")) {
                indf = read.table(infile, header = T, comment.char = "",
                                  sep = "\t", check.names = F)
              } else if(endsWith(infile, ".csv")) {
                indf = read.table(infile, header = T,  comment.char = "",
                                  sep = ",", check.names = F)
              } else {
                print("Trying as tsv")
                indf = read.table(infile, header = T,  comment.char = "",
                                  sep = "\t", check.names = F)
              }
    
              rownames(indf) = indf[,1]
              return(indf)
            }
    
            read_annotation_color_table = function(infile = NULL) {
              if(endsWith(infile, ".tsv")) {
                indf = read.table(infile, header = T,  comment.char = "",
                                  sep = "\t", check.names = F)
              } else if(endsWith(infile, ".csv")) {
                indf = read.table(infile, header = T,  comment.char = "",
                                  sep = ",", check.names = F)
              } else {
                print("Trying as tsv")
                indf = read.table(infile, header = T, comment.char = "",
                                  sep = "\t", check.names = F)
              }
    
              return(indf)
            }
    
    
            PlotHeatmap = function(expression_file = NULL,
                                   out_png = NULL,
                                   out_pdf = NULL,
                                   annotation_file = NULL,
                                   ht_color_file = NULL,
                                   column_cluster_file = NULL,
                                   row_cluster_file = NULL,
                                   max_row_name = 75,
                                   max_column_name = 50,
                                   column_ann_border_max_sample = 100,
                                   ht_annotation_border_color = "#444444",
                                   matrix_name = "Expression",
                                   missing_ann_color = "lightgrey",
                                   ht_max_percentile = 99,
                                   ht_mid_percentile = 50,
                                   ht_min_percentile = 1,
                                   ht_max_color = "red",
                                   ht_mid_color = "white",
                                   ht_min_color = "blue") {
              expression_df = read_expression_table(infile = expression_file)
              ht_colsplit = rep("Patients", ncol(expression_df))
              ht_rowsplit = rep("Genes", nrow(expression_df))
              annotation_color_df = data.frame()
              ht_top_annotation = NULL
    
              if(!is.null(column_cluster_file)) {
                column_cluster_df = read_clustering_table(infile = column_cluster_file)
                expression_df = expression_df[,colnames(expression_df) %in% rownames(column_cluster_df)]
                expression_df = expression_df[,rownames(column_cluster_df)]
    
                if(nrow(column_cluster_df) > 0) {
                  ht_colsplit = as.character(column_cluster_df[,2])
                }
              }
    
              if(!is.null(row_cluster_file)) {
                row_cluster_df = read_clustering_table(infile = row_cluster_file)
                rownames(row_cluster_df) = row_cluster_df[,1]
                expression_df = expression_df[rownames(expression_df) %in% rownames(row_cluster_df),]
                expression_df = expression_df[rownames(row_cluster_df),]
                if(nrow(row_cluster_df) > 0) {
                  ht_rowsplit = as.character(row_cluster_df[,2])
                }
              }
    
              ann_border_width = 0.5 + (100 - ncol(expression_df))/100
    
              if(ncol(expression_df) > column_ann_border_max_sample) {
                ht_annotation_border_color = NULL
              }
    
              if(!is.null(annotation_file)) {
                annotation_df = read_expression_table(infile = annotation_file)
                annotation_df = annotation_df[colnames(expression_df), , drop = FALSE]
                annotation_df = annotation_df
                if("Cluster" %in% colnames(annotation_df)) {
                  annotation_df$Cluster = as.character(annotation_df$Cluster)
                }
    
                anno_legend_param_list = list()
                break_size = 4
                for(i in 1:ncol(annotation_df)) {
                  n_breaks = ceiling(length(unique(annotation_df[,i]))/break_size)
                  anno_legend_param_list[[i]] = list(ncol = 2,
                                                     title = colnames(annotation_df)[i])#,
                                                     #labels_gp = gpar(fontsize = 18),
                                                     #title_gp = gpar(fontsize = 22, fontface = "bold"))
                }
    
                names(anno_legend_param_list) = colnames(annotation_df)
                annotation_size_par = unit(4, "mm")
                annot_name_font = 10
    
                if(!is.null(ht_color_file)) {
                  annotation_color_df = read_annotation_color_table(infile = ht_color_file)
                  color_tmp = annotation_color_df$color
                  names(color_tmp) = annotation_color_df$type
                  ht_color_list = split(color_tmp, annotation_color_df$column)
                  annotation_df = annotation_df[,colnames(annotation_df) %in% names(ht_color_list)]    
    
                  if(!is.null(ht_annotation_border_color)) {
                    ht_top_annotation = HeatmapAnnotation(df = annotation_df,
                                                          col = ht_color_list, border = TRUE,
                                                          gp = gpar(col = ht_annotation_border_color, 
                                                          lwd = 1*ann_border_width),
                                                          annotation_name_gp = gpar(fontsize = annot_name_font, 
                                                          fontface = "bold"),
                                                          na_col = missing_ann_color,
                                                          annotation_name_side = "left",
                                                          simple_anno_size = annotation_size_par,
                                                          annotation_legend_param = anno_legend_param_list)
                  } else {
                    ht_top_annotation = HeatmapAnnotation(df = annotation_df,
                                                          col = ht_color_list, border = TRUE,
                                                          annotation_name_gp = gpar(fontsize = annot_name_font, 
                                                          fontface = "bold"),
                                                          na_col = missing_ann_color,
                                                          annotation_name_side = "left",
                                                          simple_anno_size = annotation_size_par,
                                                          annotation_legend_param = anno_legend_param_list)
                  }
                } else {
                  if(!is.null(ht_annotation_border_color)) {
                    ht_top_annotation = HeatmapAnnotation(df = annotation_df,
                                                          gp = gpar(col = ht_annotation_border_color,
                                                                    lwd = 1*ann_border_width),
                                                          annotation_name_gp = gpar(fontsize = annot_name_font, 
                                                          fontface = "bold"),
                                                          na_col = missing_ann_color,
                                                          annotation_name_side = "left",
                                                          simple_anno_size = annotation_size_par,
                                                          annotation_legend_param = anno_legend_param_list)
                  } else {
                    ht_top_annotation = HeatmapAnnotation(df = annotation_df, border = TRUE,
                                                          annotation_name_gp = gpar(fontsize = annot_name_font, 
                                                          fontface = "bold"),
                                                          na_col = missing_ann_color,
                                                          annotation_name_side = "left",
                                                          simple_anno_size = annotation_size_par,
                                                          annotation_legend_param = anno_legend_param_list)
                  }
                }
              }
    
    
              row_name_font_size = 20 * (100 - nrow(expression_df))/100
                ht_color_limits = c(quantile(as.matrix(expression_df),
                                             as.numeric(ht_min_percentile)/100,
                                             na.rm = T)[1],
                      (quantile(as.matrix(expression_df),
                                             as.numeric(ht_min_percentile)/100,
                                             na.rm = T)[1] + quantile(as.matrix(expression_df),
                                             as.numeric(ht_max_percentile)/100,
                                             na.rm = T)[1])/2,
                      quantile(as.matrix(expression_df),
                                             as.numeric(ht_max_percentile)/100,
                                             na.rm = T)[1]
                      )
                                              
                  if(min(expression_df, na.rm=T) < 0) {
                    ht_color_limits[2] = 0
                  }
                  
              heatmap_color_func = colorRamp2(ht_color_limits,
                                              c(ht_min_color, ht_mid_color, ht_max_color))
              logical_value_dendro <- ifelse(dendrogram_display == "True", TRUE, FALSE)
              logical_value_dendro_gene <- ifelse(dendrogram_display_gene == "True", TRUE, FALSE)
              heatmap_result = Heatmap(as.matrix(expression_df),
                                       top_annotation = ht_top_annotation,
                                       show_row_names = ifelse(nrow(expression_df) > as.numeric(max_row_name), F, T),
                                       #show_column_names = ifelse(ncol(expression_df) > as.numeric(max_column_name), F, T),
                                       show_column_names =F,
                                       name = matrix_name,
                                       show_column_dend = logical_value_dendro,
                                       show_row_dend = logical_value_dendro_gene,
                                       column_split = ht_colsplit,
                                       use_raster = T,
                                       row_split = ht_rowsplit,
                                       row_names_gp = gpar(fontsize = row_name_font_size, fontfamily = "Arial"),
                                       column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
                                       row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
                                       col = heatmap_color_func,
                                       border_gp=grid::gpar(col = "black", lwd = 1),
                                       border = TRUE,
                                       heatmap_legend_param = list(
                                           title_gp = gpar(fontsize = 12, fontfamily = "Arial", fontface = "bold"),
                                           labels_gp = gpar(fontsize = 10, fontfamily = "Arial"),
                                           legend_direction = "vertical"
                                       ))
    
              if(!is.null(out_png)) {
                #png(out_png, width = 2750, height = 2750,
                #    units = "px", pointsize = 6,
                #    res = 175)
                png(out_png, width = 2000, height = 2000, units = "px", pointsize = 32, res = 215)
                print(heatmap_result)
                dev.off()
              }
              if(!is.null(out_pdf)) {
                CairoPDF(out_pdf, width = 14, height = 16, family = "Arial")
                print(heatmap_result)
                dev.off()
              }
            }
            
            ######################
            # Survival functions #
            ######################
    
    
            read_cluster_colors = function(infile = NULL) {
              if(endsWith(infile, ".tsv")) {
                indf = read.table(infile, header = T, sep = "\t", comment.char = "")
              } else if(endsWith(infile, ".csv")) {
                indf = read.table(infile, header = T, sep = ",", comment.char = "")
              } else {
                print("Trying as tsv")
                indf = read.table(infile, header = T, sep = "\t", comment.char = "")
              }
    
              colnames(indf)[1:2] = c("cluster", "color")
              indf = indf[complete.cases(indf[,1:2]),]
              indf$cluster = as.character(indf$cluster)
              indf = indf[order(indf$cluster, decreasing = F),]
              return(indf[,1:2])
            }
    
            read_survival_table = function(infile = NULL) {
              if(endsWith(infile, ".tsv")) {
                indf = read.table(infile, header = T, sep = "\t", comment.char = "")
              } else if(endsWith(infile, ".csv")) {
                indf = read.table(infile, header = T, sep = ",", comment.char = "")
              } else {
                print("Trying as tsv")
                indf = read.table(infile, header = T, sep = "\t", comment.char = "")
              }
    
              colnames(indf)[1:4] = c("patient", "time", "status", "cluster")
              indf = indf[complete.cases(indf[,1:4]),]
              indf$cluster = as.character(indf$cluster)
              return(indf[,1:4])
            }
    
    
            PlotSurvival = function(survival_df_file = NULL,
                                    cluster_color_df_file = NULL,
                                    sur_x_time = "Months",
                                    output_png = NULL,
                                    output_pdf = NULL) {
              surv_df = read_survival_table(infile = survival_df_file)
              cluster_color_df = data.frame()
              if(!is.null(cluster_color_df_file)) {
                cluster_color_df = read_cluster_colors(infile = cluster_color_df_file)

              }
                
              if(nrow(cluster_color_df) > 0) {
                surv_df$cluster = factor(surv_df$cluster, levels = cluster_color_df$cluster)
                fit <- survfit(Surv(time, status) ~ cluster, data=surv_df)
    
                ggsurv_plot = ggsurvplot(fit, data = surv_df, risk.table = TRUE, conf.int = FALSE,
                                         legend.title = "Cluster",
                                         ylim = c(0, 1),
                                         color = "cluster",
                                         palette = cluster_color_df$color,
                                         xlab = paste0("Time in ", sur_x_time),
                                         ylab = "Survival Probability",
                                         main = "Kaplan-Meier Plot by Cluster", dpi = 500)
              } else {
                fit <- survfit(Surv(time, status) ~ cluster, data=surv_df)
    
                ggsurv_plot = ggsurvplot(fit, data = surv_df, risk.table = TRUE, conf.int = FALSE, legend.title = "Cluster",
                                         ylim = c(0, 1),
                                         xlab = paste0("Time in ", sur_x_time),
                                         ylab = "Survival Probability",
                                         main = "Kaplan-Meier Plot by Cluster", dpi = 500)
              }
    
              plot_text = ""
    
              if(length(unique(as.character(surv_df$cluster))) == 2) {
                # Define the method based on the number of rows in surv_df
                method_type <- if (nrow(surv_df) < 25) "exact" else "efron"
                coxres = summary(coxph(Surv(time, status) ~ cluster, data=surv_df, method = method_type))
                surv_p_value = coxres$coefficients[[5]]
    
                if(surv_p_value < 0.01) {
                  surv_p_value = formatC(surv_p_value, format = "e", digits = 2)
                } else {
                  surv_p_value = format(round(surv_p_value, 2), nsmall = 2)
                }
    
                hr = format(coxres$conf.int[[1]], scientific = FALSE,digits = 2)
                hr_lower = format(round(coxres$conf.int[[3]], 2), nsmall = 2)
                hr_upper = format(round(coxres$conf.int[[4]], 2), nsmall = 2)
                plot_text = paste0("HR = ", hr, " (", hr_lower, "-", hr_upper, ")\n",
                                   "logrank p = ", surv_p_value)
              } else if (length(unique(as.character(surv_df$cluster))) > 2) {
                surv_p_value = surv_pvalue(surv_fit(Surv(time, status) ~ cluster, data=surv_df))$pval
    
                if(surv_p_value < 0.01) {
                  surv_p_value = formatC(surv_p_value, format = "e", digits = 2)
                } else {
                  surv_p_value = format(round(surv_p_value, 2), nsmall = 2)
                }
    
                plot_text = paste0("logrank p = ", surv_p_value)
              }
    
              ggsurv_plot$plot <- ggsurv_plot$plot +
                ggplot2::annotate(
                  "text",
                  x = Inf, y = Inf,
                  vjust = 1, hjust = 1,
                  label = plot_text,
                  size = 5
                )
    
    
              if(!is.null(output_png)) {
                png(output_png, width = 1800, height = 1800,
                    units = "px", pointsize = 6,
                    res = 150)
                print(ggsurv_plot)
                dev.off()  
              }
    
              if(!is.null(output_pdf)) {
                CairoPDF(out_pdf, family = "Arial")
                print(ggsurv_plot)
                dev.off()
                #ggsave(file = output_pdf, print(ggsurv_plot))
    
              }
            }
    
            #################
            # Analysis part #
            #################
    
            if(!is.null(ht_expression)) {
              if(is.null(ht_png) & is.null(ht_pdf)) {
                stop("[ERROR] no output heatmap pdf or png specified")
              }
              PlotHeatmap(expression_file = ht_expression,
                          out_png = ht_png,
                          out_pdf = ht_pdf,
                          annotation_file = ht_top_annotation,
                          ht_color_file = ht_top_annotation_color,
                          column_cluster_file = ht_column_clusters, 
                          row_cluster_file = ht_row_clusters, 
                          max_row_name = ht_max_row_name, 
                          max_column_name = ht_max_column_name, 
                          column_ann_border_max_sample = 100, 
                          ht_annotation_border_color = "black",
                          matrix_name = ht_matrix_name,
                          missing_ann_color = "lightgrey", 
                          ht_max_percentile = 99, 
                          ht_mid_percentile = 50, 
                          ht_min_percentile = 1, 
                          ht_max_color = 'red',
                          ht_mid_color = 'white', 
                          ht_min_color = 'blue')
            }
    
            """

            # Execute the R script
            robjects.r(rscript)


def gene_ont_all_gl_kegg(user_id, dataset_option, ht_row_cluster_path, df_index, gene_list, sub_ont, go_fun, max_gene,
                         corr_meth, q_value_cutoff, p_value_cutoff, title, max_category, order):
    """
    Perform KEGG enrichment analysis for all genes.

    Parameters:
        user_id (str): ID of the user.
        dataset_option (str): Selected dataset option.
        ht_row_cluster_path (str): Path to the row cluster data.
        df_index (list): List of data frame indices.
        gene_list (list): List of genes.
        sub_ont (str): Sub ontology type.
        go_fun (str): Gene Ontology function.
        max_gene (int): Maximum number of genes.
        corr_meth (str): Correction method for p-values.
        q_value_cutoff (float): Q-value cutoff.
        p_value_cutoff (float): P-value cutoff.

    Returns:
        None
    """

    with conversion.localconverter(robjects.default_converter + pandas2ri.converter):
        pandas2ri.activate()

        # Assign variables to R environment
        robjects.r.assign('sub_ont', sub_ont)
        robjects.r.assign('go_fun', go_fun)
        robjects.r(f'order <- "{order}"')
        robjects.r(f'title <- "{title}"')
        robjects.r(f'max_category <- "{max_category}"')
        robjects.r.assign('max_gene', max_gene)
        robjects.r.assign('corr_meth', corr_meth)
        robjects.r.assign('q_value_cutoff', q_value_cutoff)
        robjects.r.assign('p_value_cutoff', p_value_cutoff)
        robjects.r.assign('gene_clust_path', ht_row_cluster_path)
        robjects.r.assign('user_id_r', user_id)

        # Convert lists to R vectors
        df_index_r = robjects.StrVector(df_index.tolist())
        gene_list_r = robjects.StrVector(gene_list.tolist())
        robjects.globalenv['df_index_r'] = df_index_r
        robjects.globalenv['gene_list_r'] = gene_list_r

        # Define the R script to run
        r_script = '''
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(biomaRt)

        gene_clusters <- read.csv(gene_clust_path, header = TRUE, row.names = 1, sep = ",")
        input_genes <- gene_list_r

        if (length(input_genes) == 0) {
          stop("Input gene list is empty.")
        }
        print("Input Genes:")
        print(input_genes)

        gene_map <- bitr(geneID = input_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        if (nrow(gene_map) == 0) {
          stop("No genes were mapped to ENTREZID.")
        }
        print("Gene Map:")
        print(gene_map)
        gene_vector <- gene_map$ENTREZID

        max_gene <- as.numeric(max_gene)
        q_value_cutoff <- as.numeric(q_value_cutoff)
        p_value_cutoff <- as.numeric(p_value_cutoff)
        max_category = as.numeric(max_category)

        input_universe <- df_index_r
        if (length(input_universe) == 0) {
          stop("Universe gene list is empty.")
        }

        universe_vector <- bitr(geneID = input_universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        if (nrow(universe_vector) == 0) {
          stop("No universe genes were mapped to ENTREZID.")
        }
        print("Universe Vector:")
        print(universe_vector)
        universe_vector <- as.character(na.omit(universe_vector$ENTREZID))

        go_enrich_result <- enrichKEGG(
          gene = gene_vector,
          universe = universe_vector,
          pAdjustMethod = corr_meth,
          pvalueCutoff = p_value_cutoff,
          minGSSize = 10,
          maxGSSize = max_gene,
          qvalueCutoff = q_value_cutoff,
          organism = "hsa",
          keyType = "kegg"
        )

        if (is.null(go_enrich_result) || nrow(as.data.frame(go_enrich_result)) == 0) {
          stop("KEGG enrichment analysis returned no results.")
        }

        go_enrich_result_df <- as.data.frame(go_enrich_result)
        print("KEGG Enrichment Result:")
        print(go_enrich_result_df)

        if ("geneID" %in% colnames(go_enrich_result_df)) {
          kegg_ids <- unlist(strsplit(as.character(go_enrich_result_df$geneID), "/"))
          entrez_ids <- sub("hsa:", "", kegg_ids)
          if (length(entrez_ids) == 0) {
            stop("No KEGG IDs found in the result.")
          }
          print("Entrez IDs:")
          print(entrez_ids)
          hugo_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
          if (length(hugo_symbols) == 0) {
            stop("No HUGO symbols were mapped.")
          }
          print("HUGO Symbols:")
          print(hugo_symbols)
          hugo_symbols <- as.character(na.omit(hugo_symbols))
          go_enrich_result_df$HUGO_Symbols <- sapply(strsplit(as.character(go_enrich_result_df$geneID), "/"), function(ids) {
            paste(na.omit(hugo_symbols[entrez_ids %in% ids]), collapse = "/")
          })
        }

        output_csv_path <- paste0("result_data/result_data_", user_id_r, "/gene_ont_result_df_", user_id_r, ".csv")
        write.table(go_enrich_result_df, file = output_csv_path, row.names = FALSE, sep = ',')
        
        all_gs = dotplot(go_enrich_result, 
                 showCategory = max_category,    # Number of categories to display
                 color = "p.adjust",   # Color by adjusted p-value
                 size = "Count",       # Size by the number of genes in each category
                 orderBy = order, # Order categories by gene ratio
                 font.size = 12,       # Set font size
                 title = title)

        dot_plot_path <- paste0("result_data/result_data_", user_id_r, "/gene_ont_", user_id_r, ".png")
        png(dot_plot_path, res = 'max', height = 720, width = 640, pointsize = 32)
        print(all_gs)
        dev.off()

        dot_plot_path_pdf <- paste0("result_data/result_data_", user_id_r, "/gene_ont_", user_id_r, ".pdf")
        pdf(dot_plot_path_pdf, width = 8, height = 12, title = "gene_ontology")
        print(all_gs)
        dev.off()
        '''

        # Execute the R script
        robjects.r(r_script)


def gene_ont_all_gl_go(user_id, dataset_option, ht_row_cluster_path, df_index, gene_list, sub_ont, go_fun, max_gene,
                       corr_meth, q_value_cutoff, p_value_cutoff, title, max_category, order):
    """
    Perform Gene Ontology enrichment analysis for all genes.

    Parameters:
        user_id (str): ID of the user.
        dataset_option (str): Selected dataset option.
        ht_row_cluster_path (str): Path to the row cluster data.
        df_index (list): List of data frame indices.
        gene_list (list): List of genes.
        sub_ont (str): Sub ontology type.
        go_fun (str): Gene Ontology function.
        max_gene (int): Maximum number of genes.
        corr_meth (str): Correction method for p-values.
        q_value_cutoff (float): Q-value cutoff.
        p_value_cutoff (float): P-value cutoff.

    Returns:
        None
    """
    with openrlib.rlock:
        with conversion.localconverter(default_converter):
            pandas2ri.activate()

            # You can use robjects.r to create and assign R variables
            robjects.r(f'sub_ont <- "{sub_ont}"')
            robjects.r(f'go_fun <- "{go_fun}"')
            robjects.r(f'order <- "{order}"')
            robjects.r(f'title <- "{title}"')
            robjects.r(f'max_category <- "{max_category}"')
            robjects.r(f'max_gene <- "{max_gene}"')
            robjects.r(f'corr_meth <- "{corr_meth}"')
            robjects.r(f'q_value_cutoff <- "{q_value_cutoff}"')
            robjects.r(f'p_value_cutoff <- "{p_value_cutoff}"')
            robjects.r(f'gene_clust_path <- "{ht_row_cluster_path}"')
            robjects.r(f'user_id_r <- "{user_id}"')
            df_index = df_index.tolist()
            gene_list = gene_list.tolist()
            df_index_r = robjects.StrVector(df_index)
            robjects.globalenv['df_index_r'] = df_index_r
            gene_list_r = robjects.StrVector(gene_list)
            robjects.globalenv['gene_list_r'] = gene_list

            r_script = '''

            gene_clusters <- read.csv(gene_clust_path, header = TRUE, row.names = 1, sep = ",")

            # It will change to pd.DataFrame / pd.Serie from Python
            input_genes = gene_list_r
            gene_map = bitr(geneID = input_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
            gene_vector = gene_map$ENTREZID

            max_gene = as.numeric(max_gene)
            q_value_cutoff = as.numeric(q_value_cutoff)
            p_value_cutoff = as.numeric(p_value_cutoff)
            max_category = as.numeric(max_category)

            # Extract the universe gene list from the expression data frame
            universe_vector <- bitr(geneID = df_index_r, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

            # Perform Gene Ontology Enrichment Analysis
            go_enrich_result <- enrichGO(
              gene           = gene_vector,
              universe       = universe_vector,
              OrgDb          = org.Hs.eg.db,
              ont            = sub_ont, 
              pAdjustMethod  = corr_meth,
              qvalueCutoff   = q_value_cutoff,
              pvalueCutoff   = p_value_cutoff,
              readable       = TRUE,
              maxGSSize      = max_gene
            )
            
            all_gs = dotplot(go_enrich_result, 
                     showCategory = max_category,    # Number of categories to display
                     color = "p.adjust",   # Color by adjusted p-value
                     size = "Count",       # Size by the number of genes in each category
                     orderBy = order, # Order categories by gene ratio
                     font.size = 12,       # Set font size
                     title = title)

            # Specify the path for the CSV file
            output_csv_path <- paste0("result_data/result_data_" ,user_id_r, "/gene_ont_result_df_",user_id_r,".csv")
                        # Write the data frame to a CSV file
            write.table(go_enrich_result, file = output_csv_path, row.names = FALSE, sep=',')

            dot_plot_path <- paste0("result_data/result_data_" ,user_id_r, "/gene_ont_",user_id_r,".png")
            png(dot_plot_path,  res = 'max', height = 720, width = 640, pointsize = 32)
            print(all_gs)
            dev.off()

            dot_plot_path_pdf <- paste0("result_data/result_data_" ,user_id_r, "/gene_ont_",user_id_r,".pdf")
            pdf(dot_plot_path_pdf, width=8, height=12, title="gene_ontology")
            print(all_gs)
            dev.off()

            '''

            robjects.r(r_script)


def gene_ont_all_cluster_go(user_id, dataset_option, ht_row_cluster_path, df_index, gene_list, sub_ont, go_fun,
                            max_gene, corr_meth, q_value_cutoff, p_value_cutoff, title, max_category):
    """
    Perform Gene Ontology enrichment analysis for all gene clusters.

    Parameters:
        user_id (str): ID of the user.
        dataset_option (str): Selected dataset option.
        ht_row_cluster_path (str): Path to the row cluster data.
        df_index (list): List of data frame indices.
        gene_list (list): List of genes.
        sub_ont (str): Sub ontology type.
        go_fun (str): Gene Ontology function.
        max_gene (int): Maximum number of genes.
        corr_meth (str): Correction method for p-values.
        q_value_cutoff (float): Q-value cutoff.
        p_value_cutoff (float): P-value cutoff.

    Returns:
        None
    """
    with openrlib.rlock:
        with conversion.localconverter(default_converter):
            pandas2ri.activate()

            # You can use robjects.r to create and assign R variables
            robjects.r(f'gene_clust_path <- "{ht_row_cluster_path}"')
            robjects.r(f'sub_ont <- "{sub_ont}"')
            robjects.r(f'max_gene <- "{max_gene}"')
            robjects.r(f'corr_meth <- "{corr_meth}"')
            robjects.r(f'title <- "{title}"')
            robjects.r(f'max_category <- "{max_category}"')
            robjects.r(f'q_value_cutoff <- "{q_value_cutoff}"')
            robjects.r(f'p_value_cutoff <- "{p_value_cutoff}"')
            robjects.r(f'user_id_r <- "{user_id}"')
            df_index = df_index.tolist()
            gene_list = gene_list.tolist()
            df_index_r = robjects.StrVector(df_index)
            robjects.globalenv['df_index_r'] = df_index_r
            gene_list_r = robjects.StrVector(gene_list)
            robjects.globalenv['gene_list_r'] = gene_list

            r_script = '''

            input_genes_df = read.csv(gene_clust_path, header = TRUE, check.names = FALSE, sep = ",", row.names = 1)
            mapped_genes <- bitr(geneID = rownames(input_genes_df), fromType = "SYMBOL", toType = "ENTREZID", 
                                 OrgDb = org.Hs.eg.db)

            # Get a list of valid keys (gene symbols) from org.Hs.eg.db
            valid_keys <- keys(org.Hs.eg.db, keytype = "SYMBOL")

            max_gene = as.numeric(max_gene)
            q_value_cutoff = as.numeric(q_value_cutoff)
            p_value_cutoff = as.numeric(p_value_cutoff)
            max_category = as.numeric(max_category)

            # Map the universe vector
            # Extract the universe gene list from the expression data frame
            #universe_vector <- bitr(geneID = df_index_r, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

            # Extract unique clusters
            unique_clusters <- unique(input_genes_df$Cluster)

            # Initialize a list to hold the gene vectors for each cluster
            geneList <- vector("list", length(unique_clusters))
            names(geneList) <- unique_clusters

            # Loop over each unique cluster and fill the geneList
            for(cluster in unique_clusters) {
              cluster_genes <- rownames(input_genes_df[input_genes_df$Cluster == cluster, , drop = FALSE])
              valid_cluster_genes <- cluster_genes[cluster_genes %in% valid_keys]  # Now valid_keys is available
              geneList[[as.character(cluster)]] <- valid_cluster_genes
            }

            geneList_Entrez <- lapply(geneList, function(cluster_genes) {
              mapped_genes <- bitr(geneID = cluster_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
              mapped_genes$ENTREZID  # Return only the mapped Entrez IDs
            })

            comparison_result <- compareCluster(
              geneCluster = geneList_Entrez,
              fun = "enrichGO",
              OrgDb = org.Hs.eg.db,
              ont = sub_ont,
              pAdjustMethod = corr_meth,
              qvalueCutoff = q_value_cutoff,
              pvalueCutoff = p_value_cutoff,
              readable = TRUE,
              maxGSSize = max_gene
            )
            
            all_gs = dotplot(comparison_result, 
                 showCategory = max_category,    # Number of categories to display
                 color = "p.adjust",   # Color by adjusted p-value
                 size = "Count",       # Size by the number of genes in each category
                 font.size = 12,       # Set font size
                 title = title)

            # Specify the path for the CSV file
            output_csv_path <- paste0("result_data/result_data_" ,user_id_r, "/gene_ont_result_df_",user_id_r,".csv")

            # Write the data frame to a CSV file
            write.table(comparison_result, file = output_csv_path, row.names = FALSE, sep=',')

            dot_plot_path <- paste0("result_data/result_data_" ,user_id_r, "/gene_ont_",user_id_r,".png")
            png(dot_plot_path,  res = 'max', height = 720, width = 640, pointsize = 32)
            print(all_gs)
            dev.off()

            dot_plot_path_pdf <- paste0("result_data/result_data_" ,user_id_r, "/gene_ont_",user_id_r,".pdf")
            pdf(dot_plot_path_pdf,, width=8, height=12, title="gene_ontology")
            print(all_gs)
            dev.off()

            '''

            robjects.r(r_script)


def gene_ont_all_cluster_kegg(user_id, dataset_option, ht_row_cluster_path, df_index, gene_list, sub_ont, go_fun,
                              max_gene, corr_meth, q_value_cutoff, p_value_cutoff, title, max_category):
    """
    Perform KEGG enrichment analysis for all gene clusters.

    Parameters:
        user_id (str): ID of the user.
        dataset_option (str): Selected dataset option.
        ht_row_cluster_path (str): Path to the row cluster data.
        df_index (list): List of data frame indices.
        gene_list (list): List of genes.
        sub_ont (str): Sub ontology type.
        go_fun (str): Gene Ontology function.
        max_gene (int): Maximum number of genes.
        corr_meth (str): Correction method for p-values.
        q_value_cutoff (float): Q-value cutoff.
        p_value_cutoff (float): P-value cutoff.

    Returns:
        None
    """

    with conversion.localconverter(robjects.default_converter + pandas2ri.converter):
        pandas2ri.activate()

        # Assign variables to R environment
        robjects.r.assign('gene_clust_path', ht_row_cluster_path)
        robjects.r.assign('user_id_r', user_id)
        robjects.r(f'title <- "{title}"')
        robjects.r(f'max_category <- "{max_category}"')
        robjects.r.assign('max_gene', max_gene)
        robjects.r.assign('corr_meth', corr_meth)
        robjects.r.assign('q_value_cutoff', q_value_cutoff)
        robjects.r.assign('p_value_cutoff', p_value_cutoff)

        # Convert lists to R vectors
        df_index_r = robjects.StrVector(df_index.tolist())
        gene_list_r = robjects.StrVector(gene_list.tolist())
        robjects.globalenv['df_index_r'] = df_index_r
        robjects.globalenv['gene_list_r'] = gene_list_r

        # Define the R script to run
        r_script = '''
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(biomaRt)

        input_genes_df <- read.csv(gene_clust_path, header = TRUE, check.names = FALSE, sep = ",", row.names = 1)

        if (nrow(input_genes_df) == 0) {
          stop("The input gene cluster data frame is empty.")
        }
        print("Input Gene Clusters:")
        print(head(input_genes_df))

        mapped_genes <- bitr(geneID = rownames(input_genes_df), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        if (nrow(mapped_genes) == 0) {
          stop("No genes were mapped to ENTREZID.")
        }
        print("Mapped Genes:")
        print(head(mapped_genes))

        max_gene <- as.numeric(max_gene)
        q_value_cutoff <- as.numeric(q_value_cutoff)
        p_value_cutoff <- as.numeric(p_value_cutoff)
        max_category = as.numeric(max_category)

        unique_clusters <- unique(input_genes_df$Cluster)
        geneList <- vector("list", length(unique_clusters))
        names(geneList) <- unique_clusters

        for(cluster in unique_clusters) {
          cluster_genes <- rownames(input_genes_df[input_genes_df$Cluster == cluster, , drop = FALSE])
          valid_cluster_genes <- cluster_genes[cluster_genes %in% mapped_genes$SYMBOL]
          geneList[[as.character(cluster)]] <- valid_cluster_genes
        }

        geneList_Entrez <- lapply(geneList, function(cluster_genes) {
          mapped_genes <- bitr(geneID = cluster_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
          mapped_genes$ENTREZID
        })

        input_universe <- df_index_r
        if (length(input_universe) == 0) {
          stop("Universe gene list is empty.")
        }

        universe_vector <- bitr(geneID = input_universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        if (nrow(universe_vector) == 0) {
          stop("No universe genes were mapped to ENTREZID.")
        }
        print("Universe Vector:")
        print(head(universe_vector))
        universe_vector <- as.character(na.omit(universe_vector$ENTREZID))

        comparison_result <- compareCluster(
          geneCluster = geneList_Entrez,
          fun = "enrichKEGG",
          keyType = "kegg",
          pAdjustMethod = corr_meth,
          qvalueCutoff = q_value_cutoff,
          organism = 'hsa',
          universe = universe_vector,
          pvalueCutoff = p_value_cutoff,
          minGSSize = 10,
          maxGSSize = max_gene
        )

        if (is.null(comparison_result) || nrow(as.data.frame(comparison_result)) == 0) {
          stop("KEGG enrichment analysis returned no results.")
        }

        comparison_result_df <- as.data.frame(comparison_result)
        print("KEGG Enrichment Result:")
        print(head(comparison_result_df))

        if ("geneID" %in% colnames(comparison_result_df)) {
          kegg_ids <- unlist(strsplit(as.character(comparison_result_df$geneID), "/"))
          entrez_ids <- sub("hsa:", "", kegg_ids)
          if (length(entrez_ids) == 0) {
            stop("No KEGG IDs found in the result.")
          }
          print("Entrez IDs:")
          print(head(entrez_ids))
          hugo_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
          if (length(hugo_symbols) == 0) {
            stop("No HUGO symbols were mapped.")
          }
          print("HUGO Symbols:")
          print(head(hugo_symbols))
          hugo_symbols <- as.character(na.omit(hugo_symbols))
          comparison_result_df$HUGO_Symbols <- sapply(strsplit(as.character(comparison_result_df$geneID), "/"), function(ids) {
            paste(na.omit(hugo_symbols[entrez_ids %in% ids]), collapse = "/")
          })
        }

        output_csv_path <- paste0("result_data/result_data_", user_id_r, "/gene_ont_result_df_", user_id_r, ".csv")
        write.table(comparison_result_df, file = output_csv_path, row.names = FALSE, sep = ',')

        all_gs = dotplot(comparison_result, 
                 showCategory = max_category,    # Number of categories to display
                 color = "p.adjust",   # Color by adjusted p-value
                 size = "Count",       # Size by the number of genes in each category
                 font.size = 12,       # Set font size
                 title = title)
        
        

        dot_plot_path <- paste0("result_data/result_data_", user_id_r, "/gene_ont_", user_id_r, ".png")
        png(dot_plot_path, res = 'max', height = 720, width = 640, pointsize = 32)
        print(all_gs)
        dev.off()

        dot_plot_path_pdf <- paste0("result_data/result_data_", user_id_r, "/gene_ont_", user_id_r, ".pdf")
        pdf(dot_plot_path_pdf, width = 8, height = 12, title = "gene_ontology")
        print(all_gs)
        dev.off()
        '''

        # Execute the R script
        robjects.r(r_script)
