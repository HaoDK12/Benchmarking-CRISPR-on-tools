#!/usr/bin/env Rscript

###Input args

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_correlation_analysis.R <datafile.xlsx> <target_score> <compare_scores_csv>")
}

datafile <- args[1]  # ?????? Excel ????????????
target_score <- args[2]  # ?????? score
compare_scores <- unlist(strsplit(args[3], ","))

###Load library
suppressMessages({
  library(openxlsx)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(ggprism)
  library(psych)
})


###Load datasets

data <- openxlsx::read.xlsx(datafile, sheet = 1, rowNames = FALSE)
data$class <- paste0(data$Research, "_", data$Cell.line)


###Measure correlation
score_variables <- colnames(data)[grepl("_score", colnames(data))] 
correlation_values <- data.frame(class = character(), score_variable = character(), correlation = numeric(), Pcorrelation = numeric(), stringsAsFactors = FALSE)

for (class_value in unique(data$class)) {
  class_data <- subset(data, class == class_value)
  
  for (score_variable in score_variables) {
    correlation <- cor(class_data[[score_variable]], class_data$Actual.freq, method = "spearman")
    Pcorrelation <- cor.test(class_data[[score_variable]], class_data$Actual.freq, method = "spearman")$p.value
    correlation_values <- rbind(correlation_values, data.frame(class = class_value, score_variable = score_variable, correlation = correlation, Pcorrelation = Pcorrelation, stringsAsFactors = FALSE))
  }
}


###Plot heatmap
correlation_values$class <- factor(correlation_values$class, levels = unique(correlation_values$class))
correlation_values$score_variable <- factor(correlation_values$score_variable, levels = score_variables)

p <- ggplot(correlation_values, aes(class, score_variable, fill = correlation, label = round(correlation, digits = 4))) +
  geom_tile() +
  scale_fill_gradient(high = "#FBFCB5", low = "dodgerblue3", limit = c(0, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman rank's\nCorrelation") +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_prism(base_fontface = "plain", base_size = 12)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
  

ggsave("correlation_heatmap.png", plot = p, width = 12, height = 6, dpi = 600)


###Steiger's Test for Correlation

df_list <- data %>% group_split(class)
group_names <- data %>% group_keys(class) %>% pull(class)
names(df_list) <- group_names

result_list_CRISPRon <- list()

for (i in 1:length(df_list)) {
  df <- df_list[[i]]
  
  result_df <- data.frame(class = unique(df$class))
  
  for (col in compare_scores) {
    correlation <- cor(df$Actual.freq, df[[col]], method = "spearman", use = "pairwise.complete.obs")
    num_obs <- sum(!is.na(df$Actual.freq))
    
    if (col == target_score) {
      p_value <- "-"
    } else {
      cor1 <- cor(df$Actual.freq, df[[target_score]], method = "spearman")
      cor2 <- cor(df[[col]], df[[target_score]], method = "spearman")
      p_value <- psych::r.test(num_obs, correlation, cor1, cor2)$p
    }
    
    result_df[[paste0("p_value_", col)]] <- p_value
    result_df[[paste0("correlation_", col)]] <- correlation
  }
  
  result_list_CRISPRon[[i]] <- result_df
}

final_result_CRISPRon <- do.call(rbind, result_list_CRISPRon)
write.csv(final_result_CRISPRon, "./final_result_CRISPRon.csv", row.names = FALSE)

print("Analysis complete. Results saved as correlation_heatmap.png and final_result_CRISPRon.csv")
