#!/usr/bin/env Rscript

# Command-line argument parsing
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop('Usage: Rscript ndcg_analysis.R "input_file.xlsx" "score1 score2 score3"')
}

input_file <- args[1]                      # First argument: Excel file
score_vars <- unlist(strsplit(args[2], ","))  # Second argument: Space-separated scores


library(openxlsx)
library(dplyr)

# Read the data
df <- openxlsx::read.xlsx(input_file, sheet = 1, rowNames = FALSE)
df$class <- paste0(df$Research, "_", df$Cell.line)

# nDCG Calculation Function
calc_ndcg <- function(actual, predicted, n) {
  rank <- order(predicted, decreasing = TRUE)
  dcg <- sum(actual[rank[1:n]] / log2(1 + (1:n)))
  idcg <- sum(actual[order(actual, decreasing = TRUE)[1:n]] / log2(1 + (1:n)))
  return(ifelse(idcg == 0, 0, dcg / idcg))  # Handle division by zero
}

# Calculate nDCG for each score
ndcg_list <- list()
for (research in unique(df$Research)) {
  data <- subset(df, Research == research)
  
  gene_ndcg <- data %>%
    group_by(Targetgene2) %>%
    summarise(across(all_of(score_vars), 
                     list(ndcg10 = ~calc_ndcg(Actual.freq, ., n = 10),
                          ndcg20 = ~calc_ndcg(Actual.freq, ., n = 20)), 
                     .names = "{.col}_{.fn}"))
  
  ndcg_list[[research]] <- gene_ndcg
}

# Combine all nDCG results
gene_ndcg_df <- do.call(rbind, ndcg_list)
write.csv(gene_ndcg_df, "gene_ndcg_results.csv", row.names = FALSE)

library(ggplot2)
library(reshape2)
library(patchwork)

# Reshape for plotting
gene_ndcg_df_long <- melt(gene_ndcg_df, id.vars = "Targetgene2")
gene_ndcg_df_long$ndcg_type <- ifelse(grepl("ndcg20", gene_ndcg_df_long$variable), "nDCG@20", "nDCG@10")
gene_ndcg_df_long$tools <- sub("_ndcg[0-9]+", "", gene_ndcg_df_long$variable)
B1 <- gene_ndcg_df_long[which(gene_ndcg_df_long$ndcg_type == "nDCG@10"),] %>% 
  group_by(tools) %>% 
  mutate(upper =  quantile(value,0.75,na.rm = T),
         lower = quantile(value,0.25,na.rm = T),
         mean = mean(value,na.rm = T),
         median = median(value,na.rm = T))
B2 <- gene_ndcg_df_long[which(gene_ndcg_df_long$ndcg_type == "nDCG@20"),] %>% 
  group_by(tools) %>% 
  mutate(upper =  quantile(value,0.75,na.rm = T),
         lower = quantile(value,0.25,na.rm = T),
         mean = mean(value,na.rm = T),
         median = median(value,na.rm = T))
# Plot for nDCG@10
p1 <- ggplot(gene_ndcg_df_long[gene_ndcg_df_long$ndcg_type == "nDCG@10",], 
             aes(tools, value, Targetgene2)) +
  geom_jitter(aes(fill = Targetgene2),position = position_jitter(0.15), shape = 21, size = 3) +
  stat_summary(fun = "mean", 
               geom = "crossbar", mapping = aes(ymin=..y..,ymax=..y..),width = 0.4, size = 0.3) +
  geom_errorbar(data=B1, aes(ymin = lower, 
                             ymax = upper),width = 0.2,size=0.5)+theme_bw()+
  theme(axis.text.x=element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9,face="bold"),
        plot.background = element_blank(), # ?????????
        panel.background = element_blank())+
  ylab("nDCG@10 Values")

# Plot for nDCG@20
p2 <- ggplot(gene_ndcg_df_long[gene_ndcg_df_long$ndcg_type == "nDCG@20",], 
             aes(tools,value,Targetgene2)) +
  geom_jitter(aes(fill = Targetgene2),position = position_jitter(0.15), shape = 21, size = 3) +
  stat_summary(fun = "mean", 
               geom = "crossbar", mapping = aes(ymin=..y..,ymax=..y..),width = 0.4, size = 0.3) +
  geom_errorbar(data=B2, aes(ymin = lower, 
                             ymax = upper),width = 0.2,size=0.5)+theme_bw()+
  ylab("nDCG@20 Values") 

# Combined plot
combined_plot1 <- p1 / p2 + plot_layout(guides = "collect")
ggsave("combined_ndcg_plot.png", combined_plot1, width = 10, height = 8)

library(tsutils)

# Prepare data for test
ndcg10_data <- gene_ndcg_df[, grep("ndcg10", colnames(gene_ndcg_df))]
ndcg20_data <- gene_ndcg_df[, grep("ndcg20", colnames(gene_ndcg_df))]

# Perform Nemenyi Test
png("nemenyi_test_ndcg10.png", width = 800, height = 600)
test1 <- nemenyi(-ndcg10_data, conf.level = 0.90, plottype = "vmcb")
dev.off()

png("nemenyi_test_ndcg20.png", width = 800, height = 600)
test2 <- nemenyi(-ndcg20_data, conf.level = 0.90, plottype = "vmcb")
dev.off()
