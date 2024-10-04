setwd("***")

data <-  openxlsx::read.xlsx("Total_benchmark_result.xlsx",sheet = 1,rowNames = F)
data$class <- paste0(data$Research,"_",data$Cell.line)

###calculate correlation
library(tidyverse)
library(RColorBrewer)

correlation_values <- data.frame(class = character(), score_variable = character(), correlation = numeric(), Pcorrelation = numeric(), stringsAsFactors = FALSE)
score_variables <- c("Deepspcas9_score", "DeepHF_score", "CRISPRdict_score", "sgDesigner_socre", "CRISPRon_score", "Azimuth_score", "Baseline")


for (class_value in unique(data$class)) {
  
  class_data <- subset(data, class == class_value)
  
  for (score_variable in score_variables) {
    
    correlation <- cor(class_data[[score_variable]], class_data$Actual.freq, method = "spearman")
    Pcorrelation <- cor.test(class_data[[score_variable]], class_data$Actual.freq, method = "spearman")$p.value
    correlation_values <- rbind(correlation_values, data.frame(class = class_value, score_variable = score_variable, correlation = correlation, Pcorrelation = Pcorrelation,stringsAsFactors = FALSE))
  }
}
write.csv(correlation_values,"***/Correlation_values.csv")

###Correlation plot;add class column to correlation_values
correlation_values <- read.csv("***/Correlation_values2.csv",header = T,sep=';')
correlation_values$class <-factor(correlation_values$class,ordered=TRUE,levels=c("Chen_HEK293T","Xiang_HEK293T","Labuhn_HEL","Doench_A375","Hart_HCT116","Koike_Mouse"))
B <- correlation_values %>% 
  group_by(score_variable) %>% 
  mutate(upper =  quantile(correlation,0.75,na.rm = T),
         lower = quantile(correlation,0.25,na.rm = T),
         mean = mean(correlation,na.rm = T),
         median = median(correlation,na.rm = T))
ggplot(correlation_values, aes(score_variable, correlation))+
  geom_jitter(aes(fill = class),position = position_jitter(0.15),shape=21, size = 3)+
  stat_summary(fun = "mean",
               geom = "crossbar",
               mapping = aes(ymin=..y..,ymax=..y..),
               width=0.4,
               size=0.3)+
  geom_errorbar(data=B, aes(ymin = lower, 
                            ymax = upper),width = 0.2,size=0.5)+theme_bw()+ 
  theme(axis.text.x=element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9,face="bold"),
        plot.background = element_blank(), 
        panel.background = element_blank())+ 
  scale_fill_manual(values = colors <- c("#00205B", "#006400", "#E60000", "#FF8000", "#FFC000", "#660099"))


### Individual correlation plot for each tool
columns <- c("Azimuth_score", "Deepspcas9_score", "CRISPRdict_score", "DeepHF_score", "sgDesigner_socre", "CRISPRon_score", "Baseline")

plot_list <- list()
library(gridExtra)

for (group_name in unique(data$class)) {
  group_data <- subset(data, class == group_name)
  
  group_plot_list <- list()
  
  for (column in columns) {
    correlation <- cor(group_data[[column]], group_data$Actual.freq, method = "spearman")
    p_value <- cor.test(group_data[[column]], group_data$Actual.freq, method = "spearman")$p.value
    
    plot <- ggplot(data = group_data, aes(x = .data[[column]], y = Actual.freq)) +
      geom_point(size=0.6) +
      geom_smooth(method = "lm", se = T, color = "blue") +
      labs(x = column, y = "Actual.freq") +
      ggtitle(paste("Correlation with ", unique(group_data$class)))+
      theme_bw()+
      theme(plot.title = element_text(size = 10))
#      geom_text(x = max(group_data[[column]]), y = max(group_data$Actual.freq), 
#                label = paste("Correlation:", round(correlation, 5), "p-value:", round(p_value, 5)), 
#                hjust = 1, vjust = 1,size=2)
    
    group_plot_list[[column]] <- plot
  }
  combined_plot <- do.call(grid.arrange, c(group_plot_list, ncol = 3))
  

  plot_list[[group_name]] <- combined_plot
  

  for (i in seq_along(plot_list)) {
    file_name <- paste(group_name, "_plot", ".pdf", sep = "")
    ggsave(filename = file_name, plot = plot_list[[i]], width = 9, height = 9)
  }
}

###Steiger test
library(psych)
df_list <- data %>% group_split(class)
result_list <- list()
for (i in 1:length(df_list)) {
  df <- df_list[[i]]
  
  
  result_df <- data.frame(
    Deepspcas9_score = numeric(),
    DeepHF_score = numeric(),
    CRISPRdict_score = numeric(),
    sgDesigner_socre = numeric(),
    Azimuth_score = numeric(),
    CRISPRon_score = numeric(),
    class = character()
  )
  
  for (col in c("Deepspcas9_score","DeepHF_score","CRISPRdict_score","sgDesigner_socre","Azimuth_score" )) {
    correlation <- cor(df$Actual.freq, df[[col]], method = "spearman", use = "pairwise.complete.obs")
    num_obs <- sum(!is.na(df$Actual.freq))
    if (col == "CRISPRon_score") {
      p_value <- "-"
    } else {
      p_value <- psych::r.test(num_obs, correlation, cor(df$Actual.freq, df[["CRISPRon_score"]], method = "spearman"))$p
    }
    result_df["p_value", col] <- p_value
    result_df["correlation", col] <- correlation
    result_df["num_obs", col] <- num_obs
    result_df[,"class"] <- unique(df[["class"]])
  }
  
  result_list[[i]] <- result_df
}
final_result <- do.call(rbind, result_list)
write.csv(final_result,"***/Steiger_test_result.csv")

###NDCG calculation
df <-  openxlsx::read.xlsx("Total_benchmark_result.xlsx",sheet = 1,rowNames = F)
df$class <- paste0(df$Research,"_",df$Cell.line)
df <- df[-which(df$Research %in% c("Xiang","Doench")),]


calc_ndcg <- function(actual, predicted, n) {

  rank <- order(predicted, decreasing = TRUE)
  
  #Discounted Cumulative Gain (DCG)
  dcg <- sum(actual[rank[1:n]] / log2(1 + (1:n)))
  
  ideal_rank <- order(actual, decreasing = TRUE)
  
  #Ideal Discounted Cumulative Gain (IDCG)
  idcg <- sum(actual[ideal_rank[1:n]] / log2(1 + (1:n)))
  
  #Normalized Discounted Cumulative Gain (nDCG)
  ndcg <- dcg / idcg
  

  return(ndcg)
}

ndcg_list <- lapply(unique(df$class), function(x) {
  
  data <- subset(df, class == x)
  
  actual <- data$Actual.freq
  preds <- data[, c("Deepspcas9_score", "DeepHF_score", "CRISPRdict_score", "sgDesigner_socre", "Azimuth_score","CRISPRon_score","Baseline")]
  
  ndcg_list <- sapply(preds, function(x) {
    calc_ndcg(actual, x, n = 20)
  })
  
  result <- data.frame(nDCG = ndcg_list, class = x, source = colnames(preds))
  
  return(result)
})

result_df <- do.call(rbind, ndcg_list)
print(result_df)
result <- do.call(rbind, ndcg_list)
openxlsx::write.xlsx(result,"***/Ndcg_values.xlsx",sheet=1)

###NDCG plot 
df$class <-factor(df$class,ordered=TRUE,levels=c("Chen_HEK293T","Xiang_HEK293T","Labuhn_HEL","Doench_A375","Hart_HCT116","Koike_Mouse"))
 df %>% ggplot(aes(source,nDCG,fill=class))+
  geom_col(position = position_dodge(0.8),width = 0.6)+
   theme_bw()+ # ??????
   theme(axis.text.x=element_text(color="black",size=8),
         axis.title = element_text(color="black",size=9,face="bold"),
         plot.background = element_blank(),
         panel.background = element_blank())+
   scale_fill_manual(values = colors <- c("#00205B", "#006400", "#E60000", "#FF8000", "#FFC000", "#660099", "#FF66A1", "#00BFFF", "#00FF00", "#D3D3D3")) # ????????????

####
dfr <- reshape2::recast(df,class~source)
openxlsx::write.xlsx(dfr,"***/For_Ndcgtest.xlsx",sheet=1)


### The Friedman test and a post hoc Nemenyi test
data <- openxlsx::read.xlsx("***/For_Ndcgtest.xlsx",sheet = 3,rowNames = T)
data2 <- data.frame(avg = apply(data,2,mean))
b <- tsutils::nemenyi(data,conf.level=0.95,plottype="vmcb")


### Correlation plot from two Xgboost model toward consistent internal test set
data1 <-  openxlsx::read.xlsx("***/train/single-feature_pred.xlsx",sheet = 1,rowNames = F)
data2 <-  openxlsx::read.xlsx("***/train/multi-feature_pred.xlsx",sheet = 1,rowNames = F)
p1 <- ggplot(data = data1, aes(x = data1$Efficiency, y = data1$Pred)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  labs(x = "Real efficiency" , y = "Predicted efficiency of 'multi-feature model'") +
  theme_bw()+
  ggpubr::stat_cor(data=data1, method = "spearman")
p2 <- ggplot(data = data2, aes(x = data2$Efficiency, y = data2$Pred)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  labs(x = "Real efficiency" , y = "Predicted efficiency of 'single-feature model'") +
  theme_bw()+
  ggpubr::stat_cor(data=data2, method = "spearman")


