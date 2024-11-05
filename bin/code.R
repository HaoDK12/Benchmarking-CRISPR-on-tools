setwd("***")


###load data
data <-  openxlsx::read.xlsx("Supplementary Data2.xlsx",sheet = 1,rowNames = F)
data$class <- paste0(data$Research,"_",data$Cell.line)


###data filter
library(dplyr)
library(stringr)

Deephf <- openxlsx::read.xlsx("./Training/DeepHF_training.xlsx",sheet = 1)
Deepsp <- read.csv("./Training/DeepSpCas9 (Library).csv",header = T)
mm <- read.csv("./Training/Moreno-Mateos.csv",header = T)

Deephf$gRNA_pam <- paste0(Deephf$gRNA_Seq,Deephf$PAM)
Deephf_grna <- Deephf$gRNA_pam
Deepsp_grna <- Deepsp$sequence23
Crispron_grna <- data[which(data$Research == "Xiang"),]$Seq

library(eulerr)
euler_plot <- euler(list(Xiang_grna = Crispron_grna ,Deephf_grna2 = unique(Deephf_grna),Deepsp_grna2 = unique(Deepsp_grna)),)
p1 <- plot(euler_plot,
           quantities=list(cex=0.8,col='gray20'))

xiang_overlap <- append(intersect(Crispron_grna,Deepsp_grna),intersect(Crispron_grna,Deephf_grna))
data01 <- data[-which(data$Seq %in% xiang_overlap & data$Research == "Xiang"),]


euler_plot <- euler(list(Chen_grna = data[which(data$Research == "Chen"),]$Seq,CRISPRon_grna = Crispron_grna, Deephf_grna2 = unique(Deephf_grna),Deepsp_grna2 = unique(Deepsp_grna)),)
p2 <- plot(euler_plot,
           quantities=list(cex=0.8,col='gray20'))


euler_plot <- euler(list(Gagnon_grna = data[which(data$Research == "Gagnon"),]$Seq,Moreno_Mateos_dataset = mm$sequence23))
p3 <- plot(euler_plot,
           quantities=list(cex=0.8,col='gray20'))

euler_plot <- euler(list(Varshney_grna = data[which(data$Research == "Varshney"),]$Seq,Moreno_Mateos_dataset = mm$sequence23))
p4 <- plot(euler_plot,
           quantities=list(cex=0.8,col='gray20'))
gridExtra::grid.arrange(p1,p2,p3,p4,nrow = 2)


# Apply the z-score function to each column and add new columns
data_raw <- data
data <- data01
data <- data %>%
  mutate(
    Azimuth_score_z = Azimuth_score*100,
    Deepspcas9_score_z = Deepspcas9_score,
    CRISPRdict_score_z = CRISPRdict_score*100,
    DeepHF_score_z = DeepHF_score*100,
    sgDesigner_score_z = sgDesigner_socre,
    CRISPRon_score_z = CRISPRon_score
  )
data$Baseline <- apply(data[,14:19],1,mean)

library(tidyverse)

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


###correlation analysis
library(ggplot2)
library(reshape2)
library(ggprism)

correlation_values$class<- factor(correlation_values$class, levels=c("Doenchv2_A375","Hart_HCT116","Labuhn_HEL","Koike_Mouse","Chen_HEK293T" ,                                                                       
                                                                     "Xiang_HEK293T","Gagnon_Zebrafish","Varshney_Zebrafish","Telboul_Mouse"))
correlation_values$score_variable<- factor(correlation_values$score_variable, levels=c("Baseline","Azimuth_score","sgDesigner_socre","CRISPRdict_score","CRISPRon_score", "Deepspcas9_score", "DeepHF_score"))
ggplot(correlation_values,aes(class,score_variable, fill=correlation, label=round(correlation,digits = 4))) +
  geom_tile() +
  scale_fill_gradient(high = "#FBFCB5",low="dodgerblue3",limit=c(0,1)) + labs(x = NULL, y = NULL, fill = "Spearman rank's\nCorrelation")+
  geom_text() +
  theme_classic()+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_prism(base_fontface = "plain",base_size = 12)

###Steiger's test for correlation
library(psych)

df_list <- data %>% group_split(class)

result_list_CRISPRon <- list()
result_list_DeepHF <- list()

for (i in c(1,2,4,5,6,9)) {
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
  
  result_list_CRISPRon[[i]] <- result_df
}
final_result_CRISPRon <- do.call(rbind, result_list_CRISPRon)
write.csv(final_result_CRISPRon,"./Cor_p_crispron.csv")

for (i in c(3,7,8)) {
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
  
  
  for (col in c("Deepspcas9_score","CRISPRon_score","CRISPRdict_score","sgDesigner_socre","Azimuth_score" )) {
    correlation <- cor(df$Actual.freq, df[[col]], method = "spearman", use = "pairwise.complete.obs")
    num_obs <- sum(!is.na(df$Actual.freq))
    if (col == "DeepHF_score") {
      p_value <- "-"
    } else {
      p_value <- psych::r.test(num_obs, correlation, cor(df$Actual.freq, df[["DeepHF_score"]], method = "spearman"))$p
    }
    result_df["p_value", col] <- p_value
    result_df["correlation", col] <- correlation
    result_df["num_obs", col] <- num_obs
    result_df[,"class"] <- unique(df[["class"]])
  }
  
  result_list_DeepHF[[i]] <- result_df
}
final_result_DeepHF <- do.call(rbind, result_list_DeepHF)
write.csv(final_result_DeepHF,"./Cor_p_deephf.csv")


###nDCG analysis
df <-  openxlsx::read.xlsx("Supplementary Data3.xlsx",sheet = 1,rowNames = F)
df$class <- paste0(df$Research,"_",df$Cell.line)

library(dplyr)
df <- df %>%
  mutate(
    Azimuth_score_z = Azimuth_score*100,
    Deepspcas9_score_z = Deepspcas9_score,
    CRISPRdict_score_z = CRISPRdict_score*100,
    DeepHF_score_z = DeepHF_score*100,
    sgDesigner_score_z = sgDesigner_score,
    CRISPRon_score_z = CRISPRon_score
  )
df$Baseline <- apply(df[,16:21],1,mean)

ndcg_list <- list()
# Define nDCG calculation function
calc_ndcg <- function(actual, predicted, n) {
  rank <- order(predicted, decreasing = TRUE)
  
  # Discounted Cumulative Gain (DCG)
  dcg <- sum(actual[rank[1:n]] / log2(1 + (1:n)))
  
  # Ideal Discounted Cumulative Gain (IDCG)
  ideal_rank <- order(actual, decreasing = TRUE)
  idcg <- sum(actual[ideal_rank[1:n]] / log2(1 + (1:n)))
  
  # Normalized Discounted Cumulative Gain (nDCG)
  ndcg <- dcg / idcg
  return(ndcg)
}
for (research in unique(df$Research)) {
  # Subset data for the current research group
  data <- subset(df, Research == research)
  
  # Calculate nDCG@2, nDCG@5, and nDCG@10 for each gene
  gene_ndcg <- data %>%
    group_by(Targetgene2) %>%
    summarise(
      Deepspcas9_score_ndcg10 = calc_ndcg(Actual.freq, Deepspcas9_score, n = 10),
      Deepspcas9_score_ndcg20 = calc_ndcg(Actual.freq, Deepspcas9_score, n = 20),
      DeepHF_score_ndcg10 = calc_ndcg(Actual.freq, DeepHF_score, n = 10),
      DeepHF_score_ndcg20 = calc_ndcg(Actual.freq, DeepHF_score, n = 20),
      CRISPRedict_score_ndcg10 = calc_ndcg(Actual.freq, CRISPRdict_score, n = 10),
      CRISPRedict_score_ndcg20 = calc_ndcg(Actual.freq, CRISPRdict_score, n = 20),
      sgDesigner_score_ndcg20 = calc_ndcg(Actual.freq, sgDesigner_score, n = 20),
      sgDesigner_score_ndcg10 = calc_ndcg(Actual.freq, sgDesigner_score, n = 10),
      Azimuth_score_ndcg20 = calc_ndcg(Actual.freq, Azimuth_score, n = 20),
      Azimuth_score_ndcg10 = calc_ndcg(Actual.freq, Azimuth_score, n = 10),
      CRISPRon_score_ndcg20 = calc_ndcg(Actual.freq, CRISPRon_score, n = 20),
      CRISPRon_score_ndcg10 = calc_ndcg(Actual.freq, CRISPRon_score, n = 10),
      Baseline_ndcg20 = calc_ndcg(Actual.freq, Baseline, n = 20),
      Baseline_ndcg10 = calc_ndcg(Actual.freq, Baseline, n = 10)
    )
  
  # Combine results and add to list
  ndcg_list[[research]] <- gene_ndcg
}

# Convert list to data frame
gene_ndcg_df <- do.call(rbind, ndcg_list)
openxlsx::write.xlsx(gene_ndcg_df, "./ndcg.xlsx", sheet = 1)

test_entirendcg <- data %>%
  summarise(
    Targetgene2 = "Konstantakos_nDCG_Cal",
    Deepspcas9_score_ndcg10 = calc_ndcg(Actual.freq, Deepspcas9_score, n = 10),
    
    DeepHF_score_ndcg10 = calc_ndcg(Actual.freq, DeepHF_score, n = 10),
    
    CRISPRedict_score_ndcg10 = calc_ndcg(Actual.freq, CRISPRdict_score, n = 10),
    
    sgDesigner_score_ndcg10 = calc_ndcg(Actual.freq, sgDesigner_score, n = 10),
    
    CRISPRon_score_ndcg10 = calc_ndcg(Actual.freq, CRISPRon_score, n = 10),
    Baseline_ndcg10 = calc_ndcg(Actual.freq, Baseline, n = 10),
    
    Deepspcas9_score_ndcg20 = calc_ndcg(Actual.freq, Deepspcas9_score, n = 20),
    
    DeepHF_score_ndcg20 = calc_ndcg(Actual.freq, DeepHF_score, n = 20),
    
    CRISPRedict_score_ndcg20 = calc_ndcg(Actual.freq, CRISPRdict_score, n = 20),
    
    sgDesigner_score_ndcg20 = calc_ndcg(Actual.freq, sgDesigner_score, n = 20),
    
    CRISPRon_score_ndcg20 = calc_ndcg(Actual.freq, CRISPRon_score, n = 20),
    Baseline_ndcg20 = calc_ndcg(Actual.freq, Baseline, n = 20),
    Azimuth_score_ndcg20 = NA,
    Azimuth_score_ndcg10=NA
  )
gene_ndcg_df <- rbind(gene_ndcg_df,test_entirendcg)

gene_ndcg_df_Doench <- reshape2::melt(gene_ndcg_df,id.vars="Targetgene2")
library(stringr)
library(ggplot2)
gene_ndcg_df_Doench$ndcg <-  sapply(str_split(gene_ndcg_df_Doench$variable, "_"), function(x) {
  ndcg_elements <- grep("ndcg", x, value = TRUE)
  if (length(ndcg_elements) > 0) {
    return(ndcg_elements)
  } else {
    return(NA)
  }
})
gene_ndcg_df_Doench$tools <- sapply(strsplit(as.character(gene_ndcg_df_Doench$variable), "_"), function(x) {
  if (x[1] == 'Baseline') {
    return(x[1]) 
  }else{
    paste(x[1:2], collapse = "_")
  }
})
gene_ndcg_df_Doench <- gene_ndcg_df_Doench[-which(gene_ndcg_df_Doench$tools == "Azimuth_score"),]

B1 <- gene_ndcg_df_Doench[which(gene_ndcg_df_Doench$ndcg == "ndcg20" & gene_ndcg_df_Doench$Targetgene2 != "Konstantakos_nDCG_Cal"),] %>% 
  group_by(tools) %>% 
  mutate(upper =  quantile(value,0.75,na.rm = T),
         lower = quantile(value,0.25,na.rm = T),
         mean = mean(value,na.rm = T),
         median = median(value,na.rm = T))
p1 <- ggplot(gene_ndcg_df_Doench[which(gene_ndcg_df_Doench$ndcg == "ndcg20" & gene_ndcg_df_Doench$Targetgene2 != "Konstantakos_nDCG_Cal"),], aes(tools, value,Targetgene2))+
  geom_jitter(aes(fill = Targetgene2),position = position_jitter(0.15),shape=21, size = 3)+
  stat_summary(fun = "mean",
               data = gene_ndcg_df_Doench[which(gene_ndcg_df_Doench$ndcg == "ndcg20" & gene_ndcg_df_Doench$Targetgene2 != "Konstantakos_nDCG_Cal"),],
               geom = "crossbar",
               mapping = aes(ymin=..y..,ymax=..y..),
               width=0.4,
               size=0.3)+
  ylab("nDCG@20 values")+
  geom_errorbar(data=B1, aes(ymin = lower, 
                             ymax = upper),width = 0.2,size=0.5)+theme_bw()+ # ??????
  theme(axis.text.x=element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9,face="bold"),
        plot.background = element_blank(), # ?????????
        panel.background = element_blank())+
  scale_fill_manual(values = colorRampPalette(c("#91bfdb", "#ffffbf", "#fc8d59", "#d73027"))(8))

B2 <- gene_ndcg_df_Doench[which(gene_ndcg_df_Doench$ndcg == "ndcg10" & gene_ndcg_df_Doench$Targetgene2 != "Konstantakos_nDCG_Cal"),] %>% 
  group_by(tools) %>% 
  mutate(upper =  quantile(value,0.75,na.rm = T),
         lower = quantile(value,0.25,na.rm = T),
         mean = mean(value,na.rm = T),
         median = median(value,na.rm = T))
p2 <- ggplot(gene_ndcg_df_Doench[which(gene_ndcg_df_Doench$ndcg == "ndcg10" & gene_ndcg_df_Doench$Targetgene2 != "Konstantakos_nDCG_Cal"),], aes(tools, value))+
  geom_jitter(aes(fill = Targetgene2),position = position_jitter(0.15),shape=21, size = 3)+
  stat_summary(fun = "mean",
               data = gene_ndcg_df_Doench[which(gene_ndcg_df_Doench$ndcg == "ndcg10" & gene_ndcg_df_Doench$Targetgene2 != "Konstantakos_nDCG_Cal"),],
               geom = "crossbar",
               mapping = aes(ymin=..y..,ymax=..y..),
               width=0.4,
               size=0.3)+
  ylab("nDCG@10 values")+
  geom_errorbar(data=B2, aes(ymin = lower, 
                             ymax = upper),width = 0.2,size=0.5)+theme_bw()+ # ??????
  theme(axis.text.x=element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9,face="bold"),
        plot.background = element_blank(), # ?????????
        panel.background = element_blank())+
  scale_fill_manual(values = colorRampPalette(c("#91bfdb", "#ffffbf", "#fc8d59", "#d73027"))(8))
combined_plot1 <- (p2/p1) + patchwork::plot_layout(guides = "collect")


###calculate Konstantakos et al' nDCG
p3 <- ggplot(gene_ndcg_df_Doench[which(gene_ndcg_df_Doench$ndcg == "ndcg20"),], aes(tools, value,Targetgene2))+
  geom_jitter(aes(fill = Targetgene2, shape = Targetgene2),position = position_jitter(0.15),size = 3)+
  stat_summary(fun = "mean",
               data = gene_ndcg_df_Doench[which(gene_ndcg_df_Doench$ndcg == "ndcg20" & gene_ndcg_df_Doench$Targetgene2 != "Konstantakos_nDCG_Cal"),],
               geom = "crossbar",
               mapping = aes(ymin=..y..,ymax=..y..),
               width=0.4,
               size=0.3)+
  geom_errorbar(data=B1, aes(ymin = lower, 
                             ymax = upper),width = 0.2,size=0.5)+theme_bw()+ # ??????
  theme(axis.text.x=element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9,face="bold"),
        plot.background = element_blank(), # ?????????
        panel.background = element_blank())+
  ylab("nDCG@20 values")+
  scale_shape_manual(values = c("Konstantakos_nDCG_Cal" = 8, # Triangle for entirendcg
                                "CCDC101" = 21,    # Default circle shapes for other categories
                                "CUL3" = 21,
                                "HPRT1" = 21,
                                "MED12" = 21,
                                "NF1" = 21,
                                "NF2" = 21,
                                "TADA1" = 21,
                                "TADA2B" = 21)) + 
  scale_fill_manual(values = c("Konstantakos_nDCG_Cal" = "black","CCDC101" ="#91BFDB","CUL3" ="#C0DACF","HPRT1" ="#EFF5C3","MED12" ="#FEDEA1","NF1" ="#FCAD76",
                               "NF2" ="#F67F51","TADA1" ="#E6573C","TADA2B" ="#D73027"))
p4 <- ggplot(gene_ndcg_df_Doench[which(gene_ndcg_df_Doench$ndcg == "ndcg10"),], aes(tools, value,Targetgene2))+
  geom_jitter(aes(fill = Targetgene2, shape = Targetgene2),position = position_jitter(0.15),size = 3)+
  stat_summary(fun = "mean",
               data = gene_ndcg_df_Doench[which(gene_ndcg_df_Doench$ndcg == "ndcg10" & gene_ndcg_df_Doench$Targetgene2 != "Konstantakos_nDCG_Cal"),],
               geom = "crossbar",
               mapping = aes(ymin=..y..,ymax=..y..),
               width=0.4,
               size=0.3)+
  geom_errorbar(data=B2, aes(ymin = lower, 
                             ymax = upper),width = 0.2,size=0.5)+theme_bw()+ # ??????
  theme(axis.text.x=element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9,face="bold"),
        plot.background = element_blank(), # ?????????
        panel.background = element_blank())+
  ylab("nDCG@10 values")+
  scale_shape_manual(values = c("Konstantakos_nDCG_Cal" = 8, # Triangle for entirendcg
                                "CCDC101" = 21,    # Default circle shapes for other categories
                                "CUL3" = 21,
                                "HPRT1" = 21,
                                "MED12" = 21,
                                "NF1" = 21,
                                "NF2" = 21,
                                "TADA1" = 21,
                                "TADA2B" = 21)) + 
  scale_fill_manual(values = c("Konstantakos_nDCG_Cal" = "black","CCDC101" ="#91BFDB","CUL3" ="#C0DACF","HPRT1" ="#EFF5C3","MED12" ="#FEDEA1","NF1" ="#FCAD76",
                               "NF2" ="#F67F51","TADA1" ="#E6573C","TADA2B" ="#D73027"))
combined_plot2 <- (p4/p3) + patchwork::plot_layout(guides = "collect")

###Friedman's test for nDCG 
library(scmamp)
library("Rgraphviz")
library(dplyr)
library(tidyr)

df1 <- gene_ndcg_df[,c("Deepspcas9_score_ndcg10",   "DeepHF_score_ndcg10",     
                       "CRISPRedict_score_ndcg10", "sgDesigner_score_ndcg10",   
                       "CRISPRon_score_ndcg10",   "Baseline_ndcg10")]

df2 <- gene_ndcg_df[,c("Deepspcas9_score_ndcg20",   "DeepHF_score_ndcg20",     
                       "CRISPRedict_score_ndcg20", "sgDesigner_score_ndcg20",   
                       "CRISPRon_score_ndcg20",   "Baseline_ndcg20")]


test1 <- tsutils::nemenyi(-df1,conf.level=0.90,plottype="vmcb")

test2 <- tsutils::nemenyi(-df2,conf.level=0.90,plottype="vmcb")


