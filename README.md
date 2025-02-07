# A comparative analysis of CRISPR-Cas9 editing efficiency prediction tools

Design of guide RNA (gRNA) with high efficiency and specificity is vital for successful application of the CRISPR gene editing technology. Although many machine and deep learning-based tools have been developed to predict gRNA activities, a systematic and unbiased evaluation of their predictive performance is still needed. Here, we provide a brief overview of *in silico* tools for CRISPR design and assess the CRISPR datasets and statistical metrics used for evaluating model performance. We benchmark seven machine and deep learning-based CRISPR-Cas9 editing efficiency prediction tools across nine CRISPR datasets covering six cell types and three species. The deep-learning models *CRISPRon* and *DeepHF* outperform the other models exhibiting greater accuracy and higher Spearman correlation coefficient across multiple datasets. We compile all CRISPR datasets and *in silico* prediction tools into a GUIDEnet resource web portal, aiming to facilitate and streamline the sharing of CRISPR datasets. Furthermore, we summarize features affecting CRISPR gene editing activity, providing important insights into model performance and the further development of more accurate CRISPR prediction models.


![image](https://github.com/HaoDK12/Benchmarking-CRISPR-on-tools/blob/main/bin/Figure3.png)

## Requirements
The scripts are written in R 3.6.2 and run on Windows OS. The versions of R packages which we used are, specifically:
``` 
patchwork_1.1.1 tsutils_0.9.4   psych_2.3.6
ggprism_1.0.3   reshape2_1.4.4  ggplot2_3.3.5
dplyr_1.0.7    openxlsx_4.2.4
```

## Contents
  - ./Supplementary Data: This folder contains details of the collected datasets, tools, and predictions from different models.
  - ./bin: Custom R scripts for performing in-depth comparative analyses.
  - ./Training: The training datasets used for filtering out dataset to evaluate specific models.

## Evaluate the new datasets and tools
The scripts ```correlation_analysis.R``` and ```ndcg_analysis.R``` in the ```/bin``` folder perform Spearman correlation analysis and nDCG calculation, along with corresponding statistical tests. These scripts are designed to evaluate new CRISPR gRNA activity prediction tools and compare their ranking performance with existing models. For detailed information, please refer to our study (see Citation below).

Important: The test datasets used for benchmarking should be fully independent, even if you are using datasets from our study or adding new ones. We provide original test datasets and prediction results as examples to demonstrate how to run these analyses (```Supplementary_Data2.xlsx``` and ```Supplementary_Data3.xlsx``` in the ```/Supplementary Data``` folder). Moreover, for benchmarking new tools, please add their prediction results as new columns in the example files, using the _score suffix; Similarly, for evaluating new datasets, please add experimental records as new rows, ensuring the inclusion of the ```Actual.freq``` column to represent observed gRNA efficiency. Now you can start by running the workflow with the provided example data.
#### Run correlation_analysis.R
```
Rscript correlation_analysis.R <datafile.xlsx> <Reference_score> <compare_scores>
#E.g. Rscript correlation_analysis.R "Supplementary_Data2.xlsx" "CRISPRon_score" "Deepspcas9_score,DeepHF_score,CRISPRdict_score,sgDesigner_score,Azimuth_score"
```
#### Notes 
* Inputs 
1. An Excel file (e.g., Supplementary_Data2.xlsx) containing the gRNA data and prediction results from tools;
2. A reference model name to compare against
3. Other model name list separated by comma, which names must match the column name in Excel file (requisite for Steiger's test). 
* Outputs
1. correlation_heatmap.png: showing the Spearman correlation coefficients between scores and observed frequencies; 
2. cor_p_results.csv: a table provides Steiger’s test results, indicating whether differences between scores are statistically significant.

#### Run ndcg_analysis.R
```
Rscript ndcg_analysis.R <datafile.xlsx> "score1 score2 score3"
#E.g. Rscript ndcg_analysis.R "Supplementary_Data3.xlsx" "CRISPRon_score,Deepspcas9_score,DeepHF_score,CRISPRdict_score,sgDesigner_score,Azimuth_score"
```
#### Notes 
* Inputs
1. An Excel file (e.g., Supplementary_Data3.xlsx) containing the data that multiple gRNAs for single gene KO and various prediction results;
2. A comma-separated list of model names to compare, matching the column names in the Excel file. 
* Outputs
1. gene_ndcg_results.csv: Reports nDCG@10 and nDCG@20 values for each gene.
2. combined_ndcg_plot.png: Combined boxplot comparing nDCG@10 and nDCG@20.
3. nemenyi_test_ndcg10.png and nemenyi_test_ndcg20.png: Nemenyi test results nDCG comaprison, respectively.

## Citation
Hao Yuan, Chunping Song et al. An overview and comparative analysis of CRISPR-SpCas9 gRNA activity prediction tools. 2024 (Manuscript under revision for *The CRISPR journal*)
