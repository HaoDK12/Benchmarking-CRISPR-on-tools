# A comparative analysis of CRISPR-Cas9 editing efficiency prediction tools

*In silico* design of CRISPR guide RNA (gRNA) with high efficiency and specificity is vital for successful application of the CRISPR-based gene editing technology. Although many machine and deep learning (MDL)-based tools have been developed to predict the editing activity of gRNAs, a systematic and unbiased evaluation of their predictive performance is still needed. Here, we provide a brief overview of the in silico tools for CRISPR design and assess the CRISPR datasets and statistical metrics used for evaluating model performance. We benchmark six representative CRISPR/-Cas9 editing efficiency prediction tools across six nine CRISPR datasets covering five six cell types and two three species. The deep-learning model *CRISPRon*, and *DeepHF* outperformed other models, exhibiting greater accuracy and better overall performance across U6 and T7-promoter based expression system, respectively. We compile all CRISPR datasets and in silico prediction tools into a GUIDEnet resource web portal, aiming to facilitate and streamline CRISPR design for diverse experiment requests. Furthermore, we summarize features affecting CRISPR gene editing activity, providing importanceimportant insights into model performance and the further development of more accurate CRISPR prediction models.


![image](https://github.com/HaoDK12/Benchmarking-CRISPR-on-tools/blob/main/bin/Figure3.png)

## Requirements
The scripts are written in R 3.6.2 and run on Windows OS. The versions of R packages which we used are, specifically:
``` 
patchwork_1.1.1 tsutils_0.9.4   psych_2.3.6
ggprism_1.0.3   reshape2_1.4.4  ggplot2_3.3.5
dplyr_1.0.7    openxlsx_4.2.4
```

## Contents
  - ./Supplementary Data:  This folder cointains the detail of collected datasetsand tools, the predictions of different tools.
  - ./bin: custom R scripts perform comparative analysis in detail.
  - ./Training: The training datasets of part of evaluated models

## Evaluate the new datasets and tools
The script Run_correlation_analysis.R saved in /bin folder performs Spearman correlation analysis and Steiger’s test for gRNA predicted scores. It allows comparing multiple scoring methods, visualizing the results and evaluates statistical differences against a reference score. you could run with the example data
```
Rscript Run_correlation_analysis.R <datafile.xlsx> <Reference_score> <compare_scores>
#E.g. Rscript Run_correlation_analysis.R "Supplementary_Data2.xlsx" "CRISPRon_score" "Deepspcas9_score,DeepHF_score,CRISPRdict_score,sgDesigner_score,Azimuth_score"
```
Which should end with
```
"Analysis complete. Results saved as correlation_heatmap.png and final_result_CRISPRon.csv"
```
### Note 
* Ensure you have the required R packages installed if run.
* The inputs include three parts: an Excel file (take Supplementary_Data2.xlsx as example, which presents the collected dataset and tool predictions from our results) containing the gRNA data to be evaluated and predictions from tools; a reference score to be compared against and other prediction scores to compare, which should be consistent with the column name of excel file and are utilized for Steiger's test. 
* The outputs include a heatmap showing the Spearman correlation coefficients between scores and observed frequencies; And a table provides Steiger’s test results, showing whether different scores have statistically different correlations.
* If you want to benchmark several novel tools, please add the predictions from tools as the new column end with _score to Supplementary_Data2.xlsx file and then run this script; Similarly, if need to evaluate based on new datasets, please add experiment records as the new rows to this file and make sure the add provides Actual.freq, representing observed gRNA efficiency.
## Citation
Hao Yuan, Chunping Song et al. An overview and comparative analysis of CRISPR-SpCas9 gRNA activity prediction tools. 2024 (Manuscript in revision version of *The CRISPR journal*)
