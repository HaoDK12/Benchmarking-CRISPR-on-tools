# A comparative analysis of CRISPR-Cas9 editing efficiency prediction tools

*In silico* design of CRISPR guide RNA (gRNA) with high efficiency and specificity is vital for successful application of the CRISPR-based gene editing technology. Although many machine and deep learning (MDL)-based tools have been developed to predict the editing activity of gRNAs, a systematic and unbiased evaluation of their predictive performance is still needed. Here, we provide a brief overview of the in silico tools for CRISPR design and assess the CRISPR datasets and statistical metrics used for evaluating model performance. We benchmark six representative CRISPR/-Cas9 editing efficiency prediction tools across six nine CRISPR datasets covering five six cell types and two three species. The deep-learning model *CRISPRon*, and *DeepHF* outperformed other models, exhibiting greater accuracy and better overall performance across U6 and T7-promoter based expression system, respectively. We compile all CRISPR datasets and in silico prediction tools into a GUIDEnet resource web portal, aiming to facilitate and streamline CRISPR design for diverse experiment requests. Furthermore, we summarize features affecting CRISPR gene editing activity, providing importanceimportant insights into model performance and the further development of more accurate CRISPR prediction models.


![image](https://github.com/HaoDK12/Benchmarking-CRISPR-on-tools/blob/main/bin/Figure3.png)

## Requirements
The scripts are written in Python 3.8.17 and run on Windows OS. The versions of Python packages which we used are, specifically:
```
  Numpy version: 1.22.0
  Pandas version: 2.0.3
  Bio version: 1.78
  Scikit-learn version: 0.22
  Scipy versions: 1.10.1
  xgboost version: 1.7.3
  Shap version: 0.43.0
```

## Contents
  - ./Supplementary Data:  This folder cointains the detail of collected datasetsand tools, the predictions of different tools.
  - ./bin: custom R scripts perform comparative analysis in detail.
  - ./Training: The training datasets of part of evaluated models

## Citation
Hao Yuan, Chunping Song et al. An overview and comparative analysis of CRISPR-SpCas9 gRNA activity prediction tools. 2024 (Manuscript submitted to *The CRISPR journal*)
