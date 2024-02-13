# A comparative analysis of CRISPR-Cas9 editing efficiency prediction tools

The *in silico* design and selection of guide RNA (gRNA) with high efficiency is an essential step for successful application of the Clustered Regularly Interspaced Short Palindromic Repeats (CRISPR)/associated protein 9 (Cas9)-based gene editing technology. Several computational tools based on machine learning and/or deep learning algorithms have been developed to predict the editing activity of gRNAs. However, there is still a need of systematically and unbiasedly assessing the predictive performance and application scenarios for these CRISPR design tools, thus facilitating context-dependent selection of the right predicting tool for CRISPR-based gene editing applications. In this study, we conducted a thorough comparative analysis of six most-used in silico CRISPR/Cas9 on-target activity prediction tools using six benchmark datasets representing diverse cell types and species. The benchmarking results generate a practical guide for the pre-experimental tools selection of these prediction tools, providing valuable assistance to researchers in increasing the success of CRISPR experiments. 

![image](https://github.com/HaoDK12/Benchmarking-CRISPR-on-tools/blob/main/img/Figure01.png)

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
  - ./Supplementary Data:  Comparison of the prediction results of different tools on six test sets  in our analysis.
  - ./Scripts: custom Python scripts to train our "seq-only model" and "full-features model".

## Citation
Hao Yuan, Chunping Song et al. A comparative analysis of CRISPR-Cas9 editing efficiency prediction tools. 2024 (Manuscript submitted)
