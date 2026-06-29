# Thigh Muscle Fat Infiltration (TMFI): Multi-Omics and Systems-Level Analysis in UK Biobank

## Overview

This repository contains the analysis code for the study:

**"Integrative analysis uncovers the effects and multi-omics features of thigh muscle fat infiltration (TMFI)"**

We systematically investigated the associations between thigh muscle fat infiltration (TMFI) and human health outcomes, and explored its underlying genetic, cellular, proteomic, and metabolomic mechanisms using UK Biobank data.

Key analyses include:
- Epidemiological associations of TMFI with mortality and disease outcomes
- Mediation analysis of lifestyle factors via TMFI
- Genome-wide association study (GWAS) of TMFI
- Polygenic risk score (PRS) construction and validation
- Single-cell transcriptomic integration (scPagwas)
- SMR and TWAS analyses
- Proteomic and metabolomic profiling
- Deep learning-based biomarker signatures

---

## Study Design

We used data from the **UK Biobank (UKB)** imaging cohort (V2 visit), including:

- ~46,000 participants for GWAS of TMFI
- ~50,000+ participants with MRI-derived TMFI measures
- ~54,000 proteomics samples (Olink)
- ~280,000 metabolomics samples (NMR)
- Independent validation cohorts for PRS and omics signature analyses

TMFI was derived from Dixon MRI of thigh muscles, representing intramuscular fat infiltration across four thigh regions.

---

## Key Outcomes

We evaluated associations between TMFI and:

- All-cause mortality
- Cause-specific mortality (CVD, cancer, respiratory, etc.)
- 12 system-specific disease categories (ICD-10)
- Aging-related and inflammation-related biomarkers

---

## Main Analytical Components

### 1. Epidemiological Analysis
- Cox proportional hazards models
- Logistic regression for disease prevalence
- Restricted cubic spline (RCS) models for dose–response
- Stratified analyses by age and sex

### 2. Mediation Analysis
- Lifestyle factors (diet, smoking, physical activity, sleep, BMI, alcohol)
- TMFI as mediator
- Bootstrap (1,000 resamples)

### 3. Genome-Wide Association Study (GWAS)
- PLINK2 QC filtering
- fastGWA-GLMM (GCTA v1.94.1)
- Covariates: age, sex, education, genotype batch, 20 PCs

### 4. Post-GWAS Annotation
- FUMA SNP-to-gene mapping
- MAGMA gene-based analysis
- GTEx tissue enrichment
- GWAS Catalog enrichment

### 5. Polygenic Risk Score (PRS)
- Clumping and thresholding (C+T)
- PRS-CS Bayesian model (LD reference: 1000G EUR)
- Independent validation in non-GWAS UKB cohort

### 6. Single-Cell Analysis
- scPagwas integration with muscle single-cell RNA-seq
- Cell-type trait relevance scoring (TRS)
- Identification of myogenic cell subpopulations

### 7. SMR & TWAS
- GTEx v8 tissue expression integration
- SMR (HEIDI test)
- S-PrediXcan TWAS analysis

### 8. Proteomics & Metabolomics
- Olink Explore (2,923 proteins)
- NMR metabolomics (251 metabolites)
- Linear regression + multiple testing correction (BH-FDR)
- Mediation analysis of circulating biomarkers

### 9. Deep Learning Signature Models
- 1D Convolutional Neural Networks (CNN)
- Input: full proteome / metabolome
- Output: predicted TMFI
- Evaluation: Pearson correlation, MAE, MSE

---

## Data Sources

### UK Biobank
- https://www.ukbiobank.ac.uk/

### Single-cell dataset
- GEO accession: **GSE214544**
- Vastus lateralis muscle cell atlas

### GTEx v8
- Tissue-specific eQTL and expression data

### HPA (Human Protein Atlas)

### DGIdb
- Drug–gene interaction database

---

## Repository Structure
