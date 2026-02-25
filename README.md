# Peptide-to-protein data aggregation using Fisher’s method improves target identification in chemical proteomics 
PeptideFisher is a specialized Shiny application designed to aggregate peptide-level statistics into protein-level significance using Fisher’s Method. It is optimized for proteomics workflows, including **Differential Expression/Solubility** and **Partial Proteolysis (AFDIP/HOLSER)**.
---

## 📋 Table of Contents
1. [Prerequisites](#-prerequisites)
2. [Data Preparation](#-data-preparation)
3. [Step-by-Step Workflow](#-step-by-step-workflow)
4. [Statistical Methodology](#-statistical-methodology)

---
## 🛠 Prerequisites

To run this application locally, ensure you have **R (version 4.0+)** and the following libraries installed:

```R
install.packages(c("shiny", "dplyr", "readr", "readxl", "DT", "shinydashboard", "scales"))
```

---

## 📊 Data Preparation

Before starting the analysis, prepare your files according to the selected mode. Note that column names are **case-sensitive**.

### 1. **Expression / Solubility Mode**
Use this for standard quantitative proteomics where you compare protein abundance changes.

| Requirement | Description |
| :--- | :--- |
| **Peptide File** | Must contain an **`Accession`** column and numerical intensity columns for each replicate. |
| **Protein File** | Must contain an **`Accession`** column to match with the peptide data. |

### 2. **Partial Proteolysis (AFDIP / HOLSER) Mode**
Use this for structural change analysis.

| Requirement | Description |
| :--- | :--- |
| **Peptide File** | Must contain **`Accession`**, **`Gene name`**, and **`Protein description`** columns. |
| **Protein File** | **Not required** for this mode (directionality is calculated from peptides). |

> [!IMPORTANT]
> **Column Format**: Ensure all intensity values are **numeric**. Missing values (NA) are allowed but will be ignored during t-test calculations.

---



Overview

This application performs:

Protein-level differential analysis (Fold change + t-test)

Peptide-level differential analysis

Fisher’s method aggregation of top N peptides per protein

Integrated ranking score calculation

The method is described in:

	

Workflow

Analysis steps

Upload preprocessed protein-level data

Upload preprocessed peptide-level data

Select control and drug sample columns

Select Top N peptides used for Fisher aggregation

Run analysis

Download ranked protein table


	

Input Requirements

⚠️ The application does NOT perform preprocessing or filtering.

Users must prepare cleaned and normalized datasets before upload.


	
Protein File Requirements

The protein file must contain the following columns:


| Column name         | Description               |
| ------------------- | ------------------------- |
| Accession           | Unique protein identifier |
| Gene name           | Gene symbol               |
| Protein description | Protein annotation        |

Additionally:

At least 2 control sample columns

At least 2 drug sample columns

Sample columns must be numeric

Column names are case-sensitive


	

Peptide File Requirements

The peptide file must contain:

| Column name | Description                              |
| ----------- | ---------------------------------------- |
| Accession   | Protein identifier matching protein file |

Additionally:

At least 2 control sample columns

At least 2 drug sample columns

Sample columns must be numeric

Column names are case-sensitive



	
Output

The final output table contains:

| Column              | Description                    |
| ------------------- | ------------------------------ |
| Accession           | Protein ID                     |
| Gene name           | Gene symbol                    |
| Protein description | Annotation                     |
| log2FoldChange      | Protein-level fold change      |
| -log10fisher_p      | Fisher aggregated significance |
| score               | Combined ranking score         |
| rank_score          | Final protein ranking          |
