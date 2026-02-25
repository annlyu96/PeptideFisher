# Peptide-to-protein data aggregation using Fisher’s method improves target identification in chemical proteomics 
PeptideFisher is a specialized Shiny application designed to aggregate peptide-level statistics into protein-level significance using Fisher’s Method. It is optimized for proteomics workflows, including **Differential Expression/Solubility** and **Partial Proteolysis (AFDIP/HOLSER)**.
---

## 📋 Table of Contents
1. [Prerequisites](#-prerequisites)
2. [Data Preparation](#-data-preparation)
3. [Step-by-Step Workflow](#-step-by-step-workflow)
4. [Statistical Methodology](#-statistical-methodology)

---
## Prerequisites

To run this application locally, ensure you have **R (version 4.0+)** and the following libraries installed:

```R
install.packages(c("shiny", "dplyr", "readr", "readxl", "DT", "shinydashboard", "scales"))
```

**Input Requirements**

⚠️ The application does NOT perform preprocessing or filtering.

Users must prepare cleaned and normalized datasets before upload.

---

## Data Preparation

Before starting the analysis, prepare your files according to the selected mode. Note that column names are **case-sensitive**.

### 1. **Expression / Solubility Mode**
Use this for standard quantitative proteomics where you compare protein abundance or solubility changes.

| Requirement | Description |
| :--- | :--- |
| **Peptide File** | Must contain an **`Accession`**, **`Gene name`**, and **`Protein description`** columns and >2 intensity columns for each replicate. |
| **Protein File** | Must contain an **`Accession`** column to match with the peptide data. |

### 2. **Partial Proteolysis (AFDIP / HOLSER) Mode**
Use this for structural change analysis.

| Requirement | Description |
| :--- | :--- |
| **Peptide File** | Must contain **`Accession`**, **`Gene name`**, and **`Protein description`** columns. |
| **Protein File** | **Not required** for this mode (directionality is calculated from peptides). |

* **`Accession`**: Unique protein identifier.
* **`Gene name`**: Official gene symbol.
* **`Protein description`**: Full name or description of the protein.

> [!TIP]
> **Why 2+ replicates?** A minimum of 2 samples per group is statistically required to calculate the standard deviation for the Welch's T-test. If you select only one column, the app will display a validation error.

---

---

## Step-by-Step Workflow

Follow these steps to process your proteomics data and generate the Fisher scores.

### **Step 1: Launch the Application**
Run the R script in RStudio. Once the dashboard opens, you will see the **Sidebar** for configuration and the **Main Body** for results.

### **Step 2: Select Data Mode**
Choose the appropriate analysis mode based on your experimental design:
* **Expression / Solubility**: Use this if you have both peptide-level and protein-level abundance files.
* **Partial Proteolysis (AFDIP / HOLSER)**: Use this if you are analyzing partial proteolysis proteomics data or structural changes using only peptide-level data.

### **Step 3: Data Upload & Column Mapping**
1.  **Upload Peptide File**: Click `Browse` to upload your peptide data (`.csv`, `.xlsx`, or `.txt`).
2.  **Select Peptide Replicates**:
    * **Peptide Control**: Select at least **2 columns** representing your control samples.
    * **Peptide Drug**: Select at least **2 columns** representing your treated samples.
3.  **Handle Protein Data (Expression Mode Only)**:
    * Upload the **Protein File**.
    * **Protein Control/Drug**: Similarly, select at least **2 columns** for each group.
    * **Note**: Both files MUST contain the columns: **`Accession`**, **`Gene name`**, and **`Protein description`**.

### **Step 4: Configure Parameters & Run**
1.  **Top N peptides**: Enter the number of top-ranking peptides (by p-value) to include in the aggregation (default is **4**).
2.  **Click "Run Analysis"**: 
    * A **Progress Bar** will appear at the top of the dashboard.
    * The app will first **validate** your inputs. If you selected fewer than 2 replicates or are missing mandatory columns, a red error message will appear.

### **Step 5: Review and Download**
1.  **Live Preview**: Once the progress bar reaches **"Done!"**, the **Fisher Result** table will update automatically.
2.  **Sorting**: You can sort by **rank** or **score** directly in the table.
3.  **Download**: Click the **Download Result** button to save the full analysis as a CSV file.

**The final output table contains:**

| Column              | Description                                             |
| ------------------- | ------------------------------------------------------- |
| Accession           | Protein ID                                              |
| Gene name           | Gene symbol                                             |
| Protein description | Annotation                                              |
| log2FoldChange      | Protein or peptide-level fold change, depending on mode |
| -log10fisher_p      | Fisher aggregated significance                          |
| score               | Combined ranking score                                  |
| rank_score          | Final protein ranking                                   |


---
---

## Statistical Methodology

The application automates a multi-stage statistical pipeline to bridge the gap between peptide-level measurements and protein-level significance.

### **1. Peptide-level Statistics**
For each peptide, the app performs a **Welch’s T-test** (unpaired, two-sided) between the Control and Drug groups. 
* **Input**: At least 2 replicates per group.
* **Output**: A raw $p$-value and a $\log_2(\text{Fold Change})$.
* **Handling Errors**: If a peptide has zero variance or insufficient data, the $p$-value is returned as `NA`.

### **2. Top N Selection**
To reduce noise from low-quality or non-responsive peptides, the app:
1.  Groups all peptides by their **`Accession`**.
2.  Sorts them by significance ($-\log_{10} p\text{-value}$).
3.  Selects only the **Top N** (user-defined) most significant peptides for the next step.

### **3. Fisher’s Combination Method**
The $p$-values from the Top N peptides are aggregated into a single protein-level chi-squared ($\chi^2$) statistic using **Fisher’s Method**:

$$\chi^2 = -2 \sum_{i=1}^{n} \ln(p_i)$$

Where:
* $n$ is the number of peptides selected ($n \le N$).
* The resulting $\chi^2$ follows a Chi-squared distribution with $2n$ degrees of freedom.
* A final **Fisher p-value** is calculated from this distribution, representing the combined significance of the protein.



### **4. Integrated Scoring & Ranking**
To highlight biological relevance alongside statistical significance, a final **Score** is calculated:

$$\text{Score} = |\text{Directional Fold Change}| \times -\log_{10}(\text{Fisher } p\text{-value})$$

* **In Expression Mode**: The "Directional Fold Change" is taken directly from the Protein File's $\log_2(\text{Fold Change})$.
* **In Partial Mode**: The direction is determined by the majority vote of the Top N peptides' fold changes.

Proteins are then **Ranked** based on this score, with higher scores indicating top candidates for further investigation.

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
