# Peptide-to-protein data aggregation using Fisher’s method improves target identification in chemical proteomics 
<p align="left">
<sub><b>PeptideFisher</b> is a specialized Shiny application designed to aggregate peptide-level statistics into protein-level significance using Fisher’s Method. It is optimized for proteomics workflows, including <b>Differential Expression/Solubility</b> and <b>Partial Proteolysis (AFDIP/HOLSER)</b>.</sub>
</p>
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
| **Peptide File** | Must contain an **`Accession`**, **`Gene name`**, and **`Protein description`** columns and >2 intensity columns for each treatment. |
| **Protein File** | Must contain an **`Accession`** column to match with the peptide data, and >2 intensity columns for each treatment. |

### 2. **Partial Proteolysis (AFDIP / HOLSER) Mode**
Use this for structural change analysis.

| Requirement | Description |
| :--- | :--- |
| **Peptide File** | Must contain **`Accession`**, **`Gene name`**, and **`Protein description`** columns, and >2 intensity columns for each treatment.. |
| **Protein File** | **Not required** for this mode (directionality is calculated from peptides). |

* **`Accession`**: Unique protein identifier.
* **`Gene name`**: Official gene symbol.
* **`Protein description`**: Full name or description of the protein.

> [!TIP]
> **Why 2+ replicates?** A minimum of 2 samples per group is statistically required to calculate the standard deviation for the Welch's T-test. If you select only one column, the app will display a validation error.

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

* **In Expression/Solubility Mode**: The "Directional Fold Change" is taken directly from the Protein File's $\log_2(\text{Fold Change})$.
* **In Partial Proteolysis Mode (AFDIP / HOLSER)**:
    Since no Protein File is provided, the value is synthesized from the **Top N peptides**:
    1.  **Direction (+/-)**: Determined by a **majority vote** of the Top N peptides' $\log_2(\text{Fold Change})$ signs. For example, if 3 peptides are positive and 1 is negative, the direction is positive ($+1$).
    2.  **Magnitude (Value)**: Calculated as the **mean of the absolute $\log_2(\text{Fold Change})$ values** of those Top N peptides.
    3.  **Result**: The directional fold change = $(\text{Direction}) \times (\text{Mean Absolute } \log_2\text{FC})$.

Proteins are then **Ranked** based on this score, with higher scores indicating top candidates for further investigation.

---
## Troubleshooting

If you encounter issues while running the analysis, please check the following common scenarios:

### **1. Validation Errors (Red Messages)**
* **"Insufficient replicates"**: Ensure you have selected **at least 2 columns** for both Control and Drug groups in the dropdown menus. T-tests require variance, which cannot be calculated from a single sample ($n=1$).
* **"Required columns missing"**: Double-check your file headers. The app looks for exact matches for **`Accession`**, **`Gene name`**, and **`Protein description`**.
* **"File size limit"**: The default upload limit is **100MB**. If your file is larger, please contact the administrator to adjust the `shiny.maxRequestSize` option.

### **2. Results Issues**
* **All p-values are NA**: This happens if the intensity values across replicates are identical (zero variance) or if there are too many missing values for the t-test to perform.
* **Progress bar is slow**: For datasets with >50,000 peptides, the "Computing peptide statistics" stage involves thousands of t-tests and may take 30–60 seconds. Please be patient.

---

## 📄 Contact & Citation

### **How to Cite**
If you use this tool in your research, please cite it as follows:
> *Peptide-to-protein data aggregation using Fisher’s method improves target identification in chemical proteomics. (2026). [H Lyu, et. al], [Karolinska Institutet].*

### **Contribution & Feedback**
We welcome contributions to improve **FisherPep**!
* **Bug Reports**: Please open an **Issue** on GitHub with a description of the error and a small sample of the data causing it.
* **Feature Requests**: Feel free to suggest new scoring methods or visualization features.

---

