# Peptide-to-protein data aggregation using Fisher’s method improves target identification in chemical proteomics 
R implementation of peptide-level statistical testing in proteomics using Fisher’s method to combine peptide p-values. This approach avoids biases from deviant or missing peptides and improves detection of regulated or shifted proteins compared with traditional protein-level analyses.

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
