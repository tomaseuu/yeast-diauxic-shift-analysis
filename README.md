# Yeast Diauxic Shift Gene Expression Analysis

This project explores how gene expression in *Saccharomyces cerevisiae* (baker’s yeast) changes during the **diauxic shift**—a critical metabolic transition from glucose fermentation to ethanol respiration. Using a dataset from a 1997 study by Joseph DeRisi, I applied clustering and statistical analysis techniques to identify key genes involved in this switch and evaluated their biological significance.

## Dataset
- Raw gene expression: `data/diauxic_raw_ratios.txt`
- Authors’ top variable genes: `data/230genes_log_expression.txt`

## Objectives
- Visualize gene expression changes across 7 time points.
- Cluster genes with similar expression patterns using K-Means.
- Filter out low-variability genes to improve clarity.
- Compare custom-selected genes to the authors’ list.
- Perform GO (Gene Ontology) enrichment analysis.
- Simulate and calculate statistical significance using binomial models.

## Tools & Libraries

- **Python 3**
- `pandas`, `numpy` – data handling and math
- `matplotlib`, `seaborn` – visualizations
- `scikit-learn` – clustering and silhouette score
- `scipy.stats` – binomial probability analysis

## Key Findings

- Clear gene expression clusters appear at time points R6 and R7, aligning with the diauxic shift.
- 129 of the 230 most variable genes overlapped with the authors’ list (Jaccard = 0.396).
- GO enrichment revealed significant overrepresentation of metabolic processes like citrate metabolism and glucose import.
- Statistical simulations confirmed the robustness of enrichment results.

## Project Report

For a full breakdown of methods, results, figures, and discussion, see:  
📄 [`yeast_diauxic_project.pdf`](https://docs.google.com/document/d/1xNxNMKbp5Wmi72WcWUPU86XW0EFWy6G59CQSiiGvCJA/edit?usp=sharing)

## Author

**Thomas Le**  
Student in Computer Science & Bioinformatics
