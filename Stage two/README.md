# Single-Cell RNA-Seq Analysis: Bone Marrow Mononuclear Cells

## Overview
This repository contains a reproduction of a core single-cell RNA-seq analysis pipeline applied to a human bone marrow dataset. The workflow transforms raw count data into biologically meaningful insights, identifying distinct immune cell populations and assessing the physiological state of the tissue.

---

## Key Findings (Biological Interpretation)

### 1. Identified Cell Types
Based on the Leiden clustering (resolution 0.02) and marker annotation, we identified three major immune lineages:
* **Monocytes (Cluster 0):** The dominant population.
* **B Cells (Cluster 1):** Distinct lymphoid population.
* **T Cells (Cluster 2):** Distinct lymphoid population.

### 2. Biological Roles
* **Monocytes:** Innate immune sentinels that circulate in the blood and marrow; they differentiate into macrophages (phagocytosis) or dendritic cells (antigen presentation) upon tissue entry.
* **B cells:** The arm of humoral immunity; they originate in the marrow and differentiate into plasma cells to produce antibodies that neutralize pathogens.
* **T cells:** The effectors of cell-mediated immunity; they kill infected host cells (CD8+) or coordinate the broader immune response (CD4+).

### 3. Tissue Source Assessment: **BMMC (Ficoll-Separated)**
The sample is likely **Bone Marrow Mononuclear Cells (BMMCs)**, not whole marrow. This is evidenced by the complete absence of **Granulocytes** and **Erythroid** lineages, which are typically filtered out during density gradient (Ficoll) separation.

### 4. Patient Status: **Infected / Inflammatory State**
The profile suggests an **Active Infection** due to a significant shift in cell proportions:
* **Monocytosis:** Monocytes are the largest cluster, inverting the normal ratio where T cells should dominate.
* **Lymphopenia:** Relative depletion of B and T cell compartments compared to the myeloid (monocyte) expansion.
* **Emergency Myelopoiesis:** Evidence of Low-Density Granulocytes (LDGs) in the Monocyte cluster, a signature of acute inflammation.

---

## Analysis Pipeline (Methodology)
The analysis was performed using **Scanpy** and follows these steps:
1.  **Data Ingestion:** Loading `bone_marrow.h5ad` data.
2.  **Quality Control:** Filtering cells based on mitochondrial content, ribosomal content, and library size.
3.  **Preprocessing:** Doublet removal, Normalization, and Log1p transformation.
4.  **Dimensionality Reduction:** Feature selection (HVGs), PCA, and UMAP projection.
5.  **Clustering:** Graph-based clustering (Leiden algorithm, res=0.02).
6.  **Annotation:** Automated cell type annotation using `decoupler`.

---

## How to Run
### Prerequisites
Ensure you have Python installed with the following libraries:
```bash
pip install scanpy anndata igraph celltypist decoupler scrublet
