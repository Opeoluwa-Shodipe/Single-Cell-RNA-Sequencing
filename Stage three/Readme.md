# SARS-CoV-2 Infection Dynamics - Single-Cell RNA-seq Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)\](https://opensource.org/licenses/MIT)

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)\](https://www.python.org/downloads/)

[![Scanpy](https://img.shields.io/badge/scanpy-1.11.5-orange.svg)\](https://scanpy.readthedocs.io/)

[![Colab](https://img.shields.io/badge/Google%20Colab-Ready-green.svg)\](https://colab.research.google.com/)

**Single-cell RNA sequencing analysis of SARS-CoV-2 infection dynamics across mock, 1dpi, 2dpi, and 3dpi conditions using 10x Genomics datasets (~100k cells × 33k genes).**

## 🎯 Project Summary

Comprehensive scRNA-seq pipeline analyzing SARS-CoV-2 infection progression in respiratory epithelial cells. Identifies temporal cell type dynamics, ACE2 expression patterns, and infection-responsive pathways.

### 🧬 Key Biological Insights

- **Primary targets**: Ciliated/olfactory epithelial cells (peak ACE2 at 3dpi)

- **Immune dynamics**: Mast cell activation correlates with cytokine signaling

- **Neuroinvasion**: ENO2 upregulation in glia populations

- **Pathway enrichment**: Interferon response peaks at 2-3dpi

## 📁 Repository Contents

SARS-CoV-2-Infection-Dynamics/
├── SARS-Cov-2-Infection-Dynamics.ipynb # Main analysis (scanpy + CellTypist)
├── SARS-Cov-2-Infection-Dynamics-1.ipynb # Decoupler pathway analysis
├── README.md # 📄 This file
├── data/ # 10x Genomics (GSM5082289-92)
│ ├── Mock/
│ ├── 1dpi/
│ ├── 2dpi/
│ └── 3dpi/
└── outputs/ # UMAPs, annotations, scores

text

## 🚀 Quick Start (Google Colab)

### 1. Clone Repository

git clone

https://github.com/YOUR_USERNAME/SARS-CoV-2-Infection-Dynamics.git

cd SARS-CoV-2-Infection-Dynamics

text

### 2. Open in Colab

Runtime → Change runtime type → GPU (recommended)
!pip install scanpy==1.11.5 celltypist decoupler leidenalg igraph fa2-modified

text

### 3. Mount Data & Run

from google.colab import drive
drive.mount('/content/drive')

Load datasets (update paths to your 10x folders)
mock_adata = sc.read_10x_mtx('/content/drive/MyDrive/HackBio/Mock',
prefix='GSM5082289_mock_')

Run notebook cells sequentially
text

## 🔬 Complete Analysis Pipeline

🔄 Data Loading (4 conditions: mock/1/2/3 dpi)
→ 22k-28k cells × 33k genes each
🧹 Preprocessing
→ QC → Normalize → Log1p → HVGs (top 1000-2000)
📍 Dimensionality Reduction
→ PCA → Neighbors → UMAP → Leiden clustering
🏷️ Annotation
→ CellTypist (automated cell type classification)
🎯 Pathway Analysis
→ decoupler ULM (markers network, bsize=10000)
📊 Visualization
→ Condition-colored UMAPs + feature plots (ACE2, ENO2)
text

## 📊 Results Dashboard

| Timepoint | Cells (n) | Dominant Types | ACE2 Status | Key Pathways |

|-----------|-----------|----------------|-------------|--------------|

| **Mock** | 22,609 | Mixed epithelial | Baseline | Housekeeping |

| **1dpi** | 11,834 | Ciliated ↑ | Constitutive | Early ISG |

| **2dpi** | 14,695 | Olfactory ↑ | Stable | Interferon |

| **3dpi** | 28,530 | Mast ↑, Ciliated peak | **Highest** | Cytokine storm |

## ⚙️ Technical Specifications

### Dependencies

scanpy==1.11.5 # scRNA-seq core
celltypist==1.7.1 # Cell type annotation
decoupler==2.1.2 # Pathway enrichment
leidenalg==0.11.0 # Clustering
anndata==0.12.6 # Data format

text

### Memory Optimization

For free Colab (12GB RAM):
sc.pp.highly_variable_genes(adata, n_top_genes=1000)
dc.mt.ulm(..., bsize=10000) # Critical for ULM

Colab Pro (52GB): Full 33k genes possible
text

## 🔍 Key Outputs Generated

1. **UMAP plots**: Condition + cell type colored

2. **Cell annotations**: Hepatocytes, ciliated, olfactory, mast cells

3. **ACE2 violin plots**: Temporal expression dynamics

4. **Pathway scores**: ULM z-scores (interferon, cytokine)

5. **Marker networks**: decoupler results

## 📝 Citation

Shodipe, O. (2025). SARS-CoV-2 Infection Dynamics:
Single-Cell RNA-seq Analysis of Temporal Cell Type Changes.
GitHub:

https://github.com/YOUR_USERNAME/SARS-CoV-2-Infection-Dynamics

Data: GEO GSE165063 (10x Genomics)

text

## 🛠️ Troubleshooting

| Issue | Solution |

|-------|----------|

| **RAM Crash** | `n_top_genes=1000`, `bsize=5000` |

| **CellTypist slow** | First run downloads model (~500MB) |

| **UMAP poor** | `n_neighbors=30`, restart PCA |

| **decoupler fails** | Try `dc.mt.vision()` alternative |

## 📄 License

**MIT License** - Free for research/education/commercial use.

Copyright (c) 2025 Opeoluwa Shodipe

Permission is hereby granted, free of charge, to any person obtaining a copy...

text

## 🙌 Acknowledgments

- **Scanpy developers** - scRNA-seq gold standard

- **CellTypist team** - Automated annotation

- **Decoupler/SaezLab** - Pathway inference

- **10x Genomics** - Public SARS-CoV-2 datasets

---

**🔬 Reproducible bioinformatics research made easy**

*Last updated: December 10, 2025*
