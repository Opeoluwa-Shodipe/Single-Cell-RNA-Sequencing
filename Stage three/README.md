# SARS-CoV-2 Infection Dynamics — Single-Cell RNA-seq Analysis

Single-cell RNA sequencing analysis of SARS-CoV-2 infection dynamics across mock, 1dpi, 2dpi, and 3dpi conditions using 10x Genomics datasets (~100k cells × 33k genes).

## 🎯 Project Summary
Comprehensive scRNA-seq pipeline analyzing SARS-CoV-2 infection progression in respiratory epithelial cells. Identifies temporal cell type dynamics, ACE2 expression patterns, and infection-responsive pathways.

## 🧬 Key Biological Insights
- Primary targets: Ciliated and olfactory epithelial cells (peak ACE2 at 3dpi)
- Immune activation: Mast cell inflammatory signatures increase across dpi
- Neuroinvasion markers: ENO2 upregulation in glial-like populations
- Pathway enrichment: Interferon signaling peaks at 2–3 dpi

## 📁 Repository Structure
```
SARS-CoV-2-Infection-Dynamics/
├── SARS-Cov-2-Infection-Dynamics.ipynb
├── SARS-Cov-2-Infection-Dynamics-1.ipynb
├── README.md
├── data/
│   ├── Mock/
│   ├── 1dpi/
│   ├── 2dpi/
│   └── 3dpi/
└── outputs/
```

## 🚀 Quick Start (Google Colab)

### 1. Clone Repository
```bash
git clone https://github.com/YOUR_USERNAME/SARS-CoV-2-Infection-Dynamics.git
cd SARS-CoV-2-Infection-Dynamics
```

### 2. Install Dependencies
```python
!pip install scanpy==1.11.5 celltypist decoupler leidenalg igraph fa2-modified
```

### 3. Load Data
```python
from google.colab import drive
drive.mount('/content/drive')

mock_adata = sc.read_10x_mtx(
    '/content/drive/MyDrive/HackBio/Mock',
    prefix='GSM5082289_mock_'
)
```

## 🔬 Complete Analysis Pipeline
1. Data Loading (mock, 1dpi, 2dpi, 3dpi)
2. Preprocessing: QC → Normalize → Log1p → HVGs (1k–2k)
3. Dimensionality Reduction: PCA → Neighbors → UMAP → Leiden
4. Cell Typing using CellTypist
5. Pathway Analysis using decoupler ULM
6. Visualization: UMAPs, ACE2/ENO2 expression

## 📊 Results Overview
| Timepoint | Cells | Dominant Types | ACE2 | Key Pathways |
|----------|--------|----------------|------|--------------|
| Mock | 22,609 | Mixed epithelial | Baseline | Housekeeping |
| 1dpi | 11,834 | Ciliated ↑ | Constitutive | Early ISG |
| 2dpi | 14,695 | Olfactory ↑ | Stable | Interferon |
| 3dpi | 28,530 | Mast ↑, Ciliated peak | Highest | Cytokine storm |

## ⚙️ Technical Specifications

### Dependencies
```
scanpy==1.11.5
celltypist==1.7.1
decoupler==2.1.2
leidenalg==0.11.0
anndata==0.12.6
```

### Memory Tips
```python
sc.pp.highly_variable_genes(adata, n_top_genes=1000)
dc.mt.ulm(..., bsize=10000)
```

## 📝 Citation
Shodipe, O. (2025). *SARS-CoV-2 Infection Dynamics: Single-Cell RNA-seq Analysis of Temporal Cell Type Changes.*  
GitHub: https://github.com/YOUR_USERNAME/SARS-CoV-2-Infection-Dynamics  
Data: GEO GSE165063

## 🛠️ Troubleshooting
| Issue | Solution |
|-------|----------|
| RAM crash | Use n_top_genes=1000, bsize=5000–10000 |
| CellTypist slow | First run downloads 500MB model |
| Poor UMAP | Increase n_neighbors=30 |
| decoupler errors | Try dc.mt.vision() |

## 📄 License (MIT)
Permission is hereby granted...

## 🙌 Acknowledgments
- Scanpy team  
- CellTypist  
- SaezLab / decoupler  
- 10x Genomics  
