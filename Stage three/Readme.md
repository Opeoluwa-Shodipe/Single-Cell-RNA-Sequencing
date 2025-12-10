# SARS-CoV-2 Infection Dynamics - scRNA-seq Analysis


**Reproduction of neighborhood clustering, cell type identification, and pseudotime analysis from [PLoS Biology 2021](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001143) using GSE166766 data.**

## 🎯 Project Overview

Single-cell RNA-seq analysis of SARS-CoV-2 infected human bronchial epithelial cells across mock, 1dpi, 2dpi, and 3dpi conditions. Reproduces Figures 1G(i-iii), 3A, 3B, 4A, 4B from the reference paper.

**Data**: [GSE166766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166766) (~100k cells × 33k genes)

## 🧬 Key Findings

| Timepoint | Cells | Dominant Cell Types | ACE2 Status |

|-----------|-------|-------------------|-------------|

| Mock | 22,609| Mixed epithelial | Baseline |

| 1dpi | 11,834| Ciliated ↑ | Constitutive|

| 2dpi | 14,695| Olfactory ↑ | Stable |

| 3dpi | 28,530| **Ciliated peak** | **Highest** |

**Primary targets**: Ciliated/olfactory epithelial cells with peak ACE2 expression at 3dpi

## 🚀 Quick Start (Google Colab)

1. Clone repo
git clone

https://github.com/YOUR_USERNAME/SARS-CoV-2-Infection-Dynamics.git

cd SARS-CoV-2-Infection-Dynamics

2. Install dependencies
pip install scanpy==1.11.5 celltypist decoupler leidenalg

text

undefined

3. Load 10x data (update paths)
import scanpy as sc
mock_adata = sc.read_10x_mtx('data/Mock', prefix='GSM5082289_mock_')

Run notebook sequentially
text

## 🔬 Analysis Pipeline

1. **Data Loading** → 4 conditions (mock/1/2/3 dpi)

2. **Preprocessing** → QC → Normalize → HVGs (top 1000)

3. **Clustering** → PCA → UMAP → Leiden

4. **Annotation** → CellTypist

5. **Trajectory** → Pseudotime analysis

6. **Pathways** → decoupler ULM

**Memory fix for Colab**:

sc.pp.highly_variable_genes(adata, n_top_genes=1000)
dc.mt.ulm(..., bsize=10000)

text

## 📁 Repository Structure

├── SARS-Cov-2-Infection-Dynamics.ipynb # Main analysis
├── SARS-Cov-2-Infection-Dynamics-1.ipynb # Pathway analysis
├── data/ # GSE166766 MTX/TSV files
│ ├── Mock/, 1dpi/, 2dpi/, 3dpi/
└── README.md

text

## 📊 Outputs

- UMAPs (condition/cell type colored)

- Cell annotations (ciliated, olfactory, mast cells)

- ACE2 expression dynamics

- Pseudotime trajectories

- Pathway enrichment scores

## 🛠️ Dependencies

scanpy==1.11.5 # Core analysis
celltypist==1.7.1 # Annotation
decoupler==2.1.2 # Pathways
leidenalg==0.11.0 # Clustering

text

## 🆘 Troubleshooting

| Issue | Fix |

|-------|-----|

| RAM crash | `n_top_genes=1000`, `bsize=5000` |

| Slow CellTypist | First run downloads model |

| decoupler fails | Use `dc.mt.vision()` |

