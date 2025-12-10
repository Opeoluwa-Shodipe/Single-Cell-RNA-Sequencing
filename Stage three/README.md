# SARS-CoV-2 Infection Dynamics — Single-Cell RNA-seq Analysis

**Reference paper:** Ravindra *et al.*, *PLOS Biology* (2021) — *Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium*  
https://doi.org/10.1371/journal.pbio.3001143

**Data (10x mtx / tsv):** GEO GSE166766  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166766

---

## 🎯 Project Goal
Reproduce neighbourhood clustering and cell‑type identification from Ravindra *et al.* (PLOS Biology) and perform **pseudotime / trajectory analysis** to order cell state transitions across infection (mock, 1 dpi, 2 dpi, 3 dpi). Specifically target reproducing Figures: **1G(i–iii), 3A, 3B, 4A, 4B** from the paper.

---

## 🔬 Biological summary (from Ravindra et al.)
- The study identifies **eight epithelial cell types** in differentiated human bronchial epithelial cells (HBECs): **ciliated, basal, club, BC/club (basal→club intermediate), goblet, neuroendocrine, ionocytes, and tuft cells**. (See Fig 3A/B in Ravindra et al.). cite: PLOS Biology.  
- **Ciliated cells** are the primary early target of SARS‑CoV‑2 infection (validated by microscopy and reporter virus). Infection expands later to **basal and club**/BC‑club populations. cite: PLOS Biology.  
- **ACE2** expression is enriched at the *cell‑type level* in susceptible populations (ciliated, basal, club, BC/club) but **correlates poorly with infection on a per‑cell basis** (Spearman’s rho ≈ −0.06 in ciliated cells). ACE2 behaves as a *cell-type risk marker*, not a perfect per-cell infection tracer. cite: PLOS Biology.  

(Primary source: Ravindra *et al.* PLOS Biology 2021).  
Paper: https://doi.org/10.1371/journal.pbio.3001143  
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166766

---

## ✅ Reproduction checklist
1. Download processed 10x `mtx` and `tsv` files from GEO (GSE166766).  
2. Load each condition as `AnnData` (mock, 1dpi, 2dpi, 3dpi).  
3. Preprocess: QC → Normalize (CPM) → `log1p` → identify HVGs (1k–2k).  
4. Correct batch effects using **BB-kNN** (paper used BB‑kNN) and generate UMAP + PHATE.  
5. Cluster using **Louvain** (high resolution), merge clusters using canonical bronchial markers to match the 8 reported types.  
6. Validate tropism by mapping SARS‑CoV‑2 reads to viral genome and flagging infected cells (paper: ≥10 viral transcript counts).  
7. Run **trajectory / pseudotime** (PHATE + diffusion pseudotime `scanpy.tl.dpt` or Palantir) and visualize gene dynamics.  
8. Reproduce Figures: 1G, 3A/B, 4A/B (UMAP + expression heatmaps / violin plots).  

---

## 🧪 Key code snippets (Scanpy + BB-kNN + PHATE + DPT)
```python
# basics
import scanpy as sc
import bbknn
import phate
import anndata
import pandas as pd

# 1) Load 10x (example for mock)
adata = sc.read_10x_mtx("data/Mock/", var_names='gene_symbols', cache=True)

# 2) Basic QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 15, :]

# 3) Normalize & HVG
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')
adata = adata[:, adata.var['highly_variable']]

# 4) PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', n_comps=100)

# 5) Build BB-kNN graph (recommended per paper)
# If multiple conditions, concatenate AnnData objects first, then:
bbknn.bbknn(adata, batch_key='condition', neighbors_within_batch=3, n_pcs=50)

# 6) UMAP + Leiden/Louvain
sc.tl.umap(adata, min_dist=0.5)
sc.tl.louvain(adata, resolution=1.0)  # paper used Louvain / high resolution

# 7) PHATE (optional, paper used PHATE for trajectories)
ph = phate.PHATE(n_components=2)
adata.obsm['X_phate'] = ph.fit_transform(adata.obsm['X_pca'])

# 8) Pseudotime: diffusion pseudotime (DPT)
sc.tl.diffmap(adata)
sc.tl.dpt(adata, n_dcs=10)
# visualize pseudotime on UMAP/PHATE
sc.pl.umap(adata, color=['dpt_pseudotime','ACE2','cluster_name'], cmap='viridis')

# 9) Map viral reads (if raw counts include viral transcripts):
# infected = adata.obs['n_viral_counts'] >= 10
# adata.obs['infected'] = infected
```

---

## 🔬 Reproduction targets (figures)
- **Figure 1G(i–iii)**: UMAP / PHATE showing infection scores across time and clusters.  
- **Figure 3A, 3B**: UMAP colored by manual cell-type annotation (8 epithelial types) + violin plots of canonical markers.  
- **Figure 4A, 4B**: ACE2 expression overlay and heatmaps comparing susceptible vs non‑susceptible cell types across time.

Use BB-kNN + UMAP/PHATE + Louvain → manual annotation with canonical markers (e.g., *FOXJ1* for ciliated, *KRT5* for basal, *SCGB1A1* for club, *MUC5AC* for goblet, *PAX6/ENO2* for neuroendocrine/neuronal-like cells, *FOXI1* for ionocytes, *POU2F3* for tuft).

---

## ❓ Answers to the Project Questions (concise, evidence-backed)

### 1) **What cell types did you identify at different stages of infection?**
**Answer (from Ravindra et al.):** Eight epithelial cell types were identified across time: **ciliated, basal, club, BC/club (basal→club intermediate), goblet, neuroendocrine, ionocytes, and tuft cells**. Early infection (1 dpi) is concentrated in **ciliated** cells, while later time points (2–3 dpi) show expanded infection into **basal, club, and BC/club** populations. cite: PLOS Biology. citeturn0view0

### 2) **Why do these cell types correlate with COVID‑19 infection?**
**Short answer:** These cell types express the host entry factors and are located on the airway surface where virus exposure occurs; in particular **ciliated cells** present accessible apical surfaces and express ACE2/TMPRSS2 proteases at levels sufficient to permit entry, explaining early tropism. Basal/club populations become infected later as infection spreads or local microenvironments change (e.g., cytokine/IFN responses, cell-state transitions). Mechanistically, ACE2/TMPRSS2 (and related proteases) and cell surface accessibility determine tropism. cite: PLOS Biology. citeturn1view2

### 3) **Is ACE2 a good marker for tracking infection rate (based on this dataset)?**
**Answer:** **Partially.** ACE2 is useful at the **cell‑type level** (it is enriched in cell types that are more likely to be infected: ciliated, basal, club, BC/club). However, **on a per‑cell basis ACE2 poorly correlates with viral transcript load** (Spearman’s rho ~ −0.06 in ciliated cells reported by Ravindra *et al.*). ACE2 can highlight susceptible populations but **cannot reliably predict which individual cells are infected** in this dataset. cite: PLOS Biology. citeturn1view2

### 4) **What is the difference between ENO2 and ACE2 as biomarkers in the two studies?**
**Interpretation & evidence:**
- **ACE2** is the **viral entry receptor** — a mechanistic marker of cell susceptibility to SARS‑CoV‑2. It’s a functional surface receptor and used to explain which cell types are permissive. (Ravindra *et al.* show ACE2 enrichment in susceptible cell types but poor per‑cell correlation with viral load.) cite: PLOS Biology. citeturn1view2  
- **ENO2 (Neuron‑specific enolase)** is a **metabolic / neuronal lineage marker** (glycolytic enzyme commonly used as a neuronal/neuroendocrine marker). ENO2 upregulation in some single‑cell studies has been linked to **cellular stress, metabolic rewiring, or neuroendocrine‑like states** and can correlate with disease severity or marker expression in neuroendocrine-like epithelial cells — but it is **not** an entry receptor. ENO2 therefore reports **cell state or damage/metabolic changes**, not entry competence. (See studies showing elevated NSE/ENO2 in severe COVID‑19 and ENO2's role as a stress/neuronal marker). cite: PLOS Biology + literature on ENO2. citeturn3search11turn3search17

### 5) **Which cell cluster has the highest abundance of ACE2 expression after 3 dpi and what does that mean biologically (interpret visually)?**
**Answer (paper + interpretation):**  
Ravindra *et al.* report ACE2 enrichment in the **ciliated** cluster (and in basal/club/BC‑club) relative to non‑susceptible types. Visualizing ACE2 on a UMAP/PHATE (as in Fig 4A–4B) shows a **localized high ACE2 signal overlapping the ciliated cluster**, especially in infected conditions. Biologically, this means **ciliated epithelial cells are the primary portal of entry and early viral replication**, consistent with microscopy (immunofluorescence) showing viral protein in FOXJ1+ ciliated cells (Fig 3F). After 3 dpi, high ACE2 in ciliated (and some basal/club) clusters indicates sustained susceptibility and may reflect either intrinsic higher receptor expression in those lineages or state changes (e.g., ISG/IFN effects) that alter ACE2 levels or the sampling of infected/bystander cells. Visual interpretation: on UMAP/PHATE, expect bright ACE2 signal 'hotspots' overlapping the annotated ciliated cluster and co-localization with infected cell markers/viral read counts. cite: PLOS Biology. citeturn1view2

---

## 📦 Deliverables for submission
- `README.md` (this file) — project overview & reproduction steps.  
- `notebook_repro_clustering.ipynb` — code to: load GEO mtx, BB‑kNN, UMAP/PHATE, Louvain, manual annotation.  
- `notebook_pseudotime.ipynb` — DPT/Palantir or scVelo RNA‑velocity pipeline + pseudotime.  
- `figures/` — reproduction of Figures 1G, 3A/B, 4A/B (UMAPs, violin/heatmaps).  
- `results/` — table with per-cluster ACE2 average expression and infected cell counts (infection threshold = 10 viral transcripts as in Ravindra *et al.*).

---

## 🔁 Notes & tips
- **Viral transcript mapping:** ensure viral genome is included in the reference before quantifying reads (paper counted cells with ≥10 viral transcript counts as infected).  
- **Batch correction:** use BB‑kNN exactly as paper to match neighborhood structure.  
- **Pseudotime choices:** PHATE + DPT or Palantir are recommended; scVelo (RNA velocity) adds directionality if spliced/unspliced counts exist.  
- **ACE2 caveat:** ACE2 is an ISG in some contexts; interpret ACE2 changes carefully (cell type vs cell state). cite: PLOS Biology. citeturn1view2

---

## 📚 References
- Ravindra NG, Alfajaro MM, Gasque V, et al. *Single-cell longitudinal analysis of SARS‑CoV‑2 infection in human airway epithelium identifies target cells, alterations in gene expression, and cell state changes.* PLOS Biology (2021). https://doi.org/10.1371/journal.pbio.3001143. citeturn0view0  
- GEO accession: GSE166766 — primary dataset (10x mtx/tsv). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166766. citeturn0view1  
- ENO2 / NSE literature (marker & COVID context): Cione *et al.* 2021; Moreno Jr *et al.* 2022 (metabolic remodelling). citeturn3search11turn3search17

---

**If you want:** I can now (pick one)  
1. Generate the two analysis notebooks (`notebook_repro_clustering.ipynb` and `notebook_pseudotime.ipynb`) and save them in the project folder, or  
2. Produce the exact `README.md` file on disk for download (I already created a draft; I can overwrite it with this improved version), or  
3. Run the reproduction steps on the uploaded notebook you provided and produce figures and a short results summary.

Tell me which of the three you'd like me to do next.
