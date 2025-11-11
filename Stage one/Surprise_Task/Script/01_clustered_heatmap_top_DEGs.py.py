import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

gene_expression_url = 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv'

expr = pd.read_csv(gene_expression_url, index_col=0)
sns.clustermap(expr, cmap="Blues", metric='euclidean', linewidths=0.5, annot=False)
plt.title("a. Top DEG Clustered Heatmap")
plt.show()
