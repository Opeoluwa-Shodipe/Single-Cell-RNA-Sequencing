import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

deg_chr22_url = 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv'

deg = pd.read_csv(deg_chr22_url)
colors = {'up': 'green', 'down': 'orange', 'ns': 'grey'}
plt.figure(figsize=(6, 5))
plt.scatter(deg['log2FoldChange'], -np.log10(deg['padj']), c=deg['significance'].map(colors), alpha=0.7)
plt.axvline(1, ls='--', c='grey')
plt.axvline(-1, ls='--', c='grey')
plt.xlabel('log2FoldChange')
plt.ylabel('-log10Padj')
legend_labels = [plt.Line2D([0], [0], marker='o', color='w', label='Upregulated', markerfacecolor='green', markersize=7),
                 plt.Line2D([0], [0], marker='o', color='w', label='Downregulated', markerfacecolor='orange', markersize=7),
                 plt.Line2D([0], [0], marker='o', color='w', label='Not Significant', markerfacecolor='grey', markersize=7)]
plt.legend(handles=legend_labels)
plt.title("b. Volcano Plot")
plt.show()
