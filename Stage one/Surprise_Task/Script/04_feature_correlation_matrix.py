"""
Hackbio Internship - Step 1 Surprise Task
Team: Glycine

# Task: Reproduce Analysis
# Author: Opeoluwa Shodipe
# Github: https://github.com/Opeoluwa-Shodipe/Single-Cell-RNA-Sequencing/tree/main/Stage%20one
# Linkedin: https://www.linkedin.com/in/sopeoluwa/
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

bc_url = 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv'

bc = pd.read_csv(bc_url)
features = ['radius_mean', 'texture_mean', 'perimeter_mean', 'area_mean', 'smoothness_mean', 'compactness_mean']
corr = bc[features].corr()
sns.heatmap(corr, annot=True, cmap="Blues", fmt=".1f", linewidths=0.5)
plt.title('d. Correlation Matrix')
plt.show()
