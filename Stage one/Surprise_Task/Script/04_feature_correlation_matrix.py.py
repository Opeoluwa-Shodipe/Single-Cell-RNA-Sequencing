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
