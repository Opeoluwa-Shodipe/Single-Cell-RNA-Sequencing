import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

bc_url = 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv'

bc = pd.read_csv(bc_url)
sns.scatterplot(data=bc, x='radius_mean', y='texture_mean', hue='diagnosis', palette={'M': 'orange', 'B': 'blue'}, alpha=0.7)
plt.title('c. Texture vs Radius by Diagnosis')
plt.show()
