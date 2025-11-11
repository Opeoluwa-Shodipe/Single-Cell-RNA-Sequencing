import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

bc_url = 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv'

bc = pd.read_csv(bc_url)
sns.kdeplot(data=bc, x='area_mean', hue='diagnosis', common_norm=False, fill=True, palette={'M': 'orange', 'B': 'blue'}, alpha=0.5)
plt.xlabel('Area Mean')
plt.ylabel('Density')
plt.legend(title='Diagnosis')
plt.title('f. Area Mean Density by Diagnosis')
plt.show()
