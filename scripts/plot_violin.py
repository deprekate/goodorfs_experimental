import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


df = pd.read_csv(sys.argv[1],sep=' ',keep_default_na=False)
df.columns = ['ID','STOP','TYPE','COUNT']

sns.set(style="whitegrid", palette="pastel", color_codes=True)

fig, ax = plt.subplots(figsize=(8, 8))

ax = sns.violinplot(x="ID", y="COUNT", hue="TYPE", data=df, 
		split=True, 
		cut=0, 
		palette={True: "b", False: "y"},
		linewidth=0.1
		)
ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
#sns.despine(left=True)

fig.suptitle('Title', fontsize=18, fontweight='bold')
ax.set_xlabel("Genome",size = 16,alpha=0.7)
ax.set_ylabel("% in 1st cluster",size = 16,alpha=0.7)
plt.legend(loc='upper left')
plt.show()
