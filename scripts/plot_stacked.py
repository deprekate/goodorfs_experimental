import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
#import seaborn as sns

df = pd.read_csv("../se_km_3_len.tsv", sep=' ', header=None)
df.columns = ['ID','min','mid','max','TP','FP','FN','TN']
df['FN'] = -1 * df['FN']
df['TN' ]= -1 * df['TN']
df = df.drop(columns=['min','mid','max'])
df = df.sort_values(by=['TP'])
dat1 = df.iloc[0:10,:]
print(dat1)
exit()
d = pd.melt(dat, id_vars=['ID'], value_vars=['TP','FP','FN','TN'])
#sns.set()

#fig = plt.figure(figsize=(8, 6))
ax1 = plt.subplot2grid((1, 3), (0, 0), colspan=1)
dat.plot(x='ID', kind='bar', stacked=True, ax=ax1)

ax2 = plt.subplot2grid((1, 3), (0, 1), colspan=2)
dat.plot(x='ID', kind='bar', stacked=True, ax=ax2)

plt.show()
