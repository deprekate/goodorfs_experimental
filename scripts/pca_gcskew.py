import pandas as pd
import glob, os
import functools
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn import decomposition
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#from plotnine import *
#from matplotlib.mlab import PCA

# LOAD THE DATA
#all_files = glob.glob("../data/aa/*")

#dat = pd.read_csv(path + "/data/aa/" + args.genome_id + ".tsv.gz", compression='gzip', header=0,sep='\t' )
df = pd.concat([pd.read_csv(fp, header=0, sep='\t').assign(ID=os.path.basename(fp)) for fp in glob.glob("../data/aa/*")])

df.loc[:,'A':'V'] = df.loc[:,'A':'V'].div(df.loc[:,'A':'V'].sum(axis=1), axis=0)

dat = df[df.TYPE].groupby(['ID','TYPE']).mean()



# STRIP OUT ONLY THE COUNTS
dat_norm = dat[list("ARNDCEQGHILKMFPSTWYV")]

# DIVIDE BY THE ROW SUMS TO GET THE FREQUENCY

# SCALE THE VALUES
dat_scaled = StandardScaler().fit_transform(dat_norm)

# CLASSIFY EACH ROW USING KMEANS
#clust = KMeans(n_clusters=2).fit(dat_scaled).labels_
# CALCULATE THE PRINCIPLE COMPONENTS
pca = decomposition.PCA(n_components = 2, svd_solver='full').fit(dat_scaled)
dat_pca = pca.transform(dat_scaled)

x_vector = pca.components_[0]
y_vector = pca.components_[1]
# PLOT
colors = {'noncoding':'#F2766E', 'coding':'#3CC9CF', True:'#3CC9CF', False:'#F2766E'}
df = pd.DataFrame({'X':dat_pca[:,0],'Y':dat_pca[:,1],'TYPE':data.TYPE})

fig, ax = plt.subplots()
ax.scatter(df['X'], df['Y'], c=df['TYPE'].apply(lambda x: colors[x]), marker='.', linewidths=0.0, alpha=0.1, zorder=5)
for i in range(len(x_vector)):
    x = (1.2*x_vector[i]*max(dat_pca[:,0]))
    y = (1.2*y_vector[i]*max(dat_pca[:,0]))
    plt.arrow(0, 0, x, y,  color='black', width=0.00005, zorder=10)
    plt.text(x*1.1, y*1.1, dat.columns[i], color='black', zorder=10)
print("done")
blue_patch = mpatches.Patch(color='#3CC9CF', label='coding')
pink_patch = mpatches.Patch(color='#F2766E', label='non-coding')
# LEGEND
plt.legend(handles=[blue_patch,pink_patch])
ax.set_title('amino-acid frequency of potential ORFs from Lambda phage')
ax.set(xlabel='PC1', ylabel='PC2')
#plt.show()
fig.set_size_inches(20, 10)
fig.savefig('test.png', dpi=100)

