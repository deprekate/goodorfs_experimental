import os
import sys
import time
import argparse
import numpy as np
import pandas as pd 
import gzip
from math import sqrt
from statistics import mean
from math import log10
from math import log2

#from matplotlib.mlab import PCA
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn.preprocessing import StandardScaler

from sklearn import decomposition
from sklearn import datasets
from sklearn.mixture import GaussianMixture
from sklearn.mixture import BayesianGaussianMixture
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn_extra.cluster import KMedoids
from sklearn.cluster import OPTICS
from sklearn.cluster import SpectralClustering
from sklearn.cluster import Birch
from scipy import stats

import warnings
warnings.filterwarnings('ignore')

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
                      width=ell_radius_x * 2,
                      height=ell_radius_y * 2,
                      facecolor=facecolor,
                      **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = matplotlib.transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def argmed(a):
    if(len(a) == 2):
	    return -1
    if len(a) % 2 == 1:
        return np.where( a == np.median(a) )[0][0]
    else:
        l,r = len(a)/2 -1, len(a)/2
        left = np.partition(a, l)[l]
        right = np.partition(a, r)[r]
        return [np.where(a == left)[0][0], np.where(a==right)[0][0]]
setattr(np, 'argmed', argmed)

def xx(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

def xxx(data, axis=None):
	#print( np.median(data - np.mean(data, axis), axis) )
	return np.absolute(np.median(data - np.mean(data, axis), axis))
setattr(np, 'xxx', xxx)

def sums(data, axis=None):
	return np.sum(np.absolute(data - np.mean(data, axis)), axis)
setattr(np, 'sums', sums)


def aad(data, axis=None):
	return np.mean(np.absolute(data - np.mean(data, axis)), axis)
setattr(np, 'aad', aad)
def amd(data, axis=None):
	return np.median(np.absolute(data - np.median(data, axis)), axis)
setattr(np, 'amd', amd)
def mmd(data, axis=None):
	return np.mean(np.absolute(data - np.median(data, axis)), axis)
setattr(np, 'mmd', mmd)
def mad(data, axis=None):
	return np.median(np.absolute(data - np.mean(data, axis)), axis)
setattr(np, 'mad', mad)
def sad(data, axis=None):
	return np.sum(np.absolute(data - np.mean(data, axis)), axis)
setattr(np, 'sad', sad)

def argmid(x):
    for i, item in enumerate(x):
        if item != min(x) and item != max(x):
            return i
setattr(np, 'argmid', argmid)


def color(x):
	if(x):
		return 'red'
	else:
		return 'blue'

def get_best(variances, counts):
	index = None
	v = float("+Inf")
	C = dat.groupby('STOP').START.nunique().to_dict()
	for i, (var,count) in enumerate(zip(variances, counts)):
		c = dat[dat.CLUSTER==i].groupby('STOP').START.nunique().to_dict()
		#print(len(c))
		if var<v and var>=0 and (len(c)/len(C))>0.05: #(count/sum(counts))>0.05 and len(c)>1:
			index = i
			v = var
	return index
			




#################################
parser = argparse.ArgumentParser()
parser.add_argument('-i','--genome_id', action="store", dest="genome_id", required=True)
parser.add_argument('-d','--data_type', action="store", dest="data_type", required=True)
parser.add_argument('-t','--clust_type', action="store", dest="clust_type", required=True)
parser.add_argument('-n','--clust_num', action="store", dest="clust_num", required=True, type=int)
parser.add_argument('-a','--annotate', action="store_true", dest="annotate", required=False)
parser.add_argument('-o','--offset', action="store_true", dest="offset", required=False)
parser.add_argument('-m','--maxd', type=int, default=80, dest="maxd", required=False)

args = parser.parse_args()

#-----------------------READ IN THE DATA
path = os.path.dirname(os.path.realpath(__file__))
dat = pd.read_csv(path + "/data/aa/" + args.genome_id + ".tsv.gz", compression='gzip', header=0,sep='\t' )

with gzip.open(path + '/genomes/fna/' + args.genome_id + '.fna.gz') as f:
	first_line = f.readline().decode("utf-8")

counts = dat.groupby('STOP').START.nunique().to_dict()

#args.clust_num = int(log10(len(counts)))
#args.clust_num = 3 #round(log10(len(dat)-len(counts)))

count = dict()
mask = []
longest = dict()

coding_frame = dict()
for index,row in dat[['START', 'STOP', 'TYPE']].iterrows():
	r = list(row.values)
	count[r[1]] = count.get(r[1], 0) + 1
	#if count[r[1]] <= 10:
	#if count[r[1]] <= 4*log2(counts[r[1]]):
	#if count[r[1]] <= 10+sqrt(counts[r[1]]):
	if True: #count[r[1]] <= args.clust_num * sqrt(counts[r[1]]):

	#if count[r[1]] <= 3 or count[r[1]] <= (counts[r[1]]/3):
	#if count[r[1]] <= (counts[r[1]]/1.5):
		mask.append(False)
	else:
		mask.append(True)
	length = max(r[1]-r[0], r[0]-r[1])
	longest[r[1]] = max(longest.get(r[1],0), length) 
	if r[2] == True:
		for i in range(min(r[0], r[1]), max(r[0], r[1])):
			if (i not in coding_frame) or (r[0]>r[1]):
				coding_frame[i] = (max(r[0],r[1]) % 3) + 1
				if r[0] > r[1]:
					coding_frame[i] *= -1


dat = dat.drop(dat[mask].index)
#dat = dat[dat.CODON != 'TTG']

def sign(n):
	if n != 0:
		return abs(n)/n
	else:
		return 1.0

nums = []
offset = []
for index,row in dat[['START','STOP']].iterrows():
	r = list(row.values)
	length = max(r[1]-r[0], r[0]-r[1])
	#nums.append(r[0] / counts[r[1]])
	nums.append(length/longest[r[1]])
	#nums.append(length)
	frames = {0:0, 1:0, 2:0, 3:0, -1:0, -2:0, -3:0}
	l = (max(r[0],r[1]) % 3) + 1
	if r[0] > r[1]:
		l *= -1
	for i in range(min(r[0], r[1]), max(r[0], r[1])):
		frames[coding_frame.get(i, 0)] += 1
	#del frames[0]
	a = max(frames, key=frames.get)
	if a != 0:
		#print(r[0], r[1])
		#print(a, l)
		#print( sign(a)*sign(l), ((abs(a)+3-abs(l))%3) )
		offset.append( str(sign(a)*sign(l) * ((abs(a)+3-abs(l))%3)).replace('.0', '') )
	else:
		offset.append('IG')
dat['OFFSET'] = offset
dat.loc[dat.TYPE==True, 'OFFSET'] = True


dat.loc[dat.CODON=='ATG', 'M'] = dat.loc[dat.CODON=='ATG', 'M'] - 1
dat.loc[dat.CODON=='GTG', 'V'] = dat.loc[dat.CODON=='GTG', 'V'] - 1
dat.loc[dat.CODON=='TTG', 'L'] = dat.loc[dat.CODON=='TTG', 'L'] - 1

X = dat[list("ARNDCEQGHILKMFPSTWYV")]
X = X.div(X.sum(axis=1), axis=0)
#X.loc[:, list("ARNDCEQGHILKMFPSTWYV")]  = StandardScaler().fit_transform(X)

dat.loc[:,'A':'V'] = dat.loc[:,'A':'V'].div(dat.loc[:,'A':'V'].sum(axis=1), axis=0)

#x = X
if(args.data_type == 'se'):
	#this is the shannon entropy
	X = X * np.log(X)
	X.fillna(0, inplace=True)
	#X = np.nan_to_num(X, nan=0)
	X = X.div(X.sum(axis=1), axis=0)
	#X = X / X.sum(axis=1, keepdims=True)
	X = StandardScaler().fit_transform(X)
	title = 'shannon-entropy'
else:
	X = StandardScaler().fit_transform(X)
	title = 'aminoacid-percent'

#X = StandardScaler().fit_transform(np.concatenate((X,x),axis=1))

#X['NUM'] = dat['NUM']
#X['NUM'] = dat['NUM'].div(dat['NUM'].mean()/X.to_numpy().mean(), axis=0
#X['NUM'] = np.array(nums) * -np.log(np.array(nums))
#X['NUM'] = nums


print('len', len(X), 'clustnum', args.clust_num) 


if(args.clust_type == 'km'):
	dat['CLUSTER'] = KMeans(n_clusters=args.clust_num, n_init=50).fit(X).labels_
elif(args.clust_type == 'gm'):
	#dat['CLUSTER'] = GaussianMixture(n_components=args.clust_num,covariance_type='spherical').fit(X).predict(X)
	model = GaussianMixture(n_components=args.clust_num, n_init=10, covariance_type='spherical', reg_covar=0.00001).fit(X)
	print(model.covariances_)
	print(model.converged_)
	print(model.n_iter_)
	print(model.lower_bound_)
	dat['CLUSTER'] = model.predict(X)
elif(args.clust_type == 'bgm'):
	dat['CLUSTER'] = BayesianGaussianMixture(n_components=args.clust_num,covariance_type='spherical').fit_predict(X)
elif(args.clust_type == 'h'):
	from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
	Z = linkage(X, 'ward')
	print(args.maxd)
	dat['CLUSTER'] = fcluster(Z, args.maxd, criterion='distance') - 1
elif(args.clust_type == 'ag'):
	#dat['CLUSTER'] = AgglomerativeClustering(n_clusters=args.clust_num, affinity='cosine', linkage='complete').fit_predict(X)
	dat['CLUSTER'] = AgglomerativeClustering(n_clusters=None, distance_threshold=1.8, affinity='cosine', linkage='complete').fit_predict(X)
elif(args.clust_type == 'ward'):
	dat['CLUSTER'] = AgglomerativeClustering(n_clusters=args.clust_num, linkage='ward', compute_full_tree=True).fit_predict(X)
elif(args.clust_type == 'pam'):
	dat['CLUSTER'] = KMedoids(n_clusters=args.clust_num).fit_predict(X)
elif(args.clust_type == 'op'):
	dat['CLUSTER'] = OPTICS(min_samples=10).fit_predict(X)
elif(args.clust_type == 'sp'):
	dat['CLUSTER'] = SpectralClustering(n_clusters=args.clust_num, assign_labels="discretize").fit_predict(X)
elif(args.clust_type == 'b'):
	dat['CLUSTER'] = Birch(n_clusters=args.clust_num).fit_predict(X)

print(dat)

pca = decomposition.PCA(n_components=2).fit(X)
dat['x'],dat['y'] = pca.fit_transform(X).T

# PICK THE CLUSTER WITH THE LOWEST VARIANCE
print("----------")
v = [np.sum(np.var(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)] ; print('var', v)
v = [np.sum(np.amd(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)] ; print('amd', v)
v = [np.sum(np.mmd(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)] ; print('mmd', v)
v = [np.sum(np.aad(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)] ; print('aad', v)
v = [np.median(np.mad(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)] ; print('mmm', v)
v = [np.linalg.norm(np.mad(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)] ; print('euc', v)
v = [np.sum(np.xxx(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)] ; print('xxx', v)
v = [np.sum(np.sums(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)] ; print('sum', v)
v = [np.sum(np.mad(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)] ; print('mad', v)
l = [len(X[dat.CLUSTER==i]) for i in range(args.clust_num)]
print('   len', l)
#index_min = np.argmin(variance)
#index_min = np.argmin(v)
index_min = get_best(v, l)
print(index_min)
#index_mid = np.argmed(variance)
#index_max = np.argmax(variance)
print(dat.loc[dat.CLUSTER==index_min,].groupby('STOP').TYPE.agg(['first'])['first'].value_counts())


print('ratio:' , dat[ (dat.CLUSTER==index_min) & (dat.TYPE==True) ].groupby('STOP').START.nunique().to_dict().__len__(), dat[ (dat.TYPE==True) ].groupby('STOP').START.nunique().to_dict().__len__() , sep='\t')



title = "[" + ', '.join(str(round(e, 2)) for e in v) + "]"

dat['GOOD'] = dat.STOP.map(dat[dat.CLUSTER==index_min].STOP.value_counts())
dat['TOTAL'] = dat.STOP.map(dat.STOP.value_counts())
dat.fillna(0, inplace=True)

#print(dat.to_string())
#exit()
#########THIS MAKES THE VIOLIN DATA###########
#df = dat.drop_duplicates('STOP')
#for index, row in df.iterrows():
#	print(args.genome_id, row.STOP, row.TYPE, row.GOOD/row.TOTAL)
#exit()
#########THIS MAKES THE BAR DATA###########
#sd
#TP = len(dat[(dat.CLUSTER == index_min) & (dat.TYPE) & (dat.GOOD/dat.TOTAL > 0.5)].STOP.unique() )
#FP = len(dat[(dat.CLUSTER == index_min) & (~dat.TYPE) & (dat.GOOD/dat.TOTAL > 0.5)].STOP.unique() )
#TN = len(dat[~dat.TYPE].STOP.unique()) - FP
#FN = len(dat[dat.TYPE].STOP.unique()) - TP
#print(args.genome_id, round(np.var(x0),4), round(np.var(x1),4), round(np.var(x2),4), TP, FP, FN, TN)
#exit()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
fig, ax = plt.subplots(figsize=(3.97,3.49))
#fig.set_size_inches(4, 2.5)

'''
dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
fig.savefig(args.genome_id + '.png', dpi=300, bbox_inches='tight')
exit()
'''


# PLOT
markers = {k:v for k,v in zip({0,1,2,3}.difference({index_min}), ['o','s','H'])}

#fig.suptitle(args.genome_id + " (" + title + ")")
fig.suptitle(args.genome_id)
ax.set_xlabel('PC1', fontsize=12)
plt.ylabel('PC2')

if(args.annotate):
	colors = {True:'#3CC9CF', False:'#F2766E'}
	ax.scatter(dat['x'], dat['y'], c=dat['TYPE'].apply(lambda x: colors[x]), marker='.', linewidths=0.0, alpha=0.4, zorder=5)
	p = dat[(dat.NUM==1) & (dat.CLUSTER!=6)] #index_min]
	for i, row in p.iterrows():
		ax.annotate(str(row.STOP), (row.x, row.y), size=2, zorder=10)
		print(row.STOP)
elif(args.offset):
	colors = {True:'#3CC9CF', 'IG':'black', '-0':'b', '1':'g', '2':'r', '-1':'m', '-2':'y'}
	ax.scatter(dat['x'], dat['y'], c=dat['OFFSET'].apply(lambda x: colors[x]), marker='.', linewidths=0.0, alpha=0.4, zorder=5)
	for label in ['1', '2', '-0', '-1', '-2']:
		confidence_ellipse(dat.loc[dat.OFFSET == label, 'x'], dat.loc[dat.OFFSET == label, 'y'], ax, edgecolor=colors[label], linewidth=0.5, zorder=0)
else:
	colors = {True:'#3CC9CF', False:'#F2766E'}
	ax.scatter(dat['x'], dat['y'], c=dat['TYPE'].apply(lambda x: colors[x]), marker='.', linewidths=0.0, alpha=0.4, zorder=5)
	ax.scatter(dat[dat.CLUSTER==index_min].x, dat[dat.CLUSTER==index_min].y, facecolor='none', cmap='Spectral', s=80, linewidths=0.3, alpha=0.3, marker='d', edgecolor='black', label='Predicted', zorder=10)
	for i in range(args.clust_num):
		if i != index_min:
			ax.scatter(dat[dat.CLUSTER==i].x, dat[dat.CLUSTER==i].y, facecolor='none', cmap='Spectral', s=80, linewidths=0.3, alpha=0.5, marker=markers.get(i, '1'), edgecolor='green', label='Predicted', zorder=10)
			pass
	pass

## visualize projections
xvector = pca.components_[0] # see 'prcomp(my_data)$rotation' in R
yvector = pca.components_[1]
xs = 2.2*pca.transform(X)[:,0] # see 'prcomp(my_data)$x' in R
ys = 2.2*pca.transform(X)[:,1]
for i in range(len(xvector)):
	#arrows project features (ie columns from csv) as vectors onto PC axes
	plt.arrow(0, 0, xvector[i]*max(xs), yvector[i]*max(ys), color='black', alpha=0.2, width=0.00001, head_width=0.0025, zorder=11)
	plt.text(xvector[i]*max(xs)*1.1, yvector[i]*max(ys)*1.1, list("ARNDCEQGHILKMFPSTWYV")[i], color='black', alpha=0.5, zorder=11)
#	pass


ax.legend(handles=[mpatches.Patch(color=col, label=str(lab)) for lab,col in colors.items()])

fig.savefig(args.genome_id + '.png', dpi=300, bbox_inches='tight')
#plt.show() #block=False)
#time.sleep(4)
#plt.close("all") 
exit()







