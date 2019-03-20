import sys
import time
import numpy as np
import pandas as pd 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from matplotlib.mlab import PCA
from sklearn.preprocessing import StandardScaler

from sklearn import decomposition
from sklearn import datasets
from sklearn.mixture import GaussianMixture
from sklearn.mixture import BayesianGaussianMixture
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering

import warnings
warnings.filterwarnings('ignore')

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

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))
setattr(np, 'mad', mad)
def aad(data, axis=None):
	    return np.mean(np.absolute(data - np.mean(data, axis)), axis)
setattr(np, 'aad', aad)

dat = pd.read_csv("aa/" + sys.argv[1] + ".aa", header=0,sep=' ' )
end = pd.read_csv("prodigal/" + sys.argv[1] + ".ends", header=None)

dat['TYPE'] = dat['STOP'].isin(end[0])

X = dat[list("ARNDCEQGHILKMFPSTWYV")]
X = X.div(X.sum(axis=1), axis=0)
#this is the shannon entropy
X = X * np.log(X)
X = X.div(X.sum(axis=1), axis=0)
X.fillna(0, inplace=True)

X_std = StandardScaler().fit_transform(X)
tag = sys.argv[1]

pca = decomposition.PCA(n_components=2)
pca.fit(X_std)
P = pca.transform(X_std)

#DO THE MIXTUREMODEL
#gm = GaussianMixture(n_components=2,covariance_type='tied').fit(X_std).predict(X_std)
#gm = BayesianGaussianMixture(n_components=3,covariance_type='tied').fit_predict(X_std)
gm = KMeans(n_clusters=3).fit(X_std).labels_
dat['CLUSTER'] = gm
#gm = AgglomerativeClustering(n_clusters=2).fit_predict(X_std)

x0 = X_std[gm==0]
x1 = X_std[gm==1]
x2 = X_std[gm==2]
variance = [np.var(x0),np.var(x1),np.var(x2)]
variance = [v for v in variance if str(v) != 'nan']
index_min = np.argmin(variance)
index_mid = np.argmed(variance)
index_max = np.argmax(variance)
x0, x1, x2 = X_std[gm==index_min], X_std[gm==index_mid], X_std[gm==index_max]
p0 = P[gm==index_min]
p1 = P[gm==index_mid]
p2 = P[gm==index_max]

p = P[dat.TYPE]
plt.scatter(p[:, 0], p[:, 1], c='g', alpha=0.5, edgecolors='none', label='real')
p = P[~dat.TYPE]
#plt.scatter(p[:, 0], p[:, 1], c='r', alpha=0.5, edgecolors='none', label='fake')
#plt.scatter(p0[:, 0], p0[:, 1], facecolor='none', cmap='Spectral', alpha=0.5, marker='d', edgecolor='black', label='Predicted')
plt.suptitle(tag + " (shannon entropy,kmeans)", fontsize=20)
plt.legend()
plt.savefig('fig/' + tag + '.png', dpi=150)


