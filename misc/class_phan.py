import sys
import time
import numpy as np
import pandas as pd 

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

#np.random.seed(5)
S = []
E = []
L = np.array([])
X = []
Y = []
convert = {'fake':'r', 'real':'g', 'other':'gray'}

dat = pd.read_csv(sys.argv[1], header=0,sep=' ' )
tag = sys.argv[1].split("/")[-1].replace("aa", "ends")
end = pd.read_csv("phan/" + tag, header=None, sep=' ' )

X = dat[list("ARNDCEQGHILKMFPSTWYV")]
X = X.div(X.sum(axis=1), axis=0)

X = X * np.log(X)
X.fillna(0, inplace=True)
X = X.div(X.sum(axis=1), axis=0)
X.fillna(0, inplace=True)

X_std = StandardScaler().fit_transform(X)

#DO THE MIXTUREMODEL
#gm = GaussianMixture(n_components=3,covariance_type='tied').fit(X_std).predict(X_std)
#gm = BayesianGaussianMixture(n_components=3,covariance_type='tied').fit_predict(X_std)
gm = KMeans(n_clusters=3).fit(X_std).labels_
#gm = AgglomerativeClustering(n_clusters=2).fit_predict(X_std)

dat['gm'] = gm

x0 = X_std[gm==0]
x1 = X_std[gm==1]
x2 = X_std[gm==2]
index_min = np.argmin([np.var(x0),np.var(x1),np.var(x2)])
index_mid = np.argmed([np.var(x0),np.var(x1),np.var(x2)])
index_max = np.argmax([np.var(x0),np.var(x1),np.var(x2)])
x0, x1, x2 = X_std[gm==index_min], X_std[gm==index_mid], X_std[gm==index_max]

G = dict()
T = dict()
for index, row in dat.iterrows():
	if(row.gm == index_min):
		val = G.get(row.STOP, 0)
		G[row.STOP] = val+1
	val = T.get(row.STOP, 0)
	T[row.STOP] = val+1

for o in T:
	if(o in end.values):
		print(tag.replace(".ends","")+'_'+str(o), G.get(o,0)/T[o], "coding")
	else:
		print(tag.replace(".ends","")+'_'+str(o), G.get(o,0)/T[o], "noncoding")

exit()



