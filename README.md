# goodorfs
Experimental  code to classify open reading frames


## Quick start
This command will take the amino-acid counts of the file `data/aa/1002725.9.tsv.gz` and perform the following:
<pre><code>
1. Divide the amino-acid counts by rowsums to get frequency
2. Take the `-d` flag to take either the amino-acid percent [ap] or the Shannon Entropy [se]
3. Normalize the above by the Z-score
4. Cluster the ORFs using the type `-t` of algorithm: kmeans [km], Gaussian Mixture Model [gm], Bayesian Gaussian Mixturea [bgm], or Agglomerative Clustering [ag]
5. Find the variances of the points for each cluster, and assigning the cluster with the lowest variance as the 'good' cluster
6. Decompose the normalized scores using Principle Component Analysis from 20 dimensions into 2
7. Plot all the ORFs in 2-dimensional space using the 2 components using `Matplotlib`

</code></pre>
Command Line:
```sh
python3 main.py -i 1002725.9 -d se -t km -n 3
```

## Output
If **goodorfs** ran correctly you should get an output image titled 1002725.9.png, which should look like the image below:
![](https://github.com/deprekate/goodorfs/blob/master/1002725.9.png)
