from sklearn.datasets import make_blobs
import pandas as pd
import hdbscan

blobs, labels = make_blobs(n_samples=2000, n_features=10)

print (labels)

print (pd.DataFrame(blobs).head())

clusterer = hdbscan.HDBSCAN(metric='euclidean')

clusterer.fit(blobs)

print (clusterer.labels_)

print (clusterer.labels_.max())