#!@PYTHON_EXECUTABLE@
# 
# Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
# Florida International University, Miami, FL# 
# This program is licensed under the
# Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
# See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
# 


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