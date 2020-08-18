#!@PYTHON_EXECUTABLE@
# R5equired Imports
import numpy as np
import pandas as pd
import hdbscan
from sklearn.cluster import DBSCAN
from sklearn import metrics
from joblib import Memory
import sys
import distance
from sklearn.datasets import make_blobs

# Read all CSV files of peplen and return the matrix of spactra
def readCSV(dirname, peplen):
    filename = dirname + "/" + "ms2data" + str(peplen) + "_0.csv"

    # Read CSV as Data Frame
    df = pd.read_csv(filename)

    # Return the array of arrays
    return df.values


# The main function
if __name__ == '__main__':
    # Status variable
    status = 0

    # Read the data directory
    if len(sys.argv) > 1:
        dirname = sys.argv[1]
    else:
        status = -1
        print ("Error: Enter the data dir path")
        sys.exit(status)

    # Read the spectra lengths
    if len(sys.argv) > 2:
        peplen = int(sys.argv[2])
    else:
        status = -1
        print ("Error: Enter the peplen param\n")
        sys.exit(status)

    # Extract data from the CSV format
    data = readCSV(dirname, peplen)

#    print (data[0])
#    print (data[1])
#    print (list(set(data[0]) & set(data[1])))
#    print (distance.jaccard(data[0], data[1]))
#    print (data)
#    print (data.shape)

    # Construct the clustering model
    #model = hdbscan.HDBSCAN(algorithm='best', allow_single_cluster=False, alpha=1.0,
    #approx_min_span_tree=True, core_dist_n_jobs=8, gen_min_span_tree=True, memory=mem, 
    #metric='euclidean', min_cluster_size=3, min_samples=3, p=None)

    model = DBSCAN(eps=0.7, min_samples=6, metric='jaccard', 
	               metric_params=None, algorithm='auto', 
				   leaf_size=30, p=2, n_jobs=8).fit(data)

#algorithm='best', metric='euclidean', min_cluster_size=3)

    # Construct clusters and display the clusters list
    clist = model.labels_

    # Print the clusters
    print (clist)
    print (clist.shape)
