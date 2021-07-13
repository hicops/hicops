#!@PYTHON_EXECUTABLE@

# 
# Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
# Florida International University, Miami, FL# 
# This program is licensed under the
# Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
# See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
# 

# Required Imports
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

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

    if len(sys.argv) > 2:
        minlen = int (sys.argv[2])
    else:
        minlen = 6

    if len(sys.argv) > 3:
        maxlen = int (sys.argv[3])
    else:
        maxlen = 40

    if len(sys.argv) > 4:
        pcnt = int (sys.argv[4])/100
    else:
        pcnt = 15/100
		
    if len(sys.argv) > 5:
        mpcn = 1.0 + int (sys.argv[5])/100
    else:
        mpcn = 1.0 + 27/100

    # Figure and axis object
    fig, ax = plt.subplots()
    ax.grid(linestyle=':', linewidth=0.5)
    ax.set(xlabel='Number of Features', ylabel='Frequency', title='Feature (peak) frequency vs Number of Features')

    # Store number of features
    features = np.zeros(maxlen + 1 - minlen)
    index = np.arange(minlen, maxlen + 1)

    totalpk = 0
    saved = 0
    # Extract data from the CSV format
    for kk in np.arange(minlen, maxlen + 1):
        filename = dirname + "/" + "ms2databA" + str(kk) + "_0.csv"
#        print (filename)
        df = pd.read_csv(filename, header=None, names=['m_z', 'freq'])
#        print (df.head(10))
        cols = list(df.columns.values)
        max = df[cols[1]].max()
        totalpk += df[cols[1]].sum()
#        print (max)
        df_filtered = df.loc[df[cols[1]] >= (max * pcnt)].sort_values(by=[cols[1]], ascending=False)
        min = df_filtered[cols[1]].min()
        df_filtered = df_filtered.loc[df_filtered[cols[1]] > (min * mpcn)]
        saved += df_filtered[cols[1]].sum()
        ax.plot(np.arange(1, df_filtered[cols[1]].count() + 1), df_filtered[cols[1]], label= 'plen=' + str(kk), linewidth=1.2)
#        ax.legend()
        features[kk-minlen] = df_filtered.shape[0]
#        print (df_filtered.head(32))

    # Set axis stuff
    start, end = ax.get_xlim()
    ax.set_xticks(np.arange(1, end + 1, 2))
    print ('Total Peaks in Database =', totalpk)
    print ('Memory for Indexing =', (totalpk * 4/(1024*1024)).round(2), 'MB')
    print ('Extracted Peaks =', saved)
    print ('Memory Saved =', (saved * 4/(1024*1024)).round(2), 'MB')
    print ('Percentage Conservation =', str((saved*100/totalpk).round(2))+'%')
    # Show the plot
#    plt.show()
    plt.savefig(dirname +'/localization.jpg', dpi=350)

    #  Construct the data frame for number of features for each pep len
    fset = pd.DataFrame({'Peptide Lengths':index, 'No. of Features':features})
#    print (fset)

    # Plot the data frame
    tt = 'Features within ' + str(pcnt*100) + '% vs Peptide Lengths'
    fset = pd.DataFrame({'Peptide Lengths':index, 'No. of Features':features})
    fset.plot('Peptide Lengths', 'No. of Features', title = tt)
    plt.grid(linestyle=':', linewidth=0.5)

    # Save figure
#    plt.show()
    plt.savefig(dirname + '/features.jpg', dpi=350)