#!@PYTHON_EXECUTABLE@
# Required Imports
import sys
import os
import glob
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

    for kk in np.arange(minlen, maxlen + 1):
        # Figure and axis object
        fig, ax = plt.subplots()
        ax.grid(linestyle=':', linewidth=0.5)
        ax.set(xlabel='Number of Features', ylabel='Frequency', title='Feature Freq vs Num Features at ' + str(pcnt*100) + '% and plen= ' + str(kk))

        # Extract data from the CSV format
        for filename in glob.iglob(dirname + '/**', recursive=True):
            if os.path.isfile(filename):
                if filename[-(15 + len(str(kk))):] == "ms2databA" + str(kk) + "_0.csv":
#                    print (filename)
                    df = pd.read_csv(filename)
#                    print (df.head(10))
                    cols = list(df.columns.values)
                    max = df[cols[1]].max()
#                    print (max)
                    df_filtered = df.loc[df[cols[1]] >= (max * pcnt)].sort_values(by=[cols[1]], ascending=False)
                    ax.plot(np.arange(1, df_filtered[cols[1]].count() + 1), df_filtered[cols[1]], label= filename[len(dirname) + 1:-(15 + len(str(kk)) + 1)], linewidth=1.2)
                    ax.legend()
#                    print (df_filtered.head(32))

    # Set axis stuff
        ax.xaxis.set_major_locator(LinearLocator(4))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        start, end = ax.get_xlim()
        ax.set_xticks(np.arange(start, end + 1, 2))

    # Show the plot
#        plt.show()
        plt.savefig(dirname +'/locmods' + str(kk) + '.jpg', dpi=350)