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

    for kk in np.arange(minlen, maxlen + 1):
        # Figure and axis object
        fig, ax = plt.subplots()
        ax.grid(linestyle=':', linewidth=0.5)
        ax.set(xlabel='Spectrum Position', ylabel='Unique Values', title='Unique Values at Spectrum Position for plen= ' + str(kk))

        # Extract data from the CSV format
        for filename in glob.iglob(dirname + '/**', recursive=True):
            if os.path.isfile(filename):
                if filename[-(15 + len(str(kk))):] == "ms2datafA" + str(kk) + "_0.csv":
#                    print (filename)
                    df = pd.read_csv(filename, header=None, names=['pos', 'uniq'])
                    ax.plot(np.arange(1, df['uniq'].count() + 1), df['uniq'], label= filename[len(dirname) + 1:-(15 + len(str(kk)) + 1)], linewidth=1.2, marker='x')
                    ax.legend()

    # Set axis stuff
        start, end = ax.get_xlim()
        ax.set_xticks(np.arange(1, end + 1, int(end/15)))
    # Show the plot
#        plt.show()
        plt.savefig(dirname +'/uniqs' + str(kk) + '.jpg', dpi=350)

        plt.close()