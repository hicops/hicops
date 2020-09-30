#!@PYTHON_EXECUTABLE@
# Required Imports
import numpy as np
import random as rd
import pandas as pd
import sys
import pickle

# Read all CSV files of peplen and return the matrix of spactra
def readCSV(dirname, peplen):
    filename = dirname + "/" + str(peplen) + ".peps"

    # Read CSV as Data Frame
    df = pd.read_csv(filename, header=None)

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
        minlen = int(sys.argv[2])
    else:
        minlen = 6

    if len(sys.argv) > 3:
        maxlen = int(sys.argv[3])
    else:
        maxlen = 30

    if len(sys.argv) > 4:
        total = int(sys.argv[4])
    else:
        total = 100

    each = int(total/(maxlen-minlen+1))
    seqs = []
    

    for len in range(minlen, maxlen+1):
        # Extract data from the CSV format
        data = readCSV(dirname, len)

        data = data.flatten()

        for ee in range(0,each):
            idx = rd.randint(0,data.size)
            seqs.append(str(data[idx]))

    ll = np.array(seqs)
    np.savetxt(dirname + 'chosen.txt', ll, fmt="%s")
