#!/bin/bash

# Path to peps
PATH=$1
# Min length
MIN=$2
# Max length
MAX=$3


./SLMTransform.exe /lclhome/mhase003/Data/Database/Digested/parts $MIN $MAX
