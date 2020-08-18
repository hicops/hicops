#!@BASH_EXECUTABLE@
# 
# kth occurrence of a pattern finder
# Copyrights(C) 2019 PCDS Laboratory
# Muhammad Haseeb, and Fahad Saeed
# School of Computing and Information Sciences
# Florida International University (FIU), Miami, FL
# Email: {mhaseeb, fsaeed}@fiu.edu
#

awk '/$2/{i++}i==$3{print; exit}' $1
