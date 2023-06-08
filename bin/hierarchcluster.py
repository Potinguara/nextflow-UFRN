#!/usr/bin/env python

from Bio import SeqIO
from Bio import Entrez
from sklearn.cluster import KMeans
from shutil import copyfile
from urllib.request import urlopen
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import nan_euclidean_distances
from sklearn.ensemble import RandomForestClassifier as rf

import os
import glob
import time
import re
import pickle
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# pasta onde você pode encontrar os arquivos
folder = sys.argv[1]

samples_sflavius = pd.read_csv(sys.argv[1], sep=" ", index_col=0)
loci_sflavius=np.transpose(samples_sflavius)


#Computando matriz de distância ignorando NaNs
mydist_samples=nan_euclidean_distances(samples_sflavius,samples_sflavius)
mydist_loci=nan_euclidean_distances(loci_sflavius,loci_sflavius)

#Utilizando dados computados na matriz customizada para construir o plot do heatmap (metodo ward)
fig, axs = plt.subplots(1,2)
Z_columns = hierarchy.linkage(np.transpose(mydist_loci), metric='euclidean', method = "ward")
hierarchy.dendrogram(Z_columns, ax=axs[0])
Z_rows = hierarchy.linkage(mydist_samples, metric='euclidean', method="ward")
hierarchy.dendrogram(Z_rows, orientation='left', ax=axs[1])

g = sns.clustermap(samples_sflavius, row_linkage=Z_rows, col_linkage=Z_columns)

g.savefig('fplot.pdf')



