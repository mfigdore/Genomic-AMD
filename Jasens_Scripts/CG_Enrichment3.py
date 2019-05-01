# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 22:37:57 2019

@author: Jasen Zhang
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.decomposition import FastICA
from CG_GSEA import CG_GSEA
from CG_GSEA2 import CG_GSEA2
from statistics import median
import seaborn as sns

i = 'ap_clustered_noreduction'
temp = pd.read_csv('5_cluster_assignments_' + i + '.csv')
order = temp['sample'].values

df = pd.read_csv('processed_expression.csv')
df = df.set_index('Unnamed: 0')
df = df.reindex(order)

MGS = pd.read_csv('assignments_MGS.csv')
MGS = MGS.set_index('Genes')
MGS = MGS.reindex(order)

'''
GSEA
'''

CG_GSEA2(df, MGS['Cluster'].values-1).to_csv('Ground_Truth_sig_genes.csv')

titles = ['ap_clustered_noreduction', 'ap_clustered_pca',
          'gmm_clustered_ica', 'gmm_clustered_pca',
          'kmeans_ica', 'kmeans_noreduction', 'kmeans_pca']

for i in titles:
    temp = pd.read_csv('5_cluster_assignments_' + i + '.csv')
    CG_GSEA2(df, temp['clusters'].values).to_csv(i + '_sig_genes.csv')