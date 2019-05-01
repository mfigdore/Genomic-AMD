# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 18:41:11 2019

@author: Jasen
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CG_Filter_Nonexpression import CG_Filter_Nonexpression
from sklearn import preprocessing
from sklearn.decomposition import PCA
from CG_PCA import CG_PCA
from sklearn.cluster import KMeans
from CG_Kmeans import CG_Kmeans
from sklearn.metrics import adjusted_rand_score
from sklearn.decomposition import FastICA
from CG_GSEA import CG_GSEA




labels2 = []
first = True
f = open('EyeGEx_meta_combined_inferior_retina_summary_deidentified.txt', 'r')
for line in f:
    temp = line.rstrip()
    temp = temp.split(' ')
    if first:
        col_names = temp
        first = False
    else:
        if 'no' in temp:
            del temp[6]
            temp[5] = 'NA'
        if 'or' in temp:
            del temp[6]
            del temp[6]
        if line == '406_1 406 OS  79 M 1 8.2 15.9\n':
            del temp[3]
        labels2.append(temp)

labels2_df = pd.DataFrame(labels2, columns = col_names)

tsv = []
first = True
f = open('EyeGEx_retina_combined_genelevel_expectedcounts_byrid_nooutliers.counts.matrix.tsv', 'r')
for line in f:
    temp = line.rstrip()
    temp = temp.split('\t')
    if first:
        col_names = temp
        first = False
    else:
        tsv.append(temp)

col_names[0] = 'Gene'
tsv_df = pd.DataFrame(tsv, columns = col_names)

tsv_df = tsv_df.set_index('Gene')

index_name = labels2_df.columns[0]
labels2_df = labels2_df.set_index(index_name)

tsv_df2 = np.transpose(tsv_df)

new_index = []
for i in list(tsv_df2.index):
    new_index.append(i[1:-1])

tsv_df2['Genes'] = new_index

tsv_df2 = tsv_df2.set_index('Genes')

final = tsv_df2.join(labels2_df)

real_labels = final['mgs_level']


MGS = pd.DataFrame(real_labels.values, index = real_labels.index, columns = ['Cluster'])
MGS.to_csv('assignments_MGS.csv')

