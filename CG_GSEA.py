# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 18:10:43 2019

@author: Jasen Zhang
"""

import pandas as pd
import numpy as np
from scipy.stats import ranksums

df = pd.read_csv('final.csv')
labels = np.array(real_labels.values, dtype = int)

total_labels = []
for i in labels:
    if i not in total_labels:
        total_labels.append(i)
        
p_values = np.zeros([len(total_labels),len(total_labels)])

count = 0
for k in range(len(df.columns)):
    gene = df[df.columns[k]].values
    gene_labels = np.vstack((gene, labels))
    gene_labels_df = pd.DataFrame(np.transpose(gene_labels))
    p_value_matrix = np.zeros([4,4])
    for i in range(len(total_labels)):
        genes_i = gene_labels_df[gene_labels_df[1] == total_labels[i]][0]
        for j in range(len(total_labels)):
            genes_j = gene_labels_df[gene_labels_df[1] == total_labels[j]][0]
            p_value_matrix[i, j] = ranksums(genes_i, genes_j)[1]
    p_values = np.dstack((p_values, p_value_matrix))
    count += 1
    print(count)
    if count == 5000:
        break
    
'''
Number significant
'''
significant = sum(sum(p_values < 0.05/(p_values.shape[2]-1)))

num_significant = sum(significant > 0) - 1
        
    
    
