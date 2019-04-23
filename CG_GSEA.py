# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 18:10:43 2019

@author: Jasen Zhang
"""

import pandas as pd
import numpy as np
from scipy.stats import ranksums

def CG_GSEA(df, labels):
    '''
    df: dataframe of n samples by k genes
    labels: vector of n length of cluster assignment
    
    return:
        1) n by 3 array with significant genes, first column is gene index, second and third column describe the two clusters that are significantly different.
    '''
    total_labels = []
    for i in labels:
        if i not in total_labels:
            total_labels.append(i)
    
    alpha = 0.05*len(total_labels)
    
    p_values_columns = []
    for i in range(len(total_labels)-1):
        for j in range(i+1,len(total_labels)):
            p_values_columns.append(str(i)+'_' + str(j))
        
    total_p_values = np.zeros(len(p_values_columns))
    count = 0
    for k in range(len(df.columns)):
        p_values = np.zeros(len(p_values_columns))
        gene = df[df.columns[k]].values
        gene_labels = np.vstack((gene, labels))
        gene_labels_df = pd.DataFrame(np.transpose(gene_labels))
        index = 0
        for i in range(len(total_labels)-1):
            genes_i = gene_labels_df[gene_labels_df[1] == total_labels[i]][0]
            for j in range(i+1,len(total_labels)):
                genes_j = gene_labels_df[gene_labels_df[1] == total_labels[j]][0]
                p_values[index] = (ranksums(genes_i, genes_j)[1])
                index += 1
        total_p_values = np.vstack((total_p_values, p_values))
        count += 1
        if count % 100 == 0:
            print(count)
            
    total_p_values = np.delete(total_p_values,0,0)
    total_p_values_1 = np.reshape(total_p_values,-1)
    total_p_values_1 = np.vstack((total_p_values_1, np.arange(len(total_p_values_1))))
    total_p_values_1 = pd.DataFrame(np.transpose(total_p_values_1), columns = ['p_value', 'ID'])
    total_p_values_1 = total_p_values_1.sort_values(by=['p_value'])
    
    BH = np.arange(1, len(total_p_values_1)+1)*(alpha/len(total_p_values_1))
    total_p_values_1['BH'] = BH
    total_p_values_1['Bonferroni'] = alpha/len(total_p_values_1)
    total_p_values_1['BH_sig'] = total_p_values_1['p_value'] < total_p_values_1['BH']
    total_p_values_1['Bonferroni_sig'] = total_p_values_1['p_value'] < total_p_values_1['Bonferroni']
    
    a = sum(total_p_values_1['BH_sig'])
    b = sum(total_p_values_1['Bonferroni_sig'])
    sig_ID = np.array(total_p_values_1[total_p_values_1['BH_sig'] == True]['ID'].values, dtype = int)
    sig_genes = []
    for i in sig_ID:
        sig_genes.append([(i//len(p_values_columns)), int(p_values_columns[i%len(p_values_columns)].split('_')[0]), int(p_values_columns[i%len(p_values_columns)].split('_')[1])])
    
    new_sig_genes = []
    for i in sig_genes:
        temp = [df.columns[i[0]], i[1], i[2]]
        new_sig_genes.append(temp)
    
    return new_sig_genes
