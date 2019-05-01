# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 18:10:43 2019

@author: Jasen Zhang
"""

import pandas as pd
import numpy as np
from scipy.stats import ranksums
import matplotlib.pyplot as plt

def CG_GSEA2(df, labels):
    '''
    df: dataframe of n samples by k genes
    labels: vector of n length of cluster assignment
            0 means no MGS
            1-4 are clusters of samples with MGS
    
    return:
        1) n by 3 array with significant genes, first column is gene index, second and third column describe the two clusters that are significantly different.
        2) list of p-values in order
    '''
    
    total_labels = []
    for i in labels:
        if i not in total_labels:
            total_labels.append(i)
    
    if 0 in total_labels:
        total_labels.remove(0)
            
            
    df['Cluster'] = labels
    
    noMGS = df[df['Cluster'] == 0]
    noMGS = noMGS.drop(['Cluster'], axis = 1)
    alpha = 0.05
    num_genes = noMGS.shape[1]
    p_values = np.zeros(num_genes)
    first = True
    for i in total_labels:
        yesMGS = df[df['Cluster'] == i]
        yesMGS = yesMGS.drop(['Cluster'], axis = 1)
        for j in range(num_genes):
            no_gene = noMGS.iloc[:,j]
            yes_gene = yesMGS.iloc[:,j]
            p_values[j] = (ranksums(no_gene, yes_gene)[1])
            if j % 1000 == 0:
                print(j)
        p_values_df = pd.DataFrame(np.transpose(p_values), columns = ['p_value'])
        p_values_df = p_values_df.sort_values(by = 'p_value')
        BH = np.arange(1, len(p_values_df)+1)*(alpha/len(p_values_df))
        p_values_df['BH'] = BH
        p_values_df['BH_sig'] = p_values_df['p_value'] < p_values_df['BH']
        
        plt.scatter(np.arange(1, len(p_values_df)+1), p_values_df['p_value'], s = 1)
        plt.scatter(np.arange(1, len(p_values_df)+1), p_values_df['BH'], s = 1)
        plt.xlabel('Genes')
        plt.ylabel('p-value')
        title1 = 'Cluster ' + str(i)
        plt.title(title1)
        plt.legend('p_value', 'BH Threshold')
        plt.show()
        
        
        sig_ID = np.array(p_values_df[p_values_df['BH_sig'] == True].index, dtype = int)
        sig_ID = np.vstack((sig_ID, i*np.ones(len(sig_ID))))
        if first:
            total_sig_ID = sig_ID
            first = False
        else:
            total_sig_ID = np.hstack((total_sig_ID, sig_ID))
    
    total_sig_ID_df = pd.DataFrame(np.transpose(total_sig_ID), columns = ['Gene', 'Cluster'])
    ENSG = []
    for i in total_sig_ID_df['Gene']:
        ENSG.append(df.columns[int(i)])
    total_sig_ID_df['Gene'] = ENSG

    
    
    return total_sig_ID_df
