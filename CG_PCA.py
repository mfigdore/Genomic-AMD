# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 23:39:30 2019

@author: Jasen Zhang
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def CG_PCA(df, elbow, manual):
    '''
    input:
            df: n rows by m columns numpy array where each row corresponds to a sample and each column is a gene
            elbow: manual elbow threshold
            manual: manual number of principal components
    return: 
            1) n by k dataframe after PCA selects appropriate k by elbow rule
            2) n by k dataframe after PCA selects manual k
    '''
    
    model = PCA(n_components = len(df.columns)).fit(df)
    exp_var = model.explained_variance_ratio_
    plt.scatter(np.arange(1, len(exp_var) + 1), np.log(exp_var))
    plt.title('Percent Variance Explained')
    plt.xlabel('Principal Component')
    plt.ylabel('Log Percent')
    plt.savefig('PCA_Var_Explained.png')
    plt.show()

    cum_variances = []
    sum_variances = 0
    for i in exp_var:
        sum_variances += i
        cum_variances.append(sum_variances)

    plt.scatter(np.arange(1,len(cum_variances)+1), cum_variances, s = 1)
    plt.scatter(np.arange(1,len(cum_variances)+1), elbow*np.ones(len(cum_variances)), s = 1, alpha = 0.3)
    plt.ylim(0,1.1)
    plt.title('Selection of K Principal Components')
    plt.xlabel('Principal Component')
    plt.ylabel('Cumulative Percent')
    plt.legend(['Cumulative Percent', 'Elbow Rule Threshold'])
    plt.savefig('PCA_Cumulative_Var_Explained.png')
    plt.show()

    best_components = sum(np.array(cum_variances) < elbow)

    return PCA(n_components = best_components).fit_transform(df), PCA(n_components = manual).fit_transform(df)