# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 23:39:30 2019

@author: Jasen Zhang
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def CG_PCA(df):
    '''
    df: n rows by m columns numpy array where each row corresponds to a sample and each column is a gene
    return: n by k dataframe after PCA selects appropriate k by elbow rule
    '''
    
    model = PCA(n_components = len(df.columns)).fit(df)
    exp_var = model.explained_variance_ratio_
    plt.scatter(np.arange(1, len(exp_var) + 1), exp_var)
    plt.savefig('PCA_Var_Explained.png')
    plt.show()

    cum_variances = []
    sum_variances = 0
    for i in exp_var:
        sum_variances += i
        cum_variances.append(sum_variances)

    plt.scatter(np.arange(1,len(cum_variances)+1), cum_variances, s = 1)
    plt.scatter(np.arange(1,len(cum_variances)+1), 0.865*np.ones(len(cum_variances)), s = 1, alpha = 0.3)
    plt.ylim(0,1.1)
    plt.title('Selection of K Principal Components')
    plt.savefig('PCA_Cumulative_Var_Explained.png')
    plt.show()

    best_components = sum(np.array(cum_variances) < 0.865)

    return PCA(n_components = best_components).fit_transform(df)