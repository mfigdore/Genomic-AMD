# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 17:09:09 2019

@author: Jasen
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

def CG_Kmeans(df, ks, manual, feature_red, want_graphs):
    '''
    input:
            df: n rows by m columns numpy array where each row corresponds to a sample and each column is a gene
            ks: maximum amount of clusters considered
            manual: number of clusters you want to see
            feature_red: string. what type of feature reduction was implemented
            want_graphs: boolean for saving graphs or not
    return: 1) number of clusters selected by BIC
            2) array of n entries which maps to cluster, optimized by BIC
            3) array of n entries with maps to manually set # of clusters
    '''
    observations = df.shape[0]
    features = df.shape[1]
    SSE = []
    BIC = []
    for k in np.arange(1,ks):
        model = KMeans(n_clusters = k, init = 'random', max_iter = 10).fit(df)
        SSE.append(model.inertia_)
        
        temp = observations*np.log(SSE[k-1]/observations) + features*k*np.log(observations)
        '''
        temp = -2*np.log(SSE[k-1]) + features*k*np.log(observations)
        '''
        BIC.append(temp)
    
    
    plt.scatter(np.arange(1,ks), BIC)
    title_str = 'Selection of Number of Clusters with BIC after Applying ' + feature_red
    plt.title(title_str)
    plt.xlabel('Number of Clusters')
    plt.ylabel('BIC')
    fig_str = 'BIC_of_Kmeans_using ' + feature_red + '.png'
    if want_graphs:
        plt.savefig(fig_str)
        plt.show()
    
    best = np.argmin(BIC)+1
    
    return best, (KMeans(n_clusters = best, init = 'random', max_iter = 10).fit_predict(df)),  (KMeans(n_clusters = manual, init = 'random', max_iter = 10).fit_predict(df))