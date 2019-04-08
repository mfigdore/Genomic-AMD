# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 17:09:09 2019

@author: Jasen
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

def CG_Kmeans(df, ks, manual):
    '''
    df: n rows by m columns numpy array where each row corresponds to a sample and each column is a gene
    return: number of clusters,
            array of n entries which maps to cluster
    '''
    observations = df.shape[0]
    features = df.shape[1]
    SSE = []
    BIC = []
    for k in np.arange(1,ks):
        model = KMeans(n_clusters = k, init = 'random', max_iter = 10).fit(df)
        SSE.append(model.inertia_)
        temp = observations*np.log(SSE[k-1]/observations) + features*k*np.log(observations)
        BIC.append(temp)
    
    
    plt.scatter(np.arange(1,ks), BIC)
    plt.xlabel('Number of Clusters')
    plt.ylabel('BIC')
    plt.savefig('BIC of Kmeans')
    plt.show()
    
    best = np.argmin(BIC)+1
    
        
    '''
    ELBOW RULE
    
    SSE_diff = []
    for i in range(1,len(SSE)):
        SSE_diff.append(SSE[i] - SSE[i-1])
    SSE_diff2 = []
    for i in range(1,len(SSE_diff)):
        SSE_diff2.append(SSE_diff[i] - SSE_diff[i-1])
    best = np.argmax(SSE_diff2) + 2
    '''
    
    return best, (KMeans(n_clusters = best, init = 'random', max_iter = 10).fit_predict(df)),  (KMeans(n_clusters = manual, init = 'random', max_iter = 10).fit_predict(df))