# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 23:11:33 2019

@author: Jasen Zhang
"""

import numpy as np
import matplotlib.pyplot as plt

def CG_Filter_Nonexpression(df, high, thresh):
    '''
    input:
            df: n samples by m features dataframe
            high: how high can the number be and still be considered missing
            thresh: how much missing data is tolerated, ie 20% means you're ok with at most 20% of the feature being < high
    output:
           1) n samples by m features dataframe, after removing features
    '''
    percent_missing = []
    for i in range(len(df.columns)):
        missing_count = 0
        for j in df.iloc[:,i]:
            if float(j) < high:
                missing_count += 1
        percent_missing.append(missing_count/len(df))
        
    plt.hist(percent_missing)
    plt.title('Distribution of Genes with Missing Values')
    plt.xlabel('Percent Missing')
    plt.ylabel('Frequency')
    plt.savefig('Pre-filter_Graph.png')
    plt.show()
    
    values = np.arange(1,len(percent_missing)+1)
    percent_missing_2 = np.sort(percent_missing)

    plt.scatter(values, percent_missing_2, s =20)
    plt.title('Cumulative Distribution of Genes with Missing Values')
    plt.xlabel('Genes')
    plt.ylabel('Percent Missing')
    plt.savefig('Distribution_of_Missing_Values.png')
    plt.show()
    
    good_columns = df.columns[(np.array(percent_missing) <= thresh)]
    return df[good_columns]