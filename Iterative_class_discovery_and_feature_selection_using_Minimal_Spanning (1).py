#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy.spatial import distance
from scipy.stats import t


# In[2]:


def FS(x, u, v, m):
    n = x.shape[0]
    c = v.shape[0]

    um = u**m

    v_mean = v.mean(axis=0)

    d2 = scipy.spatial.distance.cdist(x, v)**2

    distance_v_mean_squared = np.linalg.norm(v - v_mean, axis=1, keepdims=True)**2

    return np.sum(um.T*d2) - np.sum(um*distance_v_mean_squared)


# In[6]:


def compute_mst():
    
    gene = [[]]## our graph
    amount_of_genes = 10
    edge_num  = 0
    selected_nodes = []
    fs_min = 0
    
    for i in range(amount_of_genes):
        selected_nodes.append(0)

    
    
    while (edge_num < amount_of_genes - 1):
        
        fs_array = [[]] ##all 0
        
        minimum = 9223372036854775807
        a = 0
        b = 0
        for i in range(amount_of_genes):
            if selected_nodes[i]:
                for j in range(amount_of_genes):
                    if ((not selected_nodes[j]) and gene[i][j]):  
                        if minimum > gene[i][j]:
                            minimum = gene[j][i]
                            a = i
                            b = j
        selected_nodes[b] = True
        edge_num += 1
        
        
        for i in range(len(selected_nodes)):
            if i:
                for j in range(i, len(selected_nodes)):
                    if j:
                        fs_array[i][j] = gene[i][j]
                        
        if(edge_num == 0):
            fs_min = FS(fs_array)
            fs_init = selected_nodes
        elif(fs_min > FS(fs_array)):
            fs_min = FS(fs_array)
            fs_init = selected_nodes 


# In[11]:


MaxNp = int(input())
P_thresh = 0.999
x_s_g = []
n = 100 # n - number of samples
F_set = [i for i in range(1, n+1)]
M = 1000 # is the total number of genes
N_p = 10 # Number of currently discovered partitions
T_thresh = P_thresh / 2

while N_p < MaxNp:
    F = F_set
    t = T_thresh / 2
    
    while True:
        if len(F) < 2*M*(1-P_thresh/100):
            #Not enough genes support partitions
            quit()   
        ...
        t_g = [t(i) for i in P]
        F_new = t_g
        ...
        if F_new == F and t == T_th:
            #Feature set has converged
            print(P, F)
            F_set = np.setdiff(f_set, F)
            Np = Np + 1
            break
        else:
            F = F_new
            t = t + 1


# In[ ]:




