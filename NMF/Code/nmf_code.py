# -*- coding: utf-8 -*-
"""
Created on Mon 26 Nov 20:18:47 2018

@author: Razak nart
"""

# Importing the libraries
import pandas as pd;
import numpy as np;
import matplotlib.pyplot as plt
import logging
from sklearn import metrics
import seaborn as sns
import os
from sklearn.decomposition import NMF
from sklearn.preprocessing import normalize


# Defining logger for nmf 
logger = logging.getLogger('nmf_logger')
logging.basicConfig(level=logging.INFO) 

# Load sample files        
indir = 'Subclonal_Structures/'

#Count the number of files in Subclonal_Structure
num_files = len([f for f in os.listdir(indir) if os.path.isfile(os.path.join(indir, f))])

samples = np.zeros((10, num_files), dtype=np.int32)
samples_CCF_values = np.array(np.arange(0.1, 1.1, 0.1))

""" Create Sample Matrix
    Takes a dataset containing the column (proportion/CCF) to insert into Sample matrix
    and also takes the file number from the sample to read"""
def create_matrix(ds, f_number):
    n_ssms = ds['n_ssms']
    ccf = ds['proportion']
    for ccf_row in np.arange(0, ccf.shape[0]):
        for ccf_band, samples_row in zip(np.arange(0, 1, 0.1), np.arange(0, 10)):
            if ccf_band < ccf.iloc[ccf_row] <= ccf_band+0.1:
                samples[samples_row, f_number] = n_ssms.iloc[ccf_row]

#Reading files from the folder                
file_number = 0
for f in os.listdir(indir):
    try:     
        dataset = pd.read_csv(indir + f, delimiter='\t', encoding='utf-8')
        create_matrix(dataset, file_number)
        file_number += 1
    except Exception as ex:
        logger.info (ex)
        assert ex

#Data frame is created from the sample data    
data = pd.DataFrame(data=samples)

#Method to get the explained variance to get the best n_component
def get_score(model, data, scorer=metrics.explained_variance_score):
    
    """ Estimate performance of the model on the data """
    prediction = model.inverse_transform(model.transform(data))
    
    return scorer(data, prediction)

n_cons = [2,4,6,8,10,12,14,15]
perfs_data = []

for n_con in n_cons:
    nmf = NMF(n_components=n_con, init='random').fit(data)
    perfs_data.append(get_score(nmf, data))

logger.info ('Performance Score: {}'.format(perfs_data))

fig = plt.figure(figsize=(12, 8))
plt.plot(n_cons,perfs_data, '-')
plt.title('Choosing components for the NMF Model')
plt.xlabel('Number of Components')
plt.ylabel('Explained Variance')
plt.savefig('best_component.png') # Save graph


#Method for tuning the right parameters for the model
def tune_NMF(X,n_components, init, max_iter, random_state, alpha, l1_ratio):
    nmf = NMF(n_components=n_components, init=init,  max_iter=max_iter, random_state=random_state,
              alpha=alpha, l1_ratio=l1_ratio)
    W = nmf.fit_transform(X)
    H = nmf.components_
    
    error = nmf.reconstruction_err_
    
    return error

#Parameters for the tuning
n_comps = [2,4,6,8,10]
init = ['random','nndsvd','nndsvda','nndsvdar']
max_iter = [2, 50 ,200]
random_state = [0,1]
alpha = [0, 0.5, 1]
l1_ratio = [0, 0.25, 0.75]


grid_search = []
search_parameter = [(a,b,c,d,e,f) for a in n_comps for b in init for c in max_iter 
                    for d in random_state for e in alpha for f in l1_ratio]

for paramts in search_parameter:
    n_comps, init, max_iter, random_state, alpha, l1_ratio = paramts
    
    grid_search.append((paramts,tune_NMF(samples, n_components=n_comps, init=init, 
                                         max_iter=max_iter, random_state=random_state, alpha=alpha, l1_ratio=l1_ratio)))
    
sort_grid_search = sorted(grid_search, key=lambda x:x[1], reverse= False)
logger.info ('Five Largest Error {}:'.format(sort_grid_search[-5:]))
logger.info ('Five Smallest Error {}:'.format(sort_grid_search[:5]))

ks = np.arange(1,11)

r_error= []

for k in ks:
        r_error.append(tune_NMF(samples,n_components=k,init = 'nndsvd', max_iter = 50, 
                                random_state = 1, alpha = 0, l1_ratio = 0))
plt.figure(figsize=(12, 8))
plt.title("Error Vs Number of Components")
plt.xlabel("Number of Components")
plt.ylabel("Errors")
plt.xticks(np.arange(0,11,2))
plt.plot(ks,r_error, '-')
plt.savefig('nmf_error.png') # Save graph

model = NMF(n_components=10, init='nndsvd', max_iter=50, random_state=1, alpha=0, l1_ratio=0)
W = model.fit_transform(data)
H = model.components_


# column names
t_names = ["Topic_" + str(i+1) for i in range(model.n_components)]

# index names
docnames = [str(ccf)[:3] for ccf in samples_CCF_values]

# Normalise the feature martix to get a unit value 
data_norm = normalize(W, norm='l1', axis=0)

#Features matrix
nmf_topic = pd.DataFrame(np.round(data_norm), columns=t_names, index=docnames)

#Coefficients matrix
coe_matrix = pd.DataFrame(np.round(H))

#Print out performance score
logger.info ('Performance Score : {}'.format(get_score(nmf, data)))

#Print out error
logger.info ('NMF Error: {}'.format(model.reconstruction_err_))


# Get dominant topic for each document
dominant_topic = np.argmax(nmf_topic.values, axis=1)
nmf_topic['dominant_topic'] = dominant_topic

df_topic_distribution = nmf_topic['dominant_topic'].value_counts().reset_index(name="Num Documents")

df_topic_distribution.columns = ['Topic Num', 'Num Documents']


#Plotting heatmap for the Coefficients matrix
fig, ax = plt.subplots()
fig.set_size_inches(12, 8.27)
sns.heatmap(coe_matrix.T, ax=ax)
ax.set(ylabel='Samples', xlabel='Topics')
ax.set_title('Heatmap Coefficients Matrix')
fig.savefig('coe_matrix.png')


# Distributions for each topic - CCF band
fig, ax = plt.subplots(5, 1, figsize=(12, 8), sharex=True)
for i, k in enumerate([0, 3, 5, 7, 9]):
    ax[i].stem(np.array(nmf_topic)[k, :model.n_components], linefmt='r-.', markerfmt='bo', basefmt='w-')
    ax[i].set_xlim(-1, model.n_components)
    ax[i].set_ylim(0, 1)
    ax[i].set_yticks(samples_CCF_values)
    ax[i].set_ylabel("Exp. Profile")
    ax[i].set_title("Doc {}".format(k+1))
    ax[i].grid(alpha=0.2)

ax[4].set_xlabel("Topic")
plt.tight_layout()
plt.savefig('nmf_topics_results.png') # Save graph
