# Importing the libraries
import os
import numpy as np
import pandas as pd # Import/manage datasets

# Plots
import matplotlib.pyplot as plt 

from pprint import pprint

#from sklearn.feature_extraction.text import CountVectorizer
from sklearn.decomposition import LatentDirichletAllocation # Library for LDA
from sklearn.model_selection import GridSearchCV # Library for searching best models

# Global variable for text corpora
samples = np.zeros((10, 500), dtype=np.int32)
samples_CCF_values = np.array(np.arange(0.1, 1.1, 0.1))

""" Create Sample Matrix
    Takes a dataset containing the column (proportion/CCF) to insert into Sample matrix
    and a file number to keep track of files """
def createSamplesMatrix(ds, fn):
    n_ssms = ds['n_ssms']
    ccf = ds['proportion']
    for ccf_row in np.arange(0, ccf.shape[0]):
        for ccf_band, samples_row in zip(np.arange(0, 1, 0.1), np.arange(0, 10)):
            if ccf_band < ccf.iloc[ccf_row] <= ccf_band+0.1:
                samples[samples_row, fn] = n_ssms.iloc[ccf_row]

""" Remove features that have a sparcity of non-zeros below 10% """
def removeFeatures(ds, ds_ccf):
    rows_2_del = []
    # Compute sparcity = percentage of non-zero cells
    for i in range(len(ds)):
        sparsity = np.sum(np.where(ds[i, :] > 0, 1, 0)) / 500
        if sparsity < 0.1:
            rows_2_del.extend([i])
        else:
            continue
    ds = np.delete(ds, rows_2_del, axis=0)
    ds_ccf = np.delete(ds_ccf, rows_2_del, axis=0)
    return ds, ds_ccf

""" Save the DataFrame to a CSV file"""
def saveCSV(file_path, file_name, ds):       
    try:
        ds.to_csv(path_or_buf=file_path + file_name, sep='\t')
        return True
    except Exception as ex:
        print 'Error: ' + ex
        return False


# Load sample files        
indir = '../../Data/Phase_two/Subclonal_Structure/'
file_number = 0
for f in os.listdir(indir):
    try:     
        dataset = pd.read_csv(indir + f, delimiter='\t', encoding='utf-8')
        createSamplesMatrix(dataset, file_number)
        file_number += 1
    except Exception as ex:
        print ex
        assert ex

# Set the results folder
# If the folder doesn't exist, create it
outdir = '../Results/'
if not os.path.exists(outdir):
        os.makedirs(outdir)

#samples, samples_CCF_values = removeFeatures(samples, samples_CCF_values)


# Pre-process samples CCF for plotting (cancer cell fraction)
samples_mean = np.mean(samples, axis=1) # Mean of CCFs by band (0.1, 0.2, ..., 1.0)

# mean of CCFs by band from 500 samples
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(1, 1, 1)
ax.plot(samples_CCF_values, samples_mean)
ax.set_xticks(samples_CCF_values)
ax.set_xlabel('CCF (cancer cell fraction)')
ax.set_ylabel('Mean ocurrences')
ax.set_title('Mean Subclonal Structure of 500 samples')
plt.savefig(outdir + 'mean_subclonal_500_samples.png')


# Build LDA model
lda_model = LatentDirichletAllocation(n_components=10,      # Number or topics
                                  max_iter=10,              # Max learning iterations
                                  random_state=100,         # Random state (seed)
                                  learning_method='online',
                                  batch_size=128,           # No of docs in each iter
                                  evaluate_every=-1,        # Compute perplexity every n iters
                                  n_jobs=-1)                # Use all available CPUs

lda_output = lda_model.fit_transform(samples)
print(lda_model)

# Diagnose model performance with perplexity and log-likelihood
# Log Likelyhood: Higher the better
print "Log Likelihood: ", lda_model.score(samples)

# Perplexity: Lower the better. Perplexity = exp(-1. * log-likelihood per word)
print("Perplexity: ", lda_model.perplexity(samples))

# See model parameters
pprint(lda_model.get_params())


# Perform GridSearch for the best LDA model
# Define Search Param
search_params = {
        'n_components': [6, 7, 8, 9], # take 10 topics
        'learning_decay': [0.5, 0.7, 0.9],
        'max_iter': [6, 7, 8, 9],
        'random_state': [2018]
        }

# Init the Model
lda = LatentDirichletAllocation()

# Init Grid Search Class
model = GridSearchCV(lda, param_grid=search_params)

# Do the Grid Search
model.fit(samples)

# Visualizing best topi model and its parameters
# Best Model
best_lda_model = model.best_estimator_

# Model Parameters
print "Best Model's Params: ", model.best_params_

# Log Likelihood Score
print "Best Log Likelihood Score: ", model.best_score_

# Perplexity
print "Model Perplexity: ", best_lda_model.perplexity(samples)



# Get Log Likelyhoods from Grid Search Output
n_topics = [6, 7, 8, 9]
log_likelyhoods_5 = [round(gscore.mean_validation_score) for gscore in model.grid_scores_ if gscore.parameters['learning_decay']==0.5 and gscore.parameters['max_iter'] == 9]
log_likelyhoods_7 = [round(gscore.mean_validation_score) for gscore in model.grid_scores_ if gscore.parameters['learning_decay']==0.7 and gscore.parameters['max_iter'] == 9]
log_likelyhoods_9 = [round(gscore.mean_validation_score) for gscore in model.grid_scores_ if gscore.parameters['learning_decay']==0.9 and gscore.parameters['max_iter'] == 9]

# Comparison of log likelihoods from GridSearchCV
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(1, 1, 1)
ax.plot(n_topics, log_likelyhoods_5, label='0.5')
ax.plot(n_topics, log_likelyhoods_7, label='0.7')
ax.plot(n_topics, log_likelyhoods_9, label='0.9') # Best one
plt.xlabel("Num Topics")
plt.ylabel("Log Likelyhood Scores")
plt.title("Choosing Optimal LDA Model")
plt.legend(title='Learning decay', loc='best')
plt.savefig(outdir + 'optimal_lda_model.png') # Save graph

# Create Document - Topic Matrix
lda_output = best_lda_model.transform(samples)

# Create Topic - Terms matrix
topic_terms = best_lda_model.components_

# column names
topicnames = ["Topic" + str(i+1) for i in range(best_lda_model.n_components)]

# index names
docnames = [str(ccf)[:3] for ccf in samples_CCF_values]

# Make the pandas dataframe
df_document_topic = pd.DataFrame(np.round(lda_output, 2), columns=topicnames, index=docnames)

# Distributions for each topic - CCF band
fig, ax = plt.subplots(5, 1, figsize=(12, 8), sharex=True)
for i, k in enumerate([0, 3, 5, 7, 9]):
    ax[i].stem(np.array(df_document_topic)[k, :9], linefmt='r-.', markerfmt='ro', basefmt='w-')
    ax[i].set_xlim(-1, 9)
    ax[i].set_ylim(0, 1)
    ax[i].set_xticks(np.arange(-1, 9))
    ax[i].set_yticks(samples_CCF_values)
    ax[i].set_ylabel("Prob")
    ax[i].set_title("CCF {0:.1f}".format(samples_CCF_values[k]))
    ax[i].grid(alpha=0.2)

ax[4].set_xlabel("Topic")
plt.tight_layout()
plt.savefig(outdir + 'random_topics_results.png') # Save graph

is_saved = saveCSV(outdir, 'document_topic.csv', df_document_topic)
if is_saved:
    print('Document-topic file created sucessfully!')
else:
    print('There was a problem while running the LDA model')
# bic for topic selection
