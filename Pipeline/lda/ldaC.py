# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 12:00:14 2018

@author: Jesus Vera
"""

# Importing the libraries
import os
import logging
import numpy as np
import pandas as pd # Import/manage datasets
import matplotlib.pyplot as plt # Plots
import seaborn as sns # Heatmap
from sklearn.decomposition import LatentDirichletAllocation # Library for LDA
from sklearn.model_selection import GridSearchCV # Library for searching best models

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

class lda_class():
    
    def __init__(self, path):
        self._path = path
        # Count the number of files to analyse
        num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
        self._samples = np.zeros((10, num_files), dtype=np.int32)
        self._samples_CCF_values = np.array(np.arange(0.1, 1.1, 0.1))

    """ Run the model """
    def start_lda(self):

        """ Create Sample Matrix
            Takes a dataset containing the column (proportion/CCF) to insert into Sample matrix
            and a file number to keep track of files """
        def createSamplesMatrix(ds, samples, fn):
            n_ssms = ds['n_ssms']
            ccf = ds['proportion']
            for ccf_row in np.arange(0, ccf.shape[0]):
                for ccf_band, samples_row in zip(np.arange(0, 1, 0.1), np.arange(0, 10)):
                    if ccf_band < ccf.iloc[ccf_row] <= ccf_band+0.1:
                        samples[samples_row, fn] = n_ssms.iloc[ccf_row]

        """ Save the DataFrame to a CSV file"""
        def saveCSV(file_path, file_name, ds):       
            try:
                ds.to_csv(path_or_buf=file_path + file_name, sep='\t')
                return True
            except Exception as ex:
                print(ex)
                return False
                        
        # Load sample files
        logging.info('Creating samples matrix...')
        indir = self._path
        file_number = 0
        for f in os.listdir(indir):
            try:     
                dataset = pd.read_csv(indir + f, delimiter='\t', encoding='utf-8')
                createSamplesMatrix(dataset, self._samples, file_number)
                file_number += 1
            except Exception as ex:
                logging.error(ex)
                assert ex
        
        # Set the results folder
        # If the folder doesn't exist, create it
        outdir = 'results/plot/lda'
        if not os.path.exists(outdir):
                os.makedirs(outdir)
        
        print("Plotting Mean Subclonal Structure of 500 samples...")
        # Pre-process samples CCF for plotting (cancer cell fraction)
        samples_mean = np.mean(self._samples, axis=1) # Mean of CCFs by band (0.1, 0.2, ..., 1.0)
        
        # Mean of CCFs by band from 500 samples
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self._samples_CCF_values, samples_mean)
        ax.set_xticks(self._samples_CCF_values)
        ax.set_xlabel('CCF (cancer cell fraction)')
        ax.set_ylabel('Mean ocurrences')
        ax.set_title('Mean Subclonal Structure of 500 samples')
        try:
            print("Saving plot into Results folder...")
            plt.savefig(outdir + '/' + 'mean_subclonal_500_samples.png')
        except Exception as ex:
            logging.error(ex)
            assert ex
        
        logging.info("Executing GridSearch to find best LDA model...")
        # Perform GridSearch for the best LDA model
        # Define Search Param
        search_params = {
            'n_components': [9, 10, 11, 12], # take 10 topics
            'learning_decay': [0.5, 0.7, 0.9],
            'max_iter': [8, 9, 10, 11],
            'random_state': [2018]
        }
        
        # Init the Model
        lda = LatentDirichletAllocation()
        
        # Init Grid Search Class
        model = GridSearchCV(lda, param_grid=search_params)
        
        # Do the Grid Search
        samples = self._samples.T # Work with transposed
        model.fit(samples)
        
        # Visualizing best topi model and its parameters
        # Best Model
        best_lda_model = model.best_estimator_
        
        # Model Parameters
        logging.info("Best Model's Params: {}".format(model.best_params_))
        
        # Log Likelihood Score
        logging.info("Best Log Likelihood Score: {}".format(model.best_score_))
        
        # Perplexity
        logging.info("Model Perplexity: {}".format(best_lda_model.perplexity(samples)))
        
        # Get Log Likelyhoods from Grid Search Output
        n_topics = [9, 10, 10, 12]
        log_likelyhoods_5 = [round(gscore.mean_validation_score) for gscore in model.grid_scores_ if gscore.parameters['learning_decay']==0.5 and gscore.parameters['max_iter'] == best_lda_model.max_iter]
        log_likelyhoods_7 = [round(gscore.mean_validation_score) for gscore in model.grid_scores_ if gscore.parameters['learning_decay']==0.7 and gscore.parameters['max_iter'] == best_lda_model.max_iter]
        log_likelyhoods_9 = [round(gscore.mean_validation_score) for gscore in model.grid_scores_ if gscore.parameters['learning_decay']==0.9 and gscore.parameters['max_iter'] == best_lda_model.max_iter]
        
        # Comparison of log likelihoods from GridSearchCV
        logging.info("Plotting the different log-likelihoods for choosing the optimal LDA model...")
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(n_topics, log_likelyhoods_5, label='0.5')
        ax.plot(n_topics, log_likelyhoods_7, label='0.7')
        ax.plot(n_topics, log_likelyhoods_9, label='0.9') # Best one
        plt.xlabel("Num Topics")
        plt.ylabel("Log Likelyhood Scores")
        plt.title("Choosing Optimal LDA Model")
        plt.legend(title='Learning decay', loc='best')
        try:
            logging.info("Saving plot into Results folder...")
            plt.savefig(outdir + '/' + 'optimal_lda_model.png') # Save graph
        except Exception as ex:
            logging.error(ex)
            assert ex
        
        # Create Document - Topic Matrix
        lda_output = best_lda_model.transform(samples)
        
        # Plot heatmap with lda_output   
        logging.info("Plotting heatmap of document-topic matrix...")
        fig, ax = plt.subplots()
        fig.set_size_inches(12, 8)
        sns.heatmap(lda_output, ax=ax)
        ax.set(ylabel='Samples', xlabel='Topics')
        ax.set_title('Heatmap Documents-Topics Matrix')
		
        try:
            logging.info ("Saving plot into Results folder...")
            fig.savefig(outdir + '/' + 'heatmap_doc_topics.png')
        except Exception as ex:
            logging.info ("Error: ", ex)
            assert ex

        # column names
        topicnames = ["Topic" + str(i+1) for i in range(best_lda_model.n_components)]
        
        # index names
        docnames = [str(doc) for doc in np.arange(1, len(lda_output)+1)]
        
        # Make the pandas dataframe
        df_document_topic = pd.DataFrame(np.round(lda_output, 2), columns=topicnames, index=docnames)
        
        is_saved = saveCSV(outdir+'/', 'document_topic.csv', df_document_topic)
        if is_saved:
            logging.info('Saving document-topic file into Results folder...')
        else:
            logging.error('There was a problem while saving the document-topic file.')
        
        # Distributions for each topic - CCF band
        print("Plotting five random topics results...")
        fig, ax = plt.subplots(5, 1, figsize=(12, 8), sharex=True)
        for i, k in enumerate([0, 3, 5, 7, 9]):
            ax[i].stem(np.array(df_document_topic)[k, :best_lda_model.n_components], linefmt='r-.', markerfmt='ro', basefmt='w-')
            ax[i].set_xlim(-1, best_lda_model.n_components)
            ax[i].set_ylim(0, 1)
            ax[i].set_xticks(np.arange(-1, best_lda_model.n_components))
            ax[i].set_yticks(self._samples_CCF_values)
            ax[i].set_ylabel("Prob")
            ax[i].set_title("Doc {}".format(k+1))
            ax[i].grid(alpha=0.2)
        ax[4].set_xlabel("Topic")
        plt.tight_layout()
        try:
            logging.info("Saving plot into Results folder...")
            plt.savefig(outdir + '/' + 'random_topics_results.png') # Save graph
        except Exception as ex:
            logging.error(ex)
            assert ex
        
        logging.info("LDA process finished!")
