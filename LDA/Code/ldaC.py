# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 12:00:14 2018

@author: Jesus Vera
"""

# Importing the libraries
import os
import numpy as np
import pandas as pd # Import/manage datasets
import matplotlib.pyplot as plt # Plots
from sklearn.decomposition import LatentDirichletAllocation # Library for LDA
from sklearn.model_selection import GridSearchCV # Library for searching best models

class lda_class():
    
    def __init__(self, path):
        self._path = path
        self._samples = np.zeros((10, 500), dtype=np.int32)
        self._samples_CCF_values = np.array(np.arange(0.1, 1.1, 0.1))

    """ Run the model """
    def start_lda(self):
        
        """ Create Sample Matrix
            Takes a dataset containing the column (proportion/CCF) to insert into Sample matrix
            and also takes the file number from the sample to read"""
        def createSamplesMatrix(ds, samples, sample_number):
            n_ssms = ds['n_ssms']
            ccf = ds['proportion']
            for i in np.arange(0, ccf.shape[0]):
                for j, z in zip(np.arange(0, 1, 0.1), np.arange(0, 10)):
                    if j < ccf.iloc[i] <= j+0.1:
                        samples[z, sample_number-1] = n_ssms.iloc[i]
                        
        """ Save the DataFrame to a CSV file"""
        def saveCSV(file_path, file_name, ds):       
            try:
                ds.to_csv(path_or_buf=file_path + file_name, sep='\t')
                return True
            except Exception as ex:
                print('Error: ' + ex)
                return False
                        
        # Load sample files
        print('*********************************')
        print('Creating samples matrix...')
        indir = self._path
        for f in os.listdir(indir):
            try:     
                dataset = pd.read_csv(indir + f, delimiter='\t', encoding='utf-8')
                point_pos = f[8:].find('.')
                sample_number = int(f[8:8+point_pos])
                createSamplesMatrix(dataset, self._samples, sample_number)
            except Exception as ex:
                print('Error: ' + ex)
                assert ex
        
        # Set the results folder
        # If the folder doesn't exist, create it
        outdir = '../results/lda/'
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
            print("Saving plot into Results forlder...")
            plt.savefig(outdir + 'mean_subclonal_500_samples.png')
        except Exception as ex:
            print("Error: ", ex)
            assert ex
        
        print('*********************************')
        print("Executing GridSearch to find best LDA model...")
        # Perform GridSearch for the best LDA model
        # Define Search Param
        search_params = {
                'n_components': [5, 6, 7, 8, 9], 
                'learning_decay': [0.5, 0.7, 0.9],
                'max_iter': [6, 7, 8, 9]
                }
        
        # Init the Model
        lda = LatentDirichletAllocation()
        
        # Init Grid Search Class
        model = GridSearchCV(lda, param_grid=search_params)
        
        # Do the Grid Search
        model.fit(self._samples)
        
        # Visualizing best topi model and its parameters
        # Best Model
        best_lda_model = model.best_estimator_
        
        # Model Parameters
        print("Best Model's Params: ", model.best_params_)
        
        # Log Likelihood Score
        print("Best Log Likelihood Score: ", model.best_score_)
        
        # Perplexity
        print("Model Perplexity: ", best_lda_model.perplexity(self._samples))
        
        # Get Log Likelyhoods from Grid Search Output
        n_topics = [5, 6, 7, 8, 9]
        log_likelyhoods_5 = [round(gscore.mean_validation_score) for gscore in model.grid_scores_ if gscore.parameters['learning_decay']==0.5 and gscore.parameters['max_iter'] == 9]
        log_likelyhoods_7 = [round(gscore.mean_validation_score) for gscore in model.grid_scores_ if gscore.parameters['learning_decay']==0.7 and gscore.parameters['max_iter'] == 9]
        log_likelyhoods_9 = [round(gscore.mean_validation_score) for gscore in model.grid_scores_ if gscore.parameters['learning_decay']==0.9 and gscore.parameters['max_iter'] == 9]
        
        # Comparison of log likelihoods from GridSearchCV
        print("Plotting the different log likelihoods for choosing the optimal LDA model...")
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
            print("Saving plot into Results forlder...")
            plt.savefig(outdir + 'optimal_lda_model.png') # Save graph
        except Exception as ex:
            print("Error: ", ex)
            assert ex
        
        # Create Document - Topic Matrix
        lda_output = best_lda_model.transform(self._samples)

        # column names
        topicnames = ["Topic" + str(i+1) for i in range(best_lda_model.n_components)]
        
        # index names
        docnames = [str(ccf)[:3] for ccf in self._samples_CCF_values]
        
        # Make the pandas dataframe
        df_document_topic = pd.DataFrame(np.round(lda_output, 2), columns=topicnames, index=docnames)
        
        # Get dominant topic for each document
        dominant_topic = np.argmax(df_document_topic.values, axis=1)
        df_document_topic['dominant_topic'] = dominant_topic
        
        df_topic_distribution = df_document_topic['dominant_topic'].value_counts().reset_index(name="Num Documents")
        df_topic_distribution.columns = ['Topic Num', 'Num Documents']
        
        is_saved = saveCSV(outdir, 'document_topic.csv', df_document_topic)
        if is_saved:
            print('Saving document-topic file into results folder...')
        else:
            print('There was a problem while saving the document-topic file.')
        
        # Distributions for each topic - CCF band
        print("Plotting five random topics results...")
        fig, ax = plt.subplots(5, 1, figsize=(12, 8), sharex=True)
        for i, k in enumerate([0, 3, 5, 7, 9]):
            ax[i].stem(np.array(df_document_topic)[k, :9], linefmt='r-.', markerfmt='ro', basefmt='w-')
            ax[i].set_xlim(-1, 9)
            ax[i].set_ylim(0, 1)
            ax[i].set_yticks(self._samples_CCF_values)
            ax[i].set_ylabel("Prob")
            ax[i].set_title("CCF {0:.1f}".format(self._samples_CCF_values[k]))
            ax[i].grid(alpha=0.2)
        ax[4].set_xlabel("Topic")
        plt.tight_layout()
        try:
            print("Saving plot into Results forlder...")
            plt.savefig(outdir + 'random_topics_results.png') # Save graph
        except Exception as ex:
            print("Error: ", ex)
            assert ex
        
        print("Process finished!")