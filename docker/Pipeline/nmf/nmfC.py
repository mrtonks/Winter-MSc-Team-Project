# -*- coding: utf-8 -*-
"""
Created on Mon 26 Nov 20:18:47 2018

@author: Razak nart
"""

# Importing all required libraries
import os
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from time import time
import logging
import seaborn as sns
from sklearn import metrics
from sklearn.decomposition import NMF # Library for NMF
from sklearn.preprocessing import normalize



# Main class for the NMF model 
class nmf_class():
    
    def __init__(self, path):
        self._path = path
        
        self._indir = self._path
        #Count the number of files in Subclonal_Structure
        num_files = len([f for f in os.listdir(self._indir) if os.path.isfile(os.path.join(self._indir, f))])
         
        self._samples = np.zeros((10, num_files), dtype=np.int32)
        self._samples_CCF_values = np.array(np.arange(0.1, 1.1, 0.1))
        
        
    """ Run the model """
    def start_nmf(self):
        
        # Defining logger for nmf 
        logger = logging.getLogger('nmf_logger')
        logging.basicConfig(level=logging.INFO) 

              
        """ Create Sample Matrix. Takes a dataset containing the column
        (proportion/CCF) to insert into Sample matrix and a file number to keep track of files """
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
                logger.info ('Error: ' + ex)
                return False
            
            
        # ************ Load sample files ***************
        
        logger.info ('Creating samples matrix...')
        indir = self._path
        file_number = 0
        for f in os.listdir(self._indir):
            try:     
                dataset = pd.read_csv(indir + f, delimiter='\t', encoding='utf-8')
                createSamplesMatrix(dataset, self._samples, sample_number)
                file_number += 1
            except Exception as ex:
                logger.error(ex)
                assert ex
                
        # Set the results folder
        # If the folder doesn't exist, create it
        outdir = 'results/plot/nmf/'
        if not os.path.exists(outdir):
                os.makedirs(outdir)
                
        
        print("Finding the best model for NMF............")
        
        #Method to get the score on the data
        def get_score(model, data, scorer=metrics.explained_variance_score):
            
            """ Estimate performance of the model on the data """
            prediction = model.inverse_transform(model.transform(data))
            
            return scorer(data, prediction)


        n_cons = [2,4,6,8,10,12,14,15]
        perfs_data = []

        for n_con in n_cons:
            nmf = NMF(n_components=n_con, init='random').fit(self._samples)
            perfs_data.append(get_score(nmf, self._samples))
            
        logger.info ('Performance Score: {}'.format(perfs_data))

        plt.figure(figsize=(12, 8))
        plt.plot(n_cons,perfs_data, '-')
        plt.title('Choosing components for the NMF Model')
        plt.xlabel('Number of Components')
        plt.ylabel('Cumulative Explained Variance')
        
        try:
            logger.info ("Saving plot into Results folder...")
            plt.savefig(outdir + 'best_component.png') # Save graph
        except Exception as ex:
            logger.info ("Error: ", ex)
            assert ex
            
        
        logger.info ("Tuning NMF to find the best parametres with less error............")
        
        #Method to tune the NMF model to find the best parametres
        def tune_NMF(X,n_components, init, max_iter, random_state, alpha, l1_ratio):
            nmf = NMF(n_components=n_components, init=init,  max_iter=max_iter,
                      random_state=random_state, alpha=alpha, l1_ratio=l1_ratio)
           
            W = nmf.fit_transform(X)
            H = nmf.components_
            
            error = nmf.reconstruction_err_
            
            return error
        
        # Define Search Param
        n_comps = [2,4,6,8,10]
        init = ['random','nndsvd','nndsvda','nndsvdar']
        max_iter = [2, 50 ,200]
        random_state = [0,1]
        alpha = [0, 1]
        l1_ratio = [0, 1]

        grid_search = []
        search_parameter = [(a,b,c,d,e,f) for a in n_comps for b in init for c in max_iter 
                    for d in random_state for e in alpha for f in l1_ratio]

        for paramts in search_parameter:
            n_comps, init, max_iter, random_state, alpha, l1_ratio = paramts
    
        grid_search.append((paramts,tune_NMF(self._samples, n_components=n_comps, 
                                             init=init, max_iter=max_iter, random_state=random_state, 
                                             alpha=alpha, l1_ratio=l1_ratio)))
    
        sort_grid_search = sorted(grid_search, key=lambda x:x[1], reverse= False)
        
        #Print out the five smallest and the largest error
        logger.info ('Five Largest Error {}:'.format(sort_grid_search[-5:]))
        logger.info ('Five Smallest Error {}:'.format(sort_grid_search[:5]))
        
        
        logger.info ('Plotting the error and the number of components')
        ks = np.arange(1,11)
        r_error= []


        for k in ks:
            r_error.append(tune_NMF(self._samples,n_components=k,init = 'nndsvd', max_iter = 50, 
                                random_state = 1, alpha = 0, l1_ratio = 0))

        plt.figure(figsize=(12, 8))
        plt.title("Error Vs Number of Components")
        plt.xlabel("Number of Components")
        plt.ylabel("Errors")
        plt.xticks(np.arange(0,11,2))
        plt.plot(ks,r_error, '-')
        
        try:
            logger.info ("Saving plot into Results folder...")
            plt.savefig(outdir + 'nmf_error.png') # Save graph
        except Exception as ex:
            logger.info ("Error: ", ex)
            assert ex
            
        
        logger.info ('Final model using the best parametres')
        tm = time()
        model = NMF(n_components=10, max_iter=50, alpha=0, l1_ratio=0)
        W = model.fit_transform(self._samples)
        H = model.components_

        logger.info ("Done in %0.3fs." % (time() - tm))

        # column names
        t_names = ["Topic_" + str(i+1) for i in range(model.n_components)]

        # index names
        docnames = [str(ccf)[:3] for ccf in self._samples_CCF_values]
        
        # Normalise the feature martix to get a unit value 
        data_norm = normalize(W, norm='l1', axis=0)

        # Make the pandas dataframe
        #Features matrix
        nmf_topic = pd.DataFrame(np.round(data_norm), columns=t_names, index=docnames)

        #Coefficients matrix
        coe_matrix = pd.DataFrame(np.round(H))       
               
        #Print out performance score
        logger.info ('Performance Score : {}'.format(get_score(model, self._samples)))

        #Print out error
        logger.info ('NMF Error: {}'.format(model.reconstruction_err_))
        
        # Get dominant topic for each document
        dominant_topic = np.argmax(nmf_topic.values, axis=1)
        nmf_topic['dominant_topic'] = dominant_topic

        df_topic_distribution = nmf_topic['dominant_topic'].value_counts().reset_index(name="Num Documents")

        df_topic_distribution.columns = ['Topic Num', 'Num Documents']
        
        is_saved = saveCSV(outdir, 'nmf_topic.csv', nmf_topic)
        if is_saved:
            logger.info ('Saving nmf-topic file into Results folder...')
        else:
            logger.info ('There was a problem while saving the nmf-topic file.')
            
            
        #Plotting heatmap for the Coefficients matrix
        fig, ax = plt.subplots()
        fig.set_size_inches(12, 8.27)
        sns.heatmap(coe_matrix.T, ax=ax)
        ax.set(ylabel='Samples', xlabel='Topics')
        ax.set_title('Heatmap Coefficients Matrix')
        
        try:
            logger.info ("Saving plot into Results folder...")
            fig.savefig(outdir + 'coe_matrix.png') # Save graph
        except Exception as ex:
            logger.info ("Error: ", ex)
            assert ex

        # Distributions for each topic - CCF band
        fig, ax = plt.subplots(10, 1, figsize=(12, 14), sharex=True)
        for i, k in enumerate([0, 3, 5, 7, 9]):
            ax[i].stem(np.array(nmf_topic)[k, :model.n_components], linefmt='r-.', markerfmt='bo', basefmt='w-')
            ax[i].set_xlim(-1, model.n_components)
            ax[i].set_ylim(0, 1)
            ax[i].set_yticks(self._samples_CCF_values)
            ax[i].set_ylabel("Exp. Profile")
            ax[i].set_title("Doc {}".format(k+1))
            ax[i].grid(alpha=0.2)

        ax[9].set_xlabel("Topic")
        plt.tight_layout()
        
        try:
            logger.info ("Saving plot into Results folder...")
            plt.savefig(outdir + 'nmf_topic_result.png') # Save graph
        except Exception as ex:
            logger.info ("Error: ", ex)
            assert ex
        
        
        print("NMF process finished!")
        
#if __name__ == '__main__':
    
        
        
                
        
            
            
            
                        

