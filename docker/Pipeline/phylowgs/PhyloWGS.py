# Author:
# Lukas Rubikas

from __future__ import print_function

import pandas as pd
import numpy as np

from pandas import DataFrame

import subprocess
import os

import logging
import time

import json
import gzip
import zipfile

import sys # for running separately from the pipeline.py
import argparse # for running separately from the pipeline.py

logger = logging.getLogger('phylowgs_logger')

class PhyloWGS():

	def __init__(self, num_cores, burnin_samples, mcmc_samples, mh_iterations, random_seed=None):
		self.PATH_SUBCLONAL = './results/phylowgs/subclonal_structure'
		self.PATH_MUTASS = './results/phylowgs/mutation_assignments'
		self.random_seed = random_seed
		self.num_cores = num_cores
		self.burnin_samples = burnin_samples
		self.mcmc_samples = mcmc_samples
		self.mh_iterations = mh_iterations

	def _execute_bash_script(self, sample, purity, path_to_data):
		logger.info("***************************************")
		logger.info("Starting a PhyloWGS run for sample: {} with purity: {}".format(sample, purity))

		command_line = './phylowgs/run_phylowgs.sh -s {} -p {} -c {} -b {} -m {} -i {}'.format(
							sample, 
							purity,
							self.num_cores,
							self.burnin_samples,
							self.mcmc_samples,
							self.mh_iterations
							)

		if path_to_data is not None:
			command_line += ' -d {}'.format(path_to_data)

		logger.info("Executing bash script:")
		logger.info(command_line)

		_time = time.time()
		subprocess.call(command_line, shell=True)
		total_time = time.time() - _time

		return total_time

	def _write_runtime_to_file(self, sample, time_to_complete, directory='./results/phylowgs/'):
		if not os.path.exists(directory):
			os.makedirs(directory)
		with open(directory + "runtimes.tsv",'a+') as f:
			f.write(sample + '\t' + str(time_to_complete) + '\n')

	def _obtain_subclonal_structure(self, trees, purity):
		chain_likelihoods = np.array([(i, trees[str(i)]['llh']) for i in trees], 
									 dtype=[('chain', (np.str_, 10)), ('llh', np.float128)])
		maxllh_chain = chain_likelihoods[chain_likelihoods['llh'].argmax()]
		pop_dict = trees[maxllh_chain['chain']]['populations']
		relevant_data = [[cl, pop_dict[cl]['num_ssms'], pop_dict[cl]['cellular_prevalence'][0], pop_dict[cl]['cellular_prevalence'][0] / purity] for cl in pop_dict if pop_dict[cl]['num_ssms']> 0]
		return maxllh_chain, relevant_data

	def _obtain_mutation_assignments(self, mutass, ccfs, ssms):
		mut_ccf = lambda mut : [ccfs[mut[1]]]
		mut = np.array([[mut_id, cluster] for cluster in mutass for mut_id in mutass[cluster]['ssms']], dtype=np.str_)
		mut = np.hstack((mut,map(mut_ccf, mut))) 
		mut = mut[mut[:,0].argsort()] # sort for stacking
		ssms = ssms[ssms[:,0].argsort()] #sort for stacking   
		relevant_data = np.concatenate((ssms[:,[0,3,4]], mut[:,1:]), axis=1)    
		return relevant_data


	def _convert_to_common_format(self, sample, purity, phylowgs_res_dir='./phylowgs/phylowgs.results'):

		if not os.path.exists(self.PATH_SUBCLONAL):        
			os.makedirs(self.PATH_SUBCLONAL)
			
		if not os.path.exists(self.PATH_MUTASS):
			os.makedirs(self.PATH_MUTASS)

		results_folder = "{}.results".format(sample)
		failed_sample = False

		try:
			start_time = time.time()
			
			with gzip.open("{}/{}/test_results/{}.summ.json.gz".format(phylowgs_res_dir, results_folder, sample), 
				"rb") as f:
				d = json.loads(f.read().decode("ascii"))
				maxllh_chain, subclonal_structure = self._obtain_subclonal_structure(trees=d['trees'], purity=purity)

			with zipfile.ZipFile("{}/{}/test_results/{}.mutass.zip".format(phylowgs_res_dir, results_folder, sample), 
				"r") as f:
				d = json.loads(f.read("{}.json".format(maxllh_chain['chain'])))
				ccfs = {cl[0]: cl[-1] for cl in subclonal_structure}
				chrom_pos = lambda cp: [int(cp[1].split('_')[0]), int(cp[1].split('_')[1])]

				ssms = np.loadtxt('{}/{}/ssm_data.txt'.format(phylowgs_res_dir, results_folder), 
								  delimiter='\t', skiprows=1, usecols=(0,1), dtype=(np.str_, np.str_)) 
				ssms = np.concatenate((ssms,map(chrom_pos, ssms)), axis=1)
				ssms = np.hstack((np.arange(ssms.shape[0])[:, np.newaxis], ssms))        
				mutation_assignments = self._obtain_mutation_assignments(mutass=d['mut_assignments'], ccfs=ccfs, ssms=ssms)


			ss_df = DataFrame(data=subclonal_structure, columns=['cluster', 'n_ssms', 'proportion', 'ccf'])
			ss_df['ccf'] = pd.to_numeric(ss_df['ccf']) # for allowing to calculate proportion
			ss_df['proportion'] = ss_df['ccf'] * purity

			ma_df = DataFrame(data=mutation_assignments, columns=['id', 'chr', 'pos', 'cluster', 'ccf'], 
							  index=mutation_assignments[:,0])

			ma_df['id'] = pd.to_numeric(ma_df['id']) # for sorting (to retain the original mutation ordering)
			ma_df.sort_values(by=['id'], inplace=True)
			ma_df.drop(columns=['id'], inplace=True)

			ma_df['ccf'] = pd.to_numeric(ma_df['ccf']) # for allowing to calculate proportion
			ma_df['proportion'] = ma_df['ccf'] * purity

			ss_df.to_csv(path_or_buf="{}/{}_subclonal_structure.txt".format(self.PATH_SUBCLONAL, sample), sep='\t', index=False)
			ma_df.to_csv(path_or_buf="{}/{}_mutation_assignments.txt".format(self.PATH_MUTASS, sample), sep='\t', index=False)

			running_time = time.time() - start_time
			logger.info("Converted {} to common format in {} (s)".format(sample, running_time))
		
		except IOError as e:
			logger.warning(str(e))
			logger.warning("{} possibly consists of too many polyclonal trees, not enough to report a good posterior".format(sample))
			failed_sample = True

		return failed_sample
		

	def start_phylowgs(self, sample, purity, cooldown_min=None, path_to_data=None, path_to_phylowgs = './phylowgs', path_to_cmmn_fmt = './results'):

		time_to_complete = self._execute_bash_script(sample, purity, path_to_data)
		failed_sample = self._convert_to_common_format(sample, purity)

		self._write_runtime_to_file(sample, time_to_complete)

		logger.info("Time to complete: {} (s)".format(time_to_complete))
		logger.info("Results stored in: {}/phylowgs".format(path_to_cmmn_fmt))

		if failed_sample:
			logger.warning("Sample {} failed to be converted to the common format".format(failed_sample))
		

def __parse_args():
	parser = argparse.ArgumentParser(description='PhyloWGS pipeline')

	parser.add_argument(
		'--path-to-data', 
		dest = 'path_to_data',
		help = 'Path to data folder',
		default = '../data',
		type = str
		)

	parser.add_argument('--start-with', 
		dest = 'start_with',
		help = 'Index of a sample to start with (in case of interrupted operation',
		default = 0,
		type = int
		)

	parser.add_argument('--run-for',
		dest = 'no_of_samples_to_run',
		help = 'The number of samples to run, given the starting point',
		default = -1,
		type = int
		)

	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()
	return args


if __name__ == '__main__':
	logging.basicConfig(level=logging.INFO)
	args = __parse_args()

	samples = pd.read_csv("{}/pp_table.txt".format(args.path_to_data), sep='\t')
	
	print(samples)

	path_to_data = args.path_to_data
	start_with = args.start_with
	no_of_samples = args.no_of_samples_to_run


	# samples = samples.iloc[start_with:, :] if no_of_samples < 0 else samples[start_with : start_with + no_of_samples, :]

	phylowgs = PhyloWGS(
			num_cores = 2,
			burnin_samples = 10,
			mcmc_samples = 25,
			mh_iterations = 500,
			)

	for sample in samples.iloc[start_with:start_with+1,:].itertuples():	
		phylowgs.start_phylowgs(
			sample = getattr(sample, "sample"),
			purity = getattr(sample, "purity"),
			path_to_data = '../data',
			cooldown_min = 5
			)
