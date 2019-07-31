from pycloneL.pycloneC import pyclone_class
from ccube.ccubeC import Ccube
from phylowgs.PhyloWGS import PhyloWGS
from lda.ldaC import lda_class
from nmf.nmfC import nmf_class

import os
import argparse
import numpy as np
import sys
import logging
import time

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main_pipeline_logger')

def parse_args():
	currentPath = os.getcwd()
	
	parser = argparse.ArgumentParser(description='Pipeline')
						
	parser.add_argument('--path',
		dest = 'workPlace',
		help = 'Specifies working directory for analysis. All paths in the rest of the PyClone and Ccube files are relative to this',
		default = currentPath,
		type = str)
						
	parser.add_argument('--random-samples',
		dest = 'random_samples',
		help = 'Number of randomnly selected samples to run',
		default = 0, 
		type = int)
						
	parser.add_argument('--samples-to-run',
		dest = 'selected_samples',
		help = 'A newline-seperated file of explicitly-stated tumour sample names to be tested with Stage One tools (PyClone, Ccube, PhyloWGS',
		default = None,
		type = str)
						
	parser.add_argument('--phylowgs-burnin-samples',
		dest = 'phwgs_burnin_samples',
		help = 'Number of burn-in samples for PhyloWGS (default: 1000)',
		default = 1000,
		type = int)

	parser.add_argument('--phylowgs-mcmc-samples',
		dest = 'phwgs_mcmc_samples',
		help = 'Number of true MCMC samples for PhyloWGS (default: 2500)',
		default = 2500,
		type = int)

	parser.add_argument('--phylowgs-mh-iterations',
		dest = 'phwgs_mh_iterations',
		help = 'Number of Metropolis-Hastings iterations for PhyloWGS (default: 5000)',
		default = 5000, 
		type = int)

	parser.add_argument('--pyclone-burnin-samples',
		dest = 'pyclone_burnin_samples',
		help = 'Number of burn-in samples for PyClone (10 percent of total MCMC suggested in the official documentation, 50 used here)',
		default = 50,
		type = int)
					

	parser.add_argument('--pyclone-mcmc-iterations',
		dest = 'pyclone_mcmc_iterations',
		help = 'Number of MCMC iterations for PyClone (default: 500)',
		default = 500,
		type = int)
						
						
	parser.add_argument('--ccube-max-clusters',
		dest = 'ccube_maxcluster',
		help = 'Maximum number of clusters for Ccube (default: 6)',
		default = 6, 
		type = int)

	parser.add_argument('--ccube-vbem-max-iters',
		dest = 'ccube_vbem_max_iters',
		help = 'Number of VBEM iterations for Ccube (default: 1000)',
		default = 1000,
		type = int)
						
	parser.add_argument('--ccube-repeats',
		dest = 'ccube_repeat',
		help = 'Number of repeated Ccube runs for each candidate number of clusters (default: 1)',
		default = 1, 
		type = int)

	parser.add_argument('--ccube-random-seed',
		dest = 'random_seed',
		help = 'Random seed (used by Ccube), required to run the the tool deterministically', 
		default = None,
		type = int)
						
	parser.add_argument('--num-cores',
		dest = 'num_cores',
		help = 'Number of processor cores to be employed in computations concurrently (used by Ccube and PhyloWGS) (default: 1)',
		default = 1, 
		type = int)
						
	parser.add_argument('--run-pyclone', 
		action = 'store_true',
		help = 'Flag for running PyClone (default: False)',
		default = False)

	parser.add_argument('--run-ccube', 
		action = 'store_true',
		help = 'Flag for running Ccube (default: False)', 
		default = False)

	parser.add_argument('--run-phylowgs', 
		action = 'store_true',
		help = 'Flag for running PhyloWGS (default: False)', 
		default = False)
	
	parser.add_argument('--pyclone-delete-tmp-folder',
		dest='tmpfolder',
		help='Flag to delete temporary folder ("./inputTmp") containing configuration files generated while performing PyClone runs (default: False)',
		action='store_true',default=False)
		
	parser.add_argument('--lda-nmf-input',
		dest='results_folder',
		help='Input folder for LDA and/or NMF analysis (may be one of the result folders of Stage One tools)',
		default=None)
		
	parser.add_argument('--run-nmf',
		action = 'store_true',
		help = 'Flag for running NMF (default: False)',
		default = False)

	parser.add_argument('--run-lda',
		action = 'store_true',
		help = 'Flag for running LDA (default: False)',
		default = False)
	
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()
	return args

def readPtable(num):
	if num == 0:
		logger.info('Running all available samples')
	else:
		logger.info('{:d} samples are selected'.format(num))
	with open('data/pp_table.txt','r') as f:
		lines = f.readlines()
		alldata = []
		for line in lines:
			alldata.append(line.split('\t'))
		alldata = np.array(alldata)[1:,0:2]
	randomArray = np.random.randint(0,alldata.shape[0],num)
	if num == 0:
		return alldata
	else:
		return alldata[randomArray,:]

def readSetPtable(selected_sample):
	with open(selected_sample,'r') as f:
		lines = f.readlines()
		selectedFile = []
		for line in lines:
			line = line.split('\r')[0]
			selectedFile.append(line.split('\n')[0])
	num = np.array(selectedFile).shape[0]
	logger.info(str(num)+' samples are selected')
	with open('data/pp_table.txt','r') as f:
		lines = f.readlines()
		alldata = []
		for line in lines:
			tmp = line.split('\t')
			if tmp[0] in selectedFile:
				alldata.append(tmp)
	alldata = np.array(alldata)[:,0:2]
	return alldata

if __name__ =='__main__':
	
	args = parse_args()

	workdir = args.workPlace if args.workPlace is not None else '.'

	random_samples = args.random_samples
	phwgs_burnin_samples = args.phwgs_burnin_samples
	pyclone_burnin_samples = args.pyclone_burnin_samples
	ccube_vbem_max_iters = args.ccube_vbem_max_iters
	random_seed = args.random_seed if args.random_seed is not None else np.random.randint(10000)
	num_cores = args.num_cores
	phwgs_mcmc_samples = args.phwgs_mcmc_samples
	pyclone_mcmc_iterations = args.pyclone_mcmc_iterations
	phwgs_mh_iterations = args.phwgs_mh_iterations
	ccube_repeat = args.ccube_repeat
	cluster = args.ccube_maxcluster
	tmpfolder = args.tmpfolder
	
	results_folder = args.results_folder
	#select tool
	run_pyclone = args.run_pyclone
	run_ccube = args.run_ccube
	run_phylowgs = args.run_phylowgs
	run_nmf = args.run_nmf
	run_lda = args.run_lda
	
	#action = os.path.exists(self._path),
	selected_samples = args.selected_samples

	if selected_samples is not None:
		if os.path.isfile(selected_samples):
			alldata = readSetPtable(selected_samples)
		else:
			logger.error('Selected samples file does not exist')
	else:
		alldata = readPtable(random_samples)


	phyloWGS_pipeline = PhyloWGS(num_cores, phwgs_burnin_samples, phwgs_mcmc_samples, phwgs_mh_iterations) if run_phylowgs else None
	pycloneC = pyclone_class(workdir,pyclone_burnin_samples,pyclone_mcmc_iterations) if run_pyclone else None
	ccubeC = Ccube(workdir,random_seed,num_cores,ccube_repeat,cluster,ccube_vbem_max_iters) if run_ccube else None


	for singleData in alldata:
		logger.info(singleData[0]+' is selected')
		if run_pyclone:
			logger.info("***************************************")
			logger.info("Running PyClone for sample: " + singleData[0] + " with purity " + singleData[1])
			logger.info("The parameters for pyclone is\nburnin_sample:{}\nMCMC chain iteration:{}\nsample name:{}\npurity{}".format(pyclone_burnin_samples,pyclone_mcmc_iterations,singleData[0],singleData[1]))
			
			pycloneC.start_pyclone(singleData[1],singleData[0])

		if run_ccube:
			logger.info("***************************************")
			logger.info("Running Ccube for sample: " + singleData[0] + " with purity " + singleData[1])
			logger.info("The parameters for Ccube is\nnum_cores:{}\nrepeat:{}\nmax cluster:{}\nmax iterations:{}\nsample name:{}\npurity:{}\n".format(num_cores,ccube_repeat,cluster,ccube_vbem_max_iters,singleData[0],singleData[1]))
			
			ccubeC.start_ccube(singleData[0],singleData[1])

		if run_phylowgs:
			phyloWGS_pipeline.start_phylowgs(
				sample = singleData[0], 
				purity = float(singleData[1]), 
				cooldown_min = 0, 
				path_to_data = './data'
				)

	if tmpfolder:
		pass
	elif tmpfolder:
		os.system("rm -rf inputTmp")

#stage two
	if results_folder is not None and os.path.exists(results_folder):
		if run_nmf:
			nmf_call = nmf_class(results_folder)
			nmf_call.start_nmf()
		if run_lda:
			lda_call = lda_class(results_folder)
			lda_call.start_lda()
	elif results_folder is None:
		pass
	else:
		logger.error('Input folder for Stage Two (LDA, NMF) does not exists')
