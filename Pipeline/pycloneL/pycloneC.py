# --------------------------------------------------------
# Pyclone with R
# Written by Liu Mingfeng
# --------------------------------------------------------

import os
import subprocess
import logging
import time
from pycloneL.convertInputPyclone import start_run

class pyclone_class():
    
    def __init__(self,path,burnIn,num_iter):
        self._path = path
        #self._samId = str(samId)
        self._burnIn = burnIn
        self._num_iter = num_iter
        #self._sampleName = "Sim"+self._samId
        #self._workdir = self._path + "/inputTmp/pyclone"
        #self._resultFolder = self._path + "/results/"
        self._resultFolder = "." + "/results/"
        
        assert os.path.exists(self._path), 'path does not exist: {}.'.format(self._path)
        assert burnIn < num_iter, 'burnIn should be smaller than number of iteration.'
    
    
    def start_pyclone(self,purity,prefix):
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('pyclone_logger')
        
        rPath = self._path+"/pycloneL/"
        workdir = self._path + "/inputTmp/"+ prefix +"/pyclone"
        
        start_run(prefix)
        logger.info("start generate mutations file")
        #run R script for generating mutations file
        tmp1 = subprocess.call("Rscript {:s}generateMutationsFile.R {:s} {:s} {:s} {:s} {:d} {:d} {:f}".format(rPath,self._path,prefix,prefix,self._resultFolder,self._burnIn,self._num_iter,float(purity)),shell=True)
        tmp2 = subprocess.call("PyClone build_mutations_file --in_file {:s}/pyclone_data.tsv --out_file {:s}/pyclone_mutations.yaml --prior parental_copy_number".format(workdir,workdir),shell=True)
        assert tmp1+tmp2 == 0, 'generateMutationsFile failed\n\n'
        logger.info("start generate config file")
        #run R script for generating configure file
        tmp = subprocess.call("Rscript {:s}generateConfigFile.R".format(rPath),shell=True)
        assert tmp == 0, 'generate configure file failed\n\n'
        logger.info("start analysis")
        btime = time.time()
        #run pyclone for analysing
        tmp = subprocess.call("PyClone run_analysis --config_file {:s}/pyclone_configure.yaml --seed 1234".format(workdir),shell=True)
        assert tmp ==0,'analysis process failed\n\n'
        logger.info("posting result")
        #run R script for post Process
        tmp = subprocess.call("Rscript {:s}postProcess.R".format(rPath),shell=True)
        assert tmp == 0, 'post process failed\n\n'
        atime = time.time() - btime
        logger.info("Time to complete (s): {}".format(atime))
        with open("results/pyclone/runtimes.tsv",'a+') as f:
            f.write(prefix+'\t'+str(atime)+'\n')
        logger.info("process finished,results have been store in {:s}pyclone/{:s}".format(self._resultFolder, prefix))

