import os
from ccube.ConvertInput import start_run

class Ccube:
    def __init__(self,path,seed,core,repeat,cluster,maxiter):
        self._p=path
        #self._id=sampleid
        self._s=seed
        self._c = core
        self._r = repeat
        self._cl=cluster
        self._m=maxiter
    
    def start_ccube(self,prefix,purity):
        currentPath = os.getcwd()
        print("__________*************___________")
        print("Convert the data of VCF and segments to the format of input")
        start_run(prefix,float(purity))
        print("__________*************___________")
        print("Running Ccubeoutput Rcode")
        os.system("Rscript %s/ccube/twst.R %s %s %s %s %s %s %s" %(currentPath,self._p,prefix,self._s,self._c,self._r,self._cl,self._m))
