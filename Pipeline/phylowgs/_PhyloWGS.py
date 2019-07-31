from __future__ import print_function # Only Python 2.x
import pandas as pd
import subprocess
import sys
import time


samples = pd.read_csv("./team_project_data/Data/pp_table.txt", sep='\t')
subgroup = samples.iloc[94:300,:]

def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline(), ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

for i in xrange(subgroup.shape[0]):
    sample = subgroup.iloc[i]['sample']
    purity = subgroup.iloc[i]['purity']

    cmd = "./run_phylowgs.sh -s {} -p {} -c 2 -b 10 -m 25 -i 500".format(sample, purity)
    print(cmd)
    start_time = time.time()
    # for path in execute(cmd.split()):
    #     print(path, end="")
    subprocess.call(cmd.split())

    duration = time.time() - start_time
    print("Time to complete:", duration)
    with open("./phylowgs.runtimes/" + sample + ".runtime.txt", 'w+') as f:
        _duration = str(duration)
        f.seek(0)
        f.write(_duration)
    print("Cooldown for 5 minutes:")
    time.sleep(60 * 5)
