import os
import numpy as np
def intergret(firstfilename,secondfilename):
    firsttmp = firstfilename.split('/')[-1]
    firsttmp = firsttmp.split('.')[0]
    secondtmp = secondfilename.split('/')[-1]
    secondtmp = secondtmp.split('.')[0]
    if firsttmp == secondtmp:
        outputfile = firsttmp + '.tsv'
        fileindex = int(firsttmp.split('_')[-1])
        fileS = '{0:03d}'.format(fileindex)
    with open(firstfilename,'r') as f:
        seg_lines = f.readlines()
        seg_content = []
        for line in seg_lines:
            tmp = line.split('\t')
            tmp.pop()
            seg_content.append(tmp)
        seg_content = np.array(seg_content)
        seg_content = seg_content[:,0:7]
        slabel = seg_content[0,:]
        seg_content = seg_content[seg_content[:,6]=='1']
        seg_content = np.insert(seg_content,0,values=slabel,axis=0)
    #:output of VCF
    with open(secondfilename,'r') as f:
        vcf_lines = f.readlines()
        k = 0
        for line in vcf_lines:
            if k<5:
                k = k+1
                continue
            tmp = line.split('\t')
            if k == 5:
                k=6
                vcf_content = np.array(tmp)
                continue
            vcf_content = np.vstack((vcf_content,tmp))
        #process last column
        last_line = vcf_content[:,-1]
        vcf_content = np.delete(vcf_content,-1,axis=1)
        k = 0
        for lastL in last_line:
            if k == 0:
                k = k+1
                last_column = lastL.split(';')[0:2]
                fe = last_column[0].split('=')[1]
                se = last_column[1].split('=')[1]
                lastf_column = np.array([fe,se])
                continue
            last_column = lastL.split(';')[0:2]
            fe = last_column[0].split('=')[1]
            se = last_column[1].split('=')[1]
            lastf_column = np.vstack((lastf_column,np.array([fe,se])))
        vcf_content = np.hstack((vcf_content,lastf_column))
    label = ['chr','pos','id','ref','alt','qual','filter','alt','ref']
    vcf_content = np.vstack((label,vcf_content))
    #output: rows of data
    k = 0
    for item in seg_content[1:]:
        add_item = np.array([item[3],item[4],item[5],'Sim{:s}'.format(fileS)])
        minrange = int(item[1])
        maxrange = int(item[2])
        chrindex = item[0]
        grouptmp = vcf_content[vcf_content[:,0] == chrindex]
        for subitem in grouptmp:
            final_id = 'Sim{:s}:BB:{:s}:{:s}'.format(fileS,chrindex,subitem[1])
            if int(subitem[1]) < maxrange and int(subitem[1]) > minrange:
                if k == 0:
                    variant_v = int(subitem[7]) / (int(subitem[7]) + int(subitem[8]))
                    final_result = np.hstack(([final_id,subitem[8],subitem[7]],add_item))
                    final_result = np.hstack((final_result,variant_v))
                    final_result = np.hstack((final_result,'BB'))
                    k = 1
                    continue
                variant_v = int(subitem[7]) / (int(subitem[7]) + int(subitem[8]))
                single_item = np.hstack(([final_id,subitem[8],subitem[7]],add_item))
                single_item = np.hstack((single_item,variant_v))
                single_item = np.hstack((single_item,'BB'))
                final_result = np.vstack((final_result,single_item))
            else:
                continue

    label = ['mutation_id','ref_counts','var_counts','normal_cn','minor_cn','major_cn','variant_case','variant_freq','genotype']
    final_result = np.vstack((label,final_result))
    #write file
    fshape = final_result.shape
    tmp = final_result.reshape(fshape[0]*fshape[1])
    txt = ''
    if not os.path.exists('inputTmp/'+firsttmp):
        os.makedirs('inputTmp/'+firsttmp)
    with open('inputTmp/'+firsttmp+'/'+outputfile,'w') as f:
        k = 0
        for text in tmp:
            if k == 8:
                k = 0
                txt = txt + text + '\n'
            else:
                k = k + 1
                txt = txt + text + '\t'
        f.write(txt)

def start_run(filename):
    currentPath = '.'
    firstF = currentPath + '/data/Segments/' + filename + '.segments.txt'
    secondF = currentPath + '/data/VCF/' + filename + '.no_real_info.vcf'
    intergret(firstF,secondF)
