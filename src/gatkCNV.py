#!/usr/bin/python

import re,sys,os
import pandas as pd
import random
import numpy as np

script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.append(script_path+'/config/')
sys.path.append(script_path+'/cluster/')
sys.path.append(script_path+'/util/')

from cluster import *
from util import *

random.seed(10)

def sub_lsf(configFile,s_200_df,s_rem_df,clusterID,
            file_str,tmp_bin,sb_log,workDir,launch):
    ''' Subroutine to create LSF script for 200 samples
    and remaining per cluster '''

    objC = cluster()
    objP = picard()
    objG = gatk()
    objV = cnv()
        
    #Read config Analysis XML file
    configDict = objU.getConfigDict(configFile)
 
    s_cohort_list = list(s_200_df.PATH)
    s_case_list = list(s_rem_df.PATH)

    file_str_1 =  file_str+'_'+clusterID
    projDateStrs = re.split('_',file_str_1)
    cl_type = projDateStrs[0]
    cnv_mat = '_'.join(projDateStrs[1:5])
    outDir = '/'.join([workDir,cl_type,cnv_mat,clusterID])

    #Write out-file to the outDir
    ttDir = '/nfs/research/sds/sds-dham-ukbndd/200k/analysis/cnv/cnv_mat/clusters/hpc_train_test'
    out_200_df = '/'.join([ttDir,'train_'+cnv_mat+'_'+clusterID+'.csv'])
    out_rest_df = '/'.join([ttDir,'rest_'+cnv_mat+'_'+clusterID+'.csv'])
    s_200_df.to_csv(out_200_df,sep='\t',header=True)
    s_rem_df.to_csv(out_rest_df,sep='\t',header=True)

    # Preprocess Intervals
    #fileId = 'ms.ppi.rc.flt.pld.sct.call.'+projDate
    fileId = file_str_1

    sh_type = 'LSF'
    configDict['lsf']['params1']['M'] = '70G'
    #configDict['lsf']['params1']['q'] = 'research'
    configDict['lsf']['params1']['q'] = 'research'

    sh_file, wh = objC.getClusterWriteHandle(tmp_bin,fileId,sh_type)
    wh = objC.writeClusterTop(configDict,'HPC-CNV-'+projDate,wh)
    sh_out,sh_err,wh = objC.writeClusterInit(configDict,sb_log,
                                              sh_type,fileId,wh
                                             )

    #wh.write('unset $THEANO_FLAGS\n')
    wh.write('module load singularityce/3.10.3\n')
    wh.write('echo "$HOSTNAME"\n\n')
    wh.write('mkdir -p '+outDir+'\n')

    ppIntervalList = '/hps/nobackup/dunham/ai-uk-can/ukb200k/gatk_cnv/hpc_cnv/batch_00/gdd.pp.interval_list'
    ppGCannotInterval = '/hps/nobackup/dunham/ai-uk-can/ukb200k/gatk_cnv/hpc_cnv/batch_00/gdd.pp.gc.annotated.tsv'
    scatterDir = '/hps/nobackup/dunham/ai-uk-can/ukb200k/gatk_cnv/cnv_calls/scatter_for_all'

    #wh.close()
     
    wh.write('cp '+ppGCannotInterval+' '+outDir+'\n')
    ##### Filter Intervals #######
     
    fltDir = '/'.join([workDir,projDate])
    filterIntervalList,wh = objV.getFilterIntervals(configDict,ppIntervalList,
                                                    ppGCannotInterval,s_cohort_list,
                                                    outDir,wh
                                                   )
    

    ##### Step 3: Determine Germline Contig Ploidy ######
    #pldDir = '/'.join([workDir,projDate,'ploidy'])
    pldDir = '/'.join([outDir,'ploidy'])
    wh.write('mkdir -p '+pldDir)
    wh.write('\n')

    
    ploidyDir,wh = objV.getDetermineGermlineContigPloidy(configDict,
                                                         filterIntervalList,
                                                         s_cohort_list,
                                                         pldDir,wh
                                                        )

    ''' 
    #### Step 4.1: Scatter Intervals ######
    sctDir = '/'.join([outDir,'scatter'])
    wh.write('mkdir -p '+sctDir)
    wh = objV.getIntervalListTools(configDict,filterIntervalList,sctDir,wh)
    wh.close()
    ''' 

    ##### Step 4.2: Germline CNV Calling #####

    #cnvDir = '/'.join([workDir,projDate,'cnv_call'])
    cnvDir = '/'.join([outDir,'cnv_call'])
    wh.write('mkdir -p '+cnvDir)
    wh.write('\n')

    #sctDirList = os.listdir(scatterDir)
    cohort_mode = 'COHORT'

    wh =  objV.getGermlineCNVCaller(configDict,outDir,cnvDir,
                                        s_cohort_list,
                                        cohort_mode,wh)

    wh.close()
    
    ##### Step 5: PostprocessGermlineCNVCalls for cohort #####
    '''
    #vcfDir = '/'.join([workDir,projDate,'cnv_vcf'])
    vcfDir = '/'.join([outDir,'cnv_vcf'])
    wh.write('mkdir -p '+vcfDir)
    wh.write('\n')

    for index,rows in s_200_df.iterrows():
        sample_index = index
        sample_name = str(rows['SAMPLE'])
        s_vcf_dir = '/'.join([vcfDir,sample_name])
        wh.write('mkdir -p '+s_vcf_dir)
        wh.write('\n')
        wh = objV.getPostprocessGermlineCNVCalls(configDict,outDir,
                                                 vcfDir,sample_index,
                                                 sample_name,wh)
    wh.close()
  
    if s_rem_df.shape[0] >0:
        for index,rows in s_rem_df.iterrows():
            sample_index = index
            sample_name = str(rows['SAMPLE'])
            s_dir = '/'.join([outDir,sample_name])
            wh.write('mkdir -p '+s_dir)
            wh.write('\n')
            s_case_list = [sample_name]
            pldCaseDir = '/'.join([s_dir,'ploidy'])
            ploidyDir,wh = objV.getDetermineGermlineContigPloidy(configDict,
                                                         filterIntervalList,
                                                         s_case_list,
                                                         pldCaseDir,wh
                                                        )

            #Germline CNV caller
            caseCnvDir = '/'.join([s_dir,'cnv_call'])
            wh =  objV.getGermlineCNVCaller(configDict,s_dir,
                                            caseCnvDir,s_case_list,wh)

            # PostprocessGermlineCNVCalls
            wh = objV.getPostprocessGermlineCNVCalls(configDict,s_dir,
                                                 vcfDir,sample_index,
                                                 sample_name,wh)

    
    ''' 

    ######## Launch the Bsub script ############
    if launch:
        if not re.search('hpc_cnv_mat_7_3',sh_file):
            os.system('bsub < '+sh_file)
            #print(sh_file)

 
if __name__=="__main__":

    objU = util()
    cmdDict = objU.procGatkCnvArgs()
    configFile = cmdDict['config'];workDir = cmdDict['workDir']
    scriptDir = cmdDict['scriptDir']
    pedDir = cmdDict['inpFile'];sh_type=cmdDict['shellType']
    projDate = cmdDict['projDate'];launch = cmdDict['launchFlag']

    #Read config XML file
    #configDict = objU.getConfigDict(configFile)

    # Initialize project work directory
    tmp_dict = objU.processInit(scriptDir,sys.argv,projDate)
    tmp_dir = tmp_dict['tmpDir'];tmp_bin = tmp_dict['tmpBin'];
    tmp_data = tmp_dict['tmpData']; tmp_status =tmp_dict['tmpStat'];
    tmp_log = tmp_dict['tmpLog']; sb_log = tmp_dict['sbLog'];

    objC = cluster()

    #famDict = objC.getFamRelDict(pedFile,'cc')

    objP = picard()
    objG = gatk()
    objV = cnv()
    
    #dddList = os.listdir(pedDir)
    if re.search('hpc',projDate):
        hpc_dnx_file = '/nfs/research/dunham/aak/dnx/output/manifest/sv/hpc_cnv_manifest.txt'
        clustDir = '/'.join([pedDir,'hpc'])
    elif re.search('dnx',projDate):
        hpc_dnx_file = '/nfs/research/dunham/aak/dnx/output/manifest/sv/hpc_dnx_manifest.txt'
        clustDir = '/'.join([pedDir,'dnx'])
    
    clustMatList = os.listdir(clustDir)
    
    samplesDF = pd.read_csv(hpc_dnx_file,sep='\t',header=None,names=["PATH"])
    a1 = samplesDF.PATH.str.split('/|_')
    a2 = [x[13] for x in a1]
    samplesDF['SAMPLE'] = a2
    samplesDF['SAMPLE'] = samplesDF['SAMPLE'].astype(np.int64)

    #Read each file as DF

    for _inFile in clustMatList:
        
        low_ind = re.split('_|\.',_inFile)[2]
        up_ind = re.split('_|\.',_inFile)[3]
        file_str = '_'.join([projDate,str(low_ind),str(up_ind)])

        matFile = '/'.join([clustDir,_inFile])
        
        matDF = pd.read_csv(matFile,sep='\t',names=['SAMPLE','CLUSTER'],header=None)
        for i in range(1,6):
            matDF_tmp = matDF[matDF['CLUSTER']==i]
            if matDF_tmp.shape[0] >=100:
                r_200_df = matDF_tmp.sample(100,replace=False)
                r_rem_df = pd.merge(matDF_tmp,r_200_df, how='left', indicator=True)
                r_rem_df = r_rem_df[r_rem_df['_merge'] == 'left_only']
            else:
                r_200_df = matDF_tmp
                r_rem_df = pd.DataFrame()
            
            #Query the 200 samples in sampleDF for path
            s_200_df = pd.merge(samplesDF,r_200_df,how='inner',on=['SAMPLE'])
            if r_rem_df.shape[0] >0:
                s_rem_df = pd.merge(samplesDF,r_rem_df,how='inner',on=['SAMPLE'])
            else:
                s_rem_df = pd.DataFrame()

            clusterID = 'c'+str(i)
            
            sub_lsf(configFile,s_200_df,s_rem_df,clusterID,
                    file_str,tmp_bin,sb_log,workDir,launch)
            
        #sys.exit() 
            

