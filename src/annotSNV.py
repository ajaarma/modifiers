#!/usr/bin/python

###############################################################################
#
# Description: Script to annotate UKB WES exome - 50K, 200K 
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
#
###############################################################################

import re,sys,os
import time

script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.append(script_path+'/config/')
sys.path.append(script_path+'/cluster/')
sys.path.append(script_path+'/util/')

from cluster import *
from util import *

if __name__=="__main__":

    objU = util()
    cmdDict = objU.procAnnotSNVargs()
    configFile = cmdDict['config'];workDir = cmdDict['workDir']
    scriptDir = cmdDict['scriptDir']
    jgtRecalFile = cmdDict['jgtFile'];sh_type=cmdDict['shellType']
    projDate = cmdDict['projDate'];launch = cmdDict['launchFlag']

    #Read config XML file
    configDict = objU.getConfigDict(configFile)

    # Initialize project work directory
    tmp_dict = objU.processInit(scriptDir,sys.argv,projDate)
    tmp_dir = tmp_dict['tmpDir'];tmp_bin = tmp_dict['tmpBin'];
    tmp_data = tmp_dict['tmpData']; tmp_status =tmp_dict['tmpStat'];
    tmp_log = tmp_dict['tmpLog']; sb_log = tmp_dict['sbLog'];

    objC = cluster()
    objA = annotation()

  
    if not re.search('200k',projDate):
        if not re.search('ukbb',projDate):
            chr_list = ['chr'+str(x) for x in range(1,23)]+['chrX','chrY']
        else:
            chr_list = [str(x) for x in range(1,23)]+['X','Y']


        for chrNum in chr_list:

            outChrDir = '/'.join([tmp_data,chrNum])
            fileId = chrNum+'_annotSNV'

            configDict['lsf']['params1']['M'] = '10G'
            configDict['lsf']['params1']['q'] = 'research'

            sh_file,wh = objC.getClusterWriteHandle(tmp_bin,fileId,sh_type)
            wh = objC.writeClusterTop(configDict,'NDD-WES-GATK-AnnotSNV',wh)
            sh_out,sh_err,wh = objC.writeClusterInit(configDict,sb_log,
                                                     sh_type,fileId,wh
                                                    )
            
            wh.write('echo "$HOSTNAME"\n\n')
            #wh.write('unset BCFTOOLS_PLUGINS\n')
            wh.write('mkdir -p '+outChrDir+'\n')

            ''' 
            chrSplitOutFile,wh = objA.writeSplitByChrom(jgtRecalFile,
                                                      outChrDir,chrNum,wh
                                                      )
           
            normOutFile,wh = objA.writeNormalizeVCF(configDict,chrSplitOutFile,
                                                    outChrDir,chrNum,wh
                                                   )

            vepOutFile,wh = objA.writeVEPannotation(configDict,normOutFile,
                                                    outChrDir,chrNum,wh
                                                   )
            '''

            vepOutFile = '/'.join([outChrDir,
                                   'tmpFile_'+chrNum+'.AC0.norm.vep.bcf'
                                  ])
            vepCustAnnoFile,wh = objA.writeClusterCustomAnnot(configDict,vepOutFile,
                                                            outChrDir,chrNum,wh
                                                           )

            vepAnnoICFreq,wh = objA.writeClusterICFreq(configDict,
                                                       vepCustAnnoFile,
                                                       outChrDir,chrNum,
                                                       projDate,wh
                                                      )
            wh.close()

            if launch:
                os.system('bsub < '+sh_file)
    else:
        pvcf_dir = os.path.abspath(jgtRecalFile)
        fileList = os.listdir(pvcf_dir)

        vcfDict = OrderedDict()
        vcfDict2 = OrderedDict()

        chr_list = ['c'+str(x) for x in range(1,23)]+['cX','cY']

        for cb_file in fileList:
            if re.search('.vcf.gz$',cb_file):
                strs = re.split('_',cb_file)
                chrNum = strs[1]
                if chrNum in vcfDict:
                    tmp = vcfDict[chrNum]
                    #tmp.append('/'.join([pvcf_dir,cb_file]))
                    tmp.append(cb_file)
                    vcfDict[chrNum] = tmp
                    tmp = []
                else:
                    #vcfDict[chrNum] = ['/'.join([pvcf_dir,cb_file])] 
                    vcfDict[chrNum] = [cb_file] 
        
        for ele in chr_list:
            vcfDict2[ele] = vcfDict[ele]

        for chr_num,values in vcfDict2.items():
            #print(chrNum+'\t'+','.join(values)+'\t'+str(len(values)))
            for vcf_ele in values:
          
                chrNum = 'chr'+str(re.split('^c',chr_num)[1])
                print(chrNum)
                vcf_ele_strs = re.split('.vcf.gz',vcf_ele)
                outChrDir = '/'.join([workDir,projDate,chrNum,vcf_ele_strs[0]])
                fileId = vcf_ele_strs[0]+'_annotSNV'

                configDict['lsf']['params1']['M'] = '10G'
                configDict['lsf']['params1']['q'] = 'research'

                sh_file,wh = objC.getClusterWriteHandle(tmp_bin,fileId,sh_type)
                wh = objC.writeClusterTop(configDict,'NDD-WES-GATK-AnnotSNV',wh)
                sh_out,sh_err,wh = objC.writeClusterInit(configDict,sb_log,
                                                         sh_type,fileId,wh
                                                        )
                
                wh.write('echo "$HOSTNAME"\n\n')
                #wh.write('unset BCFTOOLS_PLUGINS\n')
                wh.write('mkdir -p '+outChrDir+'\n')

                jgtRecalFile = '/'.join([pvcf_dir,vcf_ele])
                
                chrSplitOutFile,wh = objA.writeSplitByChrom(jgtRecalFile,
                                                          outChrDir,chrNum,wh
                                                          )
                #chrSplitOutFile = '/'.join([pvcf_dir,vcf_ele])
                normOutFile,wh = objA.writeNormalizeVCF(configDict,chrSplitOutFile,
                                                        outChrDir,chrNum,wh
                                                       )

                vepOutFile,wh = objA.writeVEPannotation(configDict,normOutFile,
                                                        outChrDir,chrNum,wh
                                                       )

                 
                #vepOutFile = '/'.join([outChrDir,
                #                       'tmpFile_'+chrNum+'.AC0.norm.vep.bcf'
                #                      ])
                vepCustAnnoFile,wh = objA.writeClusterCustomAnnot(configDict,vepOutFile,
                                                                outChrDir,chrNum,wh
                                                               )

                vepAnnoICFreq,wh = objA.writeClusterICFreq(configDict,
                                                           vepCustAnnoFile,
                                                           outChrDir,chrNum,
                                                           projDate,wh
                                                          )
                
                wh.write('\n\n')
                wh.write('rm '+vepOutFile+'\n')
                #wh.write('rm '+vepCustAnnoFile)
                wh.write('\n')
                wh.close()
                if launch:
                    os.system('bsub < '+sh_file)
                    time.sleep(1.5) 
                

