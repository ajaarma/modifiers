#!/usr/bin/python

import re,sys,os

script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.append(script_path+'/config/')
sys.path.append(script_path+'/cluster/')
sys.path.append(script_path+'/util/')

from cluster import *
from util import *

if __name__=="__main__":

    objU = util()
    cmdDict = objU.procGATKArgs()
    configFile = cmdDict['config'];workDir = cmdDict['workDir']
    scriptDir = cmdDict['scriptDir']
    pedFile = cmdDict['inpFile'];sh_type=cmdDict['shellType']
    projDate = cmdDict['projDate'];launch = cmdDict['launchFlag']

    #Read config XML file
    configDict = objU.getConfigDict(configFile)

    # Initialize project work directory
    tmp_dict = objU.processInit(scriptDir,sys.argv,projDate)
    tmp_dir = tmp_dict['tmpDir'];tmp_bin = tmp_dict['tmpBin'];
    tmp_data = tmp_dict['tmpData']; tmp_status =tmp_dict['tmpStat'];
    tmp_log = tmp_dict['tmpLog']; sb_log = tmp_dict['sbLog'];

    objC = cluster()
    famDict = objC.getFamRelDict(pedFile,'cc')

    objP = picard()
    objG = gatk()

    for key,value in famDict.items():
        
        pat_id = key
        cramFile = famDict[key]['path']
        tmp_data = workDir
        sampleWorkDir = '/'.join([tmp_data,pat_id])
        tmpStatFile = '/'.join([tmp_status,pat_id+'.s.txt'])
        #tmp_bin_new = '/'.join([tmp_bin,pat_id])
        #os.system('mkdir -p '+tmp_bin_new)

        #if re.search('cip',cramFile):
        if re.search('cram$',cramFile):
            #### STEP 1 ######
            fileId = pat_id
            sh_file, wh = objC.getClusterWriteHandle(tmp_bin,fileId,sh_type)
            wh = objC.writeClusterTop(configDict,'NDD-WES-GATK-SNV',wh)
            sh_out, sh_err,wh = objC.writeClusterInit(configDict,sb_log,
                                                      sh_type,fileId,wh
                                                     )
            wh.write('echo "$HOSTNAME"\n\n')

            # CRAM to BAM conversion
            bamFile, wh = objP.writeCramToBam(configDict,sampleWorkDir,pat_id,
                                              cramFile,None,wh
                                             )
            # BAM to RevertSAM; Split across RG-group
            rgIdList,rgBamList,outRGDir,wh = objP.writeRevertSam(configDict,
                                                                 sampleWorkDir,
                                                                 cramFile,
                                                                 bamFile,
                                                                 bamFile,wh
                                                                )
            
            outRGBamDirList = []
            ##### STEP 2 ######
            ''' 
            sh_file, wh = objC.getClusterWriteHandle(tmp_bin_new,fileId,
                                                     sh_type
                                                    )
            wh = objC.writeClusterTop(configDict,'NDD-WES-FQBWA-STEP2',wh)
            sh_out, sh_err,wh = objC.writeClusterInit(configDict,sb_log,
                                                      sh_type,fileId,wh
                                                     )

            wh.write('echo "$HOSTNAME"\n\n')
            outRGBamDir = '/'.join([outRGDir,rgId])
            outRGBamDirList.append(outRGBamDir)
            '''
            miaBamList,wh = objP.writeMarkIlluminaAdapters(configDict,
                                                            outRGDir,
                                                            rgBamList,
                                                            rgIdList,
                                                            tmpStatFile,
                                                            wh
                                                           )
            #Sam To Fastq
            fqGzList,wh = objP.writeSamToFastq(configDict,outRGDir,
                                               miaBamList,rgIdList,
                                               None,wh
                                              )

            #BWA Alignment - grch38
            bwaAlnList,wh = objP.writeBwaAlign(configDict,outRGDir,
                                               fqGzList,rgIdList,
                                               None,wh
                                              )
            # Merge Aligned BAM with Un-Mapped BAM per Read Group (RG)
            mergeBamAlnList,wh = objP.writeMergeBamAlignment(configDict,
                                                          rgBamList,
                                                          bwaAlnList,
                                                          outRGDir,
                                                          None,
                                                          wh
                                                         )
            #Merge Sam Files
            mergeRGSamFile,wh = objP.writeMergeRGSamFiles(configDict,mergeBamAlnList,
                                                          outRGDir,pat_id,
                                                          None,wh
                                                         )
            # MarkDuplicates
            dedupFile,wh = objP.writeMarkDuplicates(configDict,mergeRGSamFile,
                                                      None,wh
                                                    )


            ###### STEP 4 #######
            ####### BaseRecalibrator ########
            objG = gatk()
            chrList = list(range(1,23,1))
            chrList = ['chr'+str(x) for x in chrList]
            chrList.append('chrX')
            chrList.append('chrY')
   
            gvcfList = []
            chromBamList = []

            #BaseRecalibrator step
            baseRecalFile,wh = objG.writeBaseRecalibrator(configDict,dedupFile,
                                                          outRGDir,
                                                          None,
                                                          wh
                                                         )
            # ApplyBQSR
            bqsrFile,wh = objG.writeApplyBQSR(configDict,dedupFile,baseRecalFile,
                                              pat_id,outRGDir,None,wh
                                             )
            # HaplotypeCaller
            hcGvcfFile,wh = objG.writeHaplotypeCaller(configDict,bqsrFile,pat_id,
                                                      chrList,outRGDir,None,
                                                      wh
                                                     )
            hcEditGvcfFile,wh = objG.writeEditHCGvcf(configDict,pat_id,hcGvcfFile,wh)
            '''
            ############## Gather Chromosome specific VCFs ##############
            combGvcfFile,wh = objP.writeGatherVcfs(configDict,chrGvcfFileList,pat_id,
                                                    outRGDir,None,wh
                                                  )
            '''
            destPath = os.path.dirname(outRGDir)
            bqsrFileName = os.path.basename(bqsrFile)
            gvcfFileName = os.path.basename(hcEditGvcfFile)

            mvCmd = ['mv',bqsrFile,'/'.join([destPath,bqsrFileName])]
            bqsr_cmd_index = 'samtools index '+'/'.join([destPath,bqsrFileName])

            mvCmdHC = ['mv',hcEditGvcfFile,'/'.join([destPath,gvcfFileName])]
            gvcf_cmd_index = 'tabix -p vcf '+'/'.join([destPath,gvcfFileName])

            wh.write('echo "Moving hg38 aligned, dedup, recal file"\n')
            wh.write(' '.join(mvCmd))
            wh.write('\n')
            wh.write(bqsr_cmd_index)
            wh.write("\n\n")

            wh.write('echo "Moving HC gvcfFile "\n')
            wh.write(' '.join(mvCmdHC))
            wh.write('\n')
            wh.write(gvcf_cmd_index)
            wh.write("\n\n")
             
            wh.write('echo "Removing the tmp directory"\n')
            rmCmd = 'rm -r '+outRGDir
            wh.write(rmCmd)
            wh.write('\n\n')
           
            #cpStoragePath = '/hps/research1/dunham/samples/ddd/df_2017/GDD-GATK'
            #cpCmd = 'cp -r '+destPath+' '+cpStoragePath
            #wh.write(cpCmd)
            #wh.write('\n\n')
            
            wh.close()
            
            ######## Launch the Bsub script ############
            if launch:
                os.system('bsub < '+sh_file)
