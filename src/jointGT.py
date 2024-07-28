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
    cmdDict = objU.procJointGTArgs()
    configFile = cmdDict['config'];workDir = cmdDict['workDir']
    interval_dir = cmdDict['intervalDir']
    manifest = cmdDict['manifest'];sh_type=cmdDict['shellType']
    projDate = cmdDict['projDate'];launch = cmdDict['launchFlag']
    expr = cmdDict['expr']
    scriptDir = cmdDict['scriptDir']

    #Read config XML file
    configDict = objU.getConfigDict(configFile)
    
    if re.search('ukbb',projDate):
        configDict['gatk']['GenomicsDBImport']['sample-name-map'] = \
        '/nfs/research/dunham/samples/ukbb/data/23161_ren_mapFile.txt'
    
    configDict['lsf']['params1']['M'] = str(20000)

    # Initialize project work directory
    tmp_dict = objU.processInit(scriptDir,sys.argv,projDate)
    tmp_dir = tmp_dict['tmpDir'];tmp_bin = tmp_dict['tmpBin'];
    tmp_data = tmp_dict['tmpData']; tmp_status =tmp_dict['tmpStat'];
    tmp_log = tmp_dict['tmpLog']; sb_log = tmp_dict['sbLog'];

    objC = cluster()
    #famDict = objC.getFamRelDict(pedFile,'cc')

    objP = picard()
    objG = gatk()

    chrList = list(range(1,23,1))
    chrList = ['chr'+str(x) for x in chrList]
    chrList.append('chrX')
    chrList.append('chrY')


    gvcfFileList = []
    #gvcfDir = '/hps/nobackup/research/dunham/ddd/edit_gdd/tmp_data'
    
    fh = open(manifest)
    for lines in fh:
        lines = lines.strip()
        if not re.search('^proband',lines):
            #for dddId in dddList:
            strs = re.split('\t',lines)
            dddId = strs[0].strip()
            if re.search('ukbb',projDate):
                gvcfFile = strs[1]
            else:
                gvcfFile = strs[8] #'/'.join([gvcfDir,dddId,dddId+'.gvcf.edit.norm.vcf.gz'])
                #gvcfFile = '/'.join([gvcfDir,dddId,dddId+'.gvcf.edit.vcf.gz'])
            if os.path.isfile(gvcfFile):
                gvcfFileList.append(gvcfFile)
            else:
                print('File do not exist: '+gvcfFile)
                print('\n')

    fh.close()

    jointGTChrDict = OrderedDict()


    for chrNum in chrList:
        
        chrIntervalDir = '/'.join([interval_dir,chrNum])
        #print(chrIntervalDir)
        intervalFileList = os.listdir(chrIntervalDir)
        jointGTChrFileList = []
        #print(intervalFileList)

        for intervalFile in intervalFileList:

            if re.search('bed',intervalFile):
                intervalFile = re.split('.bed',intervalFile)[0]
            else:
                intervalFile = intervalFile

            inBed = '/'.join([interval_dir,chrNum,intervalFile])
            if not re.search('bed$',intervalFile):
                outBed = '/'.join([interval_dir,chrNum,intervalFile+'.bed'])
            else:
                outBed = inBed
           
            #os.system('mv '+inBed+' '+outBed)

            #Change the tmp_data directory
            tmp_data = '/'.join([workDir,expr,'tmp_data'])
            outDir = '/'.join([tmp_data,chrNum,intervalFile])
            
            #print(outDir)
            #os.system('mkdir -p '+outDir)
            
            chrIntervalFile = chrNum+'_'+intervalFile
            sh_file,wh = objC.getClusterWriteHandle(tmp_bin,
                                                    'jointGT_'+chrIntervalFile,
                                                    sh_type
                                                   )
            wh = objC.writeClusterTop(configDict,
                                      'NDD-WES-GATK-SNV-JointGT-'+chrIntervalFile,
                                      wh
                                     )
            sh_out, sh_err,wh = objC.writeClusterInit(configDict,sb_log,
                                                      sh_type,
                                                      'JGT_'+expr+'_UKBB_V1_'+chrIntervalFile,
                                                      wh
                                                     )
            wh.write('echo "$HOSTNAME"\n\n')
            wh.write('mv '+inBed+' '+outBed)
            wh.write('\n\n')
            wh.write('mkdir -p '+outDir)
            wh.write('\n\n')
            
            # Genomics DB Import step
            if expr=='comb':
                chrDbImportFile,wh = objG.writeCombineGVCFs(configDict,
                                                        gvcfFileList,
                                                        outDir,chrNum,wh
                                                       )
            elif expr=='dbi':
                chrDbImportFile,wh = objG.writeGenomicsDBImport(configDict,
                                                        gvcfFileList,
                                                        outDir,chrNum,
                                                        intervalFile,wh
                                                       )
            else:
                print('Please enter valid analysis type: comb/dbi')

            # Joint Genotyping step
            chrJointGTFile,wh = objG.writeGenotypeGVCF(configDict,chrDbImportFile,
                                                       outDir,chrNum,intervalFile,
                                                       wh
                                                      )


            # Genotype GVCF
            jointGTChrFileList.append(chrJointGTFile)
            jointGTChrDict[chrNum] = jointGTChrFileList
            
            wh.close()
            
            ######## Launch the Bsub script ############
            if launch:
                #if not re.search('chrY_',sh_file):
                os.system('bsub < '+sh_file)

    ''' 
    sh_file, wh = objC.getClusterWriteHandle(tmp_bin,'jointGT_GatherVcfs',sh_type)
    wh = objC.writeClusterTop(configDict,'NDD-WES-GATK-SNV-GatherGT',wh)

    sh_out, sh_err,wh = objC.writeClusterInit(configDict,sb_log,
                                              sh_type,'jointGT_GatherVcfs',wh
                                             )

    gatherVCFFile =  objP.writeGatherVcfs(configDict,jointGTChrFileList,
                                          'gatherVcfs',tmp_data+'/chr1',
                                          None,wh
                                         )

    wh.close()
    '''

