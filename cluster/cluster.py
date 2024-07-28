#!/usr/bin/python3

import re,os,sys
import subprocess
import datetime
from collections import OrderedDict

class cluster:

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def getClusterWriteHandle(self,bin_dir,sh_name,sh_type=[],chr_num=[]):
        '''Subroutine to genenrate SLURM/LSF/SHELL script file name '''

        if len(chr_num)!=0:
            sh_file = os.path.abspath(bin_dir+"/"+sh_name+"_chr"+chr_num+".sh")
            wh = open(sh_file,'w')
        elif re.search('slurm',sh_type,re.I):
            sh_file = os.path.abspath('/'.join([bin_dir,sh_name+'.slurm']))
            wh = open(sh_file,'w')
        elif re.search('lsf',sh_type,re.I):
            sh_file = os.path.abspath('/'.join([bin_dir,sh_name+'.lsf']))
            wh = open(sh_file,'w')

        return sh_file,wh

    def writeClusterTop(self,configDict,proj_name,wh):
        ''' Subroutine to write the top header of lsf/slurm batch script '''

        now = datetime.datetime.now()
        vers = configDict['general']['version']
        wh.write(
                '\n'.join([
                    '#!/bin/bash',
                    '#!',
                    '## date:'+'-'.join([str(now.year),str(now.month),str(now.day)]),
                    '## project name: '+proj_name,
                    '## pipeline version: '+vers,
                    '## genome version: '+configDict['general']['genomeBuild'],
                    '##\n\n'
                ])
        )

        return wh

    def processDDDLogFile(self,logFile,outFile):
        ''' Subroutine to process the DDD log file and return back the 
            Sample names as dictionary 
        '''

        sampleDict = {}
        count = 0

        fh = open(logFile,'r')

        for lines in fh:
           lines = lines.strip()
           strs = re.split('\s',lines)
           strs = [x for x in strs if x]
           if len(strs)==8 and re.search('^EGA',strs[3]):
               count = count+1
               ega_id = strs[3]
               sample_name = strs[7]
               sampleDict[count] = {'ega':ega_id,'cram':sample_name}
                
        fh.close()

        wh = open(outFile,'w')

        for keys in sampleDict:
            wh.write(' '.join([str(keys),sampleDict[keys]['ega'],sampleDict[keys]['cram'],'\n']))

        wh.close()


        return sampleDict

    def writeClusterInit(self,configDict,sb_log,sh_type,task,wh):
        ''' Subroutine to write initializtion specific details for SLURM/LSF 
            scripts 
        '''
   
        email_id = configDict['general']['email']
        script_out = sb_log+"/"+task+".o.txt"
        script_err = sb_log+"/"+task+".e.txt"
    
        if re.search('slurm',sh_type,re.I):
            slurmDict = configDict['slurm']

            for slr in slurmDict:
                dash = slurmDict['dash']
                doubDash = slurmDict['doubDash']

                if slr=='params1':
                    par1_dict = slurmDict[slr]
                    for par1,values in par1_dict.items():
                        if par1=='J':
                            wh.write(
                                ' '.join(['#SBATCH',dash+par1,values,'\n'])
                            )
                        else:
                            wh.write(
                                ' '.join(['#SBATCH',dash+par1,values,'\n'])
                            )
                if slr=='params2':
                    par2_dict = slurmDict[slr]
                    for par2,values in par2_dict.items():
                        try:
                            if par2=='mail-user':
                                wh.write(
                                    ' '.join(['#SBATCH',doubDash+par2+'='+email_id,'\n'])
                                )
                            elif par2=='output':
                                wh.write(
                                    ' '.join(['#SBATCH',doubDash+par2+'='+script_out,'\n'])
                                )
                            elif par2=='error':
                                wh.write(
                                    ' '.join(['#SBATCH',doubDash+par2+'='+script_err,'\n'])
                                )
                            else:
                                 wh.write(
                                    ' '.join(['#SBATCH',doubDash+par2+values,'\n'])
                                )
                        except:
                            wh.write(
                                ' '.join(['#SBATCH',doubDash+par2,'\n'])
                            )
                    wh.write('\n')

        elif re.search('lsf',sh_type,re.I):
            lsfDict = configDict['lsf']
            for slr in lsfDict:
                dash = lsfDict['dash']

                if slr=='params1':
                    par1_dict = lsfDict[slr]
                    for par1,values in par1_dict.items():
                        if par1=='J':
                            wh.write(
                                ' '.join(['#BSUB',dash+par1,task,'\n'])
                            )
                        elif par1=='u':
                             wh.write(
                                ' '.join(['#BSUB',dash+par1,email_id,'\n'])
                            )
                        elif par1=='o':
                             wh.write(
                                ' '.join(['#BSUB',dash+par1,script_out,'\n'])
                            )
                        elif par1=='e':
                             wh.write(
                                ' '.join(['#BSUB',dash+par1,script_err,'\n'])
                            )
                        elif par1=='N':
                             wh.write(
                                ' '.join(['#BSUB',dash+par1,'\n'])
                            )
                        else:
                            wh.write(
                                ' '.join(['#BSUB',dash+par1,values,'\n'])
                            )
                    wh.write('\n')
        return script_out, script_err, wh

    def writeClusterDownloadDDD(self,configDict,sampleDict,tmp_data,low,up,wh):
        ''' Subroutine to wrote pyega3 commands for downloading of the files 
        '''
        pyegaConfig = configDict['pyega']['cf']
        pyegaBin = configDict['pyega']['bin']
        clust_num = configDict['pyega']['c']

        for countKey in range(low,up):
            countKey = countKey+1
            print(sampleDict[int(countKey)])
            egaID = sampleDict[countKey]['ega']
            #os.system('mkdir -p '+'/'.join([tmp_data,egaID]))
            cramFile = '/'.join([tmp_data,egaID,sampleDict[countKey]['cram']])

            pyega_cmd = ' '.join([pyegaBin,'-c',clust_num,'-cf',pyegaConfig,
                                            'fetch',egaID,'--saveto',cramFile])
            wh.write('echo \' The sample counter:'+str(countKey)+'\'\n')
            wh.write('mkdir -p '+'/'.join([tmp_data,egaID])+'\n')
            wh.write(pyega_cmd+'\n')

        return wh

    def processEga(configDict):
        ''' Subroutine to process EGA datasets. Returns dictionary for all
            the respective mappings 
        '''

        egaDict = configDict['ega']
        ddd_ega_file = egaDict['ddd2ega']
        pat_pehno_file = egaDict['patPheno']
        fam_rel_file = egaDict['famRel']
        ext_id_map = egaDict['extId']


    def getDdd2EgaDict(ddd_ega_file):
        ''' Subroutine to get ddd patient to EGA id mapping dictionary '''

        dpEgaDict = OrderedDict()

        fh = open(ddd_ega_file,'r')

        for lines in fh:
            lines = lines.trip()
            strs = re.split('\t',lines);strs = [x.strip() for x in strs]
            if not re.search('person',strs[0]):
                patientId = strs[0]
                egaId = strs[1]
                dpEgaDict[patientId] = egaId

        fh.close()
        
        return dpEgaDict

    def getPatientPhenoDict(pat_pheno_file):
        ''' Subroutine to get Patient Phenotype dictionary '''

        ppDict = OrderedDict()

        fh = open(pat_pheno_file)
        for lines in fh:
            lines = lines.strip()
            if not re.search('^stable',lines):
                strs = re.split('\t',lines); strs = [x.strip() for x in strs]
                dpId = strs[0];gender = strs[1];age=str(strs[2]);
                childHpo = strs[3];childTerm = strs[4];
                patHpo = strs[5];patTerm = strs[6];
                matHpo = strs[7];matTerm = strs[8]
                
                ppDict[dpId]={'gender':gender,
                              'age':age,
                              'childHpo':childHpo,
                              'childTerm':childTerm,
                              'patHpo':patHpo,
                              'patTerm':patTerm,
                              'matHpo':matHpo,
                              'matTerm':matTerm
                             }
        fh.close()

        return ppDict

    def getFamRelDict(self,fam_rel_file,flag=[]):
        ''' Subroutine to get family relationship dictionary '''

        frDict = OrderedDict()

        fh = open(fam_rel_file)
        for lines in fh:
            lines = lines.strip()
            if not re.search('^family',lines):
                strs = re.split('\t|\s',lines); strs=[x.strip() for x in strs]
                if flag=='cc':
                    famId = strs[1]; patId = strs[0];
                else:
                    famId = strs[0]; patId = strs[1];
                dadId = strs[2];mumId = strs[3];
                gender = strs[4];affectStatus = str(strs[5])
                
                try:
                    path = strs[6]
                except:
                    path = ''
                
                frDict[patId] = {'famId':famId,
                                 'father':dadId,
                                 'mother':mumId,
                                 'gender':gender,
                                 'affected':affectStatus,
                                 'path':path
                                }

        fh.close()
        
        return frDict

    def getProbandExtDict(ext_id_map):
        ''' Subroutine to get external mapping ID for the proband '''

        extDict = OrderedDict()
        fh = open(ext_id_map)

        for lines in fh:
            lines = lines.strip()
            if not re.search('^proband',lines):
                strs = re.split('\t',lines);strs = [x.strip() for x in strs]
                probandId = strs[0]
                extId = strs[1]
                extDict[probandId] = extId

        fh.close()
        return extDict
    
    def getEgaMappingDict(self,configDict,workDir):
        ''' Subroutine to get respective DDD to EGA mappings '''

        ega_file = '/'.join([workDir,configDict['ega']['ddd2ega']])
        pheno_file = '/'.join([workDir,configDict['ega']['patPheno']])
        fam_file = '/'.join([workDir,configDict['ega']['famRel']])
        map_file = '/'.join([workDir,configDict['ega']['sampleMap']])

        # Processing ddd 2 ega mapping file
        egaDict = OrderedDict()

        fh = open(ega_file)
        for lines in fh:
            lines = lines.strip()
            if not re.search('\_person',lines):
                strs = re.split('\t',lines)
                strs = [x.strip() for x in strs]
                dddp_id = strs[0]
                ega_s_id = strs[1]
                egaDict[dddp_id] = ega_s_id
        fh.close()

        # Processing Patient Phenotypes file
        phenoDict = OrderedDict()

        fh = open(pheno_file)
        for lines in fh:
            lines = lines.strip()
            if not re.search('^stable',lines):
                strs = re.split('\t',lines)
                strs = [x.strip() for x in strs]
                dddp_id = strs[0]
                phenoDict[dddp_id] = strs[1:]
        fh.close()

        # Processing Family relationships
        famDict = OrderedDict()
        fh = open(fam_file)

        for lines in fh:
            lines = lines.strip()
            if not re.search('^family',lines):
                strs = re.split('\t',lines);strs=[x.strip() for x in strs]
                fam_id = strs[0]
                dddp_id = strs[1]
                dad_id = strs[2]
                mum_id = strs[3]
                gender = strs[4]
                affected = strs[5]
               
                print(dddp_id)
                famDict[dddp_id] = {'fam':fam_id,'dad':dad_id,'mum':mum_id,
                                    'sex':gender,'status':affected
                                   }
                if dddp_id=='DDDP139147':
                    print(lines)

                #print(famDict['DDDP139147']+'\n')
        fh.close()

        # Processing Sample Map file
        fh = open(map_file)
        mapDict = OrderedDict()
        for lines in fh:
            lines = lines.strip()
            strs = re.split('\t',lines); strs=[x.strip() for x in strs]
            ega_id = strs[1]
            bam_id = strs[2]
            egf_id = strs[3]
            mapDict[ega_id] = {'egf':egf_id,'bam':bam_id}
        fh.close()


        return egaDict,phenoDict,famDict,mapDict


    def writeClusterTrios(self,configDict,inpFile,outFile,egaDict,phenoDict,
                                          famDict,mapDict,tmp_data,low,up,wh
                         ):
        ''' Subroutine to map mild to severe phenotypes to Trios '''
        pyegaConfig = configDict['pyega']['cf']
        pyegaBin = configDict['pyega']['bin']
        clust_num = configDict['pyega']['c']


        wh1 = open(outFile,'w')
        wh1.write('\t'.join(['KG_ID','DDD_ID','PHENOTYPE\n']))

        count_key = 0
        fh = open(inpFile)
        for lines in fh:
            lines = lines.strip()
            if not re.search('^individual',lines):
                strs = re.split('\,',lines);strs=[x.strip() for x in strs]
                kg_id = strs[0]
                kg_id_strs = re.split('KG',kg_id)
                kg_id_new = 'DDDP1'+str(kg_id_strs[1])
                pheno_type = strs[1] 
                
                wh1.write('\t'.join([kg_id,kg_id_new,pheno_type]))
                wh1.write('\n')

                if kg_id_new in egaDict.keys():
                    ega_s_id = egaDict[kg_id_new]
                    egf_id = mapDict[ega_s_id]['egf']
                    #bam_id = re.split('\.cip',mapDict[ega_s_id]['bam'])[0]
                    bam_id = mapDict[ega_s_id]['bam']
                    cramFile = '/'.join([tmp_data,kg_id_new,bam_id])

                    pyega_cmd = ' '.join([pyegaBin,'-c',clust_num,'-cf',pyegaConfig,
                                                 'fetch',egf_id,'--saveto',cramFile
                                         ])
                    count_key = count_key+1
                    wh.write('echo \' The sample counter:'+str(count_key)+'\'\n')
                    wh.write('echo \'The mappings are:\''+'\t'.join([kg_id_new,
                                                                     ega_s_id,
                                                                     egf_id,
                                                                     bam_id+'\''+'\n'
                                                                    ])
                            )
#                   wh.write('mkdir -p '+'/'.join([tmp_data,kg_id_new])+'\n')
                    wh.write(pyega_cmd+'\n\n')

                else:
                    wh.write('EGA-missing: '+kg_id_new+'\n')

                '''
                if kg_id_new in famDict.keys():
                    fam_id = famDict[kg_id_new]['fam']
                    dad_id = famDict[kg_id_new]['dad']
                    mum_id = famDict[kg_id_new]['mum']
                    gender = famDict[kg_id_new]['sex']
                    affected = famDict[kg_id_new]['status']

                    wh.write('\t'.join([kg_id_new,fam_id,dad_id,mum_id,
                                                    gender,affected+'\n'
                                       ])
                            )
                else:
                    wh.write('EGA-missing-fam: '+kg_id_new+'\n')
                '''

        fh.close()
        wh1.close()
        return wh
    
    def checkErrorFile(self, out_file, wh, flag, deleteFile=None):
        ''' Subroutine to implement Error/Warning status '''
        wh.write("\nif [ ! -s "+out_file+" ]; then \n")
        wh.write("   echo -e \"ERROR: "+out_file+" doesnot exist\" \n")
        wh.write("   exit 1\n")
        wh.write("else\n")

        if deleteFile !=None:
            if isinstance(deleteFile,list):
                for ele in deleteFile:
                    wh.write('   echo \"Removing File:'+ele+'\"\n')
                    wh.write('   rm '+ele+'\n')
            else:
                wh.write('   echo \"Removing File:'+deleteFile+'\"\n')
                wh.write('   rm '+deleteFile+'\n')
        else:
            wh.write('   echo \"Removing File:'+str(deleteFile)+'\"\n')
        
        wh.write("\nfi\n\n")
        
        return wh

    def checkErrorDir(self, out_file, wh, flag, deleteFile=None):
        ''' Subroutine to implement Error/Warning status '''
        wh.write("\nif [ ! -d "+out_file+" ]; then \n")
        wh.write("   echo -e \"ERROR: "+out_file+" doesnot exist\" \n")
        wh.write("   exit 1\n")
        wh.write("else\n")

        if deleteFile !=None:
            if isinstance(deleteFile,list):
                for ele in deleteFile:
                    wh.write('   echo \"Removing File:'+ele+'\"\n')
                    wh.write('   rm '+ele+'\n')
            else:
                wh.write('   echo \"Removing File:'+deleteFile+'\"\n')
                wh.write('   rm '+deleteFile+'\n')
        else:
            wh.write('   echo \"Removing File:'+str(deleteFile)+'\"\n')
        
        wh.write("\nfi\n\n")
        
        return wh


    def checkOutputFile(self, out_file, wh, flag, deleteFile=None):
        ''' Subroutine to check if the File exist. Incorporated
        To avoid re run of the command '''

        if_cmd = "\nif [ ! -s "+out_file+" ]; then \n"
        echo_1 = "   echo -e \"ERROR: "+out_file+" doesnot exist\" \n"
        echo_2 = "   echo \"0\" | cat > \""+status_file+"\"\n"
        echo_3 = "   exit 1\n"
        echo_4 = "   else\n"
        echo_5 = "          echo \"1\" | cat > \""+status_file+"\""
        end_cmd = "\nfi\n\n"
        wh.write(if_cmd+echo_1+echo_2+echo_3+echo_4+echo_5+end_cmd)
        return wh


class picard:
    ''' Generic class for implementing GATK commands '''
    
    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def writeCramToBam(self,configDict,sampleWorkDir,pat_id,cramFile,
                       deleteFile=None,wh=None):
        ''' Subroutine to write commands for converting CRAM to BAM file '''
        
        refFile = '/'.join([configDict['general']['resourceDir'],
                            configDict['resources']['reference']['grch37']
                           ])
        bamFileBase = re.split('.cram',os.path.basename(cramFile))[0]+'.bam'
        bamFile = '/'.join([sampleWorkDir,'tmp',bamFileBase])
        cmd = ' '.join(['samtools view -T',refFile,'-@ 12',cramFile,
                        '-b >',bamFile,'\n'
                       ])
        cmdIndex = ' '.join(['samtools index',bamFile])

        #Create Sample directory
        
        wh.write('mkdir -p '+'/'.join([sampleWorkDir,'tmp']))
        wh.write('\n')

        # Write the command
        wh.write('######### CRAM to BAM Conversion ##########\n')
        wh.write('echo "CRAM to BAM conversion"\n\n')
        
        if_cmd = "\nif [ ! -s "+bamFile+" ]; then \n"
        wh.write(if_cmd)
        wh.write(cmd)
        wh.write(cmdIndex)
        wh.write('\n')
        wh.write('fi\n\n')

        # Write the Error/Warning status
        objC = cluster()
        wh = objC.checkErrorFile(bamFile,wh,'error',None)

        return bamFile, wh

    def writeRevertSam(self,configDict,sampleWorkDir,cramFile,bamFile,
                                             deleteFile=None,wh=None):
        ''' Convert the input BAM file to UnMAPPED BAM per read group '''
        
        revSamDict = configDict['picard']['RevertSam']
        outRGDir = '/'.join([sampleWorkDir,'tmp'])
        #os.system('mkdir -p '+outRGDir)
        #wh.write('mkdir -p '+outRGDir)
        #wh.write('\n\n')

        picardLib = configDict['general']['picardLib']

        cmdList = ['java -Xmx10G -jar',picardLib,'RevertSam']
        revSamDict['I'] = bamFile
        revSamDict['O'] = outRGDir
        revSamDict['TMP_DIR']=outRGDir

        for key,value in revSamDict.items():
            cmdList.append(key+'='+value)

        wh.write('\n########## RevertSam Step ##########\n')
        wh.write('echo "Convert from BAM to UnMapped"\n\n')
        wh.write(' '.join(cmdList))
        wh.write('\n')

        # Get RG BAm file id
        rgBamListCmd = "samtools view -H "+cramFile+"|grep '^@RG' |"+\
                    "awk '{print $2}'"+"| awk -F ':' '{print $2}'"

        result = subprocess.Popen(rgBamListCmd,shell=True,stdout=subprocess.PIPE)
        subprocess_return = result.stdout.read()
        rgIdList = re.split('\n',subprocess_return.decode("utf-8"))
        rgIdList = [x.strip() for x in rgIdList if x]
        
        rgBamList = []
        for rg_id in rgIdList:
            #os.system('mkdir -p '+'/'.join([outRGDir,rg_id])) 
            wh.write('mkdir -p '+'/'.join([outRGDir,rg_id]))
            wh.write('\n')
            rgBamFile = '/'.join([outRGDir,rg_id+'.bam'])
            mvrgBamFile = '/'.join([outRGDir,rg_id,rg_id+'.bam'])
            wh.write('mv '+rgBamFile+' '+mvrgBamFile)
            wh.write('\n')

            objC = cluster()
            wh = objC.checkErrorFile(mvrgBamFile,wh,'error',None)
            rgBamList.append(mvrgBamFile)

        wh.write('\n')
        #wh.write('rm '+bamFile+'\n')
            
        return rgIdList,rgBamList,outRGDir,wh


    def writeMarkIlluminaAdapters(self,configDict,outRGDir,rgBamList,
                                  rgIdList,deleteFile=None,wh=None):

        ''' Subroutine to Mark Illumina Adapter to Un-Mapped BAM file '''

        wh.write('\n########## MarkIlluminaAdapter ##########\n')
        wh.write('echo "MarK Illumina Adapter step"\n\n')
       
        outMiaBamList = []
        rgIdCmdList = []

        for rgBamFile,rgId in zip(rgBamList,rgIdList):
            miaDict = configDict['picard']['MarkIlluminaAdapters']
            picardLib = configDict['general']['picardLib']
            
            miaDict['I'] = rgBamFile
            miaDict['O'] = re.split('.bam',rgBamFile)[0]+'.mia.bam'
            outMiaBamList.append(miaDict['O'])
            miaDict['M'] = re.split('.bam',rgBamFile)[0]+'.mia_metrics.txt'
            miaDict['TMP_DIR'] = '/'.join([outRGDir,rgId])

            cmdList = []
            cmdList = ['java -Xmx6G -jar',picardLib,'MarkIlluminaAdapters']
            for key,value in miaDict.items():
                cmdList.append(key+'='+value)

            #if_cmd = "\nif [ ! -s "+miaDict['O']+" ]; then \n"
            #wh.write(if_cmd)
            rgIdCmdList.append(' '.join(cmdList))
        
        wh.write('echo " '+' | '.join(rgIdCmdList)+' "| parallel')

        objC = cluster()
        for outBamFile in outMiaBamList:        
            wh = objC.checkErrorFile(outBamFile,wh,'error',None)

        return outMiaBamList,wh

    def writeSamToFastq(self,configDict,outRGDir,miaBamList,
                        rgIdList,deleteFile=None,wh=None):

        ''' Subroutine to convert Sam to Fastq file '''
        wh.write('\n########## Sam To Fastq ##########\n')
        wh.write('echo "Sam To Fastq step"\n\n')
     
        stfDict = configDict['picard']['SamToFastq']
        picardLib = configDict['general']['picardLib']
        stfCmdList = []
        outStfList = []

        for miaBamFile,rgId in zip(miaBamList,rgIdList):
            stfDict['I'] = miaBamFile
            stfDict['F'] = re.split('.bam',miaBamFile)[0]+'.fq.gz'
            stfDict['TMP_DIR'] = '/'.join([outRGDir,rgId])
            outStfList.append(stfDict['F'])

            cmdList = []
            cmdList = ['java -Xmx6G -jar',picardLib,'SamToFastq']
            for key,value in stfDict.items():
                cmdList.append(key+'='+value)

            stfCmdList.append(' '.join(cmdList))
        
        wh.write('echo " '+' | '.join(stfCmdList)+' "| parallel')

        objC = cluster()

        for outStfFile,miaBamFile in zip(outStfList,miaBamList):
            #wh = objC.checkErrorFile(outStfFile,wh,'error',miaBamFile)
            wh = objC.checkErrorFile(outStfFile,wh,'error',None)

        return outStfList, wh

    def writeBwaAlign(self,configDict,outRGDir,fqGzList,rgIdList,
                      deleteFile=None,wh=None):
        ''' Subroutine to Mark Illumina Adapter to Un-Mapped BAM file '''
        
        wh.write('\n########## BWA Alignment ##########\n')
        wh.write('echo "BWA alignment step"\n\n')
         
        bwaDict = configDict['bwa']['mem']
        picardLib = configDict['general']['picardLib']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        bwaCmdList = []
        outBwaList = []
        for fqGzFile,rgId in zip(fqGzList,rgIdList):
            bwaDict['o'] = re.split('.gz',fqGzFile)[0]+'.aln.bam'
            cmdList = []
            cmdList = ['bwa mem']
            for key,value in bwaDict.items():
                if value !=None:
                    cmdList.append('-'+key+' '+str(value))
                else:
                    cmdList.append('-'+key)

            cmdList.insert(5,ref38)
            cmdList.insert(6,fqGzFile)

            bwaCmdList.append(' '.join(cmdList))
            outBwaList.append(bwaDict['o'])
            #if_cmd = "\nif [ ! -s "+bwaDict['o']+" ]; then \n"
            #wh.write(if_cmd)
            wh.write(' '.join(cmdList))
            wh.write('\n')
            #wh.write('fi\n\n')
        #wh.write('echo " '+' | '.join(bwaCmdList)+' "| parallel')
            objC = cluster()
            #wh = objC.checkErrorFile(bwaDict['o'],wh,'error',fqGzFile)
            wh = objC.checkErrorFile(bwaDict['o'],wh,'error',None)
            
            #wh.write('rm '+fqGzFile+'\n\n')

        return outBwaList,wh

    def writeMergeBamAlignment(self,configDict,rgBamList,bwaAlnList,
                                outRGDir,deleteFile=None,wh=None):
        ''' Subroutine to convert Sam to Fastq file '''
        
        #outRGDir = os.path.dirname(rgBamFile)
        mbDict = configDict['picard']['MergeBamAlignment']
        picardLib = configDict['general']['picardLib']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        
        wh.write('\n########## Merge Bam Alignment ##########\n')
        wh.write('echo "Merge Aligned BAM to Unmapped BAM"\n\n')
        
        cmdBamAlnList = []
        outMergeBamList = []

        for bwaAlnFile,rgBamFile in zip(bwaAlnList,rgBamList):
            mbDict['ALIGNED'] = bwaAlnFile
            mbDict['UNMAPPED'] = rgBamFile
            mbDict['O'] = re.split('.bam',bwaAlnFile)[0]+'.mumap.bam'
            mbDict['R'] = ref38 
            rgId = os.path.basename(os.path.dirname(bwaAlnFile))
            mbDict['TMP_DIR'] = '/'.join([outRGDir,rgId])

            cmdList = []
            cmdList = ['java -Xmx6G -jar',picardLib,'MergeBamAlignment']
            for key,value in mbDict.items():
                cmdList.append(key+'='+value)

            cmdBamAlnList.append(' '.join(cmdList))
            outMergeBamList.append(mbDict['O'])
        wh.write('echo " '+' | '.join(cmdBamAlnList)+' "| parallel')
            
        objC = cluster()
        for outBwaAlnUmapFile,bwaAlnFile in zip(outMergeBamList,bwaAlnList):
            #wh = objC.checkErrorFile(outBwaAlnUmapFile,wh,'error',bwaAlnFile)
            wh = objC.checkErrorFile(outBwaAlnUmapFile,wh,'error',None)
            
        return outMergeBamList,wh

    def writeMergeRGSamFiles(self,configDict,mergeBamList,outRGDir,pat_id,
                                        deleteFile=None,wh=None):
        
        ''' Subroutine to Sam files per read group '''
        
        mbDict = configDict['picard']['MergeSamFiles']
        picardLib = configDict['general']['picardLib']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        
        mbDict['I'] = mergeBamList
        mbDict['O'] = '/'.join([outRGDir,pat_id])+'.fq.hg38.bam'
        mbDict['R'] = ref38 
        mbDict['TMP_DIR'] = outRGDir

        cmdList = []
        #cmdList = ['java -Xmx30G -jar',picardLib,'MergeSamFiles']
        cmdList = ['samtools merge ',mbDict['O']]+mergeBamList+['-@ 10']
        ''' 
        for key,value in mbDict.items():
            if key=='I':
                for bam in mergeBamList:
                    cmdList.append(key+'='+bam)
            else:
                cmdList.append(key+'='+value)
        '''
        wh.write('\n########## Merge Bam Alignment ##########\n')
        wh.write('echo "Merge Aligned Bam per RG"\n\n')
        if_cmd = "\nif [ ! -s "+mbDict['O']+" ]; then \n"
        wh.write(if_cmd)
        wh.write(' '.join(cmdList))
        wh.write('\nsamtools index '+mbDict['O'])
        wh.write('\n')
        wh.write('fi\n\n')
        
        objC = cluster()
        #wh = objC.checkErrorFile(mbDict['O'],wh,'error',mergeBamList)
        wh = objC.checkErrorFile(mbDict['O'],wh,'error',None)

        return mbDict['O'],wh

    def writeMarkDuplicates(self,configDict,mergeRGSamFile,deleteFile=None,
                            wh=None):
        ''' Subroutine to Mark Duplicates '''
        
        dedupDict = configDict['picard']['MarkDuplicates']
        picardLib = configDict['general']['picardLib']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        
        dedupDict['I'] = mergeRGSamFile
        dedupDict['O'] = re.split('.bam',mergeRGSamFile)[0]+'.dedup.bam'
        dedupDict['TMP_DIR'] = os.path.dirname(mergeRGSamFile)

        cmdList = []
        cmdList = ['java -Xmx4G -jar',picardLib,'MarkDuplicates']
        for key,value in dedupDict.items():
            cmdList.append(key+'='+value)

        wh.write('\n########## Picard Mark Duplicates ##########\n')
        wh.write('echo "Picard mark duplicates "\n\n')
        if_cmd = "\nif [ ! -s "+dedupDict['O']+" ]; then \n"
        wh.write(if_cmd)
        wh.write(' '.join(cmdList))
        wh.write('\n')
        wh.write('fi\n\n')
        
        objC = cluster()
        #wh = objC.checkErrorFile(dedupDict['O'],wh,'error',mergeRGSamFile)
        wh = objC.checkErrorFile(dedupDict['O'],wh,'error',None)

        return dedupDict['O'],wh

    def writeGatherVcfs(self,configDict,chrGvcfList,pat_id,outRGDir,
                    deleteFile=None,wh=None):
        ''' Subroutine to Gather Chromosome GVCFs '''

                
        chrDict = configDict['picard']['GatherVcfs']

        picardLib = configDict['general']['picardLib']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        
        chrDict['I'] = chrGvcfList
        chrDict['O'] = '/'.join([os.path.abspath(outRGDir),
                                   pat_id+'.jgt.vcf.gz'
                                  ])
        chrDict['TMP_DIR'] = outRGDir
        chrDict['R'] = ref38

        cmdList = []
        cmdList = ['java -Xmx8G -jar',picardLib,'GatherVcfs']
        for key,value in chrDict.items():
            if key=="I":
                for ele in value:
                    cmdList.append(key+'='+ele+' \\'+'\n')
            else:
                cmdList.append(key+'='+value)

        wh.write('\n########## Picard Gather VCFs ##########\n')
        wh.write('echo "Picard Gather VCFs "\n\n')
        wh.write(' '.join(cmdList))
        wh.write('\n')
        
        objC = cluster()
        #wh = objC.checkErrorFile(chrDict['O'],wh,'error',chrGvcfList)
        wh = objC.checkErrorFile(chrDict['O'],wh,'error',None)
        #cmd_index = 'bcftools index '+chrDict['O']
        #wh.write(cmd_index)
        wh.write('\n')

        return chrDict['O'],wh

    def writeSortGatherVcf(self,configDict,gatherVcfFile,chrNum,
                        outDir,deleteFile=None,wh=None):

        ''' Subroutine to Sort the input VCF file '''

        sortDict = configDict['picard']['SortVcf']
        picardLib = configDict['general']['picardLib']
        sortDict['INPUT'] = gatherVcfFile
        sortDict['TMP_DIR'] = outDir
        sortDict['OUTPUT'] = '/'.join([outDir,chrNum+'.jgt.sorted.vcf.gz'])

        cmdList = ['java -jar',picardLib,'SortVcf']
        for key,value in sortDict.items():
            cmdList.append(key+'='+value+' \\'+'\n')

        wh.write('\n####### Picard Sort Vcf ############\n')
        wh.write('echo "Picard sort VCF by coordinates"\n\n')
        wh.write(' '.join(cmdList))
        wh.write('\n')
        
        objC = cluster()
        wh = objC.checkErrorFile(sortDict['OUTPUT'],wh,'error',None)
        wh.write('\n')

        wh.write('echo "Generating Index of the Sorted VCF"\n\n')
        wh.write('gatk IndexFeatureFile -I '+sortDict['OUTPUT'])
        wh.write('\n')

        return sortDict['OUTPUT'],wh


    def writeMergeChromBamFiles(self,configDict,chromBamList,outRGDir,pat_id,
                                                        tmpStatFile,wh):
        
        ''' Subroutine to Sam files per read group '''
        
        mbDict = configDict['picard']['MergeSamFiles']
        picardLib = configDict['general']['picardLib']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        
        mbDict['I'] = chromBamList
        mbDict['O'] = '/'.join([os.path.dirname(outRGDir),pat_id])+'.hg38.recal.bam'
        mbDict['R'] = ref38 
        mbDict['TMP_DIR'] = outRGDir

        cmdList = []
        cmdList = ['java -Xmx30G -jar',picardLib,'MergeSamFiles']
        for key,value in mbDict.items():
            if key=='I':
                for bam in chromBamList:
                    cmdList.append(key+'='+bam)
            else:
                cmdList.append(key+'='+value)

        wh.write('\n########## Merge Chromosome Bam Alignment ##########\n')
        wh.write('echo "Merge Chromosome BAM output files"\n\n')
        wh.write(' '.join(cmdList))
        wh.write('\n')
        
        objC = cluster()
        wh = objC.checkErrorFile(mbDict['O'],wh,'error',tmpStatFile)

        return mbDict['O'],wh


class gatk:
    ''' Generic class for implementing GATK commands '''


    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1


    def writeBaseRecalibrator(self,configDict,dedupFile,outRGDir,
                              deleteFile=None,wh=None):
        ''' Subroutine to implement GATK BaseRecalibrator step'''


        recalDict = configDict['gatk']['BaseRecalibrator']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        
        recalDict['input'] = dedupFile
        dedupBase = os.path.basename(dedupFile)
        recalDict['output'] = '/'.join([outRGDir,
                                        re.split('.bam$|.cram$',dedupBase)[0]+'.recal.table'
                                       ])
        recalDict['reference'] = ref38
        #recalDict['intervals'] = chrNum
        #recalDict['TMP_DIR'] = os.path.dirname(dedupFile)

        cmdList = []
        cmdList = ['gatk --java-options "-Xmx10G"','BaseRecalibrator']
        for key,value in recalDict.items():
            if re.search('^known',key):
                for ele in value:
                    ele = '/'.join([configDict['general']['resourceDir'],ele])
                    cmdList.append('--'+key+' '+ele)
            else:
                cmdList.append('--'+key+' '+value)

        wh.write('\n########## GATK BaseRecalibrator ##########\n')
        wh.write('echo "GATK BaseRecalibrator "\n\n')
        wh.write(' '.join(cmdList))
        wh.write('\n')
        
        objC = cluster()
        wh = objC.checkErrorFile(recalDict['output'],wh,'error',None)

        return recalDict['output'],wh

    def writeApplyBQSR(self,configDict,dedupFile,baseRecalFile,pat_id,
                       outRGDir,deleteFile=None,wh=None):

        ''' Subroutine to implement GATK BQSR step'''

        bqsrDict = configDict['gatk']['ApplyBQSR']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        
        bqsrDict['input'] = dedupFile
        bqsrDict['output'] ='/'.join([outRGDir,
                                      pat_id+'.hg38.dedup.recal.bam'
                                      ])
        bqsrDict['reference'] = ref38
        bqsrDict['bqsr-recal-file'] = baseRecalFile
        bqsrDict['tmp-dir'] = outRGDir
        #bqsrDict['intervals'] = chrNum

        cramFile = '/'.join([outRGDir,
                            pat_id+'.hg38.dedup.recal.cram'
                           ])
        cmdList = []
        cmdList = ['gatk --java-options "-Xmx8G"','ApplyBQSR']
        for key,value in bqsrDict.items():
            cmdList.append('--'+key+' '+value)

        wh.write('\n########## GATK ApplyBQSR ##########\n')
        wh.write('echo "GATK ApplyBQSR "\n\n')
        wh.write(' '.join(cmdList))
        wh.write('\n\n')

        objC = cluster()
        
        #wh = objC.checkErrorFile(bqsrDict['output'],wh,'error',dedupFile)
        wh = objC.checkErrorFile(bqsrDict['output'],wh,'error',None)
       
        ''' 
        wh.write('echo "Convert BAM to CRAM"\n\n')
        cmdCram = ['samtools view -T',ref38,'-C -o',cramFile,bqsrDict['output']]
        wh.write(' '.join(cmdCram))
        wh.write('\n\n')
        cmdCramIndex = 'samtools index '+cramFile
        wh.write(cmdCramIndex)
        wh.write('\n\n')
        wh = objC.checkErrorFile(cramFile,wh,'error',bqsrDict['output'])
        return cramFile,wh
        '''
        return bqsrDict['output'],wh

    def writeHaplotypeCaller(self,configDict,bqsrFile,pat_id,chrList,
                             outRGDir,deleteFile=None,wh=None):

        ''' Subroutine to implement GATK HpalotypeCaller in single 
        sample modes'''

        hcDict = configDict['gatk']['HaplotypeCaller']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        
        wh.write('\n########## GATK HaplotypeCaller - Single Mode ##########\n')
        wh.write('echo "GATK HaplotypeCaller "\n\n')
   
        outHCFileList = []

        '''  
        chrDict = {'1-5':chrList[0:5],
                   '6-10':chrList[5:10],
                   '11-15':chrList[10:15],
                   '16-20':chrList[15:20],
                   '21-Y':chrList[20:24]
                  }

        for keys,values in chrDict.items():
            hcCmdList = []
            wh.write('echo " Processing :'+' '.join(values)+' "\n')
            for chrNum in values:
        '''
        hcDict['I'] = bqsrFile
        hcDict['O']  ='/'.join([outRGDir,pat_id+'.g.vcf.gz'])
        hcDict['R'] = ref38
         
        hcDict['native-pair-hmm-threads'] = configDict['lsf']['params1']['n']
        hcDict['tmp-dir'] = outRGDir
        
        cmdList = []
        cmdList = ['gatk --java-options "-Xmx2G"','HaplotypeCaller']
        for key,value in hcDict.items():
            if key=='G':
                for ele in value:
                    cmdList.append('-'+key+' '+str(ele))
            elif key=='D':
                cmdList.append('-'+key+' '+ \
                               '/'.join([configDict['general']['resourceDir'],
                                           hcDict['D']
                                          ])
                              )
            elif key=='L':
                cmdList.append('-'+key+' '+\
                               '/'.join([configDict['general']['resourceDir'],
                                         hcDict['L']
                                        ])
                              )
            else:
                cmdList.append('-'+key+' '+str(value))

        #hcCmdList.append(' '.join(cmdList))
        #outHCFileList.append(hcDict['O'])
        #wh.write('echo "'+' | '.join(hcCmdList)+' " | parallel')
        wh.write(' '.join(cmdList))
        wh.write('\n\n')
        
        objC = cluster()
       
        #for chrGvcf in outHCFileList:
        wh = objC.checkErrorFile(hcDict['O'],wh,'error',None)

        #destPath = 
        #mvCmd = ['mv ',bqsrFile,

        return hcDict['O'],wh

    def writeGatherVcfs(self,chrGvcfList,pat_id,outRGDir,deleteFile=None,
                        wh=None):
        ''' Subroutine to gather all the Chromsome gvcf '''

        print('Inside GatherVCF\n')
        gtDict = configDict['gatk']['GatherVcfs']
        gtDict['I'] = chrGvcfList
        gtDict['O'] = '/'.join([os.path.dirname(outRGDir),
                                pat_id+'all.g.vcf.gz'
                               ])
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        gtDict['R'] = ref38
        gtDict['TMP_DIR'] = outRGDir

        cmdList = ['gatk --java-options "-Xmx8G"','GatherVcfs']
        for key,value in gtDict.items():
            if key=='I':
                for chrGvcfFile in chrGvcfList:
                    cmdList.append('-'+key+' '+chrGvcfFile)
            else:
                cmdList.append('-key'+' '+str(value))
        wh.write('\n########## Gather all HC Vcfs ##########\n')
        wh.write('echo "Gather all GC VCFs step"\n\n')
        wh.write(' '.join(cmdList))
        wh.write('\n')

        print('BeforePrint')
        print(chrGvcfList)
        objC = cluster()
        #wh = objC.checkErrorFile(gtDict['O'],wh,'error',chrGvcfList)
        wh = objC.checkErrorFile(gtDict['O'],wh,'error',None)

        return gtDict['O'],wh

    def writeGenomicsDBImport(self,configDict,gvcfFileList,outDir,chrNum,
                                                        chrIntervalFile,wh):
        ''' Subroutine to implement Genomics DB Import step 
        '''
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
 
        outDbName = 'genomicsDB_'+chrIntervalFile
        gdiDict = configDict['gatk']['GenomicsDBImport']
        gdiDict['genomicsdb-workspace-path'] = '/'.join([outDir,outDbName])
        gdiDict['tmp-dir'] = outDir

        if re.search('ddd',gvcfFileList[0]):
            gdiDict['intervals'] = '/'.join([configDict['general']['resourceDir'],
                                         'regions/grch38/split/'+\
                                         chrNum+'/'+chrIntervalFile+'.bed'
                                        ])
        elif re.search('ukbb',gvcfFileList[0]):
            if re.search('50K',gvcfFileList[0]):
                gdiDict['intervals'] =\
                '/'.join(['/nfs/research/dunham/samples/ukbb/data/bed/split',
                                             chrNum+'/'+chrIntervalFile+'.bed'
                                            ])
            elif re.search('200K|200k',gvcfFileList[0]):
                gdiDict['intervals'] =\
                '/'.join(['/nfs/research/dunham/resources/regions/grch38/ukbb/split',
                                             chrNum+'/'+chrIntervalFile+'.bed'
                                            ])
 

 
        #gdiDict['intervals'] = chrNum
        if not re.search('ukbb',gvcfFileList[0]):
            gdiDict['reference'] = ref38
        else:
            gdiDict['reference'] =\
            '/nfs/research/dunham/resources/gatk_bundle/ukbb_fe/genome.fa'
        
        #cmdList = ['gatk --java-options "-Xmx16G"','GenomicsDBImport']
        cmdList = ['gatk --java-options "-Xmx4G"','GenomicsDBImport']

        for key,value in gdiDict.items():
            if key=='variant':
                for gvcfFile in gvcfFileList:
                    cmdList.append('--'+key+' '+gvcfFile)
            else:
                cmdList.append('--'+key+' '+str(value))

        wh.write('\n########### Genomics DB Import ###############\n')
        wh.write('echo "Genomics DB Import"\n\n')
        if_cmd = "\nif [ ! -d "+gdiDict['genomicsdb-workspace-path']+" ]; then \n"
        wh.write(if_cmd)
        wh.write('\n')
        wh.write(' '.join(cmdList))
        wh.write('\n')
        wh.write('fi\n\n')

        objC = cluster()
        wh = objC.checkErrorDir(gdiDict['genomicsdb-workspace-path'],
                                 wh,'error',None 
                                )
        return gdiDict['genomicsdb-workspace-path'],wh

    def writeCombineGVCFs(self,configDict,gvcfFileList,outDir,chrNum,wh):
        ''' Subroutine to implement CombineGVCFs step 
        '''

        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        outGvcfName = chrNum+'_comb.gvcf.vcf.gz'
        gdiDict = configDict['gatk']['CombineGVCFs']
        gdiDict['output'] = '/'.join([outDir,outGvcfName])
        #gdiDict['tmp-dir'] = outDir
        gdiDict['reference'] = ref38
        #gdiDict['intervals'] = chrNum
        gdiDict['intervals'] = '/'.join([configDict['general']['resourceDir'],
                                         'regions/grch38/gt_dbi/'+\
                                         'df2017.lov.hg38.'+chrNum+'.bed'
                                        ])
         
        cmdList = ['gatk --java-options "-Xmx16G"','CombineGVCFs']


        for key,value in gdiDict.items():
            if key=='variant':
                for gvcfFile in gvcfFileList:
                    cmdList.append('--'+key+' '+gvcfFile)
            else:
                cmdList.append('--'+key+' '+str(value))

        wh.write('\n########### CombineGVCFs step ###############\n')
        wh.write('echo "CombineGVCFs step "\n\n')
        wh.write(' '.join(cmdList))
        wh.write('\n')

        objC = cluster()
        wh = objC.checkErrorFile(gdiDict['output'],wh,'error',None)

        return gdiDict['output'],wh

    def writeEditHCGvcf(self,configDict,dddId,hcGvcfFile,wh):
        ''' Subroutine to edit GVCF header; Rename Sample names '''

        gvcfDir = os.path.dirname(hcGvcfFile)

        outRenFile = re.split('.g.vcf.gz',hcGvcfFile)[0]+'.ren.g.vcf.gz'
        outEditFile = re.split('.g.vcf.gz',hcGvcfFile)[0]+'.edit.g.vcf.gz'
        
        wh.write('\n############ Sample rename & Edit the Header ###########\n')
        wh.write('echo "Sample rename"\n')
        tmp_s_file = '/'.join([gvcfDir,'sample.txt'])
        cat_cmd = 'echo "'+dddId+'" > '+tmp_s_file
        cmd_ren = 'bcftools reheader -s '+tmp_s_file+' -o '+outRenFile+\
                    ' '+hcGvcfFile
        wh.write(cat_cmd)
        wh.write('\n\n')

        wh.write(cmd_ren)
        wh.write('\n\n')

        wh.write('echo "Edit the VCF header"\n')
        cmd_header = "bcftools view "+outRenFile+\
                    " | sed 's/##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,"+\
                    "Number=A/##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=./'"+\
                    " | bgzip -c > "+outEditFile

        wh.write(cmd_header)
        wh.write('\n')

        return outEditFile,wh


    def writeGenotypeGVCF(self,configDict,chrDbImportFile,outDir,chrNum,
                                                            chrIntervalFile,wh):
        ''' Subroutine to implement joint genotyping of gvcf files 
        '''

        gtDict = configDict['gatk']['GenotypeGVCF']

        if not re.search('ukbb',outDir):
            ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        else:
            ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['ukb_grch38_fe']
                         ])
         
        gtDict['R'] = ref38
        gtDict['V'] = 'gendb://'+chrDbImportFile
        gtDict['O'] = '/'.join([outDir,'jointGT_gvcf.'+chrIntervalFile+'.vcf.gz'])
        #gtDict['L'] = chrNum
        #gtDict['L'] = '/'.join([configDict['general']['resourceDir'],
        #                        'regions/grch38/gt_dbi/'+chrNum+\
        #                        '/'+chrIntervalFile+'.bed'
        #                       ])
 
        cmdList = ['gatk --java-options "-Xmx4G"','GenotypeGVCFs']


        for key,value in gtDict.items():
            if key=='D':
                value = '/'.join([configDict['general']['resourceDir'],
                                  value
                                 ])
                cmdList.append('-'+key+' '+str(value))
            else:
                cmdList.append('-'+key+' '+str(value))

        wh.write('\n########### Genotype GVCFs step ###############\n')
        wh.write('echo "Genotype GVCFs step "\n\n')
        wh.write(' '.join(cmdList))
        wh.write('\n')

        objC = cluster()
        wh = objC.checkErrorFile(gtDict['O'],wh,'error',None)

        return gtDict['O'],wh

    def writeExtractSample(self,jointGTFile,tmp_data,s_name,dddId,wh):
        ''' Subroutine to extract samples from multi-sample GVCF File 
        '''
        
        outVcfDir = '/'.join([tmp_data,dddId])
        outVcfFile = '/'.join([outVcfDir,dddId+'.gt.vcf.gz'])

        cmd_mkdir = 'mkdir -p '+outVcfDir
        cmd_bcf = ' '.join(['bcftools view -Oz -s',s_name,jointGTFile,'>',
                            outVcfFile
                           ])
        cmd_index = 'tabix -p vcf '+outVcfFile

        wh.write('\n############## Extract Sample Step ############\n')
        wh.write('echo "Extract Sample: '+s_name+' "\n\n')
        wh.write(cmd_mkdir)
        wh.write('\n')
        wh.write(cmd_bcf)
        wh.write('\n')
        wh.write(cmd_index)
        wh.write('\n\n')
        
        objC = cluster()
        wh = objC.checkErrorFile(outVcfFile,wh,'error',None)

        return outVcfFile,wh


    def writeVariantRecal(self,configDict,sgtFile,dddId,wh):
        ''' Subroutine to implement Variant Recalibrator Step '''

        vrDict = configDict['gatk']['VariantRecalibrator']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        resDir = configDict['general']['resourceDir']

        vrDict['reference'] = ref38
        vrDict['variant'] = sgtFile
        vrDict['tranches-file'] = '/'.join([os.path.dirname(sgtFile),
                                            dddId+'.jgt.sorted.recal.tranches'
                                           ])
        '''
        vrDict['rscript-file'] = '/'.join([os.path.dirname(sgtFile),
                                           dddId+'.gt.recal.plots.AS.R'
                                          ])
        '''
        vrDict['output'] = '/'.join([os.path.dirname(sgtFile),
                                     dddId+'.jgt.sorted.recal.vcf.gz'
                                    ])

        cmdList = ['gatk VariantRecalibrator']

        for key,value in vrDict.items():
            if key=='resource':
                for ele in vrDict[key]['value']:
                    ele_strs = re.split('\:',ele)
                    ele_strs[1] = '/'.join([resDir,'gatk_bundle',ele_strs[1]])
                    ele = ' '.join(ele_strs)
                    cmdList.append('--'+key+':'+ele)
            elif key=='an':
                for ele in vrDict[key]['value']:
                    cmdList.append('-'+key+' '+ele)
            else:
                cmdList.append('--'+key+' '+value)

        wh.write('\n############## Variant Recalibrator Step ############\n')
        wh.write('echo "Variant Recalibrator Step "\n')
        wh.write(' '.join(cmdList))
        wh.write('\n\n')

        objC = cluster()
        wh = objC.checkErrorFile(vrDict['output'],wh,'error',None)

        return vrDict['output'],wh


    def writeApplyVQSR(self,configDict,sgtFile,outRecalFile,dddId,wh):
        ''' Subroutine that implement Apply VQSR file '''

        avDict = configDict['gatk']['ApplyVQSR']
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        resDir = configDict['general']['resourceDir']
 
        avDict['reference'] = ref38
        avDict['variant'] = sgtFile
        avDict['tranches-file'] = '/'.join([os.path.dirname(sgtFile),
                                            dddId+'.jgt.sorted.recal.tranches'
                                           ])
        avDict['output'] = '/'.join([os.path.dirname(sgtFile),
                                     dddId+'.jgt.sorted.recal.vqsr.vcf.gz'
                                    ])
        avDict['recal-file'] = outRecalFile

        cmdList = ['gatk ApplyVQSR']

        for key,value in avDict.items():
           cmdList.append('--'+key+' '+value)
    
        wh.write('\n############## ApplyVQSR Step ############\n')
        wh.write('echo "ApplyVQSR Step "\n')
        wh.write(' '.join(cmdList))
        wh.write('\n\n')

        return avDict['output'],wh

class annotation:
    ''' Generic class for annotating SNVs '''

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def annotSNV(self,configDict):
        pass


    def writeAddAcAnTags(self,jgtRecalFile,outChrDir,chrNum,wh):
        ''' Subroutine to split by chromosome the JointGenotyped file'''

        outFile = '/'.join([outChrDir,'tmpFile_'+'tags.bcf'])
        bcftools = '/hps/software/users/dunham/aak/anaconda3/envs/bcf/bin/bcftools'
        cmd_tags = bcftools+' +fill-tags '+ jgtRecalFile+ ' -Ob -o '+outFile+' -- -t AN,AC'
        wh.write('\n######### Add INFO/AC, INFO/AN fields  #########\n')
        wh.write('\necho "Add AC/aN tags: "\n\n')
        wh.write(cmd_tags)
        wh.write('\n')

        cmd_tags_index = bcftools+' index '+outFile
        wh.write(cmd_tags_index)
        wh.write('\n')

        objC = cluster()
        wh = objC.checkErrorFile(outFile,wh,'error',None)

        return outFile,wh


    def writeSplitByChrom(self,jgtRecalFile,outChrDir,chrNum,wh):
        ''' Subroutine to split by chromosome the JointGenotyped file'''

        outFile = '/'.join([outChrDir,'tmpFile_'+chrNum+'.bcf'])
        cmd_split_chr = 'bcftools view -r '+str(chrNum)+' '+jgtRecalFile+\
                        ' -Ob -o '+outFile
        wh.write('\n############ Split JGT file by chromosome #########\n')
        wh.write('\necho "Split the JGT file by chromosome: '+chrNum+' "\n\n')
        wh.write(cmd_split_chr)
        wh.write('\n')

        cmd_split_chr_index = 'bcftools index '+outFile
        wh.write(cmd_split_chr_index)
        wh.write('\n')

        objC = cluster()
        wh = objC.checkErrorFile(outFile,wh,'error',None)

        return outFile,wh

    def writeSplitByRegions(self,jgtRecalFile,outChrDir,regionFile,region,wh):
        ''' Subroutine to split by chromosome the JointGenotyped file'''

        outFile = '/'.join([outChrDir,'tmpFile_'+region+'.bcf'])
        cmd_split_chr = 'bcftools view -R '+regionFile+' '+jgtRecalFile+\
                        ' -Ob -o '+outFile
        wh.write('\n############ Split JGT file by chromosome #########\n')
        wh.write('\necho "Split the JGT file by chromosome: '+region+' "\n\n')
        wh.write(cmd_split_chr)
        wh.write('\n')

        cmd_split_chr_index = 'bcftools index '+outFile
        wh.write(cmd_split_chr_index)
        wh.write('\n')

        objC = cluster()
        wh = objC.checkErrorFile(outFile,wh,'error',None)

        return outFile,wh

    def writeNormalizeVCF(self,configDict,chrSplitOutFile,outChrDir,chrNum,wh):
        ''' Subroutine to Normalize indels, split multi-allelic sites by chromosome 
        the JointGenotyped file'''

        outNormFile = re.split('.bcf',chrSplitOutFile)[0]+'.AC0.norm.bcf'

        if re.search('ddd',chrSplitOutFile):
            ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])
        elif re.search('ukbb|ukb|sfari',chrSplitOutFile):
            if re.search('200k|100k',chrSplitOutFile):
                ref38 = '/'.join([configDict['general']['resourceDir'],
                              configDict['resources']['reference']['grch38']
                             ])
            else:
                ref38 = '/'.join([configDict['general']['resourceDir'],
                              configDict['resources']['reference']['ukb_grch38_fe']
                             ])

        cmd_norm = 'bcftools norm -m - -f '+ref38+' '+chrSplitOutFile+\
                        " | bcftools view -e 'AC==0' -Ob -o "+outNormFile

        wh.write('\n####### Normalize,Split Multi-allelic,purge variants AC=0 #########\n')
        wh.write('\necho "[`date` ]Normalize,Split,Multi-allelic,purge AC=0: "\n\n')
        wh.write(cmd_norm)
        wh.write('\n')

        cmd_norm_index = 'bcftools index '+outNormFile
        wh.write(cmd_norm_index)
        wh.write('\n')

        objC = cluster()
        wh = objC.checkErrorFile(outNormFile,wh,'error',None)
        #wh = objC.checkErrorFile(outNormFile,wh,'error',chrSplitOutFile+'*')

        return outNormFile,wh


    def writeVEPannotation(self,configDict,normOutFile,outChrDir,chrNum,wh):
        ''' Write subroutine to execute VEP annotation '''

        if re.search('bcf$',normOutFile):
            prefix_path = re.split('\.bcf',normOutFile)[0]
            #basename = re.split('.bcf$',os.path.basename(normOutFile))[0]
            #prefix_path = '/'.join([outChrDir,basename])
        elif re.search('.vcf.gz$',normOutFile):
            basename = re.split('.noGT.vcf.gz$',os.path.basename(normOutFile))[0]
            prefix_path = '/'.join([outChrDir,basename])

        tmp1 = prefix_path+".noGT.vcf.gz"
        tmp2 = prefix_path+".noGT.new.vcf.gz"
        vep_file_new = prefix_path+".noGT.new.vep.vcf.gz"
        stats_file = prefix_path+".noGT.vep_stats.txt"
        tmp3 = prefix_path+".vep.tmp.bcf"
        exon_anno_file_001 = prefix_path+".vep.bcf"

        wh.write("\n\n########### START: VEP Annotation #################\n")
        wh.write('\necho \"[`date`] Getting vcf without GT and pulling new'+\
                 'variants\"\n')

        ### Commented out for running VEP
        cmd_gt = 'bcftools view -G '+normOutFile+" -Oz -o "+tmp1
        cmd_gt_index = 'bcftools index '+tmp1
        wh.write(cmd_gt+'\n')
        wh.write(cmd_gt_index+'\n')
        
        vep_file = tmp1
        tmp2 = tmp1

        wh.write('\necho \"[`date`]Running VEP\"\n')
        ref38 = '/'.join([configDict['general']['resourceDir'],
                          configDict['resources']['reference']['grch38']
                         ])

        refGenome = configDict['general']['genomeBuild']
        defPar = [];allPar = []; plugIn=[]

        resDir = configDict['general']['resourceDir']
        outformat = configDict['annotation']['outFormat']
        defDict = configDict['annotation']['defPar']
        allDict = configDict['annotation']['all']
        plugInDict = configDict['annotation']['plugIn']
        doubDash = configDict["annotation"]["doubDash"]
        vep_bin = configDict['annotation']['vep-bin']

        for keys,values in defDict.items():
            if re.search('^dir\_cache',keys):
                defPar.append('--'+keys)
                defPar.append('/'.join([os.path.abspath(resDir),values,refGenome]))
            elif keys=='fasta':
                defPar.append('--'+keys)
                defPar.append(ref38)
            elif defDict[keys]==None:
                defPar.append('--'+keys)
            else:
                defPar.append('--'+keys)
                defPar.append(defDict[keys])
        
        for keys in allDict:
            if allDict[keys]==None:
                allPar.append('--'+keys)
            else:
                allPar.append(doubDash+keys)
                allPar.append(allDict[keys])


        for keys,values in plugInDict.items():
            if re.search('CADD',keys):
                plugin_str = ' '.join(['--plugin',','.join([keys,'/'.join([
                                                            resDir,values['cadd_snv']
                                                            ]),
                                                            '/'.join([resDir,
                                                                      values['cadd_indel']
                                                                     ])
                                                            ])
                                      ])
            elif re.search('SpliceAI',keys):
                plugin_str = ' '.join(['--plugin',','.join([keys,'snv='+'/'.join([
                                                            resDir,values['snv']
                                                            ]),
                                                            'indel='+'/'.join([resDir,
                                                                      values['indel']
                                                                     ])
                                                            ])
                                      ])
            elif re.search('dbNSFP',keys):
                plugin_str = ' '.join(['--plugin',','.join([keys,'/'.join([
                                                            resDir,values['txt']
                                                            ]),
                                                              values['fields']
                                                            
                                                            ])
                                      ])
            elif re.search('G2P',keys):
                plugin_str = ' '.join(['--plugin',','.join([keys,'file='+'/'.join([
                                                            resDir,values['csv']
                                                            ])
                                                            ])
                                      ])
            elif re.search('LoFtool',keys):
                plugin_str = ' '.join(['--plugin',','.join([keys,'/'.join([
                                                            resDir,values['txt']
                                                            ])
                                                            ])
                                      ])
            elif re.search('AlphaMissense',keys):
                plugin_str = ' '.join(['--plugin',','.join([keys,'/'.join(['file='+
                                                            resDir,values['txt']
                                                            ])
                                                            ])
                                      ])
  
            elif values !='.':
                plugin_str = ' '.join(['--plugin',','.join([keys,'/'.join([
                                                            resDir,values
                                                            ])
                                                        ])
                                    ])
            else:
                plugin_str = ' '.join(['--plugin',keys])

            plugIn.append(plugin_str)

 
        cmd = 'bcftools view '+tmp2+' | '+vep_bin+' --'+\
                outformat+' -o STDOUT '+' '.join(defPar)+' '+' '.join(allPar)+\
                ' '+' '.join(plugIn)
        cmd = cmd+' --stats_text --stats_file '+stats_file
        cmd = cmd+' | bgzip -c > '+vep_file_new
        cmd_index = 'bcftools index -f '+vep_file_new
        wh.write(cmd+'\n')
        wh.write(cmd_index+'\n')

        wh.write("\necho \"[`date`] Annotating bcf with VEP results\"\n")

        if re.search('.vcf.gz$',normOutFile):
            normOutFile = normOutFile
        elif re.search('.norm.bcf$',normOutFile):
            normOutFile = tmp1
        cmd_vep_res = '\nbcftools annotate -c INFO/ANN -a '+vep_file_new+" "+\
                       normOutFile+" -Ob -o "+exon_anno_file_001
        cmd_vep_res_index = 'bcftools index '+exon_anno_file_001
        wh.write(cmd_vep_res+'\n')
        wh.write(cmd_vep_res_index+'\n')

        objC = cluster()
        wh = objC.checkErrorFile(exon_anno_file_001,wh,'error',None)
        #wh = objC.checkErrorFile(exon_anno_file_001,wh,'error',tmp1+'*')
        #wh = objC.checkErrorFile(exon_anno_file_001,wh,'error',vep_file_new+'*')

        return exon_anno_file_001,wh


    def writeClusterCustomAnnot(self,configDict,vepOutFile,outChrDir,chrNum,wh):
        
        ''' Subroutine to add custom annotations : gnomAD,ExAC,HGMD, Clinvar '''

        resDir = configDict["general"]["resourceDir"]
        prefix_path = re.split("\.bcf",vepOutFile)[0]
        tmp_file = prefix_path+".tmp.bcf"
        exon_tmp_anno_file = prefix_path+".anno.bcf"
        
        cust_annot = configDict["annotation"]["cust-annot"]
        cmd_all = []

        chr_strs = re.split('chr',chrNum)
        if len(chr_strs)==2:
            chr_num = chr_strs[1]
        else:
            chr_num = chr_strs[0]
            #print(chr_num)


        first = 1
        for ann in cust_annot:
            if first:
                inp = vepOutFile
                first = 0
            else:
                inp = exon_tmp_anno_file
                
            if "tag" in cust_annot[ann]:
                c = ",".join(["INFO/" + cust_annot[ann]["tag"] + "_" + i +\
                              ":=INFO/" + i for i in cust_annot[ann]["fields"].split(',')
                             ])
            elif cust_annot[ann]["fields"] == "INFO":
                c = cust_annot[ann]["fields"]
            else:
                c = ",".join(["INFO/" + i for i in cust_annot[ann]["fields"].split(',')])

            if "vcfDir" in cust_annot[ann]:
                vcfDir = cust_annot[ann]["vcfDir"]
            else:
                vcfDir = resDir
           
            #if re.search('clinvar',ann,re.IGNORECASE):

            if ann=='GNOMADV3':
                if re.search('^chr',chrNum):
                    vcfFile = re.split('chr\%s',cust_annot[ann]['vcf']['hg38'])[0]+chrNum+'.vcf.gz'
                else:
                    vcfFile = re.split('%s',cust_annot[ann]['vcf']['grch38'])[0]+chrNum+'.vcf.bgz'
            elif ann=='GNOMADg':
                if re.search('^chr',chrNum):
                    chr_str = re.split('^chr',chrNum)[1]
                    vcfFile = re.split('%s',cust_annot[ann]['vcf']['hg38'])[0]+chr_str+'.liftover_grch38.vcf.bgz'
                else:
                    vcfFile = re.split('%s',cust_annot[ann]['vcf']['grch38'])[0]+chrNum+'.liftover_grch38.vcf.bgz'
            elif ann=='GNOMADe':
                if re.search('^chr',chrNum):
                    chr_str = re.split('^chr',chrNum)[1]
                    vcfFile = re.split('%s',cust_annot[ann]['vcf']['hg38'])[0]+chr_str+'.liftover_grch38.vcf.bgz'
                else:
                    vcfFile = re.split('%s',cust_annot[ann]['vcf']['grch38'])[0]+chrNum+'.liftover_grch38.vcf.bgz'
                            
            elif ann=='EXAC':
                if re.search('^chr',chrNum):
                    vcfFile = cust_annot[ann]['vcf']['hg38'] 
                else:
                    vcfFile = cust_annot[ann]['vcf']['grch38']
            elif ann=='GNOMADe4':
                if re.search('^chr',chrNum):
                    chr_str = re.split('^chr',chrNum)[1]
                    vcfFile = re.split('%s',cust_annot[ann]['vcf']['hg38'])[0]+chr_str+'.vcf.bgz'
                else:
                    vcfFile = re.split('%s',cust_annot[ann]['vcf']['grch38'])[0]+chrNum+'.vcf.bgz'
            elif ann=='GNOMADg4':
                if re.search('^chr',chrNum):
                    chr_str = re.split('^chr',chrNum)[1]
                    vcfFile = re.split('%s',cust_annot[ann]['vcf']['hg38'])[0]+chr_str+'.vcf.bgz'
                else:
                    vcfFile = re.split('%s',cust_annot[ann]['vcf']['grch38'])[0]+chrNum+'.vcf.bgz'
            elif ann=='CLINVAR':
                if re.search('^chr',chrNum):
                    vcfFile = cust_annot[ann]['vcf']['hg38'] 
                else:
                    vcfFile = cust_annot[ann]['vcf']['grch38']
            
            ''' 
            if "chr" in cust_annot[ann]:
                vcfFile = cust_annot[ann]["vcf"] % chr_num
                if "excl_chr" in cust_annot[ann] and chr_num in cust_annot[ann]["excl_chr"].split(","):
                    if "repl_chr" in cust_annot[ann]:
                        vcfFile = cust_annot[ann]["vcf"] % cust_annot[ann]["repl_chr"]
                    else:
                        continue 
            '''
            cmd = "bcftools annotate -c " + c + " -a " + vcfDir + "/" +\
                    vcfFile + " -r " + chrNum + " " + inp + " -Ob -o " + tmp_file
            cmd_all.append(cmd)

        wh.write("\n\n######## START: Add Custom Annotations ############# \n")
        wh.write("\necho \"[`date`] Adding custom annotations\"\n")

        for cmd in cmd_all:
            #20220526: For GDD-rerun: Clinvar annotation missing
            #if re.search('clinvar',cmd,re.IGNORECASE): 
            wh.write(cmd+'\n')
            wh.write("mv -f %s %s" % (tmp_file, exon_tmp_anno_file))
            wh.write('\n')
            wh.write("bcftools index -f %s\n" % exon_tmp_anno_file)
            wh.write('\n')

        wh.write("\n\n######## END: Add Custom Annotations ######### \n")
        
        return exon_tmp_anno_file, wh


    def writeClusterICFreq(self,configDict,vepCustAnnoFile,outChrDir,chrNum,
                           projDate,wh):
        ''' Subroutine to add PASS Flags and Internal Cohort frequency '''

        prefix_path = re.split("\.bcf",vepCustAnnoFile)[0]
        pf_in = vepCustAnnoFile
        pf_out = prefix_path+".sites.bcf"

        ic_str = []
        if re.search('gdd',projDate):
            ic_str = 'GDD'
            fr_tmp = prefix_path+".sites.gdd.freq.tmp"
            fr_out = prefix_path+".sites.gdd.freq.bcf"
        elif re.search('ukbb',projDate):
            if re.search('200k',projDate):
                ic_str = 'UKB200K'
            elif re.search('500k',projDate):
                ic_str = 'UKB500K'
            else:
                ic_str = 'UKB50K'
            fr_tmp = prefix_path+".sites.ukbb.freq.tmp"
            fr_out = prefix_path+".sites.ukbb.freq.bcf"
        elif re.search('sfari',projDate):
            ic_str = 'sfari'
            fr_tmp = prefix_path+".sites.sfari.freq.tmp"
            fr_out = prefix_path+".sites.sfari.freq.bcf"

        wh.write("\n\n###### START: Add Flags, add internal cohort frequency steps ######\n")
        wh.write('\necho \"[`date`] Adding pass flags\"\n')
        
        cmd1 = " bcftools annotate -x FILTER -Ou "+pf_in+\
                " | bcftools +fill-AN-AC | bcftools filter -e "+\
                "\'AN<(1*N_SAMPLES)\' -s \'LOWCALL\' -m 'LOWPF\' -m + -Ob > "+pf_out
        cmd1_index = "\nbcftools index "+pf_out

        wh.write(cmd1+'\n')
        wh.write(cmd1_index+'\n')

        objC = cluster()
        wh = objC.checkErrorFile(pf_out,wh,'error',None)
        #wh = objC.checkErrorFile(pf_out,wh,'error',vepCustAnnoFile+"*")

        cmd4 = "bcftools view "+pf_out
        cmd4_sed1 = "sed '/^#CHROM/i ##INFO=<ID="+ic_str+"_AF,Number=A,Type=Float,Description=\""+ic_str+" alternate allele frequencies\">'"
        cmd4_sed2 = "sed '/^#CHROM/i ##INFO=<ID="+ic_str+"_AC,Number=A,Type=Integer,Description=\""+ic_str+" allele count in genotypes\">'"
        cmd4_sed3 = "sed '/^#CHROM/i ##INFO=<ID="+ic_str+"_AN,Number=1,Type=Integer,Description=\""+ic_str+" total number of alleles in called genotypes\">'"
        cmd4_sed4 = "perl -F\\\\t -ane 'print && next if /^#/; @in = split /;/,$F[7]; for $i (@in) {$an = $1 if $i =~ /^AN=([\d]+)/; $ac = $1 if $i =~/^AC=([\d]+)/;} $af = 0; $af =$ac/$an if $an > 0;$F[7] .= \";"+ic_str+"_AN=$an;"+ic_str+"_AC=$ac;"+ic_str+"_AF=$af\"; print join(\"\\t\",@F)'"
        cmd4_sed5 = "bcftools view - -Ob -o "+fr_tmp
        
        cmd4 = " | ".join([cmd4,cmd4_sed1,cmd4_sed2,cmd4_sed3,cmd4_sed4,cmd4_sed5])
        cmd4_index = "\nbcftools index "+fr_tmp

        ### 05.09.2023: Commeted out ICF for UKBB ###
        #wh.write(cmd4+'\n')
        #wh.write(cmd4_index+'\n')
        #wh.write("\necho \"[`date`] Annotating original file with Internal Cohort AFs\"\n")

        cmd5 = "\nbcftools annotate -Ob -o "+fr_out+" -c INFO/"+ic_str+"_AC,INFO/"+ic_str+"_AN,INFO/"+ic_str+"_AF -a "+fr_tmp+" "+pf_out
        cmd5_index = "\nbcftools index "+fr_out

        #### 05.09.2023: Commenting out the ICG frequency as it is same as JGT Step #####
        
        #wh.write(cmd5+'\n')
        #wh.write(cmd5_index+'\n')
        #wh.write('rm '+fr_tmp+'*')
        #wh.write('\n')

        #wh = objC.checkErrorFile(fr_out,wh,'error',None)
        #wh = objC.checkErrorFile(fr_out,wh,'error',pf_out)

        #return fr_out,wh
        return pf_out,wh

    def writeExtractExonicRegion(self,configDict,anno_bcf,outChrDir,wh):
        ''' Subroutine to extract exonic regions fromt the annotated VCF file '''

        if re.search('gdd|ukb200k|ukb',anno_bcf):
            exonFile = '/'.join([configDict['general']['resourceDir'],
                                 configDict['resources']['ensembl']['hg38']
                                ])
        elif re.search('ukb50k',anno_bcf):
            exonFile = '/'.join([configDict['general']['resourceDir'],
                                 configDict['resources']['ensembl']['grch38']
                                ])
 
        anno_bcf_baseName = os.path.basename(anno_bcf)
        exon_bcf = '/'.join([outChrDir,
                             re.split('.bcf',anno_bcf_baseName)[0]+'.exon.bcf'
                           ])
        exon_sort_bcf = '/'.join([outChrDir,
                             re.split('.bcf',anno_bcf_baseName)[0]+'.exon.sorted.bcf'
                           ])

        wh.write("\n######## Extract Exonic Region ######\n")
        wh.write("echo 'Extract exonic region'\n\n")
        cmd_exon = ' '.join(['bcftools view','-R',exonFile,anno_bcf,'-Ob -o',exon_bcf])

        cmd_sort = ' '.join(['bcftools view '+exon_bcf+\
                             '| bcftools sort -m 1G -T '+outChrDir,
                             '| bcftools norm -D ',
                             '| bcftools view -Ob -o '+exon_sort_bcf
                            ])

        cmd_exon_index = 'bcftools index '+exon_bcf
        wh.write(cmd_exon)
        wh.write('\n')
        wh.write(cmd_sort+'\n')
        cmd_mv = 'mv '+exon_sort_bcf+' '+exon_bcf
        wh.write(cmd_mv+'\n')

        wh.write(cmd_exon_index)
        wh.write('\n')
        
        objC = cluster()
        wh = objC.checkErrorFile(exon_bcf,wh,'error',None)
        #wh = objC.checkErrorFile(exon_bcf,wh,'error',anno_bcf)

        return exon_bcf,wh

    def writeVepAnnoMergeGT(self,configDict,vepAnnoICFreq,normOutFile,outChrDir,wh):
        ''' Subroutine to merge annoated noGT with GT Bcf'''

        prefix_path = re.split('\.bcf',vepAnnoICFreq)[0]
        outMergeFile = prefix_path+'.merge.bcf'
        cmd_merge = ['bcftools annotate -a ',vepAnnoICFreq,normOutFile,
                     '-c ID,FILTER,INFO','-Ob -o ',outMergeFile]

        cmd_merge_index = 'bcftools index '+outMergeFile
        wh.write('#### Merging the Annotation without GT to with GT ####\n')
        wh.write('echo "Merging the annotation with GT"\n')
        wh.write(' '.join(cmd_merge)+'\n')
        wh.write(cmd_merge_index+'\n')

        objC = cluster()
        #wh = objC.checkErrorFile(outMergeFile,wh,'error',vepAnnoICFreq+"*")
        #wh = objC.checkErrorFile(outMergeFile,wh,'error',normOutFile+"*")
        wh = objC.checkErrorFile(outMergeFile,wh,'error',None)

        return outMergeFile,wh

    def writeFilterMAF(self,configDict,chr_exon_bcf,maf_cutoff,wh):
        ''' Subroutine to filter the annotated variants w.r.t GNOMAD,EXAC with
        MAF threshold'''
        
        maf_bcf = re.split('.bcf',chr_exon_bcf)[0]+'.0_01.bcf'

        #cmd_maf = ' '.join(["bcftools view -Ob -e 'GNOMADgV3_AF>"+maf_cutoff+"\'",
        #                    chr_exon_bcf+"|",
        #                    "bcftools view -Ob -e 'GNOMADg_AF>"+maf_cutoff+"\'",
        #                    "| bcftools view -Ob -e 'GNOMADe_AF>"+maf_cutoff+"\'",
        #                    "| bcftools view -Ob -e 'EXAC_AF>"+maf_cutoff+"\'",
        #                    ">",maf_bcf
        #                   ])

        # 07.11.2023: Added this code for filtering with gnomADv4.0
        if not re.search('chrY',chr_exon_bcf):
            cmd_maf = ' '.join(["bcftools view -Ob -e 'GNOMADg4_AF>"+maf_cutoff+"\'",
                            chr_exon_bcf+"|",
                            "bcftools view -Ob -e 'GNOMADe4_AF>"+maf_cutoff+"\'",
                            ">",maf_bcf
                           ])
        elif re.search('chrY',chr_exon_bcf):
            cmd_maf = ' '.join(["bcftools view -Ob -e 'GNOMADe4_AF>"+maf_cutoff+"\'",
                            chr_exon_bcf + " >",maf_bcf
                           ])
 
        cmd_maf_index = 'bcftools index '+maf_bcf

        wh.write("\n######## Frequency Filters ######\n")
        wh.write("echo 'MAF Filter <=0.01 or 1% '\n\n")
 
        wh.write(cmd_maf+'\n')
        wh.write(cmd_maf_index+'\n')

        objC = cluster()
        wh = objC.checkErrorFile(maf_bcf,wh,'error',None)
        #wh = objC.checkErrorFile(maf_bcf,wh,'error',chr_exon_bcf)

        return maf_bcf,wh

    def writeFilterImpact(self,configDict,maf_bcf,wh):
        ''' Subroutine for applying Impact filters '''

        imp_bcf = re.split('.bcf',maf_bcf)[0]+'.imp.bcf'
        
        cmd_imp = ' '.join(["bcftools view",maf_bcf,
                            "| grep -E '^#|^MT|MODERATE|HIGH|"+\
                            "splice_donor_5th_base_variant|"+\
                            "splice_donor_0th_base_variant|"+\
                            "splice_donor_region_variant|"+\
                            "splice_polypyrimidine_tract_variant|HGMD|CLNSIG'"+\
                            "| bcftools view -Ob -o",imp_bcf
                            ])
        cmd_imp_ind = 'bcftools index '+imp_bcf

        wh.write("\n######## Impact Filters ######\n")
        wh.write("echo 'Impact Filters' \n\n")
         
        wh.write(cmd_imp+'\n')
        wh.write(cmd_imp_ind+'\n')

        objC = cluster()
        wh = objC.checkErrorFile(imp_bcf,wh,'error',None)
        #wh = objC.checkErrorFile(imp_bcf,wh,'error',maf_bcf)

        return imp_bcf,wh

    def writeMildSevereSamples(self,configDict,chr_anno_bcf,outChrDir,ms_type,wh):
        ''' Subroutine to extract severe and mild samples '''
       
        chr_anno_base = os.path.basename(chr_anno_bcf)
        ms_anno_bcf = '/'.join([outChrDir,
                               re.split('.bcf',chr_anno_base)[0]
                              ])+'.'+ms_type+'.bcf'
        if re.search('mild',ms_type):
            mms_samples = '/nfs/research/dunham/samples/ddd/data/DDD_ID_mild_GDD.txt'
        elif re.search('severe',ms_type):
            mms_samples = '/nfs/research/dunham/samples/ddd/data/DDD_ID_severe_GDD.txt'
        elif re.search('mod',ms_type):
            mms_samples = '/nfs/research/dunham/samples/ddd/data/DDD_ID_moderate_GDD.txt'
        elif re.search('general',ms_type):
            mms_samples = '/nfs/research/dunham/samples/ddd/data/DDD_ID_general_GDD.txt'

        cmd_ms = ' '.join(['bcftools view -S',mms_samples,
                            '--force-samples',chr_anno_bcf,
                           '-Ob -o',ms_anno_bcf
                          ])
        
        cmd_index = 'bcftools index '+ms_anno_bcf
        wh.write(cmd_ms+'\n')
        wh.write(cmd_index+'\n')

        objC = cluster()
        wh = objC.checkErrorFile(ms_anno_bcf,wh,'error',None)

        return ms_anno_bcf,wh

    def writeExtractUKBB(self,configDict,maf_filter_bcf,outchrDir,batchNum,wh):
        ''' Subroutine to extract UKBB batch specific samples '''
       
        imp_batch_base = os.path.basename(maf_filter_bcf)
        imp_batch_bcf = '/'.join([outchrDir,
                                  re.split('.bcf',imp_batch_base)[0]
                                 ])+'.'+batchNum+'.vcf'

        if re.search('ukb50k',outchrDir):
            batch_samples = '/nfs/research/dunham/samples/ukbb/data/batches/'+batchNum
        elif re.search('ukb200k',outchrDir):
            batch_samples = '/nfs/research/dunham/samples/ukbb/data/200k/hpc_manifest/'+batchNum
 
        cmd_batch = ' '.join(['bcftools view -S',batch_samples,
                            '--force-samples',maf_filter_bcf,
                           '-Ov -o',imp_batch_bcf
                          ])
        cmd_index = 'bcftools index '+imp_batch_bcf
        wh.write(cmd_batch+'\n')
        wh.write(cmd_index+'\n')
        objC = cluster()
        wh = objC.checkErrorFile(imp_batch_bcf,wh,'error',None)

        return imp_batch_bcf,wh
      
    def writeConvertVcf2Tab(self,configDict,imp_filter_bcf,wh):
        ''' Subroutine to convert Vcf to Tab '''

        vcf_2_tab = re.split('.vcf',imp_filter_bcf)[0]+'.tab'

        """ 
        cmd_tab = ' '.join(['python /homes/aak/scripts/AI-UK/src/vcf2tab.v2.py',
                           '-t',imp_filter_bcf,'-o',vcf_2_tab
                           ]) 
        """
        cmd_tab = ' '.join(['python /homes/aak/scripts/AI-UK/src/vcf2tab.v2.py',
                           imp_filter_bcf,'-o',vcf_2_tab
                           ]) 


        wh.write("\n######## Convert VCF 2 TAB ######\n")
        wh.write("echo 'Convert Vcf 2 Tab' \n\n")
         
        wh.write(cmd_tab+'\n')
        objC = cluster()
        wh = objC.checkErrorFile(vcf_2_tab,wh,'error',None)
        #wh = objC.checkErrorFile(vcf_2_tab,wh,'error',imp_filter_bcf)

        return vcf_2_tab,wh

    def writeCompressBgz(self,configDict,vcf_2_tab,wh):
        ''' Subroutine to compress vcf2tab file bgz '''
        
        #vcf_2_tab_bgz = re.split('.tab',vcf_2_tab)[0]+'.tab.bgz'
        vcf_2_tab_bgz = re.split('.tab',vcf_2_tab)[0]+'.tab.bg.gz'

        """ 
        cmd_tab = ' '.join(['python /homes/aak/scripts/AI-UK/src/vcf2tab.v2.py',
                           '-t',imp_filter_bcf,'-o',vcf_2_tab
                           ]) 
        """
        cmd_tab = ' '.join(['bgzip -c ',vcf_2_tab,' > ',vcf_2_tab_bgz]) 

        wh.write("\n######## Compress VCF-TAB to BGZ ######\n")
        wh.write("echo 'Compress Vcf 2 Tab tp Bgz' \n\n")
         
        wh.write(cmd_tab+'\n')
        objC = cluster()
        #wh = objC.checkErrorFile(vcf_2_tab_bgz,wh,'error',vcf_2_tab)
        wh = objC.checkErrorFile(vcf_2_tab_bgz,wh,'error',None)

        return vcf_2_tab_bgz,wh



    def writeTranscriptPrioritization(self,configDict,vcf_2_tab,chrNum,
                                      projDate,memDict,wh,ms_type=[],isCustom=[]):
        ''' Subroutine to Prioritize Transcripts '''

        if isCustom:
            if re.search('chrY|Y',chrNum):
                configDict['lsf']['params1']['M'] = '25G'
                configDict['lsf']['params1']['q'] = 'research'
            else:
                #configDict['lsf']['params1']['M'] = memDict[chrNum]
                configDict['lsf']['params1']['M'] = memDict[chrNum]
                configDict['lsf']['params1']['q'] = 'bigmem'

            configDict['lsf']['params1']['n'] = str(20)

            tmp_dir_path =  os.path.dirname(
                                os.path.dirname(
                                    os.path.dirname(vcf_2_tab)
                                )
                            )
            tmp_bin = '/'.join([tmp_dir_path,'tmp_binaries'])
            tmp_data = '/'.join([tmp_dir_path,'tmp_data'])
            sb_log = '/'.join([tmp_dir_path,'sb_log'])

            sh_type = 'LSF'
            if not ms_type:
                fileId = chrNum+'_TransPrior'
            else:
                fileId = chrNum+'_TransPrior_'+ms_type

            objC = cluster()
            cluster_file,wh = objC.getClusterWriteHandle(tmp_bin,fileId,projDate,sh_type)
            wh = objC.writeClusterTop(configDict,'NDD-WES-GATK-TransPrior-'+chrNum,wh)
            sh_out,sh_err,wh = objC.writeClusterInit(configDict,
                                                     sb_log,sh_type,
                                                     fileId,wh
                                                    )

            #wh = objC.writeclusterSpecific(wh)
            wh.write('\nmodule load r-4.0.3-gcc-9.3.0-xiarbub\n')
            wh.write('echo "$HOSTNAME"\n\n')
        elif isCustom==False:
            pass

        trs_prior_tab = re.split('.tab',vcf_2_tab)[0]+'.prior.tab'
        transLen = '/'.join([configDict['general']['resourceDir'],
                            configDict['resources']['ensembl']['transcriptLength']
                            ])
        
        cmd_tr = ' '.join(['Rscript',
                           '/homes/aak/scripts/AI-UK/src/variantPrioritization.v3.R',
                           '-v',vcf_2_tab,'-o',trs_prior_tab,
                            '-c 20','-l',transLen
                          ])

        wh.write("\n######## Transcript Prioritization ######\n")
        wh.write("echo 'Transcript Prioritization' \n\n")
 
        wh.write(cmd_tr+'\n')
        objC = cluster()
        wh = objC.checkErrorFile(trs_prior_tab,wh,'error',None)

        wh.close()

        return trs_prior_tab

    def writeQueryGenesUkbbDv(self,tabFile,chrNum,outTabDir,wh):
        ''' Subroutine to query Genes in DV-called UKBB dataset '''
        
        batchNum = os.path.basename(outTabDir)
        wh.write('######## Querying genes per Batches ##### \n\n')
        wh.write('echo "Querying '+batchNum+'"\n')
        
        if re.search('^chr',chrNum):
            strs = re.split('chr',chrNum)
            chrNumStr = strs[1]

        r_script = '/homes/aak/scripts/AI-UK/src/dv/getUKBB_DV.R'
        cmd_query = ['Rscript',r_script,tabFile,outTabDir,chrNumStr]
        wh.write(' '.join(cmd_query))
        wh.write('\n\n')

        return wh

    def writeQueryGenesUKBB(self,configDict,vcf_2_tab,chrNum,
                                      projDate,memDict,ms_type=[]):

        ''' Subroutine to Prioritize Transcripts '''

        if re.search('ukb50k|ukb200k|ukb300k',projDate):
            chr_num = chrNum
            script = '/homes/aak/scripts/AI-UK/src/getUKBB.R'
        elif re.search('gdd',projDate):
            chr_num = re.split('chr',chrNum)[1]
            script = '/homes/aak/scripts/AI-UK/src/getGDD.R'

        if re.search('chrY|Y',chrNum):
            configDict['lsf']['params1']['M'] = '10G'
            configDict['lsf']['params1']['q'] = 'research'
        else:
            #configDict['lsf']['params1']['M'] = memDict[chrNum]
            configDict['lsf']['params1']['q'] = 'research'

        configDict['lsf']['params1']['n'] = str(10)

        tmp_dir_path =  os.path.dirname(
                            os.path.dirname(
                                os.path.dirname(vcf_2_tab)
                            )
                        )
        tmp_bin = '/'.join([tmp_dir_path,'tmp_binaries'])
        tmp_data = '/'.join([tmp_dir_path,'tmp_data'])
        sb_log = '/'.join([tmp_dir_path,'sb_log'])
        out_dir = '/'.join([tmp_dir_path,'id_genes',chrNum])

        sh_type = 'LSF'
        if not ms_type:
            fileId = 'QueryCohort_'+chrNum
        else:
            fileId = 'QueryCohort_'+ms_type+'_'+chrNum

        objC = cluster()
        cluster_file,wh = objC.getClusterWriteHandle(tmp_bin,fileId,sh_type)
        wh = objC.writeClusterTop(configDict,'NDD-WES-GATK-QueryCohort-'+chrNum,wh)
        sh_out,sh_err,wh = objC.writeClusterInit(configDict,
                                                 sb_log,sh_type,
                                                 fileId,wh
                                                )

        #wh = objC.writeclusterSpecific(wh)
        wh.write('\nmodule load r-4.0.3-gcc-9.3.0-xiarbub\n')
        wh.write('echo "$HOSTNAME"\n\n')
        wh.write('mkdir -p '+out_dir)

        cmd_tr = ' '.join(['Rscript',
                           script,
                           vcf_2_tab,chr_num,out_dir
                          ])
        
        #print(cmd_tr)

        wh.write("\n######## Query ID/High-Capabilities genes ######\n")
        wh.write("echo 'Query genes' \n\n")
 
        wh.write(cmd_tr+'\n')
        #objC = cluster()
        #wh = objC.checkErrorFile(trs_prior_tab,wh,'error',None)

        wh.close()
        #return trs_prior_tab

    def writeCombChromTab(self,configDict,inp_dir,objC,maf_cutoff,case_type,
                          projDate,sh_type):
        ''' Subroutine to Combine all transcript prioritized variants split
        across chromosome '''

        maf_cutoff = '_'.join(re.split('\.',maf_cutoff))
        tmp_dir_path = os.path.dirname(inp_dir)
        tmp_bin = '/'.join([tmp_dir_path,'tmp_binaries'])
        tmp_data = '/'.join([tmp_dir_path,'tmp_data'])
        merge_data = '/'.join([tmp_dir_path,'merge_data'])
        sb_log = '/'.join([tmp_dir_path,'sb_log'])

        os.system('mkdir -p '+merge_data)
        fileId = 'combineAllChromTab_'+case_type
        configDict['lsf']['params1']['M'] = '1100G'
        configDict['lsf']['params1']['q'] = 'bigmem'

        cluster_file,wh = objC.getClusterWriteHandle(tmp_bin,fileId,projDate,sh_type)
        wh = objC.writeClusterTop(configDict,'NDD-WES-GATK-CombChrom-'+case_type,wh)
        sh_out,sh_err,wh = objC.writeClusterInit(configDict,
                                                 sb_log,sh_type,
                                                 fileId,wh
                                                )

        #wh = objC.writeclusterSpecific(wh)
        wh.write('\nmodule load r-4.0.3-gcc-9.3.0-xiarbub\n')
        wh.write('echo "$HOSTNAME"\n\n')

        wh.write('mkdir -p '+merge_data+'\n')

        script_path = os.path.dirname(os.path.dirname(__file__))
        
        combChromBin = '/'.join([script_path,'src/combinePriorTab.R'])

        merge_vcf_file = []
        merge_tab_file = []
        
        if re.search('ukbb',projDate):
            db_str = 'ukbb'

            configDict['combChromTab']['bcfSuffix'] =\
                    '/'.join([inp_dir,"*","tmpFile_*."+db_str+"*"+case_type+"*"+\
                              maf_cutoff+".bcf"
                             ])
            configDict['combChromTab']['tabSuffix'] =\
                    '/'.join([inp_dir,"*","tmpFile_*."+db_str+"*"+case_type+"*"+\
                              maf_cutoff+"*.tab"
                             ])
            merge_vcf_file = '/'.join([merge_data,
                                      'merged.all.AC0.norm.anno.sites.'+\
                                       db_str+'.freq.'+case_type+'.exon.'+maf_cutoff+'.bcf'])
            merge_tab_file = '/'.join([merge_data,
                                      'merged.all.AC0.norm.anno.sites.'+\
                                       db_str+'.freq.'+case_type+'.exon.'+maf_cutoff+\
                                       '.imp.prior.tab'])
        else :
            db_str = 'gdd'
            configDict['combChromTab']['bcfSuffix'] =\
                    '/'.join([inp_dir,"chr*","tmpFile_*."+db_str+"*"+case_type+"*"+\
                              maf_cutoff+".bcf"
                             ])
            configDict['combChromTab']['tabSuffix'] =\
                    '/'.join([inp_dir,"chr*","tmpFile_*."+db_str+"*"+case_type+"*"+\
                              maf_cutoff+".imp.prior.tab"
                             ])
            merge_vcf_file = '/'.join([merge_data,
                                      'merged.all.AC0.norm.anno.sites.'+\
                                       db_str+'.freq.'+case_type+'.exon.'+maf_cutoff+'.bcf'
                                      ])
            merge_tab_file = '/'.join([merge_data,
                                      'merged.all.AC0.norm.anno.sites.'+\
                                       db_str+'.freq.'+case_type+'.exon.'+maf_cutoff+\
                                       '.imp.prior.tab'
                                      ])

        #print(configDict['combChromTab']['bcfSuffix'])
        cmd_vcf = "bcftools concat `ls -v "+\
                   configDict['combChromTab']['bcfSuffix']+"` -Ob -o "+\
                   merge_vcf_file 

        cmd_comb = ' '.join(['Rscript',combChromBin,
                             tmp_data,merge_tab_file,
                             case_type
                            ])
        
        wh.write('\n#### Merging Impact Tab file ####\n')
        wh.write("echo 'Merging Impact tab file'\n")
        wh.write(cmd_comb+'\n')
        wh = objC.checkErrorFile(merge_tab_file,wh,'error',None)

        wh.write('\n### Combining Frequency Variants for All Chromsomes ####\n')
        wh.write("echo 'Merging Frequency file'\n")
        wh.write(cmd_vcf+'\n')
        wh = objC.checkErrorFile(merge_vcf_file,wh,'error',None)
        wh.close()

        return merge_tab_file, merge_vcf_file

class utility:
    
    ''' Specific Class to process GVCF files '''

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def writeExonCoverage(self,jgtRecalFile,exon_coord,chrNum,outExonChrCov,wh):
        
        ''' Subroutine to extract coverage or RD per variant within given  '''

        cmd = ' '.join(['bcftools query -r',exon_coord,'-f',
                        '\'%CHROM+%REF+%ALT+%POS\\t[%DP\\t]\\n\'',jgtRecalFile,
                        '> ',outExonChrCov
                       ])
        wh.write(cmd)
        wh.write('\n')
        
        return wh



class cnv:
    ''' Specific class to process CNV calling '''

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def getPreprocessIntervals(self,configDict,outDir,wh):
        ''' Subroutine to generate preprocess intervals '''
        
        ppDict = configDict['sv']['PreprocessIntervals']
        refFile = '/'.join([configDict['general']['resourceDir'],
                            configDict['resources']['reference']['grch38']
                           ])
        intervalFile = '/'.join([configDict['general']['resourceDir'],
                                 ppDict['intervals']
                                ])

        outIntervalFile = '/'.join([outDir,'gdd.pp.interval_list'])
        ppDict['reference'] = refFile
        ppDict['output'] = outIntervalFile
        ppDict['intervals'] = intervalFile

        cmdList = ['gatk PreprocessIntervals']

        for key,value in ppDict.items():
            cmdList.append('--'+key+' '+str(value))
         
        wh.write("\n######### Preprocess Interval #########\n")
        wh.write("echo 'Preprocess Intervals'\n\n")
        wh.write(' \\\n'.join(cmdList))
        wh.write('\n')

        objC = cluster()
        wh = objC.checkErrorFile(outIntervalFile,wh,'error',None)
 
        return outIntervalFile,wh

    def getGCannotatedIntervals(self,configDict,ppIntervalList,ppiDir,wh):
        ''' Subroutine to generate preprocess intervals '''
        
        ppDict = configDict['sv']['AnnotateIntervals']
        refFile = '/'.join([configDict['general']['resourceDir'],
                            configDict['resources']['reference']['grch38']
                           ])
        intervalFile = '/'.join([configDict['general']['resourceDir'],
                                 ppDict['intervals']
                                ])

        outIntervalFile = '/'.join([ppiDir,'gdd.pp.gc.annotated.tsv'])
        ppDict['reference'] = refFile
        ppDict['output'] = outIntervalFile
        ppDict['intervals'] = ppIntervalList

        cmdList = ['gatk AnnotateIntervals']

        for key,value in ppDict.items():
            cmdList.append('--'+key+' '+str(value))
         
        wh.write("\n######### GC Annotated Intervals #########\n")
        wh.write("echo 'GC Annotated Intervals'\n\n")
        wh.write(' \\\n'.join(cmdList))
        wh.write('\n')

        objC = cluster()
        wh = objC.checkErrorFile(outIntervalFile,wh,'error',None)
 
        return outIntervalFile,wh

    def getSampleReadCounts(self,configDict,outIntervalList,bamFile,patId,outDir,wh):
        ''' Subroutine to get the read counts (RC) per bin '''

        rcDict = configDict['sv']['CollectReadCounts']
        refFile = '/'.join([configDict['general']['resourceDir'],
                            configDict['resources']['reference']['grch38']
                           ])

        outRCFile = '/'.join([outDir,patId+'.tsv'])
        
        rcDict['reference'] = refFile
        rcDict['intervals'] = outIntervalList
        rcDict['input'] = bamFile
        rcDict['output'] = outRCFile

        cmdList = ['gatk CollectReadCounts']
        for key,value in rcDict.items():
            cmdList.append('--'+key+' '+value)

        wh.write("\n######### Sample Read Count #########\n")
        wh.write("echo 'Sample Read count for sample:"+patId+"'\n\n")
        wh.write(' \\\n'.join(cmdList))
        wh.write("\n")
                
        objC = cluster()
        wh = objC.checkErrorFile(outRCFile,wh,'error',None)
        
        return outRCFile,wh


    def getFilterIntervals(self,configDict,ppIntervalList,ppGCannotInterval,
                           sampleRCList,outDir,wh):
        ''' Subroutine to Filter Intervals on GC-contenet and cohort extreme
        counts '''

        fiDict = configDict['sv']['FilterIntervals']
        fiDict['intervals'] = ppIntervalList
        fiDict['annotated-intervals'] = ppGCannotInterval
        outGCfilterFile = '/'.join([outDir,'gdd.pp.gc.filtered.interval_list'])
        
        fiDict['output'] = outGCfilterFile

        #gatk_sif = '/nfs/research/dunham/resources/softwares/gatk_docker/gatk_4.2.5.0.sif'
        gatk_sif = configDict['general']['gatk_docker']
        #cmdList = ['gatk FilterIntervals']
        cmdList = ['singularity exec']+[gatk_sif]+[' gatk FilterIntervals']
        for key,value in fiDict.items():
            if key =='input':
                for s_rc in sampleRCList:
                    cmdList.append('--'+key+' '+s_rc)
            else:
                cmdList.append('--'+key+' '+value)
        
        wh.write("\n######### Filter Intervals #########\n")
        wh.write("echo 'Filter Intervals'\n\n")
        wh.write(' \\\n'.join(cmdList))
        wh.write('\n\n')

        objC = cluster()
        wh = objC.checkErrorFile(outGCfilterFile,wh,'error',None)

        return outGCfilterFile,wh
        
    def getDetermineGermlineContigPloidy(self,configDict,filterIntervalList,
                                         sampleRCList,outDir,wh):
        ''' Subroutine to determine germline contig ploidy '''

        plDict = configDict['sv']['DetermineGermlineContigPloidy']
        plDict['intervals'] = filterIntervalList
        #outPloidyFile = '/'.join([outDir,'gdd.pp.gc.flt.ploidy.tsv'])
        
        plDict['output'] = outDir
        plDict['contig-ploidy-priors']='/'.join([configDict['general']['resourceDir'],
                                                   plDict['contig-ploidy-priors']
                                                  ])
        plDict['output-prefix'] = 'gdd.pp.gc.flt.ploidy'                                    
        
        #gatk_sif = '/nfs/research/dunham/resources/softwares/gatk_docker/gatk_4.2.5.0.sif'
        gatk_sif = configDict['general']['gatk_docker']
        cmdList = ['singularity exec ']+[gatk_sif]+[' gatk DetermineGermlineContigPloidy']
        
        for key,value in plDict.items():
            if key =='input':
                for s_rc in sampleRCList:
                    cmdList.append('--'+key+' '+s_rc)
            else:
                cmdList.append('--'+key+' '+value)
        
        wh.write("\n######### Determine Germline Contig Ploidy #########\n")
        wh.write("echo 'Determine Germline Contig Ploidy'\n\n")
        wh.write(' \\\n'.join(cmdList))
        wh.write('\n\n')

        objC = cluster()
        wh = objC.checkErrorDir(plDict['output'],wh,'error',None)

        return plDict['output'],wh

    def getIntervalListTools(self,configDict,filterIntervalList,sctDir,wh):
        ''' Subroutine to make interval list for scattering '''

        sctDict = configDict['sv']['IntervalListTools']
        sctDict['INPUT'] = filterIntervalList
        sctDict['OUTPUT'] = sctDir
        
        gatk_sif = configDict['general']['gatk_docker']
        cmdList = ['singularity exec ']+[gatk_sif]+[' gatk IntervalListTools']
        
        for key,value in sctDict.items():
            cmdList.append('--'+key+' '+value)

        wh.write("\n######### Scatter Intervals #########\n")
        wh.write("echo 'Scatter Interval List'\n\n")
        wh.write(' \\\n'.join(cmdList))
        wh.write('\n\n')

        objC = cluster()
        wh = objC.checkErrorDir(sctDict['OUTPUT'],wh,'error',None)

        return wh

    def getGermlineCNVCaller(self,configDict,outDir,cnvDir,sampleList,
                             cohort_mode,wh):
        ''' Subroutine to run GermlineCNVCaller tool in COHORT mode '''

        gcDict = configDict['sv']['GermlineCNVCaller']
        gcDict['intervals'] = '/'.join([outDir,'gdd.pp.gc.filtered.interval_list'])
        #gcDict['intervals'] = sctInterval
        gcDict['contig-ploidy-calls'] = '/'.join([outDir,
                                                  'ploidy/gdd.pp.gc.flt.ploidy-calls'
                                                 ])
        gcDict['annotated-intervals'] = '/'.join([outDir,
                                                  'gdd.pp.gc.annotated.tsv'
                                                 ])
        gcDict['output'] = cnvDir
        gcDict['tmp-dir'] = cnvDir
        #gcDict['output'] = 'cnv_call'
        gcDict['output-prefix'] = 'cnv_call'
        gcDict['run-mode'] = cohort_mode
        
        #sampleRCList = os.listdir('/'.join([outDir,'SRC']))
        sampleRCList = []
        for ele in sampleList:
            #sampleRCList.append(ele+'.tsv')
            sampleRCList.append(ele)
        
        gatk_sif = configDict['general']['gatk_docker']
        cmdList = ['singularity exec ']+[gatk_sif]+[' gatk GermlineCNVCaller']

        for key,value in gcDict.items():
            if key =='input':
                for s_rc in sampleRCList:
                    #cmdList.append('--'+key+' '+'/'.join([outDir,'SRC',s_rc]))
                    cmdList.append('--'+key+' '+s_rc)
            else:
                cmdList.append('--'+key+' '+value)
         
        wh.write("\n######### Germline CNV calling #########\n")
        #wh.write("echo 'Germline CNV Calling: "+sct+"'\n\n")
        wh.write("echo 'Germline CNV Calling:' \n\n")
        wh.write(' \\\n'.join(cmdList))
        wh.write('\n')

        objC = cluster()
        wh = objC.checkErrorFile(gcDict['output'],wh,'error',None)
 
        return wh

    def getPostprocessGermlineCNVCalls(self,configDict,outDir,vcfDir,
                                       sample_index,ddd_id,wh):
        ''' Subroutine to run the PostprocessGermlineCNVCalls.
            Generates the CNV vcfs
        '''
        refFile = '/'.join([configDict['general']['resourceDir'],
                            configDict['resources']['reference']['grch38_dict']
                           ])
        ppDict = configDict['sv']['PostprocessGermlineCNVCalls']
        ppDict['model-shard-path'] = '/'.join([outDir,'cnv_call/cnv_call-model'])
        ppDict['calls-shard-path'] = '/'.join([outDir,'cnv_call/cnv_call-calls'])
        ppDict['contig-ploidy-calls'] = '/'.join([outDir,
                                                  'ploidy/gdd.pp.gc.flt.ploidy-calls'
                                                 ])

        ppDict['sample-index'] = str(sample_index)
        ppDict['sequence-dictionary'] = refFile
        ppDict['output-genotyped-intervals'] = '/'.join([vcfDir,
                                                         ddd_id+'.gi.vcf.gz'
                                                        ])
        ppDict['output-genotyped-segments'] = '/'.join([vcfDir,
                                                         ddd_id+'.gs.vcf.gz'
                                                        ])
        ppDict['output-denoised-copy-ratios'] = '/'.join([vcfDir,
                                                         ddd_id+'.denoised.copy-ratio.tsv'
                                                        ])


        gatk_sif = configDict['general']['gatk_docker']
        cmdList = ['singularity exec ']+[gatk_sif]+[' gatk PostprocessGermlineCNVCalls']

        for key,value in ppDict.items():
            if key=='allosomal-contig':
                for ele in ppDict[key]['value']:
                    cmdList.append('--'+key+' '+ele)
            else:
                cmdList.append('--'+key+' '+str(value))

        wh.write("\n######### Postprocess Germline CNV calls #########\n")
        wh.write("echo 'Postprocess Germline CNV Calls'\n\n")
        wh.write(' \\\n'.join(cmdList))
        wh.write('\n')

        objC = cluster()
        wh = objC.checkErrorFile(ppDict['output-genotyped-segments'],wh,'error',None)
 
        return wh

    def getMantaSVcalls(self,configDict,svDir,bamList,famID,wh):
                                       
        ''' Subroutine to run the MantaSV 
        '''
        refFile = '/'.join([configDict['general']['resourceDir'],
                            configDict['resources']['reference']['grch38']
                           ])
        msvDict = configDict['sv']['mantaSV']
        msvDict['runDir'] = svDir

        msvDict['referenceFasta'] = refFile

        cmdList = ['configManta.py']

        for key,value in msvDict.items():
            if key=='bam':
                for ele in bamList:
                    if re.search('bam$',ele):
                        cmdList.append('--'+key+' '+ele)
            else:
                cmdList.append('--'+key+' '+str(value))

        cmdList.append('--exome')

        wh.write("\n######### Manta SV calls #########\n")
        wh.write("echo 'Manta SV Calls'\n\n")
        wh.write(' \\\n'.join(cmdList))
        wh.write('\n')

        outPy = '/'.join([msvDict['runDir'],'runWorkflow.py'])
        objC = cluster()
        wh = objC.checkErrorFile(outPy,wh,'error',None)

        wh.write('cd '+msvDict['runDir']+'\n')
        wh.write('./runWorkflow.py')
        wh.write('\n')
       
        ''' 
        wh.write("var=$(cat workflow.exitcode.txt)")
        wh.write('\n')

        wh.write("\nif [\"$var\" = 1]; then \n")
        wh.write("   echo -e \"ERROR: program did not finish \" \n")
        wh.write("   exit 1\n")
        ''' 
        
        return wh

    def getCnestCNVcalls(self,configDict,workDir,projName,sampleName,
                         samplePath,step,wh):
                                       
        ''' Subroutine to run the CNest CNV calling steps 
        '''
        refFile = '/'.join([configDict['general']['resourceDir'],
                            configDict['resources']['reference']['grch38']
                           ])

        
        msvDict = configDict['sv']['cnest']

        if step=='step1':
            cnDict = configDict['sv']['cnest'][step]
            bedFile = '/'.join([configDict['general']['resourceDir'],
                              cnDict['bed']
                           ])
            cnDict['bed'] = bedFile
            cnDict['project'] = projName
        elif step =='step2':
            cnDict = configDict['sv']['cnest'][step]
            cnDict['project'] = projName
            cnDict['sample'] = sampleName
            cnDict['input'] = samplePath
        elif step == 'step3':
            cnDict = {}
            cnDict = configDict['sv']['cnest'][step]
            cnDict['indextab'] = '/'.join([workDir,projName,'index_tab.txt'])
            cnDict['qc'] = '/'.join([workDir,projName,'qc_file.txt'])
            cnDict['gender'] = '/'.join([workDir,projName,'gender_file.txt'])
            cnDict['cov'] = '/'.join([workDir,projName,'coverage_file.txt'])
            cnDict['bindir'] = '/'.join([workDir,projName,'bin'])
            
        elif step =='step4':
            cnDict = {}
            cnDict = configDict['sv']['cnest'][step]
            cnDict['indextab'] = '/'.join([workDir,projName,'index_tab.txt'])
            cnDict['gender'] = '/'.join([workDir,projName,'gender_file.txt'])
            cnDict['bindir'] = '/'.join([workDir,projName,'bin'])
            cnDict['rbindir'] = '/'.join([workDir,projName,'rbin'])
            cnDict['cordir'] = '/'.join([workDir,projName,'cor'])
            cnDict['logrdir'] = '/'.join([workDir,projName,'logr'])
            cnDict['sample'] = sampleName
            cnDict['batch'] = str(1000)
        elif step == 'step5':
            cnDict={}
            cnDict = configDict['sv']['cnest'][step]
            cnDict['indextab'] = '/'.join([workDir,projName,'index_tab.txt'])
            cnDict['gender'] = '/'.join([workDir,projName,'gender_file.txt'])
            cnDict['rbindir'] = '/'.join([workDir,projName,'rbin'])
            cnDict['cnvdir'] = '/'.join([workDir,projName,'cnv'])
            cnDict['cordir'] = '/'.join([workDir,projName,'cor'])
            cnDict['cov'] = '/'.join([workDir,projName,'coverage_file.txt'])
            cnDict['sample'] = sampleName
            cnDict['batch'] = str(300)
 
        #cmdList = [' '.join(['singularity run',msvDict['docker'],step])]
        cmdList = [' '.join(['singularity run',msvDict['docker'],step])]

        for key,value in cnDict.items():
            cmdList.append('  --'+key+' '+value)

        wh.write("\n######### CNEST calls "+step+" #########\n")
        wh.write(' echo "Starting step: '+step+'"\n\n')
        wh.write('\\\n'.join(cmdList))
        wh.write('\n')

        
        objC = cluster()
        if step=='step1':
            index_tab = '/'.join([workDir,projName,'index_tab.txt'])
        if step=="step2":
            index_tab = '/'.join([workDir,projName,'bin',sampleName])

        #wh = objC.checkErrorFile(index_tab,wh,'error',None)
        
        return wh


