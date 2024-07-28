#!/usr/bin/python

import re,sys,argparse,os
from xml.etree import ElementTree
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
#import xmltodict
from collections import OrderedDict
import shutil

class util:

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    
    def display(self):
        print('Inside UTIL class. Contains Utilities related methods')

    def str2bool(self,v):
        ''' Subroutine to convert string to boolean '''
        return v.lower() in ('yes','true','t','1')

    def procDatasetArgs(self):
        ''' Subroutine to process command line arguments for fetching
        DDD dataset
        '''

        parser = argparse.ArgumentParser(
                 description='Script to parse the DDD log file and create '+\
                              'SLURM specific commands to download CRAM files'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-i','--inpFile',help='Input File (logFile)',
                            action='store',dest='inpFile',required=True)
        parser.add_argument('-o','--outFile',help='Out File',
                            action='store',dest='outFile',required=True) 
        parser.add_argument('-s','--saveTo',help='Download save directory',
                            action='store',dest='outDir',required=True)   
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   


        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        work_dir = os.path.abspath(results.workDir)
        inp_file = os.path.abspath(results.inpFile)
        out_file = os.path.abspath(results.outFile)
        out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate

        cmdDict = {'config':config_file,'inpFile':inp_file,'outFile':out_file,
                   'outDir':out_dir,'shellType':shell_type,'workDir':work_dir,
                   'projDate':project_date
                  }

        return cmdDict

    def procManifestCmdArgs(self):
        ''' Subroutine to process command line arguments related to creation of 
        Manifest file 
        '''
        parser = argparse.ArgumentParser(
                 description='Script to lauch GTAK SNV callig pipelie '+\
                              'SLURM/LSF specific commands for GATK'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-s','--sampleDir',help='Sample Directory',
                            action='store',dest='sampleDir',required=True)
        parser.add_argument('-i','--inpFile',help='Input DDD List',
                            action='store',dest='inpFile',required=True)
        parser.add_argument('-o','--outFile',help='Out Manifest File',
                            action='store',dest='outFile',required=True)


        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        sample_dir = os.path.abspath(results.sampleDir)
        inp_file = os.path.abspath(results.inpFile)
        out_file = os.path.abspath(results.outFile)

        cmdDict ={'config':config_file,'sampleDir':sample_dir,
                  'inpDDD':inp_file,'outManifest':out_file
                 }

        return cmdDict
 

    def procDecryptCmdArgs(self):
        ''' Subroutine to process command line arguments related to decryption
        of aspera downloaded files 
        '''
        parser = argparse.ArgumentParser(
                 description='Script to lauch Decryption of Aspera download files '
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-s','--sampleDir',help='Sample Directory',
                            action='store',dest='sampleDir',required=True)
        parser.add_argument('-m','--mapFile',help='DDD Map file',
                            action='store',dest='mapFile',required=True)
        parser.add_argument('-o','--outdir',help='Out Decryption Directory',
                            action='store',dest='outDir',required=True)
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)
        parser.add_argument('-l','--launch',help='Launch Flag',
                            action='store',dest='launchFlag',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   


        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        sample_dir = os.path.abspath(results.sampleDir)
        map_file = os.path.abspath(results.mapFile)
        out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)


        cmdDict ={'config':config_file,'sampleDir':sample_dir,
                  'map':map_file,'outDir':out_dir,'shellType':shell_type,
                  'projDate':project_date,'shellType':shell_type,
                  'launchFlag':launch
                 }

        return cmdDict

    def procGATKArgs(self):
        ''' Subroutine to process command line arguments related to GATK
        SNV calling pipeline
        '''

        parser = argparse.ArgumentParser(
                 description='Script to lauch GTAK SNV callig pipelie '+\
                              'SLURM/LSF specific commands for GATK'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-s','--scriptDir',help='Script Directory',
                            action='store',dest='scriptDir',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-i','--inpFile',help='Input File (logFile)',
                            action='store',dest='inpFile',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   


        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        script_dir = os.path.abspath(results.scriptDir)
        work_dir = os.path.abspath(results.workDir)
        inp_file = os.path.abspath(results.inpFile)
        #out_file = os.path.abspath(results.outFile)
        #out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)

        cmdDict = {'config':config_file,'inpFile':inp_file,
                   'shellType':shell_type,'workDir':work_dir,
                   'scriptDir':script_dir,'projDate':project_date,
                   'launchFlag':launch
                  }

        return cmdDict

    def procGatkCnvArgs(self):
        ''' Subroutine to process command line arguments related to GATK
        CNV calling pipeline
        '''

        parser = argparse.ArgumentParser(
                 description='Script to lauch GATK CNV callig pipelie '+\
                              'SLURM/LSF specific commands for GATK'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-s','--scriptDir',help='Script Directory',
                            action='store',dest='scriptDir',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-i','--inpFile',help='Input File (logFile)',
                            action='store',dest='inpFile',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   


        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        script_dir = os.path.abspath(results.scriptDir)
        work_dir = os.path.abspath(results.workDir)
        inp_file = os.path.abspath(results.inpFile)
        #out_file = os.path.abspath(results.outFile)
        #out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)

        cmdDict = {'config':config_file,'inpFile':inp_file,
                   'shellType':shell_type,'workDir':work_dir,
                   'scriptDir':script_dir,'projDate':project_date,
                   'launchFlag':launch
                  }

        return cmdDict

    def procGatkCnvVcfArgs(self):
        ''' Subroutine to process command line arguments related to GATK
        CNV calling pipeline
        '''

        parser = argparse.ArgumentParser(
                 description='Script to lauch GATK CNV callig pipelie '+\
                              'SLURM/LSF specific commands for GATK'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-s','--scriptDir',help='Script Directory',
                            action='store',dest='scriptDir',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-i','--inpFile',help='Input File (logFile)',
                            action='store',dest='inpFile',required=True),
        parser.add_argument('-f','--shardPath',help='Input Shard Path',
                            action='store',dest='shardPath',required=True),
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   

        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        script_dir = os.path.abspath(results.scriptDir)
        work_dir = os.path.abspath(results.workDir)
        inp_file = os.path.abspath(results.inpFile)
        shard_path = os.path.abspath(results.shardPath)
        #out_file = os.path.abspath(results.outFile)
        #out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)

        cmdDict = {'config':config_file,'inpFile':inp_file,
                   'shellType':shell_type,'workDir':work_dir,
                   'scriptDir':script_dir,'projDate':project_date,
                   'shardPath':shard_path,
                   'launchFlag':launch
                  }

        return cmdDict

    def procMantaSVArgs(self):
        ''' Subroutine to process command line arguments related to Manta'''
        
        parser = argparse.ArgumentParser(
                 description='Script to lauch Manta SV callig pipelie '+\
                              'SLURM/LSF specific commands for Manta'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-s','--scriptDir',help='Script Directory',
                            action='store',dest='scriptDir',required=True)
        parser.add_argument('-w','--workDir',help='Work Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-m','--manifest',help='Manifest File',
                            action='store',dest='manifest',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   

        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        script_dir = os.path.abspath(results.scriptDir)
        work_dir = os.path.abspath(results.workDir)
        inp_file = os.path.abspath(results.manifest)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)

        cmdDict = {'config':config_file,'manifest':inp_file,
                   'shellType':shell_type,'workDir':work_dir,
                   'scriptDir':script_dir,'projDate':project_date,
                   'launchFlag':launch
                  }

        return cmdDict

    def procCnestSVArgs(self):
        ''' Subroutine to process command line arguments related to CNest'''
        
        parser = argparse.ArgumentParser(
                 description='Script to lauch Manta CNV callig pipelie '+\
                              'SLURM/LSF specific commands for CNest'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-s','--scriptDir',help='Script Directory',
                            action='store',dest='scriptDir',required=True)
        parser.add_argument('-w','--workDir',help='Work Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-m','--manifest',help='Manifest File',
                            action='store',dest='manifest',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   

        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        script_dir = os.path.abspath(results.scriptDir)
        work_dir = os.path.abspath(results.workDir)
        inp_file = os.path.abspath(results.manifest)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)

        cmdDict = {'config':config_file,'manifest':inp_file,
                   'shellType':shell_type,'workDir':work_dir,
                   'scriptDir':script_dir,'projDate':project_date,
                   'launchFlag':launch
                  }

        return cmdDict

    
    
    def procJointGTArgs(self):
        ''' Subroutine to process command line arguments related to GATK
        SNV calling pipeline - Joint Genotype Step
        '''

        parser = argparse.ArgumentParser(
                 description='Script to lauch GATK-SNV-JointGT step '+\
                              'SLURM/LSF specific commands for GATK'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-s','--scriptDir',help='Script Directory',
                            action='store',dest='scriptDir',required=True)
        parser.add_argument('-t','--intervalDir',help='Interval directory Exome',
                            action='store',dest='intervalDir',required=True)
        parser.add_argument('-i','--manifest',help='Input Directory for GVCF files',
                            action='store',dest='manifestFile',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   
        parser.add_argument('-e','--expr',help='CombineGVCF/GenomicsDBImport', 
                            action='store',dest='analType',required=True)   

        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        work_dir = os.path.abspath(results.workDir)
        script_dir = os.path.abspath(results.scriptDir)
        interval_dir = os.path.abspath(results.intervalDir)
        manifest_file = os.path.abspath(results.manifestFile)
        #out_file = os.path.abspath(results.outFile)
        #out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)
        anal_type = results.analType

        cmdDict = {'config':config_file,'manifest':manifest_file,
                   'shellType':shell_type,'workDir':work_dir,
                   'projDate':project_date,'launchFlag':launch,
                   'expr':anal_type,'intervalDir':interval_dir,
                   'scriptDir':script_dir
                  }

        return cmdDict

    def procAnnotSNVargs(self):
        ''' Subroutine to process command line arguments related to Annotating
        SNVs
        '''

        parser = argparse.ArgumentParser(
                 description='Script to lauch Annotation SNV pipelie '+\
                              'SLURM/LSF specific commands for AnnotSNV'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-s','--scriptDir',help='Script Directory',
                            action='store',dest='scriptDir',required=True)
        parser.add_argument('-i','--jgtFile',help='Input jointGT Merged File',
                            action='store',dest='jgtFile',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   

        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        work_dir = os.path.abspath(results.workDir)
        scr_dir = os.path.abspath(results.scriptDir)
        jgt_file = os.path.abspath(results.jgtFile)
        #out_file = os.path.abspath(results.outFile)
        #out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)

        cmdDict = {'config':config_file,'jgtFile':jgt_file,
                   'shellType':shell_type,'workDir':work_dir,
                   'scriptDir':scr_dir,'projDate':project_date,
                   'launchFlag':launch
                  }

        return cmdDict

    def procCombineGvcfArgs(self):
        ''' Subroutine to process command line arguments related to GATK
        SNV calling pipeline - CombineGVCFs step.
        '''

        parser = argparse.ArgumentParser(
                 description='Script to lauch GATK-SNV-CombineGVCFs step '+\
                              'SLURM/LSF specific commands for GATK'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-s','--scriptDir',help='Script Directory',
                            action='store',dest='scriptDir',required=True)
        parser.add_argument('-i','--jgtDir',help='Input Directory for jointGT',
                            action='store',dest='jgtDir',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   
        parser.add_argument('-e','--expr',help='CombineGVCF/GenomicsDBImport', 
                            action='store',dest='analType',required=True)   

        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        work_dir = os.path.abspath(results.workDir)
        script_dir = os.path.abspath(results.scriptDir)
        jgt_dir = os.path.abspath(results.jgtDir)
        #out_file = os.path.abspath(results.outFile)
        #out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)
        anal_type = results.analType

        cmdDict = {'config':config_file,
                   'shellType':shell_type,'workDir':work_dir,
                   'projDate':project_date,'launchFlag':launch,
                   'expr':anal_type,'jgtDir':jgt_dir,
                   'scriptDir':script_dir
                  }

        return cmdDict

    def procCombDVArgs(self):
        ''' Subroutine to process command line arguments related to GATK
        SNV calling pipeline - CombineGVCFs step.
        '''

        parser = argparse.ArgumentParser(
                 description='Script to lauch DV-SNV-CombVcfs step '+\
                              'SLURM/LSF specific commands for GATK'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-s','--scriptDir',help='Script Directory',
                            action='store',dest='scriptDir',required=True)
        parser.add_argument('-i','--dvDir',help='Input Directory for DV dir',
                            action='store',dest='dvDir',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   

        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        work_dir = os.path.abspath(results.workDir)
        script_dir = os.path.abspath(results.scriptDir)
        dv_dir = os.path.abspath(results.dvDir)
        #out_file = os.path.abspath(results.outFile)
        #out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)

        cmdDict = {'shellType':shell_type,'workDir':work_dir,
                   'projDate':project_date,'launchFlag':launch,
                   'dvDir':dv_dir,'scriptDir':script_dir,
                   'config':config_file
                  }

        return cmdDict

    def procVQSRArgs(self):
        ''' Subroutine to process command line arguments related to GATK
        SNV calling pipeline - VQSR Step
        '''

        parser = argparse.ArgumentParser(
                 description='Script to lauch GATK-SNV-VQSR step '+\
                              'SLURM/LSF specific commands for GATK'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-i','--jgtDir',help='Input jointGT Dir',
                            action='store',dest='jgtDir',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   


        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        work_dir = os.path.abspath(results.workDir)
        jgt_dir = os.path.abspath(results.jgtDir)
        #out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)

        cmdDict = {'config':config_file,'jgtDir':jgt_dir,
                   'shellType':shell_type,'workDir':work_dir,
                   'projDate':project_date,'launchFlag':launch,
                  }

        return cmdDict

    def procGMFargs(self):
        ''' Subroutine to automate launch of gmfAnalysis code for each of the 3
        case scenarios.
        '''

        parser = argparse.ArgumentParser(
                 description='Script to launch GMF-Analysis - 3 cases')
        parser.add_argument('-a','--config',help='Config XML File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-i','--geneDir',help='filtered set of variants',
                            action='store',dest='geneDir',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-s','--scriptDir',help='Script Directory',
                            action='store',dest='scriptDir',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   


        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        gene_file = os.path.abspath(results.geneDir)
        work_dir = os.path.abspath(results.workDir)
        script_dir = os.path.abspath(results.scriptDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)

        cmdDict = {'geneDir':gene_file,'shellType':shell_type,
                   'workDir':work_dir,'scriptDir':script_dir,
                   'projDate':project_date,'launchFlag':launch,
                   'config':config_file
                  }

        return cmdDict

    def procCopyArgs(self):
        ''' Subroutine to process command line arguments for copying
        re-aligned BAM/VCF files backup drive
        '''

        parser = argparse.ArgumentParser(
                 description='Script to lauch GATK-SNV-JointGT step '+\
                              'SLURM/LSF specific commands for GATK'  
                )
        parser.add_argument('-a','--configFile',help='Analysis Config File',
                            action='store',dest='configFile',required=True)
        parser.add_argument('-w','--workDir',help='Working Directory',
                            action='store',dest='workDir',required=True)
        parser.add_argument('-o','--outDir',help='Out backup directory path',
                            action='store',dest='outDir',required=True)
        parser.add_argument('-i','--manifest',help='Input Directory for GVCF files',
                            action='store',dest='manifestFile',required=True)
        parser.add_argument('-c','--clusterType',help='Cluster type script: Slurm/Lsf',
                            action='store',dest='shellType',required=True)   
        parser.add_argument('-p','--projDate',help='Project Date',
                            action='store',dest='projDate',required=True)   
        parser.add_argument('-l','--launch',help='LSF/SLURM launch flag',
                            action='store',dest='launchFlag',required=True)   

        results = parser.parse_args()
        config_file = os.path.abspath(results.configFile)
        work_dir = os.path.abspath(results.workDir)
        out_dir = os.path.abspath(results.outDir)
        manifest_file = os.path.abspath(results.manifestFile)
        #out_file = os.path.abspath(results.outFile)
        #out_dir = os.path.abspath(results.outDir)
        shell_type = results.shellType
        project_date = results.projDate
        launch = self.str2bool(results.launchFlag)

        cmdDict = {'config':config_file,'manifest':manifest_file,
                   'shellType':shell_type,'workDir':work_dir,
                   'projDate':project_date,'launchFlag':launch,
                   'outDir':out_dir
                  }

        return cmdDict

    def getConfigDict(self,config_file,module_name=[]):
        
        ''' Process input configuration XML file and return as Dictionary 
            object  
        '''

        import xmltodict

        with open(config_file) as cf:
            doc = xmltodict.parse(cf.read())

        doc = doc['config']

        if module_name:
            configDict = doc[module_name]
        else:
            configDict = doc

        return configDict

    
    def processInit(self,work_dir,cmd_args,proj_date=[]):

        ''' Suborutine to create temporary (tmp) folders and directories '''
        
        print("\n\nCreating temporary directories inside the work directory: "+\
              proj_date+" \n")

        date_str = proj_date
        
        tmp_dir = os.path.abspath('/'.join([work_dir,date_str]))
        tmp_bin = os.path.abspath('/'.join([tmp_dir,'tmp_binaries']))
        tmp_data = os.path.abspath('/'.join([tmp_dir,'tmp_data']))
        tmp_status = os.path.abspath('/'.join([tmp_dir,'tmp_status']))
        tmp_log = os.path.abspath('/'.join([tmp_dir,'tmp_log']))
        sb_log = os.path.abspath('/'.join([tmp_dir,'sb_log']))

        tmp_dict = {'tmpDir':tmp_dir,'tmpBin':tmp_bin,
                    'tmpData':tmp_data,'tmpStat':tmp_status,
                    'tmpLog':tmp_log,'sbLog':sb_log
                   }

        if os.path.isdir(work_dir+'/'+date_str):
            print(''' -- The entered project space already exists.\ 
                       Do you want to delete it?''')
            
            if re.search('manta',proj_date):
                text = raw_input('Enter y/n: ')
            else:
                text = input('Enter y/n: ')
            print(text)
            if text.lower() == 'y':
                print(''' -- Deleting existing project. Creating fresh workspace 
                      with name: ''',date_str)
                shutil.rmtree('/'.join([work_dir,date_str]))
                
                for keys,values in tmp_dict.items():
                    os.makedirs(values)

            elif text.lower() =='n':
                pass
        else:
            for keys,values in tmp_dict.items():
                os.makedirs(values)

        wh = open('/'.join([tmp_dir,'Command_Entered.txt']),'a')
        wh.write('python '+' '.join(cmd_args)+'\n\n')
        wh.close()

        return tmp_dict


