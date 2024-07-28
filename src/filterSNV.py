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
    cmdDict = objU.procAnnotSNVargs()
    configFile = cmdDict['config'];workDir = cmdDict['workDir']
    scriptDir = cmdDict['scriptDir']
    jgtRecalDir = cmdDict['jgtFile'];sh_type=cmdDict['shellType']
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
   
    if not re.search('ukbb',projDate):
        chr_list = ['chr'+str(x) for x in range(1,23)]+['chrX','chrY']
    else:
        chr_list = [str(x) for x in range(1,23)]+['X','Y']

    #print(chr_list)
    maf_cutoff = str(0.01)
   
    if not re.search('ukbb',projDate):
        memList = [1200,1100,1000,800,600,600,600,600,600,600,
                   800,800,1200,600,600,600,800,800,900,600,600,
                   500,500,500]
    elif re.search('ukbb',projDate):
        memList = [1200,1100,1000,1000,1000,1000,800,800,800,800,
                   800,800,1200,700,700,700,800,800,900,700,700,
                   500,900,500]
 

    ms_type = re.split('_',projDate)[1]

    memDict = OrderedDict()
    for chrNum,memNum in zip(chr_list,memList):
        memDict[chrNum] =  str(memNum)+'G'

    if re.search('ddd|ukb50k',projDate):
        for chrNum in chr_list:

            outChrDir = '/'.join([tmp_data,chrNum])
            
            fileId = 'filterMAF_'+ms_type+'_'+chrNum

            configDict['lsf']['params1']['M'] = '6G'
            configDict['lsf']['params1']['q'] = 'research'

            sh_file,wh = objC.getClusterWriteHandle(tmp_bin,fileId,sh_type)
            wh = objC.writeClusterTop(configDict,'NDD-WES-GATK-UKBB-filterSNV',wh)
            sh_out,sh_err,wh = objC.writeClusterInit(configDict,sb_log,
                                                     sh_type,fileId,wh
                                                    )
           
            wh.write('\nmodule load r-4.0.3-gcc-9.3.0-xiarbub\n')
            wh.write('module load python-2.7.18-gcc-9.3.0-phwlzad\n\n')
            
            wh.write('echo "$HOSTNAME"\n')
            wh.write('unset BCFTOOLS_PLUGINS\n\n')
            wh.write('mkdir -p '+outChrDir+'\n\n')
          
            #Annotated BCF file
            if not re.search('ukb50k',projDate):
                chr_anno_bcf = '/'.join([jgtRecalDir,chrNum,
                                     'tmpFile_'+chrNum+'.AC0.norm.vep.anno.sites.gdd.freq.bcf'
                                    ])
            else:
                chr_anno_bcf = '/'.join([jgtRecalDir,chrNum,
                                     'tmpFile_'+chrNum+'.AC0.norm.vep.anno.sites.ukbb.freq.bcf'
                                    ])
            
            """
            maf_filter_bcf = '/'.join([outChrDir,'tmpFile_'+chrNum+\
                                       '.AC0.norm.vep.anno.sites.gdd.freq.exon.0_01.bcf'
                                      ])
            """
            # GDD and UKBB specific analysis
            if not re.search('ukbb',projDate):
                #mild_severe_list = ['mild','severe']
                
                #Extract Mild and Severe Samples
                #for ms_type in mild_severe_list:
                ms_type = re.split('_',projDate)[1]
                 
                ms_anno_bcf,wh = objA.writeMildSevereSamples(configDict,chr_anno_bcf,
                                                            outChrDir,ms_type,wh
                                                           )
                
                #Extract the exonic region (Ensembl)
                chr_exon_bcf,wh = objA.writeExtractExonicRegion(configDict,ms_anno_bcf,
                                                                outChrDir,wh
                                                               )
                
                #Frequency filtering
                maf_cutoff = str(0.01)
                maf_ms_bcf,wh = objA.writeFilterMAF(configDict,chr_exon_bcf,
                                                        maf_cutoff,wh
                                                       )
                
                #Impact filtering
                imp_filter_bcf,wh = objA.writeFilterImpact(configDict,maf_ms_bcf,wh)
                """
                vcf_2_tab = '/'.join([outChrDir,'tmpFile_'+chrNum+\
                                           '.AC0.norm.vep.anno.sites.gdd.freq.'+\
                                           ms_type+'.exon.0_01.imp.tab'
                                      ])
                
                #VCf-2-TAB creation
                """
                vcf_2_tab,wh = objA.writeConvertVcf2Tab(configDict,imp_filter_bcf,wh)
       
                #Transcript Prioritization
                """ 
                trs_prior_tab = objA.writeTranscriptPrioritization(configDict,vcf_2_tab,
                                                                      chrNum,projDate,memDict,
                                                                      ms_type
                                                                     )
                """
                objA.writeQueryGenesUKBB(configDict,vcf_2_tab,chrNum,projDate,memDict,ms_type)

            elif re.search('ukb50k',projDate):

                batchNum = re.split('\_',projDate)[1]
                
                #Extract samples Batch specific
                batch_ukb_bcf,wh = objA.writeExtractUKBB(configDict,chr_anno_bcf,
                                                          outChrDir,batchNum,wh
                                                      )
                #Extract the exonic region (Ensembl)
                chr_exon_bcf,wh = objA.writeExtractExonicRegion(configDict,batch_ukb_bcf,
                                                                outChrDir,wh
                                                               )

                #Frequency filtering
                maf_cutoff = str(0.01)
                maf_ukb_bcf,wh = objA.writeFilterMAF(configDict,chr_exon_bcf,
                                                      maf_cutoff,wh
                                                    )
                """
                maf_ukb_bcf = '/'.join([outChrDir,
                                        'tmpFile_'+chrNum+'.AC0.norm.vep.anno.sites.ukbb.freq.batch1.exon.0_01.bcf'
                                       ])
                """
                #Impact Filtering
                imp_ukb_bcf,wh = objA.writeFilterImpact(configDict,maf_ukb_bcf,wh)
           
                
                #VCf-2-TAB creation
                vcf_2_tab,wh = objA.writeConvertVcf2Tab(configDict,imp_ukb_bcf,wh)
               
                """
                vcf_2_tab = '/'.join([outChrDir,
                                        'tmpFile_'+chrNum+'.AC0.norm.vep.anno.sites.ukbb.freq.batch1.exon.0_01.imp.tab'
                                       ])
                """
                #Transcript Prioritization
                """
                trs_prior_tab = objA.writeTranscriptPrioritization(configDict,vcf_2_tab,
                                                                   chrNum,projDate,memDict
                                                                  )
                """
                objA.writeQueryGenesUKBB(configDict,vcf_2_tab,chrNum,projDate,memDict,batchNum)
                wh.close()
                
                if launch:
                    os.system('bsub <'+sh_file)

    elif re.search('ukb200k',projDate):

        f_batchNum = re.split('_',projDate)[2]

        # Get List of all the Freq Files from the annot directory
        chrDict = OrderedDict()
        inpDir = jgtRecalDir
        for chrNum in chr_list:
            batchDir = os.listdir('/'.join([inpDir,chrNum]))

            for batchNum in batchDir:
                chrNum_batch = chrNum+':'+batchNum
                #freqFile = '/'.join([inpDir,chrNum,batchNum,
                #                    'tmpFile_'+chrNum+\
                #                    '.AC0.norm.vep.anno.sites.ukbb.freq.bcf'
                #                   ])
                # Update: 29.10.2023; Added snippet to process annotated DV file
                freqFile = '/'.join([inpDir,chrNum,batchNum,
                                    'tmpFile_tags.AC0.norm.vep.anno.sites.merge.bcf'
                                   ])

                if chrNum_batch in chrDict:
                    tmp = chrDict[chrNum_batch]
                    tmp.append(freqFile)
                else:
                    chrDict[chrNum_batch] = [freqFile]
       
        
        # Iteratively do this for each batch
        for keys, values in chrDict.items():
            strs = re.split('\:',keys)
            chrNum = strs[0]
            batchNum = strs[1]
            
            #batch_ukb_bcf = values[0]
            chr_anno_bcf = values[0]
            outChrDir = '/'.join([workDir,'dv_fm_'+f_batchNum,chrNum,batchNum])

            fileId = 'fm_'+chrNum+'_'+batchNum

            configDict['lsf']['params1']['M'] = '8G'
            configDict['lsf']['params1']['q'] = 'research'
            configDict['lsf']['params1']['n'] = str(4)

            sh_file,wh = objC.getClusterWriteHandle(tmp_bin,fileId,sh_type)
            wh = objC.writeClusterTop(configDict,
                                      'WES-HPC-DV-UKB200K-FM-'+f_batchNum,wh)
            sh_out,sh_err,wh = objC.writeClusterInit(configDict,sb_log,
                                                     sh_type,fileId+'_'+f_batchNum,
                                                     wh
                                                    )
           
            #wh.write('\nmodule load r-4.2.2-gcc-11.2.0-oa3uudy\n')
            wh.write('module load python-2.7.18-gcc-9.3.0-phwlzad\n\n')
            
            #wh.write('echo "$HOSTNAME"\n')
            #wh.write('unset BCFTOOLS_PLUGINS\n\n')
            wh.write('mkdir -p '+outChrDir+'\n\n')
            
            #print(chr_anno_bcf)
            #Extract the exonic region (Ensembl)
            chr_exon_bcf,wh = objA.writeExtractExonicRegion(configDict,chr_anno_bcf,
                                                            outChrDir,wh
                                                           )
           

            #Frequency filtering
            maf_cutoff = str(0.01)
            maf_ukb_bcf,wh = objA.writeFilterMAF(configDict,chr_exon_bcf,
                                                  maf_cutoff,wh
                                                )

            wh.write('rm '+chr_exon_bcf+'*'+'\n\n')
            
            """
            maf_ukb_bcf = '/'.join([outChrDir,
                                    'tmpFile_'+chrNum+'.AC0.norm.vep.anno.sites.ukbb.freq.batch1.exon.0_01.bcf'
                                   ])
            """
            #Impact Filtering
            imp_ukb_bcf,wh = objA.writeFilterImpact(configDict,maf_ukb_bcf,wh)
            
            wh.write('rm '+maf_ukb_bcf+'*\n\n') 

            for s_batchNum in ['batch1','batch2','batch3']:

                #Extract samples Batch specific
                batch_ukb_bcf,wh = objA.writeExtractUKBB(configDict,imp_ukb_bcf,
                                                          outChrDir,s_batchNum,wh
                                                      )

                #VCf-2-TAB creation
                vcf_2_tab,wh = objA.writeConvertVcf2Tab(configDict,batch_ukb_bcf,wh)
                
                wh.write('rm '+batch_ukb_bcf+'*\n\n')

                vcf_2_tab_bgz,wh = objA.writeCompressBgz(configDict,vcf_2_tab,wh)
                wh.write('rm '+vcf_2_tab+'\n\n')
          
            wh.write('rm '+imp_ukb_bcf+"*"+'\n')
            
            """
            vcf_2_tab = '/'.join([outChrDir,
                                    'tmpFile_'+chrNum+'.AC0.norm.vep.anno.sites.ukbb.freq.batch1.exon.0_01.imp.tab'
                                   ])
            #Transcript Prioritization
            trs_prior_tab = objA.writeTranscriptPrioritization(configDict,
                                                               vcf_2_tab_bgz,
                                                               chrNum,projDate,memDict,
                                                               wh,batchNum,False
                                                              )
            #objA.writeQueryGenesUKBB(configDict,vcf_2_tab,chrNum,projDate,memDict,batchNum)
            """

            wh.close()

            if launch:
                os.system('bsub < '+sh_file)
    
    ''' 
    # Combine All chromosomes output

    
    if not re.search('ukbb',projDate):
        #mild_severe_list = ['mild','severe']
        #for ms_type in mild_severe_list:
        
        ms_type = re.split('_',projDate)[1]
        combChrTab,comChromFreqVcf = objA.writeCombChromTab(configDict,tmp_data,
                                                               objC,maf_cutoff,
                                                               ms_type,projDate,sh_type
                                                            )
    else:

        batchNum = re.split('_',projDate)[1]
        combChrTab,comChromFreqVcf = objA.writeCombChromTab(configDict,tmp_data,
                                                            objC,maf_cutoff,
                                                            batchNum,projDate,sh_type
                                                            )
    '''
