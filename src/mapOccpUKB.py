#!/usr/bin/python

import re,sys,os
from collections import OrderedDict
import math

if __name__=='__main__':

    inpFile = sys.argv[1]
    outFile = sys.argv[2]

    ocpDict = OrderedDict()

    fh = open(inpFile)

    for lines in fh:
        lines = lines.strip()
        strs =re.split('\t',lines)
        if not re.search('^coding',lines):
            code_id = int(strs[0])
            desc = strs[1]
            node_id = strs[2]
            par_id = strs[3]
            select = strs[4]

            ocpDict[code_id] = {'desc':desc,'node':node_id,
                                'par':par_id,'sel':select,'supar':''
                               }

    ocpDict_2 = OrderedDict()
    for keys in ocpDict:
        code_id = keys
        #if code_id > 100:
        if code_id == 1111:
            tmp_code = code_id
            while tmp_code >10:
                print(tmp_code)
                frac,whole = math.modf(int(ocpDict[tmp_code]['par'])/10)
                print(frac,whole)
                tmp_code = int(whole)
                print(tmp_code)
            ocpDict[code_id]['supar'] = tmp_code
            sys.exit()
        else:
            ocpDict[code_id]['supar'] = ocpDict[code_id]['par']

        
    for keys,values in ocpDict.items():
        print(values)


                





