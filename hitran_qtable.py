#!/usr/bin/env python
import sys
import numpy as np
from collections import OrderedDict
from subprocess import check_output

iso = OrderedDict()
iso.update({'H2O':[161,181,171,162,182,172]})
iso.update({'CO2':[626,636,628,627,638,637,828,728,727,838,837]})
iso.update({'O3':[666,668,686,667,676,886,868,678,768,786,776,767,888,887,878,778,787,777]})
iso.update({'N2O':[446,456,546,448,447]})
iso.update({'CO':[26,36,28,27,38,37]})
iso.update({'CH4':[211,311,212,312]})
iso.update({'O2':[66,68,67]})
iso.update({'NO':[46,56,48]})
iso.update({'SO2':[626,646]})
iso.update({'NO2':[646]})
iso.update({'NH3':[4111,5111]})
iso.update({'HNO3':[146]})
iso.update({'OH':[61,81,62]})
iso.update({'HF':[19]})
iso.update({'HCl':[15,17]})
iso.update({'HBr':[19,11]})
iso.update({'HI':[17]})
iso.update({'ClO':[56,76]})
iso.update({'OCS':[622,624,632,623,822]})
iso.update({'H2CO':[126,136,128]})
iso.update({'HOCl':[165,167]})
iso.update({'N2':[44]})
iso.update({'HCN':[124,134,125]})
iso.update({'CH3Cl':[215,217]})
iso.update({'H2O2':[1661]})
iso.update({'C2H2':[1221,1231,1222]})
iso.update({'C2H6':[1221,1231]})
iso.update({'PH3':[1111]})
iso.update({'COF2':[269]})
iso.update({'SF6':[29]})
iso.update({'H2S':[121,141,131]})
iso.update({'HCOOH':[126]})
iso.update({'HO2':[166]})
iso.update({'O':[6]})
iso.update({'ClONO2':[5646,7646]})
iso.update({'NO+':[46]})
iso.update({'HOBr':[169,161]})
iso.update({'C2H4':[221,231]})
iso.update({'CH3OH':[2161]})
iso.update({'CH3Br':[219,211]})
iso.update({'CH3CN':[2124,2134,3124,3134]})
iso.update({'CF4':[29]})
iso.update({'C4H2':[1221]})
iso.update({'HC3N':[12224,12234,12324,13224,12225,22224]})
iso.update({'C2N2':[4224,5225]})
iso.update({'CS':[22,24,32,23]})
iso.update({'H2':[11,12]})
iso.update({'SO':[26,46,28]})
iso.update({'C3H4':[1221]})
iso.update({'CH3':[2111]})
iso.update({'CS2':[222,224,223,232]})

flag = True
sys.stdout.write('Temp(K)')
for spec,nums in iso.iteritems():
    for i in nums:
        name = '%s_%d'%(spec,i)
        if flag:
            sys.stdout.write('%19s '%(name))
            flag = False
        else:
            sys.stdout.write('%26s '%(name))
sys.stdout.write('\n')
for t in np.arange(70.0,3000.1,1.0):
    flag = True
    sys.stdout.write('%6.1f'%(t))
    for i,nums in enumerate(iso.values()):
        for j in range(len(nums)):
            command = 'hitran_tips -n 2012 -o %s -m %d -i %d'%(t,i+1,j+1)
            out = check_output(command,shell=True)
            item = out.split()
            if len(item) != 4:
                raise ValueError('Error, len(item)=%d'%(len(item)))
            if flag:
                sys.stdout.write(' %19.6f'%(float(item[3])))
                flag = False
            else:
                sys.stdout.write(' %26.6f'%(float(item[3])))
    sys.stdout.write('\n')
