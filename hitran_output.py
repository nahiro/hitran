#!/usr/bin/env python
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import OrderedDict
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

# Constants
HITRAN2004_PARDIR = '/usr/local/HITRAN/hitran04/HITRAN2004/By-Molecule/Uncompressed-files'
HITRAN2004_MOLPARAM = '/usr/local/HITRAN/hitran04/Global_Data/molparam.txt'
HITRAN2008_PARDIR = '/usr/local/HITRAN/hitran08/HITRAN2008/By-Molecule/Uncompressed-files'
HITRAN2008_MOLPARAM = '/usr/local/HITRAN/hitran08/Global_Data/molparam.txt'
HITRAN2012_PARDIR = '/usr/local/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Uncompressed-files'
HITRAN2012_MOLPARAM = '/usr/local/HITRAN/HITRAN2012/Global-Data/molparam.txt'
PNAM = {0:'M',1:'I',2:'Wavelength',3:'Intensity',4:'A-coef.',5:'$\gamma_{air}$',
6:'$\gamma_{self}$',7:'E$_{low}$',8:'$n_{air}$',9:'$\delta_{air}$',}
UNITS = {0:'',1:'',2:' (nm)',3:' (cm$^{-1}$/(molecule$\cdot$cm$^{-2}$))',4:' (s$^{-1}$)',
5:' (cm$^{-1}$)',6:' (cm$^{-1})$',7:' (cm$^{-1})$',8:'',9:' (cm$^{-1}$/atm)'}

# Default values
HVER = 2012          # HITRAN version
SPEC = 2             # Species number
NPAR = 5             # Parameter number
XMIN = 1.599e3       # Min wavelength in nm
XMAX = 1.601e3       # Max wavelength in nm
XUNI = 'nm'          # X unit

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-D','--pardir',default=None,help='Data directory (%s)'%(HITRAN2012_PARDIR))
parser.add_option('-M','--molparam',default=None,help='Molecular parameter file (%s)'%(HITRAN2012_MOLPARAM))
parser.add_option('-H','--hver',default=HVER,type='int',help='HITRAN version (%default)')
parser.add_option('-N','--spec',default=SPEC,type='int',help='Species number (%default)')
parser.add_option('-I','--niso',default=None,type='int',help='Isotope number (%default)')
parser.add_option('-n','--npar',default=NPAR,type='int',help='Parameter number (%default)')
parser.add_option('-x','--xmin',default=XMIN,type='float',help='Min wavelength in nm (%default)')
parser.add_option('-X','--xmax',default=XMAX,type='float',help='Max wavelength in nm (%default)')
parser.add_option('-U','--xuni',default=XUNI,help='X unit (%default)')
parser.add_option('-y','--ymin',default=None,type='float',help='Min. Y (%default)')
parser.add_option('-Y','--ymax',default=None,type='float',help='Max. Y (%default)')
parser.add_option('-s','--show',default=False,action='store_true',help='Show species (%default)')
parser.add_option('-b','--batch',default=False,action='store_true',help='Batch mode (%default)')
(opts,args) = parser.parse_args()

if opts.pardir is None:
    if opts.hver == 2004:
        opts.pardir = HITRAN2004_PARDIR
    elif opts.hver == 2008:
        opts.pardir = HITRAN2008_PARDIR
    else:
        opts.pardir = HITRAN2012_PARDIR
if opts.molparam is None:
    if opts.hver == 2004:
        opts.molparam = HITRAN2004_MOLPARAM
    elif opts.hver == 2008:
        opts.molparam = HITRAN2008_MOLPARAM
    else:
        opts.molparam = HITRAN2012_MOLPARAM

# Read molparam
nam = None
mol = OrderedDict()
iso = OrderedDict()
with open(opts.molparam,'r') as fp:
    for line in fp:
        m = re.search('(\S+)\s*\(([\d\s]*)\)',line)
        if m:
            nam = m.group(1)
            num = int(m.group(2))
            mol[num] = nam
            iso[nam] = []
            if len(mol) != num:
                raise ValueError('Error, len(mol)=%d, num=%d >>> %s'%(len(mol),num,line))
        elif not re.search('Mol',line):
            item = line.split()
            if len(item) != 5:
                raise ValueError('Error, len(item)=%d >>> %s'%(len(item),line))
            iso[nam].append(item[0])
if opts.show:
    for i,spec in mol.iteritems():
        print '%2d %-6s'%(i,spec),
        for j in range(len(iso[spec])):
            print ' (%d) %-4s'%(j+1,iso[spec][j]),
        print ''
    sys.exit(0)
if not mol.has_key(opts.spec):
    raise LookupError('Invalid species number >>> %d'%(opts.spec))

spec = mol[opts.spec]
if opts.niso is not None:
    fnam = 'hitran_output_%s_%s.dat'%(spec,iso[spec][opts.niso-1])
else:
    fnam = 'hitran_output_%s.dat'%(spec)
command = 'hitran_output -n %d -x %.6e -X %.6e'%(opts.hver,opts.xmin,opts.xmax)
if opts.xuni == 'nm':
    command += ' -U 1'
    unit = 'Wavelength (nm)'
elif opts.xuni == 'cm-1':
    command += ' -U 2'
    unit = 'Wavenumber (cm$^{-1}$)'
elif opts.xuni == 'Hz':
    command += ' -U 3'
    unit = 'Frequency (Hz)'
elif opts.xuni == 'GHz':
    command += ' -U 4'
    unit = 'Frequency (GHz)'
else:
    raise ValueError('Error, XUNI should be nm, cm-1, Hz, or GHz')
title = '%s'%(re.sub('(\d+)',r'$_{\1}$',spec))
if opts.niso is not None:
    command += ' -N %d'%(opts.niso)
    title += '_%s'%(iso[spec][opts.niso-1])
command += ' <%s/%02d_hit%02d.par >%s'%(opts.pardir,opts.spec,opts.hver%100,fnam)
call(command,shell=True)
x,y = np.loadtxt(fnam,unpack=True,usecols=(2,opts.npar))
if not opts.batch:
    plt.interactive(True)
fig = plt.figure(1,facecolor='w')
fig.clear()
plt.subplots_adjust(top=0.86,bottom=0.12,left=0.15)
ax1 = plt.subplot(111)
if opts.npar == 3:
    ax1.set_yscale('log')
ax1.minorticks_on()
ax1.grid(True)
ax1.plot(x,y)
ax1.set_xlim(opts.xmin,opts.xmax)
if opts.ymin is not None:
    ax1.set_ylim(bottom=opts.ymin)
if opts.ymax is not None:
    ax1.set_ylim(top=opts.ymax)
ax1.set_title(title,y=1.06)
ax1.set_xlabel(unit)
ax1.set_ylabel('%s%s'%(PNAM[opts.npar],UNITS[opts.npar]))
ax1.xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
ax1.xaxis.set_tick_params(pad=7)
if opts.niso is not None:
    plt.savefig('hitran_output_%s_%s.pdf'%(spec,iso[spec][opts.niso-1]))
else:
    plt.savefig('hitran_output_%s.pdf'%(spec))
if not opts.batch:
    plt.draw()
