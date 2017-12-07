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

# Default values
HVER = 2012          # HITRAN version
SPEC = 2             # Species number
RISO = 1.0           # Isotope abundance ratio
XMIN = 1.599e3       # Min wavelength in nm
XMAX = 1.601e3       # Max wavelength in nm
XSTP = 0.01          # Step wavelength in nm
XUNI = 'nm'          # X unit
WGAM = 1.0e5         # Lorentz width in gamma
TEMP = 300.0         # Temperature in K
PRES = 1.0           # Pressure in atm
MIXR = 390.0e-6      # Mixing ratio in mol/mol

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-D','--pardir',default=None,help='Data directory (%s)'%(HITRAN2012_PARDIR))
parser.add_option('-M','--molparam',default=None,help='Molecular parameter file (%s)'%(HITRAN2012_MOLPARAM))
parser.add_option('-H','--hver',default=HVER,type='int',help='HITRAN version (%default)')
parser.add_option('-N','--spec',default=SPEC,type='int',help='Species number (%default)')
parser.add_option('-I','--niso',default=None,type='int',help='Isotope number (%default)')
parser.add_option('-R','--riso',default=RISO,type='float',help='Isotope abundance ratio (%default)')
parser.add_option('-x','--xmin',default=XMIN,type='float',help='Min. X in XUNI (%default)')
parser.add_option('-X','--xmax',default=XMAX,type='float',help='Max. X in XUNI (%default)')
parser.add_option('-S','--xstp',default=XSTP,type='float',help='Step X in XUNI (%default)')
parser.add_option('-U','--xuni',default=XUNI,help='X unit (%default)')
parser.add_option('-y','--ymin',default=None,type='float',help='Min. Y (%default)')
parser.add_option('-Y','--ymax',default=None,type='float',help='Max. Y (%default)')
parser.add_option('-W','--wgam',default=WGAM,type='float',help='Lorentz width in gamma (%default)')
parser.add_option('-t','--temp',default=TEMP,type='float',help='Temperature in K (%default)')
parser.add_option('-p','--pres',default=PRES,type='float',help='Pressure in atm (%default)')
parser.add_option('-r','--mixr',default=MIXR,type='float',help='Mixing ratio in mol/mol (%default)')
parser.add_option('-l','--lorentz',default=False,action='store_true',help='Lorentz profile (pressure broadening only) (%default)')
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
rat = OrderedDict()
mas = OrderedDict()
with open(opts.molparam,'r') as fp:
    for line in fp:
        m = re.search('(\S+)\s*\(([\d\s]*)\)',line)
        if m:
            nam = m.group(1)
            num = int(m.group(2))
            mol[num] = nam
            iso[nam] = []
            rat[nam] = []
            mas[nam] = []
            if len(mol) != num:
                raise ValueError('Error, len(mol)=%d, num=%d >>> %s'%(len(mol),num,line))
        elif not re.search('Mol',line):
            item = line.split()
            if len(item) != 5:
                raise ValueError('Error, len(item)=%d >>> %s'%(len(item),line))
            iso[nam].append(item[0])
            rat[nam].append(float(item[1]))
            mas[nam].append(float(item[4]))
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
    fnam = 'hitran_absorption_%s_%s.dat'%(spec,iso[spec][opts.niso-1])
else:
    fnam = 'hitran_absorption_%s.dat'%(spec)
command = 'hitran_absorption -n %d -x %.6e -X %.6e'%(opts.hver,opts.xmin,opts.xmax)
command += ' -D %.6e -W %.6e'%(opts.xstp,opts.wgam)
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
    if opts.riso > 0.0:
        riso = rat[spec][opts.niso-1]/opts.riso
        command += ' -r %.6e'%(riso)
        title += ' (%.2e)'%(opts.riso)
    else:
        title += ' (%.2e)'%(rat[spec][opts.niso-1])
if not opts.lorentz:
    if opts.niso is None:
        sr = 0.0
        sm = 0.0
        for j in range(len(iso[spec])):
            sr += rat[spec][j]
            sm += rat[spec][j]*mas[spec][j]
        mass = sm/sr
    else:
        mass = mas[spec][opts.niso-1]
    command += ' -M %.4f'%(mass)
command += ' -p %.6e -o %.6e -R %.6e'%(opts.pres,opts.temp,opts.mixr)
command += ' <%s/%02d_hit%02d.par >%s'%(opts.pardir,opts.spec,opts.hver%100,fnam)
call(command,shell=True)
x,y = np.loadtxt(fnam,unpack=True)
if not opts.batch:
    plt.interactive(True)
fig = plt.figure(1,facecolor='w')
fig.clear()
plt.subplots_adjust(top=0.86,bottom=0.12)
ax1 = plt.subplot(111)
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
ax1.text(0.5,1.01,'T: %.1f (K), P: %.2e (atm), VMR: %.2e'%(opts.temp,opts.pres,opts.mixr),
         transform=ax1.transAxes,ha='center',va='bottom',size=14)
ax1.set_xlabel(unit)
ax1.set_ylabel('Absorption Cross-Section (cm$^{2}$/molecule)')
ax1.xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
ax1.xaxis.set_tick_params(pad=7)
if opts.niso is not None:
    plt.savefig('hitran_absorption_%s_%s.pdf'%(spec,iso[spec][opts.niso-1]))
else:
    plt.savefig('hitran_absorption_%s.pdf'%(spec))
if not opts.batch:
    plt.draw()
