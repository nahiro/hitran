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
HITRAN2016_PARDIR = '/usr/local/HITRAN/HITRAN2016/HITRAN2016/By-Molecule/Uncompressed-files'
HITRAN2016_MOLPARAM = '/usr/local/HITRAN/HITRAN2016/Global-Data/molparam.txt'

# Default values
HVER = 2012          # HITRAN version
SPEC = 2             # Species number
RISO = 1.0           # Isotope abundance ratio
XMIN = 1.599e3       # Min wavelength in nm
XMAX = 1.601e3       # Max wavelength in nm
XUNI = 'nm'          # X unit
TMIN = 300.0         # Min temperature in K
TMAX = 300.0         # Min temperature in K
TSTP = 10.0          # Step temperature in K

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-D','--pardir',default=None,help='Data directory ({})'.format(HITRAN2012_PARDIR))
parser.add_option('-M','--molparam',default=None,help='Molecular parameter file ({})'.format(HITRAN2012_MOLPARAM))
parser.add_option('-H','--hver',default=HVER,type='int',help='HITRAN version (%default)')
parser.add_option('-N','--spec',default=SPEC,type='int',help='Species number (%default)')
parser.add_option('-I','--niso',default=None,type='int',help='Isotope number (%default)')
parser.add_option('-R','--riso',default=RISO,type='float',help='Isotope abundance ratio (%default)')
parser.add_option('-x','--xmin',default=XMIN,type='float',help='Min. X in XUNI (%default)')
parser.add_option('-X','--xmax',default=XMAX,type='float',help='Max. X in XUNI (%default)')
parser.add_option('-U','--xuni',default=XUNI,help='X unit (%default)')
parser.add_option('-y','--ymin',default=None,type='float',help='Min. Y (%default)')
parser.add_option('-Y','--ymax',default=None,type='float',help='Max. Y (%default)')
parser.add_option('-t','--tmin',default=TMIN,type='float',help='Min. temperature in K (%default)')
parser.add_option('-T','--tmax',default=TMAX,type='float',help='Max. temperature in K (%default)')
parser.add_option('-S','--tstp',default=TSTP,type='float',help='Step temperature in K (%default)')
parser.add_option('-s','--show',default=False,action='store_true',help='Show species (%default)')
parser.add_option('-l','--lmod',default=False,action='store_true',help='Long mode (%default)')
parser.add_option('-b','--batch',default=False,action='store_true',help='Batch mode (%default)')
(opts,args) = parser.parse_args()

if opts.pardir is None:
    if opts.hver == 2004:
        opts.pardir = HITRAN2004_PARDIR
    elif opts.hver == 2008:
        opts.pardir = HITRAN2008_PARDIR
    elif opts.hver == 2016:
        opts.pardir = HITRAN2016_PARDIR
    else:
        opts.pardir = HITRAN2012_PARDIR
if opts.molparam is None:
    if opts.hver == 2004:
        opts.molparam = HITRAN2004_MOLPARAM
    elif opts.hver == 2008:
        opts.molparam = HITRAN2008_MOLPARAM
    elif opts.hver == 2016:
        opts.molparam = HITRAN2016_MOLPARAM
    else:
        opts.molparam = HITRAN2012_MOLPARAM

# Read molparam
nam = None
mol = OrderedDict()
iso = OrderedDict()
rat = OrderedDict()
with open(opts.molparam,'r') as fp:
    for line in fp:
        m = re.search('(\S+)\s*\(([\d\s]*)\)',line)
        if m:
            nam = m.group(1)
            num = int(m.group(2))
            mol[num] = nam
            iso[nam] = []
            rat[nam] = []
            if len(mol) != num:
                raise ValueError('Error, len(mol)={}, num={} >>> {}'.format(len(mol),num,line))
        elif not re.search('Mol',line):
            item = line.split()
            if len(item) != 5:
                raise ValueError('Error, len(item)={} >>> {}'.format(len(item),line))
            iso[nam].append(item[0])
            rat[nam].append(float(item[1]))
if opts.show:
    for i,spec in mol.iteritems():
        sys.stdout.write('{:2d} {:6s}'.format(i,spec))
        for j in range(len(iso[spec])):
            sys.stdout.write(' ({:2d}) {:4s}'.format(j+1,iso[spec][j]))
        sys.stdout.write('\n')
    sys.exit(0)
if not mol.has_key(opts.spec):
    raise LookupError('Invalid species number >>> {}'.format(opts.spec))

spec = mol[opts.spec]
if opts.niso is not None:
    fnam = 'hitran_intensity_{}_{}.dat'.format(spec,iso[spec][opts.niso-1])
else:
    fnam = 'hitran_intensity_{}.dat'.format(spec)
command = 'hitran_intensity -n {} -x {:.6e} -X {:.6e}'.format(opts.hver,opts.xmin,opts.xmax)
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
title = '{}'.format(re.sub('(\d+)',r'$_{\1}$',spec))
if opts.niso is not None:
    command += ' -M {}'.format(opts.niso)
    title += '_{}'.format(iso[spec][opts.niso-1])
    if opts.riso > 0.0:
        riso = rat[spec][opts.niso-1]/opts.riso
        command += ' -r {:.6e}'.format(riso)
        title += ' ({:.2e})'.format(opts.riso)
    else:
        title += ' ({:.2e})'.format(rat[spec][opts.niso-1])
command += ' -t {:.6e} -T {:.6e} -S {:.6e}'.format(opts.tmin,opts.tmax,opts.tstp)
if opts.lmod:
    command += ' -l'
command += ' <{} >{}'.format(os.path.join(opts.pardir,'{:02d}_hit{:02d}.par'.format(opts.spec,opts.hver%100)),fnam)
call(command,shell=True)
if opts.tmax < opts.tmin:
    x,y = np.loadtxt(fnam,usecols=(0,1),unpack=True)
    t = np.repeat(opts.tmin,x.size)
    ts = np.array([opts.tmin])
else:
    x,t,y = np.loadtxt(fnam,usecols=(0,1,2),unpack=True)
    ts = np.arange(opts.tmin,opts.tmax+1.0e-5*opts.tstp,opts.tstp)
if not opts.batch:
    plt.interactive(True)
fig = plt.figure(1,facecolor='w')
fig.clear()
plt.subplots_adjust(top=0.86,bottom=0.12)
ax1 = plt.subplot(111)
ax1.minorticks_on()
ax1.grid(True)
ax1.set_yscale('log')
tmin = np.nanmin(ts)
tmax = np.nanmax(ts)
if tmax > tmin:
    tdif = tmax-tmin
else:
    tdif = 1.0
if x.size == 1:
    for tmp in reversed(ts):
        if np.fabs(t-tmp)<1.0e-5:
            yc = y
            ym = y*0.1
            ax1.vlines(x,ym,yc,color=cm.jet((tmp-tmin)/tdif),label='{:.0f}'.format(tmp))
elif x.size > 1:
    for tmp in reversed(ts):
        cnd = (np.fabs(t-tmp)<1.0e-5)
        if cnd.sum() < 1:
            continue
        yc = y[cnd]
        ym = np.nanmin(yc)
        ax1.vlines(x[cnd],ym,yc,color=cm.jet((tmp-tmin)/tdif),label='{:.0f}'.format(tmp))
ax1.set_xlim(opts.xmin,opts.xmax)
if opts.ymin is not None:
    ax1.set_ylim(bottom=opts.ymin)
if opts.ymax is not None:
    ax1.set_ylim(top=opts.ymax)
ax1.set_title(title,y=1.06)
ax1.set_xlabel(unit)
ax1.set_ylabel('Intensity (cm$^{-1}\cdot$cm$^{2}$/molecule)')
ax1.xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
ax1.xaxis.set_tick_params(pad=7)
lg = plt.legend(title='Temperature (K)',ncol=4)
lg.get_frame().set_facecolor('None')
plt.savefig('hitran_intensity_{}.pdf'.format(spec))
if not opts.batch:
    plt.draw()
