#!/usr/bin/env python
"""

Purpose:
  Run Emperical Orthogonal Function Analyisis on Arbitrary EcoFOCI Timeseries data

Input:
    EPIC files:
    TODO: CF files

Output:
    EPIC EOF files (netcdf)
    Summary file (txt)
    TODO: Plots of EOF vs original Timeseries

Additional requirements for this example:

    * netCDF4 (http://unidata.github.io/netcdf4-python/)
    * matplotlib (http://matplotlib.org/)
    * eofs

  eofs: package developed and available:
    - https://github.com/ajdawson/eofs
    - http://ajdawson.github.io/eofs/ 
    - http://doi.org/10.5334/jors.122 (journal article)

    - developed using V 1.3.0 (updated to 1.4.0)

Addtional Notes:
    Tested on python=3.8


"""
from __future__ import print_function

import argparse
import datetime
import sys

import matplotlib.pyplot as plt
import numpy as np
from eofs.standard import Eof
from netCDF4 import Dataset

from calc.EPIC2Datetime import Datetime2EPIC, EPIC2Datetime, get_UDUNITS
from io_utils.ConfigParserLocal import get_config
from io_utils.EcoFOCI_netCDF_write import NetCDF_Create_Timeseries

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2018, 6, 13)
__modified__ = datetime.datetime(2018, 6, 13)
__version__  = "0.1.0"
__status__   = "Development"

"""----------------------------------MAIN-----------------------------------"""

# parse incoming command line options
parser = argparse.ArgumentParser(description='Analyze timeseries EOF')
parser.add_argument('pfile', 
    metavar='pfile', 
    type=str,
    help='pointer file with full paths on each line')
parser.add_argument('varname',
    metavar='varname',
    type=str,
    help='name of variable, may be EPIC name')
parser.add_argument('config_file_name', 
    metavar='config_file_name', 
    type=str, 
    help='full path to config file - eof_config.yaml')
parser.add_argument('-o','--outfile', 
    type=str,
    default='data/eof_results',
    help='name of output file data/{run-name} defaults to data/eof_results')
parser.add_argument('-s', '--start_date',
    type=str, 
    help='yyyymmddhhmmss')
parser.add_argument('-e', '--end_date',
    type=str, 
    help='yyyymmddhhmmss')
parser.add_argument('--eof_num',
    type=int,
    default=10000,
    help='number of eofs. default is all (arbitrarily large number)')
parser.add_argument('--epic', 
    action="store_true", 
    help='assume EPIC time format')
parser.add_argument('--plots', 
    action="store_true", 
    help='output some basic plots - TODO')
parser.add_argument('--summary', 
    action="store_true", 
    help='output summary only')
parser.add_argument('--normalize', 
    action="store_true", 
    help='normalize each timeseries by dividing by the standard deviation')

args = parser.parse_args()

### Watch for missing dates
if not args.start_date:
    sys.exit("Exiting: No start data was specified.")
else:
    start_date=datetime.datetime.strptime(args.start_date,'%Y%m%d%H%M%S')

if not args.end_date:
    print("No end date was specified.  Setting end date to now")
    end_date=datetime.datetime.today()
else:
    end_date=datetime.datetime.strptime(args.end_date,'%Y%m%d%H%M%S')

### Associate key variables names
if args.varname in ['u_1205']:
    altvarname = 'U_320'
elif args.varname in ['U_320']:
    altvarname = 'u_1205'
elif args.varname in ['wu']:
    args.varname = 'U_320'
    altvarname = 'u_1205'
elif args.varname in ['v_1206']:
    altvarname = 'V_321'
elif args.varname in ['V_321']:
    altvarname = 'v_1206'
elif args.varname in ['wv']:
    args.varname = 'V_321'
    altvarname = 'v_1206'
else:
    altvarname = ''

fcount = 0
files = {}
with open(args.pfile) as fp:
    for line in fp:
        files.update({fcount:line.strip()})
        fcount +=1


### EPIC Flavored time word
if args.epic:
    for key,filename in (files.items()):
        print("Reading file for {}".format(filename))

        ncin = Dataset(str(filename), 'r')
        try:
            data = ncin.variables[args.varname][:,0,0,0]
        except KeyError:
            data = ncin.variables[altvarname][:,0,0,0]
        ncdata = {'time':ncin.variables['time'][:],
                  'time2':ncin.variables['time2'][:]}
        ncin.close()

        #Convert two word EPIC time to python datetime.datetime 
        # representation and then format for CF standards
        dt_from_epic =  np.array(EPIC2Datetime(ncdata['time'], ncdata['time2']))

        #Find relevant chuck of times to keep based on arguments
        dt_index = np.where((dt_from_epic >= start_date) & 
                            (dt_from_epic <= end_date) )

        if key == 0:
            eof_data = data[dt_index]
        else:
            try:
                eof_data = np.vstack((eof_data,data[dt_index]))
            except ValueError:
                sys.exit("Exiting: timeseries have different lengths")

    if args.normalize:
        eof_data_std = np.std(eof_data, axis=1)
        eof_data = eof_data.T / np.std(eof_data, axis=1)
    else:
        #transpose so time is first dimension
        eof_data = eof_data.T

# Crete an EOF solver to do the EOF analysis.  No weights
# First dimension is assumed time by program... not true if timseries is of interest, 
print("Solving for n={} modes".format(args.eof_num))
solver = Eof(eof_data, center=False)
pcs = solver.pcs(npcs=args.eof_num)
eigval = solver.eigenvalues(neigs=args.eof_num)
varfrac = solver.varianceFraction(neigs=args.eof_num)
eofs = solver.eofs(neofs=args.eof_num)
eofcorr = solver.eofsAsCorrelation(neofs=args.eof_num)
eofcov = solver.eofsAsCovariance(neofs=args.eof_num)


"""---------------------------------Report-----------------------------------"""
### Print Select Results to file
outfile = args.outfile+'.txt'
print("EOF Results:", file=open(outfile,"w"))
print("------------", file=open(outfile,"a"))

print("File path: {}".format("/".join(filename.split('/')[:-1])),
    file=open(outfile,"a"))
for key,filename in (files.items()):
    print("Files input: {}".format(filename.split('/')[-1]),
        file=open(outfile,"a"))

print("\n\n", file=open(outfile,"a"))
print("Variables used: ", file=open(outfile,"a"))
print("----------------", file=open(outfile,"a"))
print("{}, {}".format(args.varname,altvarname), file=open(outfile,"a"))
print("\n\n", file=open(outfile,"a"))

print("Start/Stop Date: ", file=open(outfile,"a"))
print("-----------------", file=open(outfile,"a"))
print("{}, {}".format(start_date,end_date), file=open(outfile,"a"))
print("\n\n", file=open(outfile,"a"))

print("mode file names:", file=open(outfile,"a"))
print("---------------", file=open(outfile,"a"))
print("\n",            file=open(outfile,"a"))

for index in range(0,eofs.shape[0],1):
    print("File output: {0}_eof{1}.nc".format(args.outfile,str(index+1).zfill(3)),
        file=open(outfile,"a"))
print("\n\n",            file=open(outfile,"a"))


print("Timeseries Normalized by:", file=open(outfile,"a"))
print("-------------------------", file=open(outfile,"a"))
if args.normalize:
    print(",\t".join([str(x) for x in eof_data_std]), file=open(outfile,"a"))
else:
    print("***Timeseries not normalized***", file=open(outfile,"a"))
print("\n\n",            file=open(outfile,"a"))


print("EOFs :", file=open(outfile,"a"))
print("-----", file=open(outfile,"a"))
print('\n'.join(["Mode {}:".format(y+1)+' '.join(['{:6.2f}'.format(item) for item in row]) 
      for y,row in enumerate(eofs)]), file=open(outfile,"a"))
print("\n\n", file=open(outfile,"a"))
 
print("EOF as Covariance :", file=open(outfile,"a"))
print("------------------", file=open(outfile,"a"))
print('\n'.join(["Mode {}:".format(y+1)+' '.join(['{:6.2f}'.format(item) for item in row]) 
      for y,row in enumerate(eofcov)]), file=open(outfile,"a"))
print("\n\n", file=open(outfile,"a"))

print("EOFs as Correlation :", file=open(outfile,"a"))
print("--------------------", file=open(outfile,"a"))
print('\n'.join(["Mode {}:".format(y+1)+' '.join(['{:6.2f}'.format(item) for item in row]) 
      for y,row in enumerate(eofcorr)]), file=open(outfile,"a"))
print("\n\n", file=open(outfile,"a"))

print("EigenValues: (largest to smallest)", file=open(outfile,"a"))
print("---------------------------------", file=open(outfile,"a"))
print(",\t".join([str(x) for x in eigval]), file=open(outfile,"a"))
print("\n\n", file=open(outfile,"a"))
 
print("Total Variance fraction for each EOF mode: (0 to 1)", file=open(outfile,"a"))
print("---------------------------------------------------", file=open(outfile,"a"))
print("\n".join(["Mode {y}:{x}".format(y=y,x=x) for y,x in enumerate(varfrac)]),
                                           file=open(outfile,"a"))
print("\n\n", file=open(outfile,"a"))

"""---------------------------------NetCDF-----------------------------------"""
### Create EOF nc files
if args.epic and not args.summary:

    # From config file, get variable attribute definitions
    EPIC_VARS_dict = get_config(args.config_file_name,'yaml')

    for index in range(0,eofs.shape[0],1):
        print("Creating EPIC file for {0}_eof{1}.nc".format(args.outfile,str(index+1).zfill(3)))

        # Link data to a dictionary to match variable names
        data_dic = {'PCS_6001':pcs.T[index]}

        time1,time2 = np.array(Datetime2EPIC(list(dt_from_epic[dt_index])), dtype='f8')
        ncinstance = NetCDF_Create_Timeseries(
            savefile="{0}_eof{1}.nc".format(args.outfile,str(index+1).zfill(3)))
        ncinstance.file_create()
        ncinstance.sbeglobal_atts(History='Run name: {}'.format(args.outfile), 
                                  Software='EcoFOCI_eof.py ' + __version__)
        ncinstance.dimension_init(time_len=len(time1))
        ncinstance.variable_init(EPIC_VARS_dict)
        ncinstance.add_coord_data(time1=time1, time2=time2)
        ncinstance.add_data(EPIC_VARS_dict,data_dic=data_dic)
        ncinstance.close()


if args.plots:
    if args.eof_num < 5:
        nmax = args.eof_num
    else:
        nmax = 5
    
    #plot PCs
    fig = plt.figure()
    ax = plt.subplot(111)

    for index in range(0,nmax,1):
        plt.plot(dt_from_epic[dt_index],pcs.T[index],label='PC mode:{}'.format(index+1))

    plt.legend()
    fig.set_size_inches( (22, 8.5) )
    plt.savefig("{0}_pcs{1}".format(args.outfile,str(index+1).zfill(3))+'.png',bbox_inches='tight', dpi=(300))

    #plot eigenvectors / EOF maps as corr
    fig = plt.figure()
    ax = plt.subplot(111)

    for index in range(0,nmax,1):
        plt.plot(eofcorr[index],range(0,len(eof_data[1])),
            label='EOF Corr. mode:{}'.format(index+1))
    
    ax.invert_yaxis()
    plt.legend()
    fig.set_size_inches( (4.25, 11) )
    plt.savefig("{0}_eofcor{1}".format(args.outfile,str(index+1).zfill(3))+'.png',bbox_inches='tight', dpi=(300))

    #plot eigenvectors / EOF maps as cov
    fig = plt.figure()
    ax = plt.subplot(111)

    for index in range(0,nmax,1):
        plt.plot(eofcov[index],range(0,len(eof_data[1])),
            label='EOF Cov. mode:{}'.format(index+1))
    
    ax.invert_yaxis()
    plt.legend()
    fig.set_size_inches( (4.25, 11) )
    plt.savefig("{0}_eofcov{1}".format(args.outfile,str(index+1).zfill(3))+'.png',bbox_inches='tight', dpi=(300))

    #plot eigenvectors / EOF maps 
    fig = plt.figure()
    ax = plt.subplot(111)

    for index in range(0,nmax,1):
        plt.plot(eofs[index],range(0,len(eof_data[1])),
            label='EOFs mode:{}'.format(index+1))
    
    ax.invert_yaxis()
    plt.legend()
    fig.set_size_inches( (4.25, 11) )
    plt.savefig("{0}_eofs{1}".format(args.outfile,str(index+1).zfill(3))+'.png',bbox_inches='tight', dpi=(300))
