#!/usr/bin/env python
"""

Purpose:
  Run Emperical Orthogonal Function Analyisis on Arbitrary EcoFOCI Timeseries data

EPIC files:

CF files:

Output:


Additional requirements for this example:

    * netCDF4 (http://unidata.github.io/netcdf4-python/)
    * matplotlib (http://matplotlib.org/)

  eofs: package developed and available:
    - https://github.com/ajdawson/eofs
    - http://ajdawson.github.io/eofs/ 

"""

#System Stack
import argparse
import datetime
import sys

#Science Stack
from eofs.standard import Eof
import numpy

#User defined Stack
from calc.EPIC2Datetime import EPIC2Datetime, get_UDUNITS, Datetime2EPIC

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
parser.add_argument('-o','--outfile', 
    type=str,
    default='eof_results.txt',
    help='name of output file ')
parser.add_argument('-s', '--start_date',
    type=str, 
    help='yyyymmddhhmmss')
parser.add_argument('-e', '--end_date',
    type=str, 
    help='yyyymmddhhmmss')
parser.add_argument('--eof_num',
    type=int,
    default=1e35,
    help='number of eofs. default is all (arbitrarily large number)')
parser.add_argument('--epic', 
    action="store_true", 
    help='assume EPIC time format')

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
elif args.varname in ['v_1206']:
    altvarname = 'V_321'
elif args.varname in ['V_321']:
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

    #transpose so time is first dimension
    #eof_data = eof_data.T

# Crete an EOF solver to do the EOF analysis.  No weights
# First dimension is assumed time by program... not true if timseries is of interest, 
solver = Eof(eof_data, center=False)
eigval = solver.eigenvalues(neigs=args.eof_num)
varfrac = solver.varianceFraction(neigs=args.eof_num)
eofs = solver.eofs(neofs=args.eof_num)


"""---------------------------------Report-----------------------------------"""
### Print Select Results to file

print("EOF Results:", file=open(args.outfile, "w"))
print("------------", file=open(args.outfile, "a"))

for key,filename in (files.items()):
    print("Files input: {}".format(filename.split('/')[-1]),
        file=open(args.outfile, "a"))

print("\n\n", file=open(args.outfile, "a"))
print("Variables used: ", file=open(args.outfile, "a"))
print("----------------", file=open(args.outfile, "a"))
print("{}, {}".format(args.varname,altvarname), file=open(args.outfile, "a"))
print("\n\n", file=open(args.outfile, "a"))

print("eof file names:", file=open(args.outfile, "a"))
print("---------------", file=open(args.outfile, "a"))
print("\n",            file=open(args.outfile, "a"))

print("File path: {}".format("/".join(filename.split('/')[:-1])),
    file=open(args.outfile, "a"))
for index in range(0,eofs.shape[0]+1,1):
    print("File output: {1}eof_{2}.nc".format(args.outfile,str(index).zfill(3)),
        file=open(args.outfile, "a"))
print("\n\n",            file=open(args.outfile, "a"))

print("EigenValues: (largest to smallest)", file=open(args.outfile, "a"))
print("---------------------------------", file=open(args.outfile, "a"))
print("\t".join([str(x) for x in eigval]), file=open(args.outfile, "a"))
print("\n\n", file=open(args.outfile, "a"))
 
print("Total Variance explained by each EOF mode: (0 to 1)", file=open(args.outfile, "a"))
print("---------------------------------------------------", file=open(args.outfile, "a"))
print("\t".join([str(x) for x in varfrac]), file=open(args.outfile, "a"))
print("\n\n", file=open(args.outfile, "a"))

### Create EOF nc files










