#!/usr/bin/env python

"""
 EcoFOCI_netCDF_write.py
 
 class for building netcdf files from specified instruments
 
 
  History:
 --------
 2018-06-14: Stripped down version of EcoFOCI_Utilitis/io_utils
"""

import datetime
import os

import numpy as np
from netCDF4 import Dataset

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 1, 13)
__modified__ = datetime.datetime(2014, 12, 2)
__version__  = "0.1.0"
__status__   = "Development"


"""-------------------------------NCFile Creation--------------------------------------"""

        
class NetCDF_Create_Timeseries(object):
    """ Class instance to generate a NetCDF file.  

    Standards
    ---------
    EPICNetCDF (PMEL) Standards  


    Usage
    -----
    
    Order of routines matters and no error checking currently exists
    ToDo: Error Checking
    
    Use this to create a nc file with all default values
        ncinstance = NetCDF_Create_Timeseries()
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.dimension_init()
        ncinstance.variable_init()
        ncinstance.add_coord_data()
        ncinstance.add_data()
        ncinstance.close()
    """ 
    
    
    nc_format = 'NETCDF3_CLASSIC'
    nc_read   = 'w'

    def __init__(self, savefile='data/test.nc'):
        """initialize output file path"""
        
        self.savefile = savefile
    
    def file_create(self):
            rootgrpID = Dataset(self.savefile, NetCDF_Create_Timeseries.nc_read, 
                                format=NetCDF_Create_Timeseries.nc_format)
            self.rootgrpID = rootgrpID
            return ( rootgrpID )
        
    def sbeglobal_atts(self, History='', Software=''):
        """
        
        """
        
        self.rootgrpID.CREATION_DATE = datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")
        self.rootgrpID.History = History
        self.rootgrpID.Software = Software
        
    def dimension_init(self, time_len=1):
        """
        Assumes
        -------
        Dimensions will be 'time', 'depth', 'lat', 'lon'
        
        Todo
        ----
        User defined dimensions
        """

        self.dim_vars = ['time']
        
        self.rootgrpID.createDimension( self.dim_vars[0], time_len ) #time
        
        
    def variable_init(self, EPIC_VARS_dict):
        """
        EPIC keys:
            passed in as a dictionary (similar syntax as json data file)
            The dictionary keys are what defines the variable names.
        """
        #exit if the variable dictionary is not passed
        if not bool(EPIC_VARS_dict):
            raise RuntimeError('Empty EPIC Dictionary is passed to variable_init.')

        #build record variable attributes
        rec_vars, rec_var_name, rec_var_longname = [], [], []
        rec_var_generic_name, rec_var_FORTRAN, rec_var_units, rec_var_epic = [], [], [], []

        #cycle through epic dictionary and create nc parameters
        for evar in EPIC_VARS_dict.keys():
            rec_vars.append(evar)
            rec_var_name.append( EPIC_VARS_dict[evar]['name'] )
            rec_var_longname.append( EPIC_VARS_dict[evar]['longname'] )
            rec_var_generic_name.append( EPIC_VARS_dict[evar]['generic_name'] )
            rec_var_units.append( EPIC_VARS_dict[evar]['units'] )
            rec_var_FORTRAN.append( EPIC_VARS_dict[evar]['fortran'] )
            rec_var_epic.append( EPIC_VARS_dict[evar]['EPIC_KEY'] )
        
        rec_vars = ['time','time2'] + rec_vars

        rec_var_name = ['', ''] + rec_var_name
        rec_var_longname = ['', ''] + rec_var_longname
        rec_var_generic_name = ['', ''] + rec_var_generic_name
        rec_var_FORTRAN = ['', ''] + rec_var_FORTRAN
        rec_var_units = ['True Julian Day', 'msec since 0:00 GMT'] + rec_var_units
        rec_var_type= ['i4', 'i4'] + ['f4' for spot in rec_vars[2:]]
        rec_var_strtype= ['EVEN', 'EVEN'] + ['' for spot in rec_vars[2:]]
        rec_epic_code = [624, 624] + rec_var_epic
        
        var_class = []
        var_class.append(self.rootgrpID.createVariable(rec_vars[0], rec_var_type[0], self.dim_vars[0]))#time1
        var_class.append(self.rootgrpID.createVariable(rec_vars[1], rec_var_type[1], self.dim_vars[0]))#time2

        for i, v in enumerate(rec_vars[2:]):  #1D coordinate variables
            var_class.append(self.rootgrpID.createVariable(rec_vars[i+2], rec_var_type[i+2], self.dim_vars))
             
        ### add variable attributes
        for i, v in enumerate(var_class): #4dimensional for all vars
            print("Adding Variable {0}".format(v))
            v.setncattr('name',rec_var_name[i])
            v.long_name = rec_var_longname[i]
            v.generic_name = rec_var_generic_name[i]
            v.FORTRAN_format = rec_var_FORTRAN[i]
            v.units = rec_var_units[i]
            v.type = rec_var_strtype[i]
            v.epic_code = rec_epic_code[i]
            
        self.var_class = var_class
        self.rec_vars = rec_vars

        
    def add_coord_data(self, time1=None, time2=None):
        """ """
        self.var_class[0][:] = time1
        self.var_class[1][:] = time2


    def add_data(self, EPIC_VARS_dict, data_dic=None, missing_values=1e35):
        """
            using the same dictionary to define the variables, and a new dictionary
                that associates each data array with an epic key, cycle through and populate
                the desired variables.  If a variable is defined in the epic keys but not passed
                to the add_data routine, it should be populated with missing data
        """
        #exit if the variable dictionary is not passed
        if not bool(EPIC_VARS_dict):
            raise RuntimeError('Empty EPIC Dictionary is passed to add_data.')
        
        #cycle through EPIC_Vars and populate with data - this is a comprehensive list of 
        # all variables expected
        # if no data is passed but an epic dictionary is, complete routine leaving variables
        #  with missing data if not found

        for EPICdic_key in EPIC_VARS_dict.keys():
            di = self.rec_vars.index(EPICdic_key)
            try:
                self.var_class[di][:] = data_dic[EPICdic_key]
            except KeyError:
                self.var_class[di][:] = missing_values
        
        
    def add_history(self, new_history):
        """Adds timestamp (UTC time) and history to existing information"""
        self.rootgrpID.History = self.rootgrpID.History + '\n'\
                                + datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")\
                                + ' ' + new_history
                    
    def close(self):
        self.rootgrpID.close()

class CF_NC(object):


    """ Class instance to generate a NetCDF file.  
    Assumes data format and information ingested is a dataframe object from ctd.py 

    Standards
    ---------
    EPICNetCDF (PMEL) Standards  


    Usage
    -----
    
    Order of routines matters and no error checking currently exists
    ToDo: Error Checking
    
    Use this to create a nc file with all default values
        ncinstance = CF_NC()
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.dimension_init()
        ncinstance.variable_init()
        ncinstance.add_coord_data()
        ncinstance.add_data()
        ncinstance.close()
    """ 
    
    
    nc_format = 'NETCDF3_CLASSIC'
    nc_read   = 'w'
    def __init__(self, savefile='ncfiles/test.nc'):
        """data is a numpy array of temperature values"""
        
        self.savefile = savefile
    
    def file_create(self):
            rootgrpID = Dataset(self.savefile, CF_NC.nc_read, format=CF_NC.nc_format)
            self.rootgrpID = rootgrpID
            return ( rootgrpID )
        
    def sbeglobal_atts(self, raw_data_file='', Water_Mass='B', Water_Depth=9999, Prog_Cmnt='',\
                        Experiment='', Edit_Cmnt='', Station_Name='', Inst_Type='', Project='', History=''):
        """
        Assumptions
        -----------
        
        Format of DataFrame.name = 'dy1309l1_ctd001'
        
        seabird related global attributes found in DataFrame.header list
        
        """
        
        self.rootgrpID.CREATION_DATE = datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")
        self.rootgrpID.COMPOSITE = 1
        self.rootgrpID.INST_TYPE = Inst_Type
        self.rootgrpID.DATA_CMNT = raw_data_file
        self.rootgrpID.EPIC_FILE_GENERATOR = 'nc_epic2udunits_time.py V' + __version__ 
        self.rootgrpID.PROG_CMNT01 = Prog_Cmnt
        self.rootgrpID.EDIT_CMNT01 = Edit_Cmnt
        self.rootgrpID.WATER_DEPTH = Water_Depth
        self.rootgrpID.MOORING = Station_Name
        self.rootgrpID.WATER_MASS = Water_Mass
        self.rootgrpID.EXPERIMENT = Experiment
        self.rootgrpID.PROJECT = Experiment
        self.rootgrpID.History = History
                        
        
    def dimension_init(self, time_len=1):
        """
        Assumes
        -------
        Dimensions will be 'time', 'depth', 'lat', 'lon'
        
        Todo
        ----
        User defined dimensions
        """

        self.dim_vars = ['time', 'depth', 'lat', 'lon']
        
        self.rootgrpID.createDimension( self.dim_vars[0], time_len ) #time
        self.rootgrpID.createDimension( self.dim_vars[1], 1 ) #depth
        self.rootgrpID.createDimension( self.dim_vars[2], 1 ) #lat
        self.rootgrpID.createDimension( self.dim_vars[3], 1 ) #lon
        
        
    def variable_init(self, nchandle, udunits_time_str='days since 1900-1-1' ):
        """
        built from knowledge about previous file
        """
        
        #build record variable attributes
        rec_vars, rec_var_name, rec_var_longname = [], [], []
        rec_var_generic_name, rec_var_FORTRAN, rec_var_units, rec_var_epic = [], [], [], []
        
        for v_name in nchandle.variables.keys():
            print(v_name)
            if not v_name in ['time','time2','depth','lat','lon','latitude','longitude']:
                print("Copying attributes for {0}".format(v_name))
                rec_vars.append( v_name )
                rec_var_name.append( nchandle.variables[v_name].name )
                rec_var_longname.append( nchandle.variables[v_name].long_name )
                rec_var_generic_name.append( nchandle.variables[v_name].generic_name )
                rec_var_units.append( nchandle.variables[v_name].units )
                rec_var_FORTRAN.append( nchandle.variables[v_name].FORTRAN_format )
                rec_var_epic.append( nchandle.variables[v_name].epic_code )

        
        rec_vars = ['time','depth','lat','lon'] + rec_vars

        rec_var_name = ['', '', '', ''] + rec_var_name
        rec_var_longname = ['', '', '', ''] + rec_var_longname
        rec_var_generic_name = ['', '', '', ''] + rec_var_generic_name
        rec_var_FORTRAN = ['', '', '', ''] + rec_var_FORTRAN
        rec_var_units = [udunits_time_str,'dbar','degree_north','degree_west'] + rec_var_units
        rec_var_type= ['f8'] + ['f4' for spot in rec_vars[1:]]
        rec_var_strtype= ['EVEN', 'EVEN', 'EVEN', 'EVEN'] + ['' for spot in rec_vars[4:]]
        rec_epic_code = [624,1,500,501] + rec_var_epic
        
        var_class = []
        var_class.append(self.rootgrpID.createVariable(rec_vars[0], rec_var_type[0], self.dim_vars[0]))#time1
        var_class.append(self.rootgrpID.createVariable(rec_vars[1], rec_var_type[1], self.dim_vars[1]))#depth
        var_class.append(self.rootgrpID.createVariable(rec_vars[2], rec_var_type[2], self.dim_vars[2]))#lat
        var_class.append(self.rootgrpID.createVariable(rec_vars[3], rec_var_type[3], self.dim_vars[3]))#lon
        
        for i, v in enumerate(rec_vars[4:]):  #1D coordinate variables
            var_class.append(self.rootgrpID.createVariable(rec_vars[i+4], rec_var_type[i+4], self.dim_vars))
            
        ### add variable attributes
        for i, v in enumerate(var_class): #4dimensional for all vars
            print("Adding Variable {0}".format(v))
            v.setncattr('name',rec_var_name[i])
            v.long_name = rec_var_longname[i]
            v.generic_name = rec_var_generic_name[i]
            v.FORTRAN_format = rec_var_FORTRAN[i]
            v.units = rec_var_units[i]
            v.type = rec_var_strtype[i]
            v.epic_code = rec_epic_code[i]
            
            
        self.var_class = var_class
        self.rec_vars = rec_vars

        
    def add_coord_data(self, depth=None, latitude=None, longitude=None, time=None, CastLog=False):
        """ """
        self.var_class[0][:] = time
        self.var_class[1][:] = depth
        self.var_class[2][:] = latitude
        self.var_class[3][:] = longitude #PMEL standard direction

    def add_data(self, data=None):
        """ """
        
        for ind, varname in enumerate(data.keys()):
            if not varname in ['time','time2','lat','lon','depth','latitude','longitude']:
                di = self.rec_vars.index(varname)
                self.var_class[di][:] = data[varname][:]
        

        
    def add_history(self, new_history):
        """Adds timestamp (UTC time) and history to existing information"""
        self.rootgrpID.History = self.rootgrpID.History + ' ' + datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")\
                    + ' ' + new_history + '\n'
                    
    def close(self):
        self.rootgrpID.close()
