
#+---------------------------------------------------------------------
#+ Python Script that creates CMAQ_ADJ v4.5 forcing files
#+ Check the 'CHANGE' comments to define your directories for the run
#+ Author: Camilo Moreno
#+ Email = cama9709@gmail.com
#+---------------------------------------------------------------------
#%%
import numpy as np
import glob 
from netCDF4 import Dataset
from datetime import datetime, date, timedelta
from os import listdir, scandir, getcwd

def create_forced_uniarray(coords, dx, dy, lons, lats):
    """Create array where cells that wanted to be forced are equal to 1
    and the rest are equal to 0.

    Keyword arguments:
    coords -- list of latlon coordinates of lower left and upper right corners
    dx -- float representing halft the width of a grid cell
    dy -- float representing halft the height of a grid cell
    lons -- array containing the cell-centered longitud values of the grid
    lats -- array containing the cell-centered latitude values of the grid
    """
    lat_min = coords[0] - dy
    lon_min = coords[1] - dx
    lat_max = coords[2] + dy
    lon_max = coords[3] + dx

    forced_uniarr = lons

    forced_uniarr[lons < lon_min] = 0
    forced_uniarr[lons > lon_max] = 0
    forced_uniarr[lats < lat_min] = 0
    forced_uniarr[lats > lat_max] = 0
    
    forced_uniarr[forced_uniarr != 0] = 1

    return forced_uniarr

def create_forced_var(ds_latlon, dic_coords):
    """Create array where cells that wanted to be forced are equal to 1
    and the rest are equal to 0 for multiple locations.

    Keyword arguments:
    ds_latlon -- Dataset of the latlon netCDF file of the grid
    dic_coords -- dictionary of latlon coordinates for lower left and 
                  upper right corners of each specified kay-value location
    """
    lats = ds_latlon.variables['LAT'][0][0][:][:]
    lons = ds_latlon.variables['LON'][0][0][:][:]
    lats = np.ma.getdata(lats)
    lons = np.ma.getdata(lons)

    meters_per_degree = 110574.61
    dx = ds_latlon.XCELL/meters_per_degree/2
    dy = ds_latlon.YCELL/meters_per_degree/2

    ltime = 25
    llay = 1
    lrow, lcol = lats.shape

    forced_arr = np.zeros([ltime, llay, lrow, lcol])

    for coords in dic_coords.values():
        forced_uniarr = create_forced_uniarray(list(coords), dx, dy, lons, lats)
        forced_uniarr = np.resize(forced_uniarr, [ltime, llay, lrow, lcol])

        forced_arr = forced_arr + forced_uniarr

    return forced_arr

def monthly_date(day):
    """convert date from YYYYDDD to YYYYMMDD.

    Keyword arguments:
    day -- int of the day in format YYYYDDD
    """
    day_str = str(day)
    year = int(day_str[0:4])
    day_y = int(day_str[4:])

    date = datetime(year, 1, 1) + timedelta(day_y - 1)
    daym = date.timetuple().tm_mday
    str_daym = '0'*(2 - len(str(daym))) + str(daym)
    month = date.timetuple().tm_mon
    str_month = '0'*(2 - len(str(month))) + str(month)

    return str(year) + str_month + str_daym

def create_ncfile(save_dir, ds_latlon, spc_name, forced_arr, ds_conc):
    """Create Final NETCDF file.

    Keyword arguments:
    save_dir -- string of the location for saving the netCDF files
    ds_latlon -- Dataset of the latlon netCDF file of the grid
    spc_name -- string of the forced species name
    forced_arr -- array containing the forced values by location for the species
    ds_conc -- conc file dataset
    """
    hr = 0
    num_vars = 1
    lays = 1
    ltime = 25
    cols = len(ds_latlon.dimensions['COL'])
    rows = len(ds_latlon.dimensions['ROW'])
    datetimes = len(ds_latlon.dimensions['DATE-TIME'])
    day = ds_conc.SDATE

    #* Create new netCDF
    day_monthly = monthly_date(day)
    new_cmaq_file = f'{save_dir}/ADJ_FORCE.{day_monthly}'
    ds_new_cmaq = Dataset(new_cmaq_file, open = True, mode = 'w', format=  "NETCDF3_64BIT")

    #* Create dimenssions
    TSTEP = ds_new_cmaq.createDimension("TSTEP", None)
    DATE_TIME = ds_new_cmaq.createDimension("DATE-TIME", datetimes)
    LAY = ds_new_cmaq.createDimension("LAY", lays)
    VAR = ds_new_cmaq.createDimension("VAR", num_vars)
    ROW = ds_new_cmaq.createDimension("ROW", rows)
    COL = ds_new_cmaq.createDimension("COL", cols)

    ds_new_cmaq.sync()

    #* Creatae attributes
    attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", 
           "SDATE", "STIME", "TSTEP", "NTHIK", "NCOLS", "NROWS", "GDTYP", "P_ALP", 
           "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG", "YORIG", "XCELL", "YCELL", 
           "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]
    
    for attr in attrs:
        if hasattr(ds_conc, attr): 
            attrVal = getattr(ds_conc, attr)
            setattr(ds_new_cmaq, attr, attrVal)

    varlist = spc_name + ' '*(16 - len(spc_name))
    cmaq_attrs = {'NLAYS': np.int32(lays),
                  'NVARS': np.int32(num_vars),
                  'UPNAM': "RD_FORCE_FILE",
                  'VAR-LIST': varlist,
                  'FILEDESC': "Adjoint forcing file. Forcing file of specified geolocations"
                  }

    for attr in cmaq_attrs:
        ds_new_cmaq.setncattr(attr, cmaq_attrs[attr])

    #* Create variables
    tflag = ds_new_cmaq.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    fill_attrs(ds_conc, tflag)

    var_temp = ds_new_cmaq.createVariable(spc_name,"f4",("TSTEP", "LAY", "ROW", "COL"))
    fill_attrs(ds_conc, var_temp)

    #* Fill variables
    ds_new_cmaq.variables[spc_name][:, :, :] = forced_arr
    dattim = np.squeeze(ds_conc.variables['TFLAG'][:][:])
    dattim = dattim[0:ltime,1,:]
    tflag[:] = np.zeros([ltime,lays,datetimes])
    tflag[:] = np.reshape(dattim,[ltime,lays,datetimes])  
    
    #* Close new netcdf file
    ds_new_cmaq.close()

    print(f"{day_monthly} Forcing file DONE")

def fill_attrs(ds_conc, nc_var):
    """Fill atribute values for variables as appear in conc file.

    Keyword arguments:
    ds_conc -- conc file dataset
    nc_var -- netCDF variable of the file
    var_name -- variable name
    """
    varattrs=["long_name","units","var_desc"]
    for varattr in varattrs:
        if hasattr(ds_conc.variables[nc_var.name], varattr): 
            varattrVal = getattr(ds_conc.variables[nc_var.name], varattr)
            setattr(nc_var, varattr, varattrVal)

def get_concfiels(conc_file_dir):
    """Get all files that begins whith 'CONC' from the given directory.

    Keyword arguments:
    conc_file_dir -- string of the directory whre conc day files
    """
    all_files = [f for f in glob.glob(f'{conc_file_dir}/CONC.*')]
    all_files.sort()
    return all_files

# %%
if __name__ == "__main__":
    #CHANGE: save_dir: directory path where the forced files will be saved
    #        latlon_file: path of the latlon netCDF of the run grid
    #        conc_file_dir: path where CONC files are located
    #        spc_name: cbo5_aero5_aq name of the one species to be forced
    #        dic_coord: dictionary containing a list of latlon coords of lower left 
    #                   and upper right corners for each key-value location.
    #                   Format fo dic_coord must be:
    #                   dic_coord = {
    #                               'name_of_location_1': [lower_left_lat_1, lower_left_lon_1, upper_left_lat_1, upper_left_lon_1]
    #                               'name_of_location_2': [lower_left_lat_2, lower_left_lon_2, upper_left_lat_2, upper_left_lon_2]
    #                                 }
    save_dir = '/Volumes/Avispa/ADJ_FORCE_files'
    latlon_file = '/Volumes/Avispa/latlon.conc'
    conc_file_dir = '//Volumes/Avispa/Conc_files'
    spc_name = 'ASO4I'
    dic_coords = {
                  'Bogota' : [4.461864, -74.223421 , 4.833805, -74.007853]
                 }
    
    ds_latlon = Dataset(latlon_file, mode = 'r',  open = True)
    all_files = get_concfiels(conc_file_dir)  
    for file in all_files:
        ds_conc = Dataset(file, mode = 'r',  open = True)
        forced_arr = create_forced_var(ds_latlon, dic_coords)
        create_ncfile(save_dir, ds_latlon, spc_name, forced_arr, ds_conc)
        ds_conc.close()

    ds_latlon.close()
    
# %%

