
#+---------------------------------------------------------------------
#+ Python Script that creates CMAQ_ADJ v4.5 forcing files
#+ Check the 'CHANGE' comments to define your directories for the run
#+ Author: Camilo Moreno
#+ Email = cama9709@gmail.com
#+---------------------------------------------------------------------
#%%
import numpy as np
from netCDF4 import Dataset

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

def create_tflag_arr(day, hr): 
    """Create the TFLAG array value for one day to fill the netCDF 
    file variable.

    Keyword arguments:
    day -- int of the day in format YYYYDDD
    hr -- initial hour of the day. May be 0
    """
    for i in range(0,25):
        if hr == 0 and i == 0:
            tflag_arr = np.array([np.tile([day, hr], (1, 1))])
        else:
            tflag_temp = np.array([np.tile([day, hr], (1, 1))])
            tflag_arr = np.concatenate((tflag_arr, tflag_temp), axis = 0)
        day = day + int((hr / 230000))
        hr = (hr + 10000) % 240000

    return tflag_arr

def create_ncfile(save_dir, day, ds_latlon, spc_name, units, forced_arr):
    """Create Final NETCDF file.

    Keyword arguments:
    save_dir -- string of the location for saving the netCDF files
    day -- int of the day in format YYYYDDD
    ds_latlon -- Dataset of the latlon netCDF file of the grid
    spc_name -- string of the forced species name
    units -- string of the units for the force species
    forced_arr -- array containing the forced values by location for the species
    """
    hr = 0
    num_vars = 1
    lays = 1
    cols = len(ds_latlon.dimensions['COL'])
    rows = len(ds_latlon.dimensions['ROW'])
    datetimes = len(ds_latlon.dimensions['DATE-TIME'])

    #* Create new netCDF
    new_cmaq_file = f'{save_dir}/ADJ_FORCE.{day}'
    ds_new_cmaq = Dataset(new_cmaq_file, open = True, mode = 'w', format=  "NETCDF3_64BIT")

    #* Create dimenssions
    TSTEP = ds_new_cmaq.createDimension("TSTEP", None)
    DATE_TIME = ds_new_cmaq.createDimension("DATE-TIME", 2)
    LAY = ds_new_cmaq.createDimension("LAY", lays)
    VAR = ds_new_cmaq.createDimension("VAR", num_vars)
    ROW = ds_new_cmaq.createDimension("ROW", rows)
    COL = ds_new_cmaq.createDimension("COL", cols)

    #* Create variables
    tflag = ds_new_cmaq.createVariable("TFLAG","i4",("TSTEP","VAR", "DATE-TIME"))
    tflag.units = '<YYYYDDD,HHMMSS>'
    tflag.long_name = 'TFLAG'
    tflag.var_desc = 'Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS'

    var_temp = ds_new_cmaq.createVariable(spc_name,"f4",("TSTEP", "LAY", "ROW", "COL"))
    var_temp.units = units
    var_temp.long_name = spc_name
    var_temp.var_desc = f'Forced species {spc_name}'

    #* Fill variables
    ds_new_cmaq.variables[spc_name][:, :, :] = forced_arr
    ds_new_cmaq.variables['TFLAG'][:, :, :] = create_tflag_arr(day, hr)

    #* Creatae attributes
    varlist = spc_name + ' '*(16 - len(spc_name))
    cmaq_attrs = {'IOAPI_VERSION': '$Id: @(#) ioapi library version 3.1 $',
                'EXEC_ID': '????????????????',
                'FTYPE': np.int32(1),
                'CDATE': np.int32(2021100),
                'CTIME': np.int32(185404),
                'WDATE': np.int32(2021100),
                'WTIME': np.int32(185404),
                'SDATE': np.int32(day),
                'STIME': np.int32(hr),
                'TSTEP': np.int32(10000),
                'NTHIK': np.int32(1),
                'NCOLS': np.int32(cols),
                'NROWS': np.int32(rows),
                'NLAYS': np.int32(lays),
                'NVARS': np.int32(num_vars),
                'GDTYP': np.int32(ds_latlon.GDTYP), # 2 is Lambert Conformal Conic
                'P_ALP': np.float64(ds_latlon.P_ALP),
                'P_BET': np.float64(ds_latlon.P_BET),
                'P_GAM': np.float64(ds_latlon.P_GAM),
                'XCENT': np.float64(ds_latlon.XCENT),
                'YCENT': np.float64(ds_latlon.YCENT),
                'XORIG': np.float64(ds_latlon.XORIG),
                'YORIG': np.float64(ds_latlon.YORIG),
                'XCELL': np.float64(ds_latlon.XCELL),
                'YCELL': np.float64(ds_latlon.YCELL),
                'VGTYP': np.int32(ds_latlon.VGTYP),
                'VGTOP': np.float32(0.0),
                'VGLVLS': np.array([np.float32(0.), np.float32(0.)]),
                'GDNAM': ds_latlon.GDNAM,
                'UPNAM': "M3WNDW",
                'VAR-LIST': varlist,
                'FILEDESC': "Forcing file of specified geolocations/ Camilo Moreno",
                'HISTORY': 'La historia es historia',}

    for attr in cmaq_attrs:
        ds_new_cmaq.setncattr(attr, cmaq_attrs[attr])

    #* Close new netcdf file
    ds_new_cmaq.close()

    (f"{day} Forcing file DONE")

# %%
if __name__ == "__main__":
    #CHANGE: save_dir: directory path where the forced files will be saved
    #        latlon_file: path of the latlon netCDF of the run grid
    #        spc_name: cbo5_aero5_aq name of the one species to be forced
    #        units: the units of the one species to be forced
    #        day_0: first day of forcing files (one file per day will be created)
    #        day_end: last day of forcing files
    #        dic_coord: dictionary containing a list of latlon coords of lower left 
    #                   and upper right corners for each key-value location.
    #                   Format fo dic_coord must be:
    #                   dic_coord = {
    #                               'name_of_location_1': [lower_left_lat_1, lower_left_lon_1, upper_left_lat_1, upper_left_lon_1]
    #                               'name_of_location_2': [lower_left_lat_2, lower_left_lon_2, upper_left_lat_2, upper_left_lon_2]
    #                                 }
    save_dir = '/Volumes/Avispa'
    latlon_file = '/Volumes/Avispa/Emissions/CMAQ_emis/latlon.emis'
    spc_name = 'ASOJK'
    units = 'g/s'

    day_0 = 2018032
    day_end = 2018033

    dic_coords = {
                  'Bogota' : [4.461864, -74.223421 , 4.833805, -74.007853]
                }
    
    ds_latlon = Dataset(latlon_file, mode = 'r',  open = True)

    for day in range(day_0, day_end + 1):
        forced_arr = create_forced_var(ds_latlon, dic_coords)
        create_ncfile(save_dir, day, ds_latlon, spc_name, units, forced_arr)

    ds_latlon.close()
# %%
