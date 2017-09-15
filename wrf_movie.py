from __future__ import print_function

from glob import glob
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt
from matplotlib import animation as manim
import netCDF4 as nc
import numpy as np
import os.path

def read_wrf_var(wrf_file, wrf_variable, wrf_level):
    wrf = nc.Dataset(wrf_file, 'r')
    try:
        wrf_val = wrf.variables[wrf_variable][0, wrf_level, :, :]
    finally:
        wrf.close()
    return wrf_val

def read_wrf_coordinates(wrf_file):
    wrf = nc.Dataset(wrf_file, 'r')
    try:
        lon = wrf.variables['XLONG'][0, :, :]
        lat = wrf.variables['XLAT'][0, :, :]
    finally:
        wrf.close()
    return lon, lat

def init_map(lon_lim, lat_lim, wrf_file_0, wrf_variable):
    wrf_lon, wrf_lat = read_wrf_coordinates(wrf_file_0)
    wrf_val_0 = read_wrf_var(wrf_file_0, wrf_variable, 0)

    wrf_map = Basemap(resolution='l', llcrnrlon=lon_lim[0], llcrnrlat=lat_lim[0], urcrnrlon=lon_lim[1], urcrnrlat=lat_lim[1])
    wrf_map.drawcoastlines()
    val_limits = [np.min(wrf_val_0), np.max(wrf_val_0)]
    pcolor_instance = wrf_map.pcolormesh(wrf_lon, wrf_lat, wrf_val_0)
    pcolor_instance.set_clim(val_limits[0], val_limits[1])
    return pcolor_instance

def update_map(pcolor_instance, wrf_file, wrf_variable):
    wrf_val = read_wrf_var(wrf_file, wrf_variable, 0)
    # According to https://stackoverflow.com/questions/29009743/using-set-array-with-pyplot-pcolormesh-ruins-figure, the
    # initial call to pcolormesh basically expects one more point in each dimension than does set_array()
    pcolor_instance.set_array(wrf_val[:-1, :-1].ravel())

def write_movie(movie_file, wrf_files, wrf_variable):
    writer = manim.FFMpegWriter(fps=8)

    # If running in IPython, sometimes previous colorbars linger
    plt.close()

    val_pcolor = init_map([-125.0, -65.0], [25.0, 50.0], wrf_files[0], wrf_variable)
    fig = plt.gcf()
    plt.savefig('init.png')
    plt.set_cmap('jet')
    #plt.colorbar()

    with writer.saving(fig, movie_file, 100):
        for f_wrf in wrf_files:
            print('Writing frame {} of {} ({})'.format(wrf_files.index(f_wrf)+1, len(wrf_files), f_wrf))
            update_map(val_pcolor, f_wrf, wrf_variable)
            #savename = 'Frames/frame{:03d}.png'.format(wrf_files.index(f_wrf))
            #plt.savefig(savename)
            writer.grab_frame()

def test_me():
    wrf_path = '/Volumes/share-wrf2/Outputs/2005/01/'
    wrf_files = sorted(glob(os.path.join(wrf_path, 'wrfout_d01_2005-01-*')))
    write_movie('test.mp4', wrf_files, 'no2')

if __name__ == '__main__':
    test_me()
