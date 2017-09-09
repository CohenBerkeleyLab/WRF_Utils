function [ sz ] = get_wrf_array_size( wrffile )
%GET_WRF_ARRAY_SIZE Gets unstaggered array sizes in WRF output file
%   SZ = GET_WRF_ARRAY_SIZE( WRFFILE ) Returns a 1x4 vector SZ which
%   contains the length of unstagger WRF arrays in the given file along the
%   west_east, south_north, bottom_top, and Time dimensions. These will be
%   the sizes for most arrays of physical variables; although some
%   variables will be staggered. This function is intended for use when you
%   want to get the size of the arrays without actually loading the arrays
%   themselves.

ncid = netcdf.open(wrffile, 'NC_NOWRITE');
cleanupObj = onCleanup(@() thisCleanup(ncid));

dim_we = netcdf.inqDimID(ncid, 'west_east');
dim_sn = netcdf.inqDimID(ncid, 'south_north');
dim_bt = netcdf.inqDimID(ncid, 'bottom_top');
dim_time = netcdf.inqDimID(ncid, 'Time');

[~,sz_we] = netcdf.inqDim(ncid, dim_we);
[~,sz_sn] = netcdf.inqDim(ncid, dim_sn);
[~,sz_bt] = netcdf.inqDim(ncid, dim_bt);
[~,sz_time] = netcdf.inqDim(ncid, dim_time);



sz = [sz_we, sz_sn, sz_bt, sz_time];

end

function thisCleanup(ncid)
% Ensure that the netCDF file ID is closed even if this function errors.
%fprintf('Closing ncid %d\n', ncid);
netcdf.close(ncid);
end