function [ wrf_z ] = calculate_wrf_altitude( wrf_filename, varargin )
%CALCULATE_WRF_ALTITUDE Calculate grid cell altitude from geopotential height
%   WRF_Z = CALCULATE_WRF_ALTITUDE( WRF_FILENAME ) Reads variables PH and
%   PHB from WRF_FILENAME and calculates elevation above sea level from
%   them.
%
%   WRF_Z = CALCULATE_WRF_ALTITUDE( ___, START, COUNT )
%   WRF_Z = CALCULATE_WRF_ALTITUDE( ___, START, COUNT, STRIDE ) 
%       Either syntax will pass START, COUNT, and optionally STRIDE to the
%       calls to NCREAD(), allowing you to read a subset of the netCDF
%       data.

geopot_base_height = ncread(wrf_filename, 'PHB', varargin{:});
geopot_pert_height = ncread(wrf_filename, 'PH', varargin{:});
% Geopotential height is height above sea level times gravitational
% acceleration, it relates to potential energy available
wrf_z = (geopot_base_height + geopot_pert_height)/9.81;

end

