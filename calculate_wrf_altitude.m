function [ wrf_z ] = calculate_wrf_altitude( wrf_filename )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

geopot_base_height = ncread(wrf_filename, 'PHB');
geopot_pert_height = ncread(wrf_filename, 'PH');
% Geopotential height is height above sea level times gravitational
% acceleration, it relates to potential energy available
wrf_z = (geopot_base_height + geopot_pert_height)/9.81;

end

