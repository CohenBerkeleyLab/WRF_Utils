function [ grid_area ] = wrf_grid_area( xlon, xlat, varargin )
%WRF_GRID_AREA Estimates the area of each WRF grid cell
%   GRID_AREA = WRF_GRID_AREA( XLON, XLAT ) takes the grid center longitude
%   (XLON) and latitude (XLAT) and returns the areas of the grid cells in
%   kilometers (GRID_AREA). This uses wrf_grid_corners to estimate the
%   corner points, so it's not as accurate as it would be if it calculated
%   them with knowledge of the map projection. The area is computed
%   assuming the WGS84 ellipsoid.

E = JLLErrors;

%%%% INPUT CHECKING %%%%%
if ~ismatrix(xlon) || ~ismatrix(xlat)
    E.badinput('XLON and XLAT should be 2D')
elseif ~isequal(size(xlon), size(xlat))
    E.badinput('XLON and XLAT must be the same size');
end

write_nc = false;

if nargin > 2
    write_nc = true;
    ncfilename = varargin{1};
    if ~ischar(ncfilename)
        E.badinput('NCFILENAME must be a string');
    end
    
    if exist(ncfilename, 'file')
        if nargin < 4 || ~varargin{2}
            warning('%s exists, will not overwrite. Give true as the fourth argument to overwrite\n', ncfilename);
            write_nc = false;
        end
    end
end
%%%%% MAIN FUNCTION %%%%%

grid_area = nan(size(xlon));
earth_ellip = referenceEllipsoid('wgs84','kilometer');

[xloncorn, xlatcorn] = wrf_grid_corners(xlon, xlat);

for a = 1:numel(xlon)
    xall = xloncorn(:,a);
    yall = xlatcorn(:,a);
    [xall, yall] = poly2cw(xall, yall);
    grid_area(a) = areaint(yall,xall,earth_ellip);
end

if write_nc
    write_netcdf(ncfilename, xlon, xlat, grid_area);
end
end

function write_netcdf(file_name, wrf_lon, wrf_lat, wrf_area)
if exist(file_name, 'file')
    delete(file_name);
end

nccreate(file_name, 'XLONG', 'Datatype', 'single','Dimensions',{'west_east',size(wrf_lon,1),'south_north',size(wrf_lon,2)});
ncwrite(file_name, 'XLONG', wrf_lon);
ncwriteatt(file_name, 'XLONG', 'description', 'LONGITUDE, WEST IS NEGATIVE');
ncwriteatt(file_name, 'XLONG', 'units', 'degree_east');
ncwriteatt(file_name, 'XLONG', 'stagger', '');

nccreate(file_name, 'XLAT', 'Datatype', 'single','Dimensions',{'west_east',size(wrf_lon,1),'south_north',size(wrf_lon,2)});
ncwrite(file_name, 'XLAT', wrf_lat);
ncwriteatt(file_name, 'XLAT', 'description', 'LATITUDE, SOUTH IS NEGATIVE');
ncwriteatt(file_name, 'XLAT', 'units', 'degree_north');
ncwriteatt(file_name, 'XLAT', 'stagger', '');

nccreate(file_name, 'AREA', 'Datatype', 'single','Dimensions',{'west_east',size(wrf_lon,1),'south_north',size(wrf_lon,2)});
ncwrite(file_name, 'AREA', wrf_area);
ncwriteatt(file_name, 'AREA', 'description', 'GRID CELL AREA, CALCULATED BY WRF_GRID_AREA.m');
ncwriteatt(file_name, 'AREA', 'units', 'km^2');
ncwriteatt(file_name, 'AREA', 'stagger', '');
end