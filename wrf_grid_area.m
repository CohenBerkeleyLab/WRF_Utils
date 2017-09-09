function [ grid_area ] = wrf_grid_area( xlon, xlat )
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


end

