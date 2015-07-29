function OMI = add2grid_BEHR_winds(Data,OMI,reslat,reslon,lonbdy,latbdy)

% add2grid_BEHR Updated version of add2grid_5km_2014 intended to be a
% usable for gridding BEHR data.
%
%   This is an intermediate function between the calling script and
%   hdf_quadrangle_general that prepares latitude and longitude corner
%   data. It requires 6 inputs:
%       1) A single, top level element (i.e. a single satellite swath) of a
%       Data structure.
%       2) A single top level element of the OMI structure it will be
%       outputting - this keeps the fields the same and in the same order.
%       3) Resolution of the grid in the latitudinal direction
%       4) Resolution of the grid in the longitudinal direction
%       5) The longitude boundaries as a 1x2 vector, i.e. [min, max]
%       6) The latitude boundaries as a 1x2 vector, i.e. [min, max]
%
%   Note that in order for this to work, the Data structure must have
%   fields "Latitude" and "Longitude" defining the center lat/lon of the
%   pixels, as well as fields "Latcorn" and "Loncorn" defining the latitude
%   and longitude corners of the pixel.  These fields must be matrices
%   arranged such that the first dimension has length 4 and corresponds to
%   the corners of the pixels.  These can be 4 x number of pixels or 4 x
%   (x) x (y), where (x) and (y) are numbers of pixels in each direction,
%   but the first dimension must have length 4.



narginchk(6,6);

if numel(Data) > 1;
    error('add2grid:DataIn','Pass only one top-level element of Data to this function');
end

if ~isfield(Data,'Latcorn') || ~isfield(Data,'Loncorn');
    error('add2grid:DataIn','Data must contain fields "Latcorn" and "Loncorn"');
elseif size(Data.Latcorn,1) > 4 || size(Data.Loncorn,1) > 4
    error('add2grid:DataIn', 'Latcorn and loncorn must have the first dimension of length 4 represent the corners of the pixel')
end
    

Dimensions=size(Data.ColumnAmountNO2Trop);
%swath=d;

x=1:1:Dimensions(1)*Dimensions(2); 
y=1;

Lon1=Data.Loncorn(1,x);          Lat1=Data.Latcorn(1,x);
if isrow(Lon1); Lon1 = Lon1'; end;  if isrow(Lat1); Lat1 = Lat1'; end;

Lon2=Data.Loncorn(2,x);          Lat2=Data.Latcorn(2,x);
if isrow(Lon2); Lon2 = Lon2'; end;  if isrow(Lat2); Lat2 = Lat2'; end;

Lon3=Data.Loncorn(3,x);          Lat3=Data.Latcorn(3,x);
if isrow(Lon3); Lon3 = Lon3'; end;  if isrow(Lat3); Lat3 = Lat3'; end;

Lon4=Data.Loncorn(4,x);          Lat4=Data.Latcorn(4,x);
if isrow(Lon4); Lon4 = Lon4'; end;  if isrow(Lat4); Lat4 = Lat4'; end;

Lon5=Data.Longitude(x);          Lat5=Data.Latitude(x);
if isrow(Lon5); Lon5 = Lon5'; end;  if isrow(Lat5); Lat5 = Lat5'; end;

CoordLon=cat(3,Lon1,Lon2,Lon3,Lon4,Lon5); %JLL 18 Mar 2014: cat concatenates along the dimension specified as the first argument
CoordLat=cat(3,Lat1,Lat2,Lat3,Lat4,Lat5);

lon1=min(lonbdy); lon2=max(lonbdy);  maxy=(abs(lon1-lon2))/reslon; miny=1; maxy=single(maxy);
lat1=min(latbdy); lat2=max(latbdy);  maxx=(abs(lat1-lat2))/reslat; minx=1; maxx=single(maxx);

lCoordLon=zeros(Dimensions(1)*Dimensions(2),y,5);
lCoordLat=zeros(Dimensions(1)*Dimensions(2),y,5); %JLL 2-14-2014: This line changed from "lCoordLat = zeros(Dimensions(1)*Dimensions,y,5);" as *Dimensions is not a scalar

for x=1:1:Dimensions(1)*Dimensions(2);
    for c=1:5;
    lCoordLon(x,y,c)=(CoordLon(x,y,c)-lon1)/reslon;
    lCoordLat(x,y,c)=(CoordLat(x,y,c)-lat1)/reslat;
    end
end

OMI = hdf_quadrangle_OMI_winds(Data, OMI, maxx, minx, maxy, miny, lCoordLon, lCoordLat, Lon1, Lon2, Lon4, Lat1, Lat2, Lat4);
OMI.MapData.LatBdy = latbdy;
OMI.MapData.LatRes = reslat;
OMI.MapData.LonBdy = lonbdy;
OMI.MapData.LonRes = reslon;

Latitude=(min(latbdy)+reslat/2):reslat:(max(latbdy)-0.025); Latitude=Latitude'; nlat = numel(Latitude);
Longitude=(min(lonbdy)+reslon/2):reslon:(max(lonbdy)-reslon/2); nlon = numel(Longitude);
Latitude=repmat(Latitude,1,nlon); Longitude=repmat(Longitude,nlat,1);

OMI.Latitude = Latitude;
OMI.Longitude = Longitude;

end