%%add2grid_5km_new
%%arr 11/19/2009

Dimensions=size(Data(d).ColumnAmountNO2);
swath=d;

x=1:1:Dimensions(1)*Dimensions(2);
y=1;

Lon1=Data(d).Loncorn(1,x)';   Lat1=Data(d).Latcorn(1,x)';
Lon2=Data(d).Loncorn(2,x)';   Lat2=Data(d).Latcorn(2,x)';
Lon3=Data(d).Loncorn(3,x)';   Lat3=Data(d).Latcorn(3,x)';
Lon4=Data(d).Loncorn(4,x)';   Lat4=Data(d).Latcorn(4,x)';
Lon5=Data(d).Longitude(x)';   Lat5=Data(d).Latitude(x)';

CoordLon=cat(3,Lon1,Lon2,Lon3,Lon4,Lon5);
CoordLat=cat(3,Lat1,Lat2,Lat3,Lat4,Lat5);

reslat=resolution;   resy=resolution2;
reslon=resolution2;   resx=resolution;
lon1=lonmin; lon2=lonmax;  maxy=(abs(lonmin-lonmax))/resy; miny=1; maxy=single(maxy);
lat1=latmin; lat2=latmax;  maxx=(abs(latmax-latmin))/resx; minx=1; maxx=single(maxx);

lCoordLon=zeros(Dimensions(1)*Dimensions(2),y,5);
lCoordLat=zeros(Dimensions(1)*Dimensions(2),y,5); %JLL 2-14-2014: This line changed from "lCoordLat = zeros(Dimensions(1)*Dimensions,y,5);" as *Dimensions is not a scalar

for x=1:1:Dimensions(1)*Dimensions(2);
    for c=1:5;
    lCoordLon(x,y,c)=(CoordLon(x,y,c)-lon1)/reslon;
    lCoordLat(x,y,c)=(CoordLat(x,y,c)-lat1)/reslat;
    end
end

hdf_quadrangle_5km_new
