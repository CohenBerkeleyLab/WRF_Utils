function [Data, ann_avgs, std_devs] = modis_avg_aod(in_range_files, year)
%Average MODIS AOD

lats = [35 35.5;36.5 37;36.2 36.5];
lons = [-119.4 -118.8;-119.5 -119.9;-119.1 -119.5];
city_names = {'Bakersfield','Fresno','Visalia'};
%lats = [35 37]; lons = [-119 -120]; city_names={'all three'};
directory = ['/Volumes/share/GROUP/SAT/MODIS/MOD04_L2/Cont_US_',num2str(year)];

olddir = cd(directory);
n = length(in_range_files);



for i=1:n %Read in each MODIS files found to be applicable
    filename = in_range_files{i};
    hdfi = hdfinfo(filename);
    modis_lon = hdfread(hdfi.Vgroup(1).Vgroup(1).SDS(1)); %Load the latitude, longitude, and best-quality aerosol optical depth data
    modis_lat = hdfread(hdfi.Vgroup(1).Vgroup(1).SDS(2));
    modis_aod = hdfread(hdfi.Vgroup(1).Vgroup(2).SDS(8));
    
    for c = 1:size(lats,1) %For each city, find the applicable data and save it to the structure "Data"
        lat_bdy = [lats(c,1), lats(c,1), lats(c,2), lats(c,2), lats(c,1)];
        lon_bdy = [lons(c,1), lons(c,2), lons(c,2), lons(c,1), lons(c,1)];
        xx = inpolygon(modis_lon, modis_lat,lon_bdy,lat_bdy);
        
        pixels_in_city = modis_aod(xx); %MODIS AOD data seems to enter a -9999 if the data is invalid; this finds the pixels inside the city limit and deletes those that = -9999
        useful_pixels = pixels_in_city(pixels_in_city~=-9999); %Further, under 'Attributes' in the original hdf file, the scale factor is noted to be 1e-3
        city_aod = mean(useful_pixels)*1e-3;
        
        Data(i).city(c).city_name = city_names{c};
        Data(i).city(c).year = year;
        Data(i).city(c).mean_aod = city_aod;
        Data(i).city(c).city_lats = modis_lat(xx);
        Data(i).city(c).city_lons = modis_lon(xx);
        
    end
    %clear modis_lon modis_lat modis_aod hdfi
end



ann_avgs = zeros(size(lats,1),1);
std_devs = zeros(size(lats,1),1);
for c = 1:size(lats,1) %Find the average AOD for each city.
    city_avg_vec = zeros(length(Data),1);
    for j=1:length(Data)
        if ~isnan(Data(j).city(c).mean_aod)
            city_avg_vec(j) = Data(j).city(c).mean_aod;
        end
    end
    ann_avgs(c) = mean(city_avg_vec);
    std_devs(c) = std(city_avg_vec);
end

cd(olddir);
%Returns Data, ann_avgs, and std_devs
%end