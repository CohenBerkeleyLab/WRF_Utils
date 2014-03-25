%function [Data, ann_avgs, std_devs] = modis_map_aod(in_range_files, year)
%Average MODIS AOD
in_range_files = useful_files; year = 2013;
%lats = [35 35.5;36.5 37;36.2 36.5];
%lons = [-119.4 -118.8;-119.5 -119.9;-119.1 -119.5];
city_names = {'Bakersfield','Fresno','Visalia'}; city_latlon = [35 -119; 36.75 -119.1; 36.35 -119.3];
lats = [34 38]; lons = [-117 -122]; %city_names={'all three'};
directory = ['/Volumes/share/GROUP/SAT/MODIS/MOD04_L2/Cont_US_',num2str(year)];

olddir = cd(directory);
n = length(in_range_files);

aod_data = struct('lats',0,'lons',0,'aod',0);
bin_lats = (min(lats):0.1:max(lats))'; bin_lons = min(lons):0.1:max(lons); %Assume that the aerosol product is gridded at 0.1 degrees
bin_aods = zeros(length(bin_lons), length(bin_lats), length(in_range_files));

for i=1:n %Read in each MODIS files found to be applicable
    filename = in_range_files{i};
    hdfi = hdfinfo(filename);
    modis_lon = hdfread(hdfi.Vgroup(1).Vgroup(1).SDS(1)); %Load the latitude, longitude, and best-quality aerosol optical depth data
    modis_lat = hdfread(hdfi.Vgroup(1).Vgroup(1).SDS(2));
    modis_aod = hdfread(hdfi.Vgroup(1).Vgroup(2).SDS(8));
    
    
    lat_bdy = [min(lats) max(lats)];
    lon_bdy = [min(lons) max(lons)];
    
    %Remove all rows without at least one pixel inside our area of interest 
    [row, ~] = find(min(modis_lat,[],2) <= lat_bdy(2) & max(modis_lat,[],2) >= lat_bdy(1));
    modis_lat = modis_lat(row,:); modis_lon = modis_lon(row,:); modis_aod = modis_aod(row,:);
    
    %Remove all columns without at least one pixel inside our area of interest
    [~, col] = find(min(modis_lon,[],1) <= lon_bdy(2) & max(modis_lon,[],1) >= lon_bdy(1));
    modis_lat = modis_lat(:,col); modis_lon = modis_lon(:,col); modis_aod = modis_aod(:,col);
    aod = double(modis_aod); aod(aod<-1000) = NaN;
    %aod_data(i).lats = modis_lat; aod_data(i).lons = modis_lon; aod_data(i).aod = modis_aod;
    %MODIS AOD data seems to enter a -9999 if the data is invalid
    %Further, under 'Attributes' in the original hdf file, the scale factor is noted to be 1e-3
    
    for y=1:(length(bin_lats)-1)
        yes = [bin_lats(y), bin_lats(y), bin_lats(y+1), bin_lats(y+1), bin_lats(y)];
        for x=1:(length(bin_lons)-1)
            xes = [bin_lons(x), bin_lons(x+1), bin_lons(x+1), bin_lons(x), bin_lons(x)];
            bin = inpolygon(modis_lon, modis_lat, xes, yes);
            bin_aods(x,y,i) = nanmean(aod(bin));
        end
    end
    
    
    %clear modis_lon modis_lat modis_aod hdfi
end
bin_lats_mat = repmat(bin_lats,1,length(bin_lons)); bin_lons_mat = repmat(bin_lons,length(bin_lats),1);
mean_aods = nanmean(bin_aods,3);
mean_aods(mean_aods==0) = NaN;

close all
m_proj('Albers Equal-Area Conic','lon',[min(lons), max(lons)], 'lat', [min(lats), max(lats)]);
ax = m_pcolor(bin_lons_mat, bin_lats_mat, mean_aods');
shading flat
alpha(ax,0.5);
m_grid;

for j = 1:length(city_names)
    m_line(city_latlon(j,2),city_latlon(j,1)+.1,'marker','.','markersize',8,'color',[1 0 1],'linestyle','none');
    m_text(city_latlon(j,2),city_latlon(j,1)+.1,['  ',city_names{j}],'vertical','top','BackgroundColor',[.7 .7 .7]);
end

%m_proj('Albers Equal-Area Conic','lon',[min(lons) max(lons)],'lat',[min(lats) max(lats)]);
%figure;
%m_pcolor(modis_lon, modis_lat, modis_aod, 'shading','flat');

cd(olddir);
%end
