% Function to find modis AOD files that pass through the area of interest.
% In particular, want to find ones that are appoximately centered on the
% area of interest (i.e. the center of the swath goes through the center
% 1/4th or so of the area) to avoid pixel distortions.  This is intended only to
% give a modis data set viable for qualitative comparison of AOD values.

% Josh Laughner <joshlaugh5@gmail.com>

%%%%%%%%%%%%%%%%%%%%%%%%
latmin = 34; latmax = 38;
lonmin = -122; lonmax = -117;
center_fraction = 1; %Change this to a smaller value to require the center of the modis swath to pass through a smaller "corridor"
check_lat = 1; %Set to 0 if you only want to check if the middle of the swath is within the latitudes (i.e. if the swath path is primarily E-W)
check_lon = 1;
lat_criterion = 0.25; %The fraction of the center swath we want to be inside our box (range from 0 to 1)
lon_criterion = 0.001;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
file = 'MOD04_L2.A*.hdf';
directory = '/Volumes/share/GROUP/SAT/MODIS/MOD04_L2/Cont_US_Jan2014/500805565'; %This will need changed if your files are elsewhere, clearly.
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Some boolean variables we will need %%%%%
lats_in = 0; lons_in = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Calculate the "corridor" if center_fraction not equal to 1 %%%%%
if center_fraction ~= 1
    delta_lat = abs(latmax - latmin); lat_median = (latmax + latmin)/2;
    lat_lower_bdy = lat_median - 0.5*delta_lat; lat_upper_bdy = lat_median + 0.5*delta_lat;
    
    delta_lon = abs(lonmax - lonmin); lon_median = (lonmax + lonmin)/2;
    lon_lower_bdy = lon_median - 0.5*delta_lon; lon_upper_bdy = lon_median + 0.5*delta_lat;
else
    lat_lower_bdy = latmin; lat_upper_bdy = latmax;
    lon_lower_bdy = lonmin; lon_upper_bdy = lonmax;
end

cd(directory);
files = dir(fullfile(directory,file)); %Find all files matching the pattern "file" in the given directory
n = length(files);

useful_files = {}; %Create empty cell array to track the file names we want
cell_count = 1; %Counter to keep track of what row we should be adding to
for a=1:n
    fprintf('Checking file %u of %u', a, n);disp(' ');
    filename = files(a).name;
    hdfi = hdfinfo(filename);
    modis_lon = hdfread(hdfi.Vgroup(1).Vgroup(1).SDS(1));
    modis_lat = hdfread(hdfi.Vgroup(1).Vgroup(1).SDS(2));
    
    %MODIS files are (usually) laid out such that each row of an SDS is a swath.  So
    %we want to look at the middle column of the SDS's.
    
    swath_width = size(modis_lat, 1); swath_center = fix(swath_width/2);
    rows_in_box = 0; %A counter.
    if check_lat == 1
        swath_length = size(modis_lat, 2);
        for b=1:swath_length
            if modis_lat(b, swath_center) > lat_lower_bdy && modis_lat(b , swath_center) < lat_upper_bdy
                rows_in_box = rows_in_box + 1;
                if rows_in_box > fix(lat_criterion*swath_length); 
                    lats_in = 1;
                    %disp('Found latitudes in the box');
                    break
                end
            end
        end
    end
    
    if check_lon == 1 && lats_in == 1 %Now check the longitudes
        swath_length = size(modis_lon,2);
        for b=1:swath_length
            if modis_lon(b, swath_center) > lon_lower_bdy && modis_lon(b, swath_center) < lon_upper_bdy
                rows_in_box = rows_in_box + 1;
                if rows_in_box > fix(lon_criterion*swath_length);
                    lons_in = 1;
                    %disp('Found longitudes in the box');
                    break
                end
            end
        end
    end
    
    if lats_in == 1 && lons_in == 1
        disp('Saving filename...');
        useful_files{cell_count,1} = filename; %#ok<SAGROW>
        cell_count = cell_count+1;
    end
    
    lats_in = 0; lons_in = 0;
end
