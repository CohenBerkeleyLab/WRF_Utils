% Function to find modis AOD files that pass through the area of interest.
% In particular, want to find ones that are appoximately centered on the
% area of interest (i.e. the center of the swath goes through the center
% 1/4th or so of the area) to avoid pixel distortions.  This is intended only to
% give a modis data set viable for qualitative comparison of AOD values.

% Josh Laughner <joshlaugh5@gmail.com> 20 Feb 2014

%%%%%%%%%%%%%%%%%%%%%%%%
latmin = 34; latmax = 38;
lonmin = -122; lonmax = -117;
center_fraction = .3; %Change this to a smaller value to require the center of the modis swath to pass through a smaller "corridor"
count_criterion = 1; %The minimum number of swath centers that have to fall within the center corridor.
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
file = 'MOD04_L2.A*.hdf';
directory = '/Volumes/share/GROUP/SAT/MODIS/MOD04_L2/Cont_US_2013'; %This will need changed if your files are elsewhere, clearly.
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Calculate the "corridor" if center_fraction not equal to 1 %%%%%
if center_fraction ~= 1
    delta_lat = abs(latmax - latmin); lat_median = (latmax + latmin)/2;
    lat_lower_bdy = lat_median - 0.5*center_fraction*delta_lat; lat_upper_bdy = lat_median + 0.5*center_fraction*delta_lat;
    
    delta_lon = abs(lonmax - lonmin); lon_median = (lonmax + lonmin)/2;
    lon_lower_bdy = lon_median - 0.5*center_fraction*delta_lon; lon_upper_bdy = lon_median + 0.5*center_fraction*delta_lat;
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
    if mod(a,10)==0; fprintf('Checking file %u of %u', a, n);disp(' '); end
    filename = files(a).name;
    hdfi = hdfinfo(filename);
    modis_lon = hdfread(hdfi.Vgroup(1).Vgroup(1).SDS(1));
    modis_lat = hdfread(hdfi.Vgroup(1).Vgroup(1).SDS(2));
    
    %MODIS files are (usually) laid out such that each row of an SDS is a swath.  So
    %we want to look at the middle column of the SDS's.
    
    swath_width = size(modis_lat, 2); swath_center = fix(swath_width/2);
    
    %Creates lists of points that define the polygon we want to check if
    %our swath falls inside.
    test_lons = [lon_upper_bdy, lon_upper_bdy, lon_lower_bdy, lon_lower_bdy, lon_upper_bdy];
    test_lats = [lat_lower_bdy, lat_upper_bdy, lat_upper_bdy, lat_lower_bdy, lat_lower_bdy];
    
    swaths_in = inpolygon(modis_lon(:,swath_center), modis_lat(:,swath_center), test_lons, test_lats);
    num_swaths_in = sum(swaths_in);
    
    if num_swaths_in >= count_criterion
        disp('Saving filename...');
        useful_files{cell_count,1} = filename; %#ok<SAGROW>
        cell_count = cell_count+1;
    end
    
end
