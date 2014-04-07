%function [Data, ann_avgs, std_devs] = modis_map_aod(in_range_files, year)
%Average MODIS AOD
in_range_files = useful_files; 
DEBUG_LEVEL = 0;
%%%%% Latitude and longitude boundaries of the mapped area. %%%%%
lats = [32 36]; lons = [-120 -115]; 

%%%%% Names and center coordinates of any cities you wish to mark on the map %%%%
city_names = {'Bakersfield','Fresno','Visalia'}; city_latlon = [35 -119; 36.75 -119.1; 36.35 -119.3];

%%%%% Time controls %%%%%
    % A cell array that contains all the date ranges that you wish to
    % average.  Each range should be a row in the cell array with dates in
    % 'mm/dd/yyyy' or 'dd-MMM-yyyy' format, where 'MMM' is the 3-letter
    % month abbreviation. Non-contiguous ranges should each be their own
    % row, i.e. to average over the summer months of 2012 and 2013, enter
    % {'01-Jun-2012', '31-Aug-2012'; '01-Jun-2013', '31-Aug-2013'}
date_ranges = {'05/01/2012','09/30/2012'; '05/01/2013','09/30/2013'};  
    % Allows for weekend effect maps.
    % Set to 1 for all days, 2 for weekdays only, 3 for weekend only
week_mode = 3; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       Begin main program - end input variables                %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert the input date ranges to date numbers that are easier to compare
date_range_nums = zeros(size(date_ranges));
for d=1:numel(date_ranges);
    date_range_nums(d) = datenum(date_ranges{d});
end

% The directory will automatically be filled in with the year from the 

olddir = pwd; %Save the current directory so we can go back there at the end
n = length(in_range_files);
ncount = 0; %A count of how many MODIS files were averaged

bin_lats = (min(lats):0.05:max(lats))'; bin_lons = min(lons):0.05:max(lons); %Create a 0.5 degree grid
bin_lats = repmat(bin_lats,1,size(bin_lons,2)); bin_lons = repmat(bin_lons,size(bin_lats,1),1);
bin_aods = zeros(size(bin_lons,1), size(bin_lats,2), length(in_range_files)); bin_aods(:) = NaN;
bin_areaweight = zeros(size(bin_lons,1), size(bin_lats,2), length(in_range_files)); bin_areaweight(:) = NaN;

for i=1:n %Read in each MODIS files found to be applicable
    filename = in_range_files{i};
    year = str2double(filename(11:14));
    weekdays = busdays(sprintf('01-Jan-%u',year),sprintf('31-Dec-%u',year),1,NaN); %Find all weekdays in the year
    fprintf('Checking %s\n',filename);
    
    modis_day = datenum(modis_day_to_date(str2double(filename(15:17)),year));
    if DEBUG_LEVEL>0; fprintf('  MODIS num = %u; Ranges = [%u, %u] and [%u %u] \n', modis_day, date_range_nums(1,1), date_range_nums(1,2), date_range_nums(2,1), date_range_nums(2,2)); end
    if ~any(modis_day >= date_range_nums(:,1) & modis_day <= date_range_nums(:,2)) %First, check if the day of the swath is within the date ranges we want. If not, skip that file.
    elseif week_mode == 2 && ~any(modis_day==weekdays); %If looking for weekdays (mode 2) skip any days not part of the week
    elseif week_mode == 3 && any(modis_day==weekdays); %If looking for weekends (mode 3) skip any days that are regular weekdays
    else
        fprintf('\t Adding to map...\n')
        
        %Move to the directory for the year that the file comes from, if
        %not already there
        directory = ['/Volumes/share/GROUP/SAT/MODIS/MOD04_L2/Cont_US_',num2str(year)];
        if ~strcmp(pwd,directory); cd(directory); end
                        
        ncount = ncount + 1;
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
        
        [modis_latcorn, modis_loncorn] = simple_corner_calc(modis_lat, modis_lon);
        
        %MODIS AOD data seems to enter a -9999 if the data is invalid
        %Further, under 'Attributes' in the original hdf file, the scale factor is noted to be 1e-3
        aod = double(modis_aod); aod(aod<-1000) = NaN;
        
        debug = struct('xx',0,'fx',[],'aod',[],'aw',[],'count',[],'lat',[],'lon',[]);
        tmp_aod = zeros(numel(bin_aods(:,:,i)),1); tmp_aod(:)=NaN;
        tmp_areaweight = zeros(numel(bin_aods(:,:,i)),1); tmp_areaweight(:)=NaN;
        tmp_count = zeros(numel(bin_aods(:,:,i)),1);
        for a=1:numel(modis_lat)
            as = []; fxes = []; bs = [];
            area = m_lldist([modis_loncorn(1,a) modis_loncorn(2,a)], [modis_latcorn(1,a), modis_latcorn(2,a)]) * m_lldist([modis_loncorn(1,a) modis_loncorn(4,a)], [modis_latcorn(1,a), modis_latcorn(4,a)]);
            
            x1 = modis_loncorn(1,a); x2 = modis_loncorn(2,a); x3 = modis_loncorn(3,a); x4 = modis_loncorn(4,a);
            y1 = modis_latcorn(1,a); y2 = modis_latcorn(2,a); y3 = modis_latcorn(3,a); y4 = modis_latcorn(4,a);
            xv = [x1, x2, x3, x4, x1]; yv = [y1, y2, y3, y4, y1];
            debug(a).lon = xv; debug(a).lat = yv;
            xx = inpolygon(bin_lons, bin_lats, xv,yv);
            fx = find(xx); if ~isrow(fx); fx=fx'; end
            if ~isempty(fx); fxes=[fxes, fx];end
            if ~isnan(aod(a)); as=[as, a]; end
            
            for b=fx
                tmp_aod(b) = nanmean([tmp_aod(b)*tmp_count(b), aod(a)]);
                tmp_areaweight(b) = nanmean([tmp_areaweight(b)*tmp_count(b), 1/area]);
                tmp_count(b) = tmp_count(b)+1; 
            end
            debug(a).xx = sum(sum(xx));
            debug(a).fx = fxes; debug(a).aod = tmp_aod; debug(a).aw = tmp_areaweight; debug(a).count = tmp_count;
            bin_aods(:,:,i) = reshape(tmp_aod,size(bin_aods(:,:,i)));
            bin_areaweight(:,:,i) = reshape(tmp_areaweight,size(bin_areaweight(:,:,1)));
        end
    end
end
mean_aods = nansum(bin_aods .* bin_areaweight,3)./nansum(bin_areaweight,3);

close all
m_proj('Albers Equal-Area Conic','lon',[min(lons), max(lons)], 'lat', [min(lats), max(lats)]);
ax = m_pcolor(bin_lons, bin_lats, mean_aods);
shading flat
%alpha(ax,0.5);
m_grid;

for j = 1:length(city_names)
    m_line(city_latlon(j,2),city_latlon(j,1)+.1,'marker','x','markersize',16,'color',[1 1 1],'linestyle','none');
    
    txt = m_text(city_latlon(j,2)+0.1,city_latlon(j,1),['  ',city_names{j}],'vertical','top','BackgroundColor',[.7 .7 .7]);
end
week_title = {'','Weekday','Weekend'};
caxis([0 200]); cb = colorbar;

%Find unique years in filenames for map title
file_years = zeros(1,numel(in_range_files));
for y=1:numel(in_range_files); file_years(y) = str2double(in_range_files{y}(11:14)); end
file_years = unique(file_years);
title(sprintf('Avg. %s MODIS AOD over SJV for %u--%u',week_title{week_mode},min(file_years),max(file_years)),'fontsize',16,'fontweight','bold');
xlabel(sprintf('n = %u', ncount),'fontsize',14);
ylabel(cb,'optical depth (unitless \times 10^{-3})','fontsize',12,'fontweight','bold');

cd(olddir);

