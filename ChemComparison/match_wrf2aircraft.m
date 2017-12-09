function [ Match ] = match_wrf2aircraft( Raw, wrf_dirs, wrf_prof_mode )
%MATCH_WRF2AIRCRAFT Generates a structure matching WRF output to aircraft data
%   MATCH = MATCH_WRF2AIRCRAFT( RAW, WRF_DIRS ) Matches the aircraft data
%   contained in RAW to the WRF data in wrfout_* files contained in
%   WRF_DIRS. RAW must be a structure that contains the fields "lon",
%   "lat", "pres", "dvec", and "campaign". "lon" and "lat" are the
%   longitude and latitude of the aircraft data (longitude must be in
%   degrees west is < 0). "pres" is the pressure of the aircraft
%   measurements in hPa. "dvec" is a vector of Matlab date numbers giving
%   the date and time of each aircraft measurement. "campaign" is a string
%   that gives the campaign name (solely for record keeping - does not need
%   to be any specific format). The other fields in RAW must have names
%   that match the names of WRF-Chem variables you want to match with the
%   aircraft data. For example, if you want to compare aircraft NO, NO2,
%   and HNO3 data with WRF files that contain that data as variables no,
%   no2, and hno3, RAW should have the aircraft data stored in fields no,
%   no2, and hno3, respectively. WRF_DIRS must be a string giving single
%   directory where all necessary WRF data can be found, or a cell array of
%   those directories.
%
%   MATCH = MATCH_WRF2AIRCRAFT( RAW, WRF_DIRS, 'monthly' ) tells this
%   function to look for monthly averaged WRF files (following the naming
%   convention WRF_BEHR_monthly_XX.nc) instead of time resolved output
%   files, where XX is the 2-digit month number.
%
%   MATCH = MATCH_WRF2AIRCRAFT( RAW, WRF_DIRS, 'v2' ) tells this to use
%   BEHR v2 NO2 profiles that are stored in .mat files with names
%   mXX_NO2_profile.mat Note that these files only contain NO2 data.
%
%   The output structure has three primary fields: "data", "wrf", and
%   "indicies". It is constructed such that all the subfields in "data" and
%   most of the subfields in "wrf" are just vectors. In "data", the fields
%   are just a concatenation in time of fields from the aircraft data for
%   the campaign. In "wrf", the vector fields are just the WRF data that
%   corresponds to the aircraft measurement; a whole series of aircraft
%   measurements that falls in the same grid cell will correspond to the
%   same WRF value.
%
%   The "indicies" field allows you to map the data in the previous two
%   fields to the WRF grid. The subfields give the index in the
%   corresponding WRF dimension. The WRF xlon and xlat fields are provided
%   as arrays, the WRF pressure for each grid cell is provided as a vector,
%   the bottom_top index should just be taken as-is, that is a model layer
%   index.
%
%   The "time" index bears mention; since a large number of WRF files do
%   not have any corresponding aircraft data, they are skipped entirely.
%   The "datevec" subfield under "wrf" has the times of the files that had
%   corresponding aircraft data, the time index is the index in this
%   vector.

% History: 
%   27 Oct 2017 - modified to (a) split out the matching code from the code
%   that loads the aircraft data and (b) to let it read monthly and BEHR
%   version 2 NO2 profiles. Verified that old functionality didn't change
%   by running it with match_wrf2campaigns.dc3 using wrf_dirs =
%   '/Volumes/share2/USERS/LaughnerJ/WRF/DC3/lnox_off-fixed_BCs' and
%   comparing vs. '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/WRF/Chem 
%   Comparison/DC3-Comparison-no_lnox-fixedBC-using_wrfgridcorners.mat'
E = JLLErrors;

%%%%%%%%%%%%%%%%%%
% INPUT CHECKING %
%%%%%%%%%%%%%%%%%%

req_fields = {'lon', 'lat', 'pres', 'dvec', 'campaign'};
if ~isstruct(Raw)
    E.badinput('RAW must be a structure')
else
    xx_fields = ~ismember(req_fields, fieldnames(Raw));
    if any(xx_fields)
        E.badinput('RAW must contain the fields %s at a minimum (missing: %s)', strjoin(req_fields, ', '), strjoin(req_fields(xx_fields), ', '));
    end
end

if ischar(wrf_dirs)
    wrf_dirs = {wrf_dirs};
elseif ~iscellstr(wrf_dirs)
    E.badinput('WRF_DIRS must be a string or cell array or strings');
end

xx_exist = cellfun(@(d) exist(d, 'dir'), wrf_dirs);
if ~all(xx_exist)
    E.badinput('The following given WRF_DIRS do not exist: \n\t%s', strjoin(wrf_dirs(~xx_exist), '\n\t'));
end

allowed_prof_modes = {'daily', 'monthly', 'v2'};
if ~exist('wrf_prof_mode', 'var')
    wrf_prof_mode = 'daily';
elseif ~ismember(wrf_prof_mode, allowed_prof_modes)
    E.badinput('WRF_PROF_MODE "%s" is not one of the allowed values (%s)', wrf_prof_mode, strjoin(allowed_prof_modes, ', '));
end

% Get a unified date number vector containing the exact time of each
% measurement. UTC is in seconds after midnight for the given day
air_dvec = Raw.dvec;

[W, wrf_dvec, time_inds] = match_in_time(wrf_dirs, air_dvec, wrf_prof_mode);

% Next we need to bin the data in space. We'll do this in two steps. First,
% by lat/lon. Then, we'll go back through and assign each value to a
% pressure bin by finding the WRF pressure in that lat/lon coordinate
% closest to the aircraft pressure

% Try to find XLONG_U and XLAT_V in the wrfout files. If there, we can use
% them as grid cell edges. Otherwise, we'll have to estimate them from the
% regular XLONG and XLAT
[wrf_edge_lon, wrf_edge_lat, wrf_lon, wrf_lat] = wrf_edge_latlon(W(1).name, wrf_prof_mode);



% Loop through each cell in the west-east and south-north directions, find
% the data points that fall in the column of grid cells, then assign them a
% vertical position by which pressure level they are closest to.
%
% Need to load wrf pressures here, may as well load the other chemical
% fields too. Match any fields from the aircraft data possible, remember;
% the fields in Raw should match WRF variables if they are desired to be
% matched.
if strcmpi(wrf_prof_mode, 'v2')
    wrf_vars = {'no2', 'pres'}; % not used in loading, but is used in assigning to the match structure 
    convert_pres = false; 
else
    wi = ncinfo(W(1).name);
    all_wrf_vars = {wi.Variables.Name};
    air_fns = fieldnames(Raw);
    wrf_vars = air_fns(ismember(air_fns, all_wrf_vars));
    convert_pres = false;
    if ~ismember('pres', wrf_vars)
        % "pres" is a variable calculated by some of my NCO scripts when
        % subsetting WRF output. If it's not there (i.e. reading direct WRF
        % output) we need to compute the grid cell pressure.
        if ismember('pres', all_wrf_vars)
            wrf_vars{end+1} = 'pres';
        else
            wrf_vars{end+1} = 'P';
            wrf_vars{end+1} = 'PB';
            convert_pres = true;
        end
    end
    if ismember('lnox_total', all_wrf_vars)
        wrf_vars{end+1} = 'lnox_total';
    end
end

Match.campaign = Raw.campaign;
Match.data = Raw;
Match.data.datevec = air_dvec;
Match.wrf = make_empty_struct_from_cell(wrf_vars);
Match.wrf.xlon = wrf_lon;
Match.wrf.xlat = wrf_lat;
Match.wrf.time = wrf_dvec;

blank_vec = nan(size(Raw.lon));
Match.indicies = struct('west_east', blank_vec, 'south_north', blank_vec, 'bottom_top', blank_vec, 'time', time_inds);
if strcmpi(wrf_prof_mode, 'v2')
    % Loading version 2 .mat files, can't use the netCDF tools. There are
    % always 28 pressures, although that dimension isn't needed in the
    % current implementation.
    wrf_sz = [size(wrf_lon), 28];
else
    wrf_sz = get_wrf_array_size(wi.Filename);
end

fprintf('Matching lat/lon\n');
for a=1:wrf_sz(1)
    for b=1:wrf_sz(2)
        % First check all data, if nothing falls in this cell, move on
        xx = Raw.lon >= wrf_edge_lon(a,b) & Raw.lon <= wrf_edge_lon(a+1,b);
        yy = Raw.lat >= wrf_edge_lat(a,b) & Raw.lat <= wrf_edge_lat(a,b+1);
        if sum(xx & yy) > 0;
            Match.indicies.west_east(xx & yy) = a;
            Match.indicies.south_north(xx & yy) = b;
        end
    end
end

for d=1:numel(W)
    fprintf('Loading WRF file %d of %d\n', d, numel(W));
    [wrf_var_struct, wrf_pres] = load_one_wrf_file(W(d).name, wrf_vars, wrf_prof_mode, convert_pres);
    
    tt = find(Match.indicies.time == d);
    i_time = d;
    
    % Now go through each data point from the aircraft and assign the vertical coordinate
    for i=1:numel(tt)
        i_we = Match.indicies.west_east(tt(i));
        i_sn = Match.indicies.south_north(tt(i));
        
        if any(isnan([i_we, i_sn, i_time]))
            continue
        end
        
        wrf_pres_vec = squeeze(wrf_pres(i_we, i_sn, :));
        [~, i_bt] = min(abs(Match.data.pres(tt(i)) - wrf_pres_vec));
        Match.indicies.bottom_top(tt(i)) = i_bt;
        
        % We can match the WRF data at the same time
        for j=1:numel(wrf_vars)
            Match.wrf.(wrf_vars{j})(tt(i)) = wrf_var_struct.(wrf_vars{j})(i_we, i_sn, i_bt);
        end
    end
end

Match.generation_date = datestr(now);

% Switch pressure back to total
if convert_pres
    Match.wrf.pres = (Match.wrf.P + Match.wrf.PB)/100;
    Match.wrf = rmfield(Match.wrf, 'P');
    Match.wrf = rmfield(Match.wrf, 'PB');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME MATCHING SUBFUNCTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W, wrf_dvec, time_inds] = match_in_time(wrf_dirs, air_dvec, wrf_prof_mode)
% Get the list of wrfout files available across all given directories, find
% the difference in time between each. The first step will be to bin the
% aircraft data by which wrfout file is closest in time.
if strcmpi(wrf_prof_mode, 'monthly') || strcmpi(wrf_prof_mode, 'v2')
    [W, wrf_dvec, time_inds] = match_monthly_wrf_files(wrf_dirs, air_dvec, wrf_prof_mode);
else
    [W, wrf_dvec, time_inds] = match_daily_wrf_files(wrf_dirs, air_dvec);
end
end

function [W, wrf_months, time_inds] = match_monthly_wrf_files(wrf_dirs, air_dvec, wrf_prof_mode)
E = JLLErrors;
W = [];
for a=1:numel(wrf_dirs)
    if strcmpi(wrf_prof_mode, 'monthly')
        W = veccat(dirff(fullfile(wrf_dirs{a}, 'WRF_BEHR_monthly*.nc')));
    elseif strcmpi(wrf_prof_mode, 'v2')
        W = veccat(dirff(fullfile(wrf_dirs{a}, 'm*_NO2_profile.mat')));
    else
        E.notimplemented('File search for wrf_prof_mode == %s not implemented', wrf_prof_mode);
    end
end

wrf_months = nan(size(W));
for a=1:numel(W)
    [~,filename] = fileparts(W(a).name);
    wrf_month_str = regexp(filename, '\d\d', 'match', 'once');
    if isempty(wrf_month_str)
        E.callError('no_wrf_month', 'Cannot identify month in filename %s', filename);
    end
    wrf_months(a) = str2double(wrf_month_str);
end

% Now bin the aircraft data in time. We won't directly bin it, rather we
% assign each data point a time index that indicates which WRF file it goes
% with. We'll also keep track
time_inds = nan(size(air_dvec));
wrf_to_keep = false(size(W));

i_t = 1;
for a=1:numel(wrf_months)
    xx = month(air_dvec) == wrf_months(a);
    if sum(xx) > 0
        wrf_to_keep(a) = true;
        time_inds(xx) = i_t;
        i_t = i_t + 1;
    end
end

W(~wrf_to_keep) = [];
wrf_months(~wrf_to_keep) = [];

end

function [W, wrf_dvec, time_inds] = match_daily_wrf_files(wrf_dirs, air_dvec)
E = JLLErrors;
W = [];
for a=1:numel(wrf_dirs)
    W = veccat(dirff(fullfile(wrf_dirs{a},'wrfout*')));
end

wrf_dvec = date_from_wrf_filenames(W);


% Compute the difference between files.
wrf_timediff = diff(wrf_dvec);
if numel(unique(wrf_timediff)) == 1
    wrf_timediff = unique(wrf_timediff);
else
    warning('Difference in time between output files varies, using average (std dev of time differences = %g)', std(wrf_timediff))
    wrf_timediff = nanmean(wrf_timediff);
end

% Now bin the aircraft data in time. We won't directly bin it, rather we
% assign each data point a time index that indicates which WRF file it goes
% with.
time_inds = nan(size(air_dvec));
wrf_to_keep = false(size(W));

i_t = 1;
for a=1:numel(wrf_dvec)
    min_time = wrf_dvec(a) - wrf_timediff/2;
    max_time = wrf_dvec(a) + wrf_timediff/2;
    xx = air_dvec >= min_time & air_dvec < max_time;
    if sum(xx) > 0
        % Double check that we haven't assigned this point yet
        if any(~isnan(time_inds(xx)))
            if sum(~isnan(time_inds(xx))) == 1
                wi = ~isnan(time_inds(xx));
                warning('One of the aircraft data points (%s) is being assigned to multiple WRF files. This is probably on the edge between two WRF files.', datestr(air_dvec(wi)));
            else
                E.callError('time_assignment', 'Multiple aircraft data points are being assigned to multiple WRF files.');
            end
        end
        % Only if there are any data points at this time will we keep it.
        % WRF files with no corresponding aircraft data will be removed.
        time_inds(xx) = i_t;
        i_t = i_t + 1;
        wrf_to_keep(a) = true;
    end
end
W(~wrf_to_keep) = [];
wrf_dvec(~wrf_to_keep) = [];
end

%%%%%%%%%%%%%%%%%%%%%
% SPATIAL FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%

function [wrf_edge_lon, wrf_edge_lat, xlon, xlat] = wrf_edge_latlon(filename, wrf_prof_mode)
% wi should be a structure for one representative wrfout file returned by
% ncinfo
if strcmpi(wrf_prof_mode, 'v2')
    % If loading version 2 profiles, they're in .mat files, so we can read
    % any variables.
    wrf_vars = {};
else
    wi = ncinfo(filename);
    wrf_vars = {wi.Variables.Name};
end
if ismember('XLONG_U', wrf_vars) && ismember('XLAT_V', wrf_vars)
    wrf_edge_lon = ncread(wi.Filename, 'XLONG_U');
    wrf_edge_lat = ncread(wi.Filename, 'XLAT_V');
    xlon = ncread(wi.Filename, 'XLONG');
    xlat = ncread(wi.Filename, 'XLAT');
else
    if strcmpi(wrf_prof_mode, 'v2')
        P = load(filename);
        % Must transpose these b/c in normal WRF files the first dimension
        % is west_east, but in the version 2 files, it's set up so that the
        % first dimension is longitude. This matters when we're figuring
        % out which grid cell each aircraft point is in.
        xlon = P.PROFILE.Longitude';
        xlat = P.PROFILE.Latitude';
    else
        warning('XLONG_U or XLAT_V not available, using wrf_grid_corners() instead')
        xlon = ncread(wi.Filename, 'XLONG');
        xlat = ncread(wi.Filename, 'XLAT');
    end
    [xloncorn, xlatcorn] = wrf_grid_corners(xlon, xlat);
    
    % Average the north and south west corners to get the west midpoint; do
    % a similar process for each of the other sides. In the end,
    % wrf_edge_lon should be staggered in the first dimension and wrf_edge
    % lat in the second.
    wloncorn = squeeze(nanmean(xloncorn([1 4], :, :)));
    eloncorn = squeeze(nanmean(xloncorn([2 3], :, :)));
    slatcorn = squeeze(nanmean(xlatcorn(1:2, :, :)));
    nlatcorn = squeeze(nanmean(xlatcorn(3:4, :, :)));
    
    weloncorn = (eloncorn(2:end,:)+wloncorn(1:end-1,:))/2;
    snlatcorn = (slatcorn(:,2:end)+nlatcorn(:,1:end-1))/2;
    
    wrf_edge_lon = cat(1, wloncorn(1,:), weloncorn, eloncorn(end, :));
    wrf_edge_lat = cat(2, slatcorn(:,1), snlatcorn, nlatcorn(:, end));
end
end

%%%%%%%%%%%%%%%%
% LOADING DATA %
%%%%%%%%%%%%%%%%

function [wrf_var_struct, wrf_pres] = load_one_wrf_file(wrf_file, wrf_vars, wrf_prof_mode, convert_pres)
if strcmpi(wrf_prof_mode, 'v2')
    P = load(wrf_file);
    PROFILE = P.PROFILE;
    pres = reshape(PROFILE.Pressure,1,1,[]);
    % As in wrf_edge_latlon, we need to reorder these so that the first
    % dimension is west_east and the second is south_north
    wrf_pres = repmat(pres, [size(PROFILE.Longitude'), 1]);
    wrf_var_struct.no2 = permute(PROFILE.NO2_profile, [3 2 1]);
    wrf_var_struct.pres = wrf_pres;
else
    wrf_var_cell = cell(size(wrf_vars));
    [wrf_var_cell{:}] = read_wrf_vars('', wrf_file, wrf_vars);
    wrf_var_struct = cell2struct(wrf_var_cell(:), wrf_vars(:), 1);
    if convert_pres
        wrf_pres = (wrf_var_struct.P + wrf_var_struct.PB)/100; % convert from Pa to hPa
    else
        wrf_pres = wrf_var_struct.pres;
    end
end
end

