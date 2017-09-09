function [ Match ] = match_wrf_noaa( wrf_dir, noaa_file )
%MATCH_WRF_NOAA Matches NOAA integrate surface data obs. to WRF grid cells
%   MATCH = MATCH_WRF_NOAA( WRF_DIR, NOAA_FILE ) Loads each WRF files in
%   WRF_DIR and matches the observations in NOAA_FILE to the WRF grid.
%   Returns the structure MATCH, which contains WRF and NOAA winds as U and
%   V components (m/s), and temperature (K). If multiple sites are present
%   in a grid cell, their values are averaged. Only values within 15 min of
%   the time of the WRF file are used, and the closest observation in time
%   is used.
%
%   WRF_DIR should point to a directory that contains WRF_BEHR_*.nc files.
%   They must include the variables XLONG, XLAT, XLONG_U, XLAT_U, Times, U,
%   V, COSALPHA, SINALPHA, and TT.
%
%   NOAA_FILE may either be a path to the .mat file containing the output
%   from load_noaa_isd or the cell array contained in that file.
%
%   In MATCH, the output fields wrf_U, wrf_V, wrf_T, noaa_U, noaa_V, noaa_T
%   are arrays with dimensions nsites x ntimes per day x ndays.

E = JLLErrors;
DEBUG_LEVEL = 1;

% Input checking
if ~ischar(wrf_dir)
    E.badinput('WRF_DIR must be a string');
elseif ~exist(wrf_dir,'dir')
    E.badinput('%s is not a directory', wrf_dir)
end
if iscell(noaa_file) && ~all(iscellcontents(noaa_sites, @(x) isa(x, 'NOAAISDSite')))
    E.badinput('If passing the NOAA sites directly, it must be a cell array containing NOAAISDSite objects')
elseif ischar(noaa_file) 
    [~,~,ext] = fileparts(noaa_file);
    if ~exist(noaa_file, 'file') || ~strcmp(ext,'.mat')
        E.badinput('%s is not a .mat file', noaa_file)
    end
elseif ~iscell(noaa_file) && ~ischar(noaa_file)
    E.badinput('NOAA_FILE must be a cell array or string')
end
    

% First get latitude and longitude out of the first WRF file, plus the
% times (to figure out how many time per day there are)
F = dir(fullfile(wrf_dir, 'WRF_BEHR*.nc'));
wi = ncinfo(fullfile(wrf_dir, F(1).name));
wrf_lon = ncread(wi.Filename, 'XLONG');
wrf_lon = wrf_lon(:,:,1);
wrf_lat = ncread(wi.Filename, 'XLAT');
wrf_lat = wrf_lat(:,:,1);
wrf_lon_u = ncread(wi.Filename, 'XLONG_U');
wrf_lon_u = wrf_lon_u(:,:,1);
wrf_lat_v = ncread(wi.Filename, 'XLAT_V');
wrf_lat_v = wrf_lat_v(:,:,1);
times = ncread(wi.Filename, 'Times')'; % must be transposed to have the proper orientation
ntimes = size(times,1);

% Load the NOAA file, which will be a cell array of NOAAISDSite objects
if ~iscell(noaa_file)
    tmp = load(noaa_file);
    noaa_sites = tmp.noaa_sites;
else
    noaa_sites = noaa_file;
end

% Set all sites to use U and V wind definition
for a=1:numel(noaa_sites)
    noaa_sites{a}.set_wind_def('uv');
end

% Record the x,y indicies of the WRF grid that each site falls in
site_wrf_inds = nan(numel(noaa_sites),2);
for a=1:numel(noaa_sites)
    xx = noaa_sites{a}.lon >= wrf_lon_u(1:end-1,:) & noaa_sites{a}.lon < wrf_lon_u(2:end,:);
    yy = noaa_sites{a}.lat >= wrf_lat_v(:,1:end-1) & noaa_sites{a}.lat < wrf_lat_v(:,2:end);
    
    if sum(xx(:) & yy(:)) > 1
        E.callError('site_assignment','Site %s falls in multiple WRF grid cells',noaa_sites{a}.usaf_id);
    elseif sum(xx(:) & yy(:)) < 1 && DEBUG_LEVEL > 0
        fprintf('Site %s cannot be assigned to a grid cell\n', noaa_sites{a}.usaf_id);
        continue
    end
    
    [x,y] = find(xx & yy);
    site_wrf_inds(a,:) = [x,y];
end

% Handle the possibility of multiple sites falling into one WRF grid cell.
[site_wrf_inds, sortvec] = sortrows(site_wrf_inds);
noaa_sites = noaa_sites(sortvec);
nans = all(isnan(site_wrf_inds),2);
site_wrf_inds(nans,:) = [];
noaa_sites(nans) = [];

[site_blocks, block_wrf_inds] = findBlock(site_wrf_inds,1);

% Now match up the NOAA observations with WRF output. We will load one day
% at a time to keep the memory load down.

blank_mat = nan(size(site_blocks,1), ntimes, numel(F));

Match.wrf_U = blank_mat;
Match.wrf_V = blank_mat;
Match.wrf_T = blank_mat;
Match.noaa_U = blank_mat;
Match.noaa_V = blank_mat;
Match.noaa_T = blank_mat;

% Save the WRF grid cell lat/lon center
Match.lon = nan(size(site_blocks,1),1);
Match.lat = nan(size(site_blocks,1),1);
% Save the site IDs as USAF-WBAN
Match.site = cell(size(site_blocks,1),1);
for a=1:numel(Match.lon)
    wrf_i = block_wrf_inds(a,1);
    wrf_j = block_wrf_inds(a,2);
    Match.lon(a) = wrf_lon(wrf_i, wrf_j);
    Match.lat(a) = wrf_lat(wrf_i, wrf_j);
    
    sites = site_blocks(a,1):site_blocks(a,2);
    c = cell(size(sites));
    for b=1:numel(sites)
        usaf = noaa_sites{sites(b)}.usaf_id;
        wban = noaa_sites{sites(b)}.wban_id;
        if isnan(wban); wban = 99999; end
        c{b} = sprintf('%s-%05d', usaf, wban);
    end
    Match.site{a} = c;
end

% Ready an array to hold datenums
Match.dnums = nan(ntimes, numel(F));

% This defines how much time we allow between WRF output time and the
% observation time
max_time_sep = 15;
time_sep_units = 'minutes'; % can be allowed units for NOAAISDSite private method nearest_obs_in_time
% This defines what quality data is allowed (1 = only best)
quality_level = 1;

for a=1:numel(F)
    if DEBUG_LEVEL > 0; fprintf('Matching obs to WRF file %s\n', F(a).name); end
    [U_wrf, V_wrf, COSALPHA, SINALPHA, T_wrf, times] = read_wrf_vars(wrf_dir, F(a), {'U','V','COSALPHA','SINALPHA','TT', 'Times'});
    % Only care about surface
    U_wrf = double(squeeze(U_wrf(:,:,1,:)));
    V_wrf = double(squeeze(V_wrf(:,:,1,:)));
    T_wrf = double(squeeze(T_wrf(:,:,1,:)));
    
    % Does not have vertical dimension and doesn't change day-to-day
    COSALPHA = squeeze(COSALPHA(:,:,1));
    SINALPHA = squeeze(SINALPHA(:,:,1));
    
    % Transpose to make a valid character array for datenum
    times = times';
    
    % Unstagger and transform to earth-relative
    [U_wrf, V_wrf] = wrf_winds_transform(U_wrf,V_wrf,COSALPHA,SINALPHA);
    
    for b=1:ntimes
        tnum = datenum(times(b,:));
        Match.dnums(b,a) = tnum;
        if DEBUG_LEVEL > 0; fprintf('\tWorking on %s\n', datestr(tnum,'HH:MM')); end
        for c=1:size(site_blocks,1)
            wrf_i = block_wrf_inds(c,1);
            wrf_j = block_wrf_inds(c,2);
            Match.wrf_U(c, b, a) = U_wrf(wrf_i, wrf_j, b);
            Match.wrf_V(c, b, a) = V_wrf(wrf_i, wrf_j, b);
            Match.wrf_T(c, b, a) = T_wrf(wrf_i, wrf_j, b);
            
            
            [u_obs, v_obs, t_obs] = avg_site_data(site_blocks(c,:), tnum);
            Match.noaa_U(c, b, a) = u_obs;
            Match.noaa_V(c, b, a) = v_obs;
            Match.noaa_T(c, b, a) = t_obs;
        end
    end
end

% Convert NOAA temperatures to Kelvin
Match.noaa_T = Match.noaa_T + 273;

    function [u,v,t] = avg_site_data(inds, dnum)
        inds_vec = inds(1):inds(2);
        nsites = numel(inds_vec);
        u = nan(1, nsites);
        v = nan(1, nsites);
        t = nan(1, nsites);
        for i=1:nsites
            t_i = noaa_sites{inds_vec(i)}.get_temperature_nearest_in_time(dnum, max_time_sep, time_sep_units, quality_level);
            if ~isempty(t_i)
                % An empty array indicates no observation within the
                % max_time_sep of the WRF time represented by dnum
                t(i) = t_i;
                [u(i), v(i)] = noaa_sites{inds_vec(i)}.get_wind_nearest_in_time(dnum, max_time_sep, time_sep_units, quality_level);
            end
        end
        u = nanmean(u);
        v = nanmean(v);
        t = nanmean(t);
    end

end

