function [ Match ] = match_wrf2aircraft( campaign_name, wrf_dir, varargin )
%MATCH_WRF2AIRCRAFT Generates a structure matching WRF output to aircraft data
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
%
%   For DC3 the following aircraft fields are:
%       no - NO from the NO_ESRL (chemiluminescence?) field
%       no2 - NO2 from our TD-LIF
%       mpn - MPN from our TD-LIF
%       hno3 - HNO3 as the average of the SAGA (mist chamber) and CIT
%       (chemical ionization time of flight mass spec).
%
%   Although this currently has inputs, they are not used in this version.
E = JLLErrors;
p = inputParser;
p.addParameter('data_model',[]);
p.addParameter('map_fxn', []);
p.parse(varargin{:});
pout = p.Results;

data_model = pout.data_model;

fprintf('Matching %s campaign data to WRF data in %s\n', campaign_name, wrf_dir);

if isempty(data_model)
    data_model = struct('no', 'NO_ESRL',...
                        'no2', 'no2_lif',...
                        'mpn', 'MPN_TDLIF',...
                        'hno3', 'HNO3_SAGA',...
                        'PHOTR_NO2', 'JNO2NOO3P',...
                        'TT', 'TEMPERATURE',...
                        'QVAPOR', 'MixingRatio');
end

% Add the location fields if necessary
if ~isfield(data_model, 'lon')
    data_model.lon = 'LONGITUDE';
end
if ~isfield(data_model, 'lat')
    data_model.lat = 'LATITUDE';
end
if ~isfield(data_model, 'pres')
    data_model.pres = 'PRESSURE';
end

% Make a cell array of the field to concatenate from the values of the
% fields in data_model
merge_fields = struct2cell(data_model);

% TODO: modify campaign_wide_ops to handle multiple requested fields? Is
% this still a todo?
% Output to structure raw; anything in it will be binned
Out = campaign_wide_ops(campaign_name, merge_fields, 'cat', 'datefmt','datenum');

% Monday JLL: to handle arbitrary fields, this conversion needs to be unit
% aware. One option is to have campaign_wide_opts return unit information
% too, then read WRF-chem, and use convert_units. Another option is to
% require a conversion function handle to be passed, though that's hard for
% a first time user. Third option is just to split out the actual matching
% code into its own function and pass it the Raw structure.

% Convert the output chemical species here to the Raw structure, also
% convert to ppm since that's how WRF outputs concentrations. Field names must
% match the variable names in WRF-Chem.
Raw.no = Out.data.NO_ESRL .* 1e-9 .* 1e6;
Raw.no2 = Out.data.no2_lif .* 1e-12 .* 1e6;
Raw.mpn = Out.data.MPN_TDLIF .* 1e-12 .* 1e6;
Raw.hno3 = (Out.data.HNO3_SAGA + Out.data.HNO3_CIT) .* 1e-12/2 .* 1e6;
Raw.PHOTR_NO2 = Out.data.JNO2NOO3P * 60; % WRF outputs in per minute
Raw.lon = Out.data.LONGITUDE; % the correction to negative is west is handled in read_merge_fields
Raw.lat = Out.data.LATITUDE;
Raw.pres = Out.data.PRESSURE;
Raw.TT = Out.data.TEMPERATURE;
Raw.QVAPOR = Out.data.MixingRatio * 1e-3; % 1) I'm assuming that "MixingRatio" is of water, since it comes from a hygrometer
                                          % 2) The units for this mixing ratio are g/kg, in WRF it is kg/kg.


Raw.co = Out.data.CO_DACOM .* 1e-9 .* 1e6;
Raw.o3 = Out.data.O3_ESRL .* 1e-9 .* 1e6;

% Get a unified date number vector containing the exact time of each
% measurement. UTC is in seconds after midnight for the given day
air_dvec = Out.dates + Out.utcs ./ (60*60*24);

% Get the list of wrfout files available, find the difference in time
% between each. The first step will be to bin the aircraft data by which
% wrfout file is closest in time.
W = dir(fullfile(wrf_dir,'wrfout*'));
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

% Next we need to bin the data in space. We'll do this in two steps. First,
% by lat/lon. Then, we'll go back through and assign each value to a
% pressure bin by finding the WRF pressure in that lat/lon coordinate
% closest to the aircraft pressure

% Try to find XLONG_U and XLAT_V in the wrfout files. If there, we can use
% them as grid cell edges. Otherwise, we'll have to estimate them from the
% regular XLONG and XLAT
wi = ncinfo(fullfile(wrf_dir, W(1).name));
[wrf_edge_lon, wrf_edge_lat] = wrf_edge_latlon(wi);



% Loop through each cell in the west-east and south-north directions, find
% the data points that fall in the column of grid cells, then assign them a
% vertical position by which pressure level they are closest to.
%
% Need to load wrf pressures here, may as well load the other chemical
% fields too. Match any fields from the aircraft data possible, remember;
% the fields in Raw should match WRF variables if they are desired to be
% matched.
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


Match.campaign = campaign_name;
Match.data = Raw;
Match.data.datevec = air_dvec;
Match.wrf = make_empty_struct_from_cell(wrf_vars);
Match.wrf.xlon = ncread(wi.Filename, 'XLONG');
Match.wrf.xlat = ncread(wi.Filename, 'XLAT');
Match.wrf.time = wrf_dvec;
 
blank_vec = nan(size(Raw.lon));
Match.indicies = struct('west_east', blank_vec, 'south_north', blank_vec, 'bottom_top', blank_vec, 'time', time_inds);
wrf_sz = get_wrf_array_size(wi.Filename);

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
    wrf_var_cell = cell(size(wrf_vars));
    [wrf_var_cell{:}] = read_wrf_vars(wrf_dir, W(d), wrf_vars);
    wrf_var_struct = cell2struct(wrf_var_cell, wrf_vars, 1);
    if convert_pres
        wrf_pres = (wrf_var_struct.P + wrf_var_struct.PB)/100; % convert from Pa to hPa
    else
        wrf_pres = wrf_var_struct.pres;
    end
   
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



end

function [wrf_edge_lon, wrf_edge_lat] = wrf_edge_latlon(wi)
% wi should be a structure for one representative wrfout file returned by
% ncinfo
wrf_vars = {wi.Variables.Name};
if ismember('XLONG_U', wrf_vars) && ismember('XLAT_V', wrf_vars)
    wrf_edge_lon = ncread(wi.Filename, 'XLONG_U');
    wrf_edge_lat = ncread(wi.Filename, 'XLAT_V');
else
    warning('XLONG_U or XLAT_V not available, using wrf_grid_corners() instead')
    xlon = ncread(wi.Filename, 'XLONG');
    xlat = ncread(wi.Filename, 'XLAT');
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

