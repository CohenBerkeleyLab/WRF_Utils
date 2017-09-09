function [ Site ] = parse_noaa_isd( file, lonlim, latlim, datelim )
%PARSE_NOAA_ISD Parse a NOAA integrate surface database file
%   PARSE_NOAA_ISD( FILE ) parses the NOAA ISD file at path FILE and
%   returns a NOAAISDSite object.
%
%   PARSE_NOAA_ISD( FILE, LONLIM, LATLIM ) will parse the file if the site
%   falls within the longitude and latitude limits. Otherwise it returns
%   the default NOAAISDSite object. You can test for this case if the
%   usaf_id field is the same as the object returned by constructing a
%   NOAAISDSite object with no arguments.

E = JLLErrors;
if ~ischar(file)
    E.badinput('FILE must be a string')
elseif ~exist(file,'file')
    E.filenotfound(file);
end

if ~exist('lonlim','var')
    lonlim = [-180, 180];
elseif ~isnumeric(lonlim) ||  numel(lonlim) ~= 2 
    E.badinput('LONLIM must be a 2 element numeric vector');
end

if ~exist('latlim','var')
    latlim = [-90, 90];
elseif ~isnumeric(latlim) || numel(latlim) ~= 2
    E.badinput('LATLIM must be a 2 element numeric vector');
end

if ~exist('datelim','var')
    datelim = [];
elseif ~isnumeric(datelim) || numel(datelim) ~= 2 || any(datelim < 0)
    E.badinput('DATELIM must be a 2 element numeric vector with both elements >0');
end
    

% Best guess for the number of observations
init_obs_vec = nan(1,1e6);
init_obs_cell = cell(1,1e6);

dvec = init_obs_vec;
wind_dir = init_obs_vec;
wind_dir_qual = init_obs_cell;
wind_vel = init_obs_vec;
wind_vel_qual = init_obs_cell;
temperature = init_obs_vec;
temperature_qual = init_obs_cell;

fid = fopen(file);
cleanup_fxn = onCleanup(@() cleanup(fid));
tline = fgetl(fid);

% Initial checking; get the static fields and make sure that we are in in
% lat/lon domain

[usaf_id, wban_id, lon, lat, elev] = parse_line(tline);
if lon < min(lonlim) || lon > max(lonlim) || lat < min(latlim) || lat > max(latlim)
    fprintf('Site outside specified longitude and latitude bounds\n')
    Site = NOAAISDSite();
    return
else
    Site = NOAAISDSite(usaf_id, wban_id, lon, lat, elev);
end

i = 0;
while ischar(tline)
    i = i + 1;
    [~, ~, ~, ~, ~, dvec(i), wind_dir(i), wind_dir_qual{i}, wind_vel(i), wind_vel_qual{i}, temperature(i), temperature_qual{i}] = parse_line( tline );
    tline = fgetl(fid);
end

dvec = dvec(1:i);
wind_dir = wind_dir(1:i);
wind_dir_qual = wind_dir_qual(1:i);
wind_vel = wind_vel(1:i);
wind_vel_qual = wind_vel_qual(1:i);
temperature = temperature(1:i);
temperature_qual = temperature_qual(1:i);

if ~isempty(datelim)
    xx = dvec >= min(datelim) & dvec <= max(datelim);
end

if sum(xx) > 0
    Site.set_dates(dvec(xx));
    Site.set_wind_dir(wind_dir(xx), wind_dir_qual(xx));
    Site.set_wind_vel(wind_vel(xx), wind_vel_qual(xx));
    Site.set_temperature(temperature(xx), temperature_qual(xx));
else
    Site = NOAAISDSite();
    fprintf('Site %s has no data within dates specified\n', usaf_id);
end

end

function [usaf_id, wban_id, lon, lat, elev, dnum, wind_dir, wind_dir_quality, wind_vel, wind_vel_quality, temperature, temperature_quality] = parse_line( line_in )
E = JLLErrors;
if ~ischar(line_in)
    E.badinput('LINE_IN must be a string')
end

wban_fill = 99999;
lon_fill = 999999;
lon_scale = 1000;
lat_fill = 99999;
lat_scale = 1000;
elevation_fill = 9999;
elevation_scale = 1;

wind_dir_fill = 999;
wind_dir_scale = 1;

wind_vel_fill = 9999;
wind_vel_scale = 10;

temperature_fill = 9999;
temperature_scale = 10;

usaf_id = line_in(5:10);
wban_id = convert_str_to_val(line_in(11:15), wban_fill, 1);
dnum = datenum(line_in(16:27),'yyyymmddHHMM');

lon = convert_str_to_val(line_in(35:41), lon_fill, lon_scale);
lat = convert_str_to_val(line_in(29:34), lat_fill, lat_scale);
elev = convert_str_to_val(line_in(47:51), elevation_fill, elevation_scale);
wind_dir = convert_str_to_val(line_in(61:63), wind_dir_fill, wind_dir_scale);
wind_dir_quality = line_in(64);
wind_vel = convert_str_to_val(line_in(66:69), wind_vel_fill, wind_vel_scale);
wind_vel_quality = line_in(70);
temperature = convert_str_to_val(line_in(88:92), temperature_fill, temperature_scale);
temperature_quality = line_in(93);
end

function val = convert_str_to_val(substr, fill, scale)
val = str2double(substr);
if val == fill
    val = nan;
else
    val = val ./ scale;
end


end

function cleanup(fid)
fprintf('Closing fid %d\n', fid);
fclose(fid);
end