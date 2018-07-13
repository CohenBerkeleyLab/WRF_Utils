function average_data = wrf_day_weighted_average(wrf_lon, target_local_hour, utc_hours, wrf_data)
%WRF_DAY_WEIGHTED_AVERAGE Create a longitude/time-weighted average of WRF data
%
%   Sometimes you may want an average of WRF data that reflects the
%   atmosphere at a given local time, and so this function does an average
%   of WRF data weighted by how close to that local time it is at each
%   location. 
%
%   AVG_DATA = WRF_DAY_WEIGHTED_AVERAGE( WRF_LON, TARGET_LOCAL_HOUR,
%   UTC_HOURS, WRF_DATA ) Averages WRF_DATA, weighting it most heavily for
%   data points where the local time is near TARGET_LOCAL_HOUR. This uses a
%   generalization of the monthly NO2 profile weighting scheme described in
%   Laughner et al. 2018 (doi: 10.5194/essd-2018-66, Eq. 11). This
%   calculates the local time from UTC time assuming that the offset is
%   given by WRF_LON/15. The inputs required are:
%
%       WRF_LON - the 2D longitude array read in (usually) from the WRF
%       variable XLONG.
%
%       TARGET_LOCAL_HOUR - a scalar number giving the local hour, e.g.
%       13.5 is 13:30 local standard time (this function does not correct
%       for daylight savings time at all).
%
%       UTC_HOURS - the UTC hour corresponding to each 3D slice in
%       WRF_DATA.
%
%       WRF_DATA - a single 4D numeric array, or cell array of such arrays,
%       that represent 3D WRF variables at multiple timesteps. The first
%       two dimensions of each array must be the same as the size of
%       WRF_LON and the fourth and final dimension must be the same length
%       as UTC_HOURS. For the ith array, WRF_DATA{i}(:,:,:,1) must
%       correspond to UTC_HOURS(1), WRF_DATA{i}(:,:,:,2) to UTC_HOURS(2)
%       and so on.
E = JLLErrors;

if ~isnumeric(wrf_lon) || ~ismatrix(wrf_lon)
    E.badinput('WRF_LON must be a 2D numeric array');
end

if ~isnumeric(target_local_hour) || ~isscalar(target_local_hour)
    E.badinput('TARGET_LOCAL_HOUR must be a scalar number')
end

if ~isnumeric(utc_hours) || ~isvector(utc_hours)
    E.badinput('UTC_HOURS must be a numeric vector')
end

if isnumeric(wrf_data)
    wrf_data = {wrf_data};
elseif ~iscell(wrf_data) || any(~cellfun(@isnumeric, wrf_data))
    E.badinput('WRF_DATA must be a numeric array or cell array of numeric arrays');
end

data_sizes = cellfun(@size, wrf_data, 'uniform', false);
lon_size = size(wrf_lon);
lon_size_check = cellfun(@(x) isequal(lon_size, x(1:2)), data_sizes);
if any(~lon_size_check)
    E.badinput('The first two dimensions of all data arrays must have the same size as WRF_LON')
end

n_data = numel(wrf_data);

n_hours = numel(utc_hours);
hr_size_check = cellfun(@(x) size(x,4) == n_hours, wrf_data);
if any(~hr_size_check)
    E.badinput('The fourth dimension of all data arrays must have the same length as utc_hours');
end

data_ndims = cellfun(@(x) ndims(x) > 4, wrf_data);
if any(data_ndims)
    E.notimplemented('WRF_DAY_WEIGHTED_AVERAGE not set up to take data arrays with >4 dimensions');
end


running_average_data = cell(size(wrf_data));
for i_dat = 1:n_data
    running_average_data{i_dat} = RunningAverage;
end

for i_hr = 1:n_hours
    lonwt_2d = calc_lonweight(wrf_lon, target_local_hour, utc_hours(i_hr));
    for i_dat = 1:n_data
        this_data = wrf_data{i_dat}(:,:,:,i_hr);
        lonwt = repmat(lonwt_2d, 1, 1, size(this_data, 3));
        running_average_data{i_dat}.addData(this_data, lonwt);
    end
end

average_data = cell(size(wrf_data));
for i_dat = 1:n_data
    average_data{i_dat} = running_average_data{i_dat}.getWeightedAverage();
end
end

function lonwt = calc_lonweight(wrf_lon, target_local_hour, utc_hour)
E = JLLErrors;
if ~isscalar(target_local_hour) || ~isnumeric(target_local_hour)
    E.badinput('TARGET_LOCAL_HOUR must be a scalar number')
elseif ~isscalar(utc_hour) || ~isnumeric(utc_hour)
    E.badinput('UTC_HOUR must be a scalar number')
end

lonwt = 1 - abs( target_local_hour - (wrf_lon/15) - utc_hour );
lonwt = clipmat(lonwt, 0, 1);

end