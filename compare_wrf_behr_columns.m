function [ col_diff, wrf_no2, behr_no2, longrid, latgrid ] = compare_wrf_behr_columns( wrf_info, start_date, end_date )
%COMPARE_WRF_BEHR_COLUMNS Calculate the difference in BEHR and WRF NO2 columns
%   Function used to compare WRF output NO2 columns with observed NO2
%   columns in the BEHR product. Takes 3 arguments:
%
%       1) A structure output from ncinfo with pointing to the netCDF files
%       output from running (slurm)run_wrf_output.sh. This can be an array,
%       where each element points to a day's worth of files.
%
%       2) The starting date to average the BEHR columns from.
%
%       3) The ending date to average the BEHR columns to.
%
%   Outputs 3 quantities: the absolute difference between the columns (in
%   molec./cm^2), the WRF columns interpolated to the 0.05 x 0.05 deg grid
%   also used when oversampling BEHR, and the BEHR columns averaged over
%   the period of time specified.
%
%   This function expects that the netCDF files containing the WRF output
%   have already been averaged over time such that there is no hourly
%   dependence. This takes advantage of the averaging code already written
%   in (slurm)run_wrf_output.sh.
%
%   Josh Laughner <joshlaugh5@gmail.com> 20 Jul 2015

DEBUG_LEVEL = 1;
E = JLLErrors;

behr_dir = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_2014';
behr_prefix = 'OMI_BEHR_omiCloudAMF_';




%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% WRF PROCESSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if DEBUG_LEVEL > 0
    fprintf('--- Calculating WRF columns ---\n');
end

% First we'll handle the WRF columns. We'll need to calculate the column up
% to the tropopause first, then average over all days.

wrf_vars = {wrf_info.Variables.Name};
ww = strcmp('no2_ndens',wrf_vars);
sz = wrf_info.Variables(ww).Size;

if numel(sz) > 3
    E.callError('wrf_unaveraged', 'This function expects that the time dimension has been averaged over in the WRF input, see documentation')
elseif numel(sz) < 3
    E.badvar('no2_ndens','expected to be 3-D')
end

wrf_no2 = nan(sz(1), sz(2), numel(wrf_info));

if DEBUG_LEVEL > 0
    fprintf('\tImporting netCDF variables... ');
end

for a=1:numel(wrf_info)
    no2 = ncread(wrf_info(a).Filename, 'no2_ndens');
    zlev = ncread(wrf_info(a).Filename, 'zlev');
    zlev = zlev * 100; % convert from m to cm 
    
    % Set all values above the tropopause to NaN so they don't add into the
    % column
    tp_lev = find_wrf_tropopause(wrf_info(a), true);
    for x=1:sz(1)
        for y=1:sz(2)
            no2(x,y,tp_lev(x,y):end) = nan;
        end
    end
    
    columns = nansum(no2 .* zlev, 3);
    wrf_no2(:,:,a) = columns;
end

if DEBUG_LEVEL > 0;
    fprintf('Done\n');
end

wrf_no2 = nanmean(wrf_no2,3);


% ================================================================== %

% Next, we need to calculate the BEHR columns over the time period.

if DEBUG_LEVEL > 0
    fprintf('--- Averaging BEHR columns ---\n');
end

[~, behr_no2, longrid, latgrid] = no2_column_map_2014(start_date, end_date, [-125, -65], [25 50],...
    'behrdir', behr_dir, 'fileprefix', behr_prefix, 'makefig', false);


% ================================================================== %

% Finally we'll interpolate the WRF columns to match the 0.05 x 0.05 deg grid
% BEHR is oversampled to. I chose to use interpolation rather than
% rewriting the regridding algorithm used in BEHR because since the WRF
% grid is the same each day, oversampling would do nothing. Interpolation
% at least accounts for the linear variation between grid points.

if DEBUG_LEVEL > 0
    fprintf('--- Interpolating WRF columns to BEHR grid ---\n');
end


if DEBUG_LEVEL > 0
    fprintf('\tLoading WRF lat/lon... ')
end

wrf_lat = ncread(wrf_info(1).Filename, 'XLAT');
wrf_lon = ncread(wrf_info(1).Filename, 'XLONG');

if DEBUG_LEVEL > 0;
    fprintf('Done\n');
end

if DEBUG_LEVEL > 0
    fprintf('\tInterpolating now... ');
end

wrf_no2 = griddata(double(wrf_lon(:)), double(wrf_lat(:)), double(wrf_no2(:)), longrid, latgrid);

if DEBUG_LEVEL > 0;
    fprintf('Done\n');
end
    
col_diff = wrf_no2 - behr_no2;

end

