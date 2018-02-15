function [ trop_no2 ] = compute_wrf_trop_columns( wrf_filename, integration_mode, tropopause )
%COMPUTE_WRF_TROP_COLUMNS Calculate WRF-Chem tropospheric NO2 columns
%   TROP_NO2 = COMPUTE_WRF_TROP_COLUMNS( WRF_FILENAME ) calculates NO2
%   tropospheric columns in molec./cm^3 for WRF_FILENAME (a path to a WRF
%   file that has the variables no2_ndens and zlev - currently this
%   restricts it to files preprocessed with slurmrun/run_wrf_output.sh in
%   the WRF-nco-tools repo). Integrates up to the tropopause calculated by
%   FIND_WRF_TROPOPAUSE().
%
%   TROP_NO2 = COMPUTE_WRF_TROP_COLUMNS( WRF_FILENAME, INTEGRATION_MODE )
%   controls whether the integration is done by assuming constant number
%   density over each box ('box', default) or by mixing ratio over pressure
%   using integPr2 from BEHR-core-utils ('mixing_ratio').
%
%   TROP_NO2 = COMPUTE_WRF_TROP_COLUMNS( WRF_FILENAME, INTEGRATION_MODE, TROPOPAUSE )
%   integrates up to TROPOPAUSE (given in hPa) instead.

E = JLLErrors;
allowed_int_modes = {'box', 'mixing_ratio'};
if ~exist('integration_mode','var') || isempty(integration_mode)
    integration_mode = allowed_int_modes{1};
elseif ~ischar(integration_mode) || ~ismember(integration_mode, allowed_int_modes)
    E.badinput('INTEGRATION_MODE must be one of: %s', strjoin(allowed_int_modes, ', '));
end

if ~exist('tropopause','var')
    tropopause = NaN;
elseif ~isnumeric(tropopause) || ~isscalar(tropopause) || tropopause < 0
    E.badinput('TROPOPAUSE must be a positive, scalar number')
end

wi = ncinfo(wrf_filename);
% Rely on error messages from NCREAD() if missing a variable.
wrf_vars = read_file_vars(wrf_filename, integration_mode);

if isnan(tropopause)
    daily_tplev = find_wrf_tropopause(wi);
else
    daily_tplev = find_fixed_tropopause(wrf_vars.pres, tropopause);
end

for a=1:size(wrf_vars.no2,1)
    for b=1:size(wrf_vars.no2,2)
        for t=1:size(wrf_vars.no2,4)
            tp = daily_tplev(a,b,t);
            if tp > 0 % tp is given -1 if the tropopause algorithm cannot find a tropopause
                wrf_vars.no2(a,b,tp:end,t) = nan;
            end
        end
    end
end

if strcmpi(integration_mode,'box')
    trop_no2 = squeeze(nansum2(wrf_vars.no2 .* (wrf_vars.zlev*100), 3));
elseif strcmpi(integration_mode, 'mixing_ratio');
    sz = size(wrf_vars.no2);
    trop_no2 = nan(sz(1:2));
    for a=1:sz(1)
        for b=1:sz(2)
            trop_no2(a,b) = integPr2(squeeze(wrf_vars.no2(a,b,:)), squeeze(wrf_vars.pres(a,b,:)), wrf_vars.pres(a,b,1));
        end
    end
end

end

function tplev = find_fixed_tropopause(pres, tropopause)
sz = size(pres);
tplev = nan(sz(1:2));
for a=1:sz(1)
    for b=1:sz(2)
        tplev(a,b) = find(pres(a,b,:) > tropopause, 1, 'last');
    end
end
end

function vars = read_file_vars(wrf_file, int_mode)
if strcmpi(int_mode, 'box')
    vars.no2 = read_wrf_preproc(wrf_file, 'no2_ndens');
    [~, vars.zlev] = read_wrf_preproc(wrf_file, 'z');
else
    vars.no2 = ncread(wrf_file, 'no2');
end
% Used if finding fixed tropopause, so always include this
vars.pres = read_wrf_preproc(wrf_file, 'pressure');
end
