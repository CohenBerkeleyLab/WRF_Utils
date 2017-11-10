function [ trop_no2 ] = compute_wrf_trop_columns( wrf_filename, tropopause )
%COMPUTE_WRF_TROP_COLUMNS Calculate WRF-Chem tropospheric NO2 columns
%   TROP_NO2 = COMPUTE_WRF_TROP_COLUMNS( WRF_FILENAME ) calculates NO2
%   tropospheric columns in molec./cm^3 for WRF_FILENAME (a path to a WRF
%   file that has the variables no2_ndens and zlev - currently this
%   restricts it to files preprocessed with slurmrun/run_wrf_output.sh in
%   the WRF-nco-tools repo). Integrates up to the tropopause calculated by
%   FIND_WRF_TROPOPAUSE().
%
%   TROP_NO2 = COMPUTE_WRF_TROP_COLUMNS( WRF_FILENAME, TROPOPAUSE )
%   integrates up to TROPOPAUSE (given in hPa) instead.

E = JLLErrors;
if ~exist('tropopause','var')
    tropopause = NaN;
elseif ~isnumeric(tropopause) || ~isscalar(tropopause) || tropopause < 0
    E.badinput('TROPOPAUSE must be a positive, scalar number')
end

wi = ncinfo(wrf_filename);
if any(~ismember({'no2', 'T', 'PH', 'PHB', 'P', 'PB'},{wi.Variables.Name}))
    E.badinput('WRF file must contain the variables no2_ndens and zlev. These are created by (slurm)run_wrf_output.sh');
end

wrf_vars = {wi.Variables.Name};
if ismember('no2_ndens', wrf_vars)
    daily_no2 = ncread(wrf_filename, 'no2_ndens'); %[NO2 in number density]
else
    wrf_ndens = calculate_wrf_air_ndens(wrf_filename);
    wrf_no2_ppm = ncread(wrf_filename, 'no2');
    daily_no2 = wrf_no2_ppm .* 1e-6 .* wrf_ndens;
end

if ismember('zlev', wrf_vars)
    daily_zlev = ncread(wrf_filename, 'zlev'); % Thickness of each layer in meters
else
    daily_z = calculate_wrf_altitude(wrf_filename);
    daily_zlev = daily_z(:,:,2:end,:) - daily_z(:,:,1:end-1,:);
end

if isnan(tropopause)
    daily_tplev = find_wrf_tropopause(wi);
else
    daily_tplev = find_fixed_tropopause(wi, tropopause);
end

for a=1:size(daily_no2,1)
    for b=1:size(daily_no2,2)
        for t=1:size(daily_no2,4)
            tp = daily_tplev(a,b,t);
            if tp > 0 % tp is given -1 if the tropopause algorithm cannot find a tropopause
                daily_no2(a,b,tp:end,t) = nan;
            end
        end
    end
end

trop_no2 = squeeze(nansum2(daily_no2 .* (daily_zlev*100), 3));

end

function tplev = find_fixed_tropopause(wi, tropopause)
pres = ncread(wi.Filename, 'P') + ncread(wi.Filename, 'PB');
p_units = ncreadatt(wi.Filename, 'P', 'units');
pb_units = ncreadatt(wi.Filename, 'P', 'units');
if ~strcmpi(p_units, pb_units)
    E.notimplemented('P and PB are in different units');
end
pres = convert_units(pres, p_units, 'hPa');
sz = get_wrf_array_size(wi.Filename);
tplev = nan(sz(1:2));
for a=1:sz(1)
    for b=1:sz(2)
        tplev = find(pres(a,b,:) > tropopause, 1, 'last');
    end
end
end

