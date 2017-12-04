function [ trop_no2 ] = compute_wrf_trop_columns( wrf_filename )
%COMPUTE_WRF_TROP_COLUMNS Calculate WRF-Chem tropospheric NO2 columns
%   TROP_NO2 = COMPUTE_WRF_TROP_COLUMNS( WRF_FILENAME ) calculates NO2
%   tropospheric columns in molec./cm^3 for WRF_FILENAME (a path to a WRF
%   file that has the variables no2_ndens and zlev - currently this
%   restricts it to files preprocessed with slurmrun/run_wrf_output.sh in
%   the WRF-nco-tools repo). Integrates up to the tropopause calculated by
%   FIND_WRF_TROPOPAUSE().

E = JLLErrors;
wi = ncinfo(wrf_filename);
if any(~ismember({'no2_ndens','zlev'},{wi.Variables.Name}))
    E.badinput('WRF file must contain the variables no2_ndens and zlev. These are created by (slurm)run_wrf_output.sh');
end

daily_no2 = ncread(wrf_filename, 'no2_ndens'); %[NO2 in number density]
daily_zlev = ncread(wrf_filename, 'zlev'); % Thickness of each layer in meters
[daily_tplev] = find_wrf_tropopause(wi);

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

