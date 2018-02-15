classdef match_wrf2campaigns
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static = true)
        function Match = discover_md(prof_mode, wrf_dirs)
            campaign_name = 'discover_md';
            if ~exist('prof_mode', 'var')
                prof_mode = 'daily';
            end
            if ~exist('wrf_dirs', 'var')
                wrf_dirs = {'/Volumes/share-wrf1/Outputs/us/2011/07'};
            end
            
            Out = campaign_wide_ops(campaign_name, {'NO2_LIF', 'LONGITUDE', 'LATITUDE', 'PRESSURE'}, 'cat', 'datefmt','datenum');
            Raw.no2 = Out.data.NO2_LIF .* 1e-12 .* 1e6;
            Raw.lon = Out.data.LONGITUDE;
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.PRESSURE;
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            Raw.campaign = campaign_name;
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
        
        function Match = discover_ca(prof_mode, wrf_dirs)
            campaign_name = 'discover_ca';
            if ~exist('prof_mode', 'var')
                prof_mode = 'daily';
            end
            if ~exist('wrf_dirs', 'var')
                wrf_dirs = '';
            end
            
            Out = campaign_wide_ops(campaign_name, {'NO2_MixingRatio_LIF', 'LONGITUDE', 'LATITUDE', 'PRESSURE'}, 'cat', 'datefmt','datenum');
            Raw.no2 = Out.data.NO2_MixingRatio_LIF .* 1e-12 .* 1e6;
            Raw.lon = Out.data.LONGITUDE;
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.PRESSURE;
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            Raw.campaign = campaign_name;
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
        
        function Match = discover_tx(prof_mode, wrf_dirs)
            campaign_name = 'discover_tx';
            if ~exist('prof_mode', 'var')
                prof_mode = 'daily';
            end
            if ~exist('wrf_dirs', 'var')
                wrf_dirs = '';
            end
            
            Out = campaign_wide_ops(campaign_name, {'NO2_MixingRatio_LIF', 'LONGITUDE', 'LATITUDE', 'PRESSURE'}, 'cat', 'datefmt','datenum');
            Raw.no2 = Out.data.NO2_MixingRatio_LIF .* 1e-12 .* 1e6;
            Raw.lon = Out.data.LONGITUDE;
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.PRESSURE;
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            Raw.campaign = campaign_name;
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
        
        function Match = discover_co(prof_mode, wrf_dirs)
            campaign_name = 'discover_co';
            if ~exist('prof_mode', 'var')
                prof_mode = 'daily';
            end
            if ~exist('wrf_dirs', 'var')
                wrf_dirs = '';
            end
            
            Out = campaign_wide_ops(campaign_name, {'NO2_LIF', 'LONGITUDE', 'LATITUDE', 'PRESSURE'}, 'cat', 'datefmt','datenum');
            Raw.no2 = Out.data.NO2_LIF .* 1e-12 .* 1e6;
            Raw.lon = Out.data.LONGITUDE;
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.PRESSURE;
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            Raw.campaign = campaign_name;
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
        
        function Match = dc3(prof_mode, wrf_dirs)
            campaign_name = 'dc3';
            if ~exist('wrf_dirs', 'var')
                wrf_dirs = '/Volumes/share2/USERS/LaughnerJ/WRF/DC3/lnox_off-fixed_BCs';
            end
            if ~exist('prof_mode', 'var')
                prof_mode = 'hourly';
            end
            
            % TODO: modify campaign_wide_ops to handle multiple requested fields
            % Output to structure raw; anything in it will be binned
            Out = campaign_wide_ops(campaign_name, {'no2_lif', 'MPN_TDLIF', 'NO_ESRL', 'HNO3_SAGA', 'HNO3_CIT', 'CO_DACOM', 'O3_ESRL', 'JNO2NOO3P', 'LONGITUDE', 'LATITUDE', 'PRESSURE', 'TEMPERATURE','MixingRatio'}, 'cat', 'datefmt','datenum');
            
            
            % Convert the output chemical species here to the Raw structure, also
            % convert to ppm since that's how WRF outputs concentrations. Field names must
            % match the variable names in WRF-Chem.
            Raw.campaign = campaign_name;
            Raw.no = Out.data.NO_ESRL .* 1e-9 .* 1e6;
            Raw.no2 = Out.data.no2_lif .* 1e-12 .* 1e6;
            Raw.mpn = Out.data.MPN_TDLIF .* 1e-12 .* 1e6;
            Raw.hno3 = (Out.data.HNO3_SAGA + Out.data.HNO3_CIT) .* 1e-12/2 .* 1e6;
            Raw.co = Out.data.CO_DACOM .* 1e-9 .* 1e6;
            Raw.o3 = Out.data.O3_ESRL .* 1e-9 .* 1e6;
            Raw.PHOTR_NO2 = Out.data.JNO2NOO3P * 60; % WRF outputs in per minute
            Raw.lon = Out.data.LONGITUDE; % the correction to negative is west is handled in remove_merge_fills
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.PRESSURE;
            Raw.TT = Out.data.TEMPERATURE;
            Raw.QVAPOR = Out.data.MixingRatio * 1e-3; % 1) I'm assuming that "MixingRatio" is of water, since it comes from a hygrometer
                                                      % 2) The units for this mixing ratio are g/kg, in WRF it is kg/kg.
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
    end
    
    methods(Access = private, Static = true)
        function dvec = convert_aircraft_times(Out)
            dvec = Out.dates + Out.utcs ./ (60*60*24);
        end
    end
    
end

