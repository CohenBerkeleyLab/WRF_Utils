classdef match_wrf2campaigns
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)
        
    end
    
    methods(Static = true)
        function Match = discover_md(prof_mode, wrf_dirs, varargin)
            [utc_range, lst_range] = match_wrf2campaigns.parse_common_args(varargin{:});
            
            campaign_name = 'discover_md';
            if ~exist('prof_mode', 'var') || isempty(prof_mode)
                prof_mode = 'daily';
            end
            if ~exist('wrf_dirs', 'var') || isempty(wrf_dirs)
                wrf_dirs = {'/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2011/07'};
            end
            
            Out = campaign_wide_ops(campaign_name, {'NO2_LIF', 'LONGITUDE', 'LATITUDE', 'PRESSURE'}, 'cat', 'datefmt','datenum');
            Raw.no2 = Out.data.NO2_LIF .* 1e-12 .* 1e6;
            Raw.lon = Out.data.LONGITUDE;
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.PRESSURE;
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            Raw.campaign = campaign_name;
            
            Raw = match_wrf2campaigns.restrict_by_time(Raw, Out.utcs, utc_range, lst_range);
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
        
        function Match = discover_ca(prof_mode, wrf_dirs, varargin)
            [utc_range, lst_range] = match_wrf2campaigns.parse_common_args(varargin{:});
            
            campaign_name = 'discover_ca';
            if ~exist('prof_mode', 'var') || isempty(prof_mode)
                prof_mode = 'daily';
            end
            if ~exist('wrf_dirs', 'var') || isempty(wrf_dirs)
                wrf_dirs = {'/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2013/01', '/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2013/02'};
            end
            
            Out = campaign_wide_ops(campaign_name, {'NO2_MixingRatio_LIF', 'LONGITUDE', 'LATITUDE', 'PRESSURE'}, 'cat', 'datefmt','datenum');
            Raw.no2 = Out.data.NO2_MixingRatio_LIF .* 1e-12 .* 1e6;
            Raw.lon = Out.data.LONGITUDE;
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.PRESSURE;
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            Raw.campaign = campaign_name;
            
            Raw = match_wrf2campaigns.restrict_by_time(Raw, Out.utcs, utc_range, lst_range);
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
        
        function Match = discover_tx(prof_mode, wrf_dirs, varargin)
            [utc_range, lst_range] = match_wrf2campaigns.parse_common_args(varargin{:});
            
            campaign_name = 'discover_tx';
            if ~exist('prof_mode', 'var') || isempty(prof_mode)
                prof_mode = 'daily';
            end
            if ~exist('wrf_dirs', 'var') || isempty(wrf_dirs)
                wrf_dirs = {'/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2013/09'};
            end
            
            Out = campaign_wide_ops(campaign_name, {'NO2_MixingRatio_LIF', 'LONGITUDE', 'LATITUDE', 'PRESSURE'}, 'cat', 'datefmt','datenum');
            Raw.no2 = Out.data.NO2_MixingRatio_LIF .* 1e-12 .* 1e6;
            Raw.lon = Out.data.LONGITUDE;
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.PRESSURE;
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            Raw.campaign = campaign_name;
            
            Raw = match_wrf2campaigns.restrict_by_time(Raw, Out.utcs, utc_range, lst_range);
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
        
        function Match = discover_co(prof_mode, wrf_dirs, varargin)
            [utc_range, lst_range] = match_wrf2campaigns.parse_common_args(varargin{:});
            
            campaign_name = 'discover_co';
            if ~exist('prof_mode', 'var') || isempty(prof_mode)
                prof_mode = 'daily';
            end
            if ~exist('wrf_dirs', 'var') || isempty(wrf_dirs)
                wrf_dirs = {'/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2014/07','/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2014/08'};
            end
            
            Out = campaign_wide_ops(campaign_name, {'NO2_LIF', 'LONGITUDE', 'LATITUDE', 'PRESSURE'}, 'cat', 'datefmt','datenum');
            Raw.no2 = Out.data.NO2_LIF .* 1e-12 .* 1e6;
            Raw.lon = Out.data.LONGITUDE;
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.PRESSURE;
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            Raw.campaign = campaign_name;
            
            Raw = match_wrf2campaigns.restrict_by_time(Raw, Out.utcs, utc_range, lst_range);
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
        
        function Match = dc3(prof_mode, wrf_dirs, varargin)
            [utc_range, lst_range] = match_wrf2campaigns.parse_common_args(varargin{:});
            
            campaign_name = 'dc3';
            if ~exist('wrf_dirs', 'var') || isempty(wrf_dirs)
                %wrf_dirs = '/Volumes/share2/USERS/LaughnerJ/WRF/DC3/lnox_off-fixed_BCs';
                wrf_dirs = {'/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2012/05','/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2012/06'};
            end
            if ~exist('prof_mode', 'var') || isempty(prof_mode)
                prof_mode = 'daily';
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
            
            Raw = match_wrf2campaigns.restrict_by_time(Raw, Out.utcs, utc_range, lst_range);
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
        
        function Match = soas(prof_mode, wrf_dirs, varargin)
            [utc_range, lst_range] = match_wrf2campaigns.parse_common_args(varargin{:});
            
            campaign_name = 'soas';
            if ~exist('prof_mode', 'var') || isempty(prof_mode)
                prof_mode = 'daily';
            end
            if ~exist('wrf_dirs', 'var') || isempty(wrf_dirs)
                wrf_dirs = {'/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2013/05','/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2013/05','/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2013/07'};
            end
            
            Out = campaign_wide_ops(campaign_name, {'NO2_ppbv', 'LONGITUDE', 'LATITUDE', 'StaticPrs'}, 'cat', 'datefmt','datenum');
            
            Raw.campaign = campaign_name;
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            Raw.no2 = Out.data.NO2_ppbv * 1e-9 * 1e6;
            Raw.lon = Out.data.LONGITUDE;
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.StaticPrs;
            
            Raw = match_wrf2campaigns.restrict_by_time(Raw, Out.utcs, utc_range, lst_range);
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
        
        function Match = seac4rs(prof_mode, wrf_dirs, varargin)
            [utc_range, lst_range] = match_wrf2campaigns.parse_common_args(varargin{:});
            
            campaign_name = 'seac4rs';
            if ~exist('prof_mode', 'var') || isempty(prof_mode)
                prof_mode = 'daily';
            end
            if ~exist('wrf_dirs', 'var') || isempty(wrf_dirs)
                wrf_dirs = {'/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2013/08','/Volumes/share-wrf1/BEHR-WRF/Outputs/us/2013/09'};
            end
            
            Out = campaign_wide_ops(campaign_name, {'NO2_TDLIF', 'LONGITUDE', 'LATITUDE', 'PRESSURE'}, 'cat', 'datefmt','datenum');
            
            Raw.campaign = campaign_name;
            Raw.dvec = match_wrf2campaigns.convert_aircraft_times(Out);
            Raw.no2 = Out.data.NO2_TDLIF * 1e-12 * 1e6;
            Raw.lon = Out.data.LONGITUDE;
            Raw.lat = Out.data.LATITUDE;
            Raw.pres = Out.data.PRESSURE;
            
            Raw = match_wrf2campaigns.restrict_by_time(Raw, Out.utcs, utc_range, lst_range);
            
            Match = match_wrf2aircraft(Raw, wrf_dirs, prof_mode);
        end
    end
    
    methods(Access = private, Static = true)
        function dvec = convert_aircraft_times(Out)
            dvec = Out.dates + Out.utcs ./ (60*60*24);
        end
        
        function varargout = parse_common_args(varargin)
            p = inputParser;
            p.addParameter('utc_range', []);
            p.addParameter('lst_range', []);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            varargout{1} = pout.utc_range;
            varargout{2} = pout.lst_range;
        end
        
        function Raw = restrict_by_time(Raw, utcs, utc_range, lst_range)
            xx_in_range = true(size(utcs));
            if ~isempty(utc_range)
                xx_in_range = xx_in_range & match_wrf2campaigns.time_range_logical(utcs, utc_range);
            end
            if ~isempty(lst_range)
                xx_in_range = xx_in_range & match_wrf2campaigns.time_range_logical(utc2local_sec(utcs, round(Raw.lon/15)), lst_range);
            end
            
            fns = fieldnames(Raw);
            for f = 1:numel(fns)
                if isnumeric(Raw.(fns{f}))
                    Raw.(fns{f}) = Raw.(fns{f})(xx_in_range);
                end
            end
        end
        
        function xx = time_range_logical(sec_after_midnight, time_range)
            % Assume the range is given in hours after midnight. Convert to
            % seconds.
            range_start = time_range(1) * 60^2;
            range_end = time_range(2) * 60^2;
            xx = sec_after_midnight >= range_start & sec_after_midnight <= range_end;
        end
    end
    
end

