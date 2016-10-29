classdef NOAAISDSite < handle
    %NOAAISDSite Represents observations from a NOAA integrated surface database site
    %   NOAA ISD sites represent data in an ASCII format specified at
    %   http://www1.ncdc.noaa.gov/pub/data/ish/ish-format-document.pdf.
    %   
    %   Constructors:
    %       NOAAISDSite() - constructs a default object. Only use to
    %       indicate that the site is invalid or unused for whatever
    %       reason.
    %
    %       NOAAISDSite(usaf_id, wban_id, lon, lat, elev) - Constructs the
    %       object with the specified properties. usaf_id should be a
    %       string, the others numeric values.
    %
    %   Public properties: usaf_id, wban_id, lon, lat, elev. Correspond to
    %   those set in the constructor. Cannot be changed later.
    %
    %   Setting data:
    %       set_wind_dir(wind_dir, wind_dir_quality)
    %       set_wind_vel(wind_vel, wind_vel_quality)
    %       set_temperature(temperature, temperature_quality)
    %   In all cases, the quality vectors should be cell arrays, the value
    %   vectors numeric arrays. The data should already be scaled and have
    %   fill values handled.
    %
    %   Getting data:
    %       get_wind_dir(quality_level)
    %       get_wind_vel(quality_level)
    %       get_temperature(quality_level)
    %   Returns the data vectors with the data points failing the quality
    %   check replaced with NaNs.
    
    
    properties (SetAccess = immutable, GetAccess = public)
        lon = nan;
        lat = nan;
        elevation = nan;
        
        usaf_id = nan;
        wban_id = nan;
        
    end
    properties (SetAccess = private, GetAccess = public)
        % These should be set by the relevant methods to ensure they are
        % not modified unexpectedly
        obs_datenums = [];
        
    end
    properties (SetAccess = private, GetAccess = public)
        % This should be set by the relevant methods to ensure they are not
        % modified unexpectedly and access by the relevant methods for
        % quality filtering
        wind_dir = [];
        wind_dir_quality = {};
        wind_dir_def = 1;
        
        wind_vel = [];
        wind_vel_quality = {};
        
        temperature = [];
        temperature_quality = {};
    
    end
    
    methods
        %%%% Constructor %%%%
        function obj = NOAAISDSite(usaf_id, wban_id, lon, lat, elev)
            if nargin == 0
                usaf_id = 'none';
                wban_id = nan;
                lon = nan;
                lat = nan;
                elev = nan;
            end
            if ~ischar(usaf_id)
                error('NOAAISDSite:bad_input','USAF_ID must be a string');
            end
            if ~isnumeric(wban_id) || ~isscalar(wban_id)
                error('NOAAISDSite:bad_input', 'WBAN_ID must be a scalar number')
            end
            if ~isnumeric(lon) || ~isscalar(lon)
                error('NOAAISDSite:bad_input', 'RAW_LON must be a scalar number')
            end
            if ~isnumeric(lat) || ~isscalar(lat)
                error('NOAAISDSite:bad_input', 'RAW_LAT must be a scalar number')
            end
            if ~isnumeric(elev) || ~isscalar(elev)
                error('NOAAISDSite:bad_input', 'RAW_ELEV must be a scalar number')
            end
            obj.usaf_id = usaf_id;
            obj.wban_id = wban_id;
            obj.lon = lon;
            obj.lat = lat;
            obj.elevation = elev;
        end
        
        %%%% SETTERS %%%%
        function set_wind_dir(obj, wind_dir_in, wind_dir_quality_in)
            if ~isnumeric(wind_dir_in)
                error('NOAAISDSite:bad_input', 'WIND_DIR_IN must be a numeric array');
            end
            if ~iscellstr(wind_dir_quality_in)
                error('NOAAISDSite:bad_input', 'WIND_DIR_QUALITY_IN must be a cell array of strings')
            end
            obj.wind_dir = wind_dir_in;
            obj.wind_dir_quality = wind_dir_quality_in;
        end
        function set_wind_vel(obj, wind_vel_in, wind_vel_quality_in)
            if ~isnumeric(wind_vel_in)
                error('NOAAISDSite:bad_input', 'WIND_VEL_IN must be a numeric array');
            end
            if ~iscellstr(wind_vel_quality_in)
                error('NOAAISDSite:bad_input', 'WIND_VEL_QUALITY_IN must be a cell array of strings')
            end
            obj.wind_vel = wind_vel_in;
            obj.wind_vel_quality = wind_vel_quality_in;
        end
        function set_temperature(obj, temperature_in, temperature_quality_in)
            if ~isnumeric(temperature_in)
                error('NOAAISDSite:bad_input', 'TEMPERATURE_IN must be a numeric array');
            end
            if ~iscellstr(temperature_quality_in)
                error('NOAAISDSite:bad_input', 'TEMPERATURE_QUALITY_IN must be a cell array of strings')
            end
            obj.temperature = temperature_in;
            obj.temperature_quality = temperature_quality_in;
        end
        function set_dates(obj, dates_in)
            obj.obs_datenums = dates_in;
        end
        
        %%%% GETTERS %%%%
        function wind_dir = get_wind_dir(obj, quality_level)
            wind_dir = obj.apply_quality_filtering(obj.wind_dir, obj.wind_dir_quality, quality_level);
        end
        function wind_vel = get_wind_vel(obj, quality_level)
            wind_vel = obj.apply_quality_filtering(obj.wind_vel, obj.wind_vel_quality, quality_level);
        end
        function temperature = get_temperature(obj, quality_level)
            temperature = obj.apply_quality_filtering(obj.temperature, obj.temperature_quality, quality_level);
        end
    end
    
    methods ( Access = private, Static )
        function vals = apply_quality_filtering(vals, qual_cell, filter_type)
            % Until we know how to filter this, just return the value array
            return
        end
    end
    
end
