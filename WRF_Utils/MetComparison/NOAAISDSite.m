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
        % not modified unexpectedly or the necessary validation can be done
        obs_datenums = [];
        wind_dir_def = 1;
    end
    properties (SetAccess = private, GetAccess = public)
        % This should be set by the relevant methods to ensure they are not
        % modified unexpectedly and access by the relevant methods for
        % quality filtering
        wind_dir = [];
        wind_dir_quality = {};
        
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
        
        function set_wind_def(obj, def_in)
            allowed_defs = {'met','vector','uv'};
            if ischar(def_in) 
                if ~ismember(lower(def_in), allowed_defs);
                    error('NOAAISDSite:bad_input', 'If giving the wind direction definition as a string, must be one of %s', strjoin(allowed_defs, ', '));
                end
            elseif isnumeric(def_in) 
                if def_in < 0 || def_in > numel(allowed_defs)
                    error('NOAAISDSite:bad_input', 'If giving the wind direction definition as a number, must be between 1 and %d', numel(allowed_defs));
                end
            elseif ~isnumeric(def_in) && ~ischar(def_in)
                error('NOAAISDSite:bad_input', 'The wind direction definition must be number or string');
            end
            if ischar(def_in)
                def = find(strcmpi(def_in, allowed_defs));
            else
                def = def_in;
            end
            obj.wind_dir_def = def;
        end
        
        %%%% GETTERS %%%%
        function [varargout] = get_wind_nearest_in_time(obj, dnum, max_time_sep, time_sep_units, quality_level)
            xx = obj.nearest_obs_in_time(dnum, max_time_sep, time_sep_units);
            if sum(xx) < 1
                varargout(1:2) = {[]};
                return
            end
            [w1, w2] = obj.get_wind(quality_level);
            varargout{1} = w1(xx);
            varargout{2} = w2(xx);
        end
        
        function t_out = get_temperature_nearest_in_time(obj, dnum, max_time_sep, time_sep_units, quality_level)
            xx = obj.nearest_obs_in_time(dnum, max_time_sep, time_sep_units);
            if sum(xx) < 1
                t_out = [];
                return
            end
            t_out = obj.get_temperature(quality_level);
            t_out = t_out(xx);
        end
        
        function [varargout] = get_wind(obj, quality_level)
            wind_dir_filt = obj.apply_quality_filtering(obj.wind_dir, obj.wind_dir_quality, quality_level);
            wind_vel_filt = obj.apply_quality_filtering(obj.wind_vel, obj.wind_vel_quality, quality_level);
            
            % There are 3 ways to define winds. Two can use a direction
            % given as an angle and a velocity, the other defines the
            % vector components. 
            switch obj.wind_dir_def
                case 1 
                    % native meteorology definition: 0 deg means the wind
                    % comes out of the north, 90 deg means the wind comes
                    % out of the east.
                    varargout{1} = wind_dir_filt;
                    varargout{2} = wind_vel_filt;
                    return
                case 2 
                    % definition of direction consitent with vector math: 0
                    % deg means the wind blows TO the east, 90 deg means
                    % the wind blows TO the north. This is like if a
                    % cartesian coordinate plane was laid over the Earth
                    % with east as +x and north as +y.
                    %
                    % The transformation between the two is x = -y
                    xmet = cosd(wind_dir_filt);
                    ymet = sind(wind_dir_filt);
                    xvec = -ymet;
                    yvec = -xmet;
                    
                    varargout{1} = atan2d(yvec, xvec);
                    varargout{2} = wind_vel_filt;
                    
                case 3 
                    % u and v components. To use the typical mathematical
                    % expressions for projection onto the x and y axis, we
                    % need to apply the x = -y transform to the
                    % meteorological wind direction, hence why U
                    % (x-component) = -sin(theta) and V = -cos(theta).
                    if nargout < 2
                        warning('Returning U and V wind components but only 1 output variable');
                    end
                    U = -sind(wind_dir_filt) .* wind_vel_filt;
                    V = -cosd(wind_dir_filt) .* wind_vel_filt;
                    varargout{1} = U;
                    varargout{2} = V;
            end
        end
        
        function temperature = get_temperature(obj, quality_level)
            temperature = obj.apply_quality_filtering(obj.temperature, obj.temperature_quality, quality_level);
        end
    end
    
    methods ( Access = private )
        function xx = nearest_obs_in_time(obj, dnum, max_time_sep, time_sep_units)
            % Convert maxdiff into days
            allowed_units = {'seconds', 'minutes', 'hours', 'days'};
            if ~ismember(time_sep_units, allowed_units)
                error('NOAAISDSite:bad_input', 'TIME_SEP_UNITS must be one of %s', strjoin(time_sep_units, ', '));
            end
            switch time_sep_units
                case 'seconds'
                    conv = 60*60*24;
                case 'minutes'
                    conv = 60*24;
                case 'hours'
                    conv = 24;
                case 'days'
                    conv = 1;
            end
            max_time_sep = max_time_sep / conv;
            
            [~,xx] = min(abs(obj.obs_datenums - dnum));
            if abs(obj.obs_datenums(xx) - dnum) > max_time_sep
                xx = false(size(xx));
            end
        end
    end
    
    methods ( Access = private, Static )
        function vals = apply_quality_filtering(vals, qual_cell, filter_type)
            switch filter_type
                case 1
                    xx = cellfun(@(x) ~strcmp(x,'1') & ~strcmp(x,'5'), qual_cell);
                otherwise
                    xx = false(size(vals));
            end
            vals(xx) = nan;
        end
    end
    
end
