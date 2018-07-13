function [ varargout ] = read_wrf_preproc( wrf_file, quantity, varargin )
%READ_WRF_PREPROC Helper function to read sometimes preprocessed WRF variables
%   
%   Certain common physical quantities (like temperature, pressure, and
%   elevation) can be computed from WRF output, but are not given directly.
%   I have previously written code in the NCL programming language that can
%   compute these quantities when subsetting files, but not all WRF files
%   that we need to read have those quantities. As a result, I ended up
%   with a situation that I have to write lots of checking code every time
%   I want to read in these quantities from an arbitrary WRF file.
%
%   This function will correct that by automatically reading the
%   preprocessed variable for a given quantity if it exists, or reading the
%   precursor variables and converting them to the desired quantity.
%
%   VALUE = READ_WRF_PREPROC( WRF_FILE, QUANTITY ) will read QUANTITY from
%   the WRF_FILE and return the VALUE requested. Special strings for
%   QUANTITY are (case insensitive):
%       - 't', 'temp', or 'temperature': gives temperature in Kelvin
%       - 'p', 'pres', or 'pressure': gives pressure in hPa
%       - 'z', 'elev', 'elevation': gives elevation above sea level in
%       meters. Will also return box height in meters as the second output.
%       - 'ndens', 'number density': gives number density of air in
%       molec./cm^3.
%       - 'no2_ndens': gives NO2 concentration in molec./cm^3
%   Any other string for QUANTITY will be assumed to be a variable in the
%   WRF file and read directly.
%
%   VALUE = READ_WRF_PREPROC( ___, START, COUNT )
%   VALUE = READ_WRF_PREPROC( ___, START, COUNT, STRIDE ) will pass START,
%   COUNT, and (optionally) STRIDE to the various NCREAD() calls, allowing
%   you to only read a subset of the data.

wi = ncinfo(wrf_file);
wrf_vars = {wi.Variables.Name};
if any(strcmpi({'t','temp','temperature'}, quantity))
    if ismember('TT', wrf_vars)
        varargout{1} = ncread(wrf_file, 'TT', varargin{:}); % assume T is always in Kelvin
    else
        wrf_T = ncread(wrf_file, 'T', varargin{:});
        wrf_P = ncread(wrf_file, 'P', varargin{:});
        wrf_PB = ncread(wrf_file, 'PB', varargin{:});
        varargout{1} = convert_wrf_temperature(wrf_T, wrf_P, wrf_PB);
    end
elseif any(strcmpi({'p','pres','pressure'}, quantity))
    if ismember('pres', wrf_vars)
        pres = ncread(wrf_file, 'pres');
        pres_units = ncreadatt(wrf_file, 'pres', 'units');
        varargout{1} = convert_units(pres, pres_units, 'hPa');
    else
        wrf_P = ncread(wrf_file, 'P', varargin{:});
        wrf_PB = ncread(wrf_file, 'PB', varargin{:});
        p_units = ncreadatt(wrf_file, 'P', 'units');
        pb_units = ncreadatt(wrf_file, 'PB', 'units');
        wrf_P = convert_units(wrf_P, p_units, 'hPa');
        wrf_PB = convert_units(wrf_PB, pb_units, 'hPa');
        varargout{1} = wrf_P + wrf_PB;
    end
elseif any(strcmpi({'z','elev','elevation'}, quantity))
    if ismember('z', wrf_vars)
        varargout{1} = ncread(wrf_file, 'z', varargin{:});
        varargout{2} = ncread(wrf_file, 'zlev', varargin{:});
    else
        z = calculate_wrf_altitude(wrf_file, varargin{:});
        zlev = diff(z,1,3);
        varargout = {z, zlev};
    end
elseif any(strcmpi({'ndens','number density'}, quantity))
    if ismember('ndens', wrf_vars)
        varargout{1} = ncread(wrf_file, 'ndens', varargin{:});
    else
        wrf_T = ncread(wrf_file, 'T', varargin{:});
        wrf_P = ncread(wrf_file, 'P', varargin{:});
        wrf_PB = ncread(wrf_file, 'PB', varargin{:});
        varargout{1} = calculate_wrf_air_ndens(wrf_T, wrf_P, wrf_PB);
    end
elseif any(strcmpi({'no2_ndens'}, quantity))
    if ismember('no2_ndens', wrf_vars)
        varargout{1} = ncread(wrf_file, 'no2_ndens', varargin{:});
    else
        wrf_T = ncread(wrf_file, 'T', varargin{:});
        wrf_P = ncread(wrf_file, 'P', varargin{:});
        wrf_PB = ncread(wrf_file, 'PB', varargin{:});
        ndens_air = calculate_wrf_air_ndens(wrf_T, wrf_P, wrf_PB);
        
        no2_mixing_ratio = ncread(wrf_file, 'no2', varargin{:});
        no2_units = ncreadatt(wrf_file, 'no2', 'units');
        no2_mixing_ratio = convert_units(no2_mixing_ratio, no2_units, 'ppp');
        varargout{1} = no2_mixing_ratio .* ndens_air;
    end
else
    varargout{1} = ncread(wrf_file, quantity, varargin{:});
end


end

