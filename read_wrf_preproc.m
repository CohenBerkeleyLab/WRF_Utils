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

p = advInputParser;
p.addParameter('error_if_missing_units', true);
p.KeepUnmatched = true;

p.parse(varargin{:});
pout = p.Results;

error_if_missing_units = pout.error_if_missing_units;
read_args = p.Unmatched;

wi = ncinfo(wrf_file);
wrf_vars = {wi.Variables.Name};
if any(strcmpi({'t','temp','temperature'}, quantity))
    if ismember('TT', wrf_vars)
        varargout{1} = ncread(wrf_file, 'TT', read_args{:}); % assume T is always in Kelvin
        % My convert_units function isn't set up to handle additive
        % offsets, so it's difficult at the moment to convert between other
        % units of temperature and Kelvin, so just verify that temperature
        % is given in Kelvin
        check_units('TT', 'K');
    else
        wrf_T = ncread(wrf_file, 'T', read_args{:});
        wrf_P = ncread(wrf_file, 'P', read_args{:});
        wrf_PB = ncread(wrf_file, 'PB', read_args{:});
        % convert_wrf_temperature assumes units of K, Pa, and Pa
        % respectively. Again, converting temperature isn't possible right
        % now. We could convert pressure to Pa, but if this is a regular
        % WRF file, there's no reason to expect that pressure and
        % base-state pressure won't be in Pa. This is just a double check.
        check_units('T','K');
        check_units('P','Pa');
        check_units('PB','Pa');
        varargout{1} = convert_wrf_temperature(wrf_T, wrf_P, wrf_PB);
    end
elseif any(strcmpi({'p','pres','pressure'}, quantity))
    % We want to return pressure in hPa, so since we need to convert the
    % regular P+PB value, we'll just always convert the units here.
    if ismember('pres', wrf_vars)
        pres = ncread(wrf_file, 'pres');
        pres_units = ncreadatt_default(wrf_file, 'pres', 'units', 'hPa', 'fatal_if_missing', error_if_missing_units);
        varargout{1} = convert_units(pres, pres_units, 'hPa');
    else
        wrf_P = ncread(wrf_file, 'P', read_args{:});
        wrf_PB = ncread(wrf_file, 'PB', read_args{:});
        p_units = ncreadatt_default(wrf_file, 'P', 'units','Pa','fatal_if_missing',error_if_missing_units);
        pb_units = ncreadatt_default(wrf_file, 'PB', 'units','Pa','fatal_if_missing',error_if_missing_units);
        wrf_P = convert_units(wrf_P, p_units, 'hPa');
        wrf_PB = convert_units(wrf_PB, pb_units, 'hPa');
        varargout{1} = wrf_P + wrf_PB;
    end
elseif any(strcmpi({'z','elev','elevation'}, quantity))
    % We could convert these from whatever to meters, but again there's no
    % reason to assume that they will be anything other than meters.
    if ismember('z', wrf_vars)
        check_units('z','m');
        check_units('zlev','m');
        varargout{1} = ncread(wrf_file, 'z', read_args{:});
        varargout{2} = ncread(wrf_file, 'zlev', read_args{:});
    else
        z = calculate_wrf_altitude(wrf_file, read_args{:});
        zlev = diff(z,1,3);
        varargout = {z, zlev};
    end
elseif any(strcmpi({'z_center'}, quantity))
    if ismember('z_center', wrf_vars)
        check_units('z_center', 'm');
        varargout{1} = ncread(wrf_file, 'z_center', read_args{:});
    else
        z = read_wrf_preproc(wrf_file, 'z', varargin{:});
        % Verify that the third dimension of the z array is
        % "bottom_top_stag", otherwise this calculation is wrong
        bottom_top_stag = get_ncdf_dim(wrf_file, 'bottom_top_stag');
        if isempty(bottom_top_stag)
            E.callError('wrong_dim_length', 'There is no bottom_top_stag dimension in the file "%s", therefore cannot verify that z is the proper dimensions to calculate z_center');
        elseif size(z,3) ~= bottom_top_stag.Length
            E.callError('wrong_dim_length', 'z was not returned with the expected third dimension (bottom_top_stag), cannot compute %s', quantity);
        end
        
        varargout{1} = (z(:,:,1:end-1,:)+z(:,:,2:end,:))/2;
    end
elseif any(strcmpi({'ndens','number density'}, quantity))
    if ismember('ndens', wrf_vars)
        check_units('ndens','molec./cm^3');
        varargout{1} = ncread(wrf_file, 'ndens', read_args{:});
    else
        varargout{1} = calculate_wrf_air_ndens(wrf_file);
    end
elseif any(strcmpi({'no2_ndens'}, quantity))
    if ismember('no2_ndens', wrf_vars)
        check_units('no2_ndens', 'molec./cm^3');
        varargout{1} = ncread(wrf_file, 'no2_ndens', read_args{:});
    else
        ndens_air = calculate_wrf_air_ndens(wrf_file);
        
        no2_mixing_ratio = ncread(wrf_file, 'no2', varargin{:});
        no2_units = ncreadatt(wrf_file, 'no2', 'units');
        no2_mixing_ratio = convert_units(no2_mixing_ratio, no2_units, 'ppp');
        varargout{1} = no2_mixing_ratio .* ndens_air;
    end
else
    varargout{1} = ncread(wrf_file, quantity, varargin{:});
end


    function check_units(wrf_var, expected_unit)
        try
            unit_in_file = ncreadatt(wrf_file, wrf_var, 'units');
        catch err
            if strcmp(err.identifier, 'MATLAB:imagesci:netcdf:libraryFailure')
                if error_if_missing_units
                    error('read_wrf_preproc:check_units:no_units', 'No "units" attribute for variable "%s" in %s', wrf_var, wrf_file);
                else
                    unit_in_file = expected_unit;
                end
            else
                rethrow(err)
            end
        end
        
        if ~strcmp(unit_in_file, expected_unit)
            error('read_wrf_preproc:check_units:wrong_unit', 'Expected units of "%s" for "%s" in %s; got "%s" instead', expected_unit, wrf_var, wrf_file, unit_in_file);
        end
        
    end

end

function dim = get_ncdf_dim(filename, dim_name)
info = ncinfo(filename);
all_dims = {info.Dimensions.Name};
xx = strcmp(all_dims, dim_name);
if sum(xx) == 1
    dim = info.Dimensions(xx);
else
    dim = [];
end
end
