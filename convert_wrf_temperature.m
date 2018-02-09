function [ T ] = convert_wrf_temperature( T, varargin )
%T = CONVERT_WRF_TEMPERATURE( WRFFILE )
%   WRF outputs temperature as perturbation to potential temperature. This
%   is not very convenient, so it must be converted to get the absolute
%   temperature. This form takes the path to a WRF output file as its sole
%   input and returns the absolute temperature in K.
%
%   T = CONVERT_WRF_TEMPERATURE( T, P, PB ) takes the arrays containing the
%   WRF variables T, P, and PB as inputs directly. Must be in units of
%   Kelvin, Pascals, and Pascals, respectively.
%
%   See http://mailman.ucar.edu/pipermail/wrf-users/2010/001896.html

p = inputParser;
p.addOptional('P',[]);
p.addOptional('PB', []);
p.addParameter('err_if_missing_units', true);

p.parse(varargin{:});
pout = p.Results;

P = pout.P;
PB = pout.PB;
err_if_missing_units = pout.err_if_missing_units;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

% First check if a path to a WRF file is given. If so, read the data from
% the WRF file.
if ischar(T)
    if ~exist(T,'file')
        E.badinput('%s does not exist', T);
    end
    try 
        wi = ncinfo(T);
    catch err
        if ismember(err.identifier,{'MATLAB:imagesci:netcdf:unableToOpenFileforRead','MATLAB:imagesci:netcdf:libraryFailure'})
            E.badinput('%s does not appear to be a valid netCDF file',T)
        else
            rethrow(err);
        end
    end
    T = ncread(wi.Filename,'T');
    P = ncread(wi.Filename,'P');
    PB = ncread(wi.Filename,'PB');
    
    T_units = ncreadatt_default(wi.Filename,'T','units','K','fatal_if_missing', err_if_missing_units);
    P_units = ncreadatt_default(wi.Filename,'P','units','Pa','fatal_if_missing', err_if_missing_units);
    PB_units = ncreadatt_default(wi.Filename,'PB','units','Pa','fatal_if_missing', err_if_missing_units);
    
    if ~strcmp(T_units, 'K')
        E.notimplemented('Cannot handle temperature not in units of Kelvin');
    end
    P = convert_units(P, P_units, 'Pa');
    PB = convert_units(PB, PB_units, 'Pa');
    
% If T is numeric, then all three of P, PB, and T must have been given as
% inputs.
elseif isnumeric(T)
    if isempty(P) || isempty(PB)
        E.badinput('All three of T, P, and PB must be passed if giving the arrays directly')
    elseif ndims(T) ~= ndims(P) || ndims(T) ~= ndims(PB) || any(size(T) ~= size(P)) || any(size(T) ~= size(PB))
        E.badinput('T, P, and PB must be the same size')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Add an offset of 300 K to the value of T since it is defined as a
% perturbation to potential temperature. See
% http://mailman.ucar.edu/pipermail/wrf-users/2013/003117.html for why we
% add 300 even though the base state temperature is always 290.

T = T+300;

% Now convert from potential temperature to absolute temperature
T = T .* ((P+PB) ./ 1e5) .^ 0.2865;
end

