function [ ndens ] = calculate_wrf_air_ndens( T, P, PB )
%NDENS = CALCULATE_WRF_AIR_NDENS( WRFFILE )
%   Calculate the number density of air from WRF outputs.  This form takes
%   the path to a WRF output file as its sole input and returns the
%   number density of air in molec./cm^3.
%
%   NDENSE = CONVERT_WRF_TEMPERATURE( T, P, PB ) takes the arrays containing the
%   WRF variables T, P, and PB as inputs directly.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

% First check if a path to a WRF file is given
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
elseif isnumeric(T)
    if nargin < 3
        E.badinput('All three of T, P, and PB must be passed if giving the arrays directly')
    elseif ndims(T) ~= ndims(P) || ndims(T) ~= ndims(PB) || any(size(T) ~= size(P)) || any(size(T) ~= size(PB))
        E.badinput('T, P, and PB must be the same size')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

Av = 6.022e23;

% Need absolute temperature
T = convert_wrf_temperature(T, P, PB);

% Now apply ideal gas law
ndens = ((P+PB) .* Av) ./ (8.314 .* T);
% convert from molec./m^3 to molec./cm^3
ndens = ndens*1e-6;
end

