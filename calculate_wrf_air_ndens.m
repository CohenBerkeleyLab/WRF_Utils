function [ ndens ] = calculate_wrf_air_ndens( T, P, PB )
%NDENS = CALCULATE_WRF_AIR_NDENS( WRFFILE )
%   Calculate the number density of air from WRF outputs.  This form takes
%   the path to a WRF output file as its sole input and returns the
%   number density of air in molec./cm^3.
%
%   NDENSE = CONVERT_WRF_TEMPERATURE( TEMP, PRES ) takes arrays of temperature in K
%   and pressure in hPa. THIS ASSUMES THAT TEMP HAS ALREADY BEEN CONVERTED FROM 
%   PERTURBATION POTENTIAL TEMPERATURE.
%
%   NDENSE = CONVERT_WRF_TEMPERATURE( T, P, PB ) takes the arrays containing the
%   WRF variables T, P, and PB as inputs directly. NOTE: in this form T is assumed
%   to be the perturbation potential temperature and P and PB assumed to be in Pa.
%   THE UNITS FOR P AND PB ARE DIFFERENT THAN IN THE TWO ARGUMENT FORM.

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
    temp = read_wrf_preproc(wi.Filename, 'temperature');
    % read_wrf_preproc returns pressure in hPa, need Pa
    pres = 100*read_wrf_preproc(wi.Filename, 'pressure');
elseif isnumeric(T)
    if nargin == 2
        if ~isequal(size(T), size(P))
            E.badinput('T and P must be the same size')
        end
        % Assume temperature is given in Kelvin and pressure in hPa
        temp = T;
        pres = P*100;
    elseif nargin == 3
        if ~isequal(size(T), size(P)) || ~isequal(size(T), size(PB))
            E.badinput('T, P, and PB must be the same size')
        end
        % For the 3 input form, assume T is the WRF perturbation potential
        % temperature and P and PB are the perturbation and base pressure,
        % in Pa
        temp = convert_wrf_temperature(T, P, PB);
        pres = P + PB;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% avogadro's constant, molecules / mole.
Av = 6.022e23;

% Apply ideal gas law. 8.314 = gas constant in m^3 Pa K^-1 mol^-1
ndens = (pres .* Av) ./ (8.314 .* temp);
% convert from molec./m^3 to molec./cm^3
ndens = ndens*1e-6;
end

