function [ T ] = convert_wrf_temperature( T, P, PB )
%T = CONVERT_WRF_TEMPERATURE( WRFFILE )
%   WRF outputs temperature as perturbation to potential temperature. This
%   is not very convenient, so it must be converted to get the absolute
%   temperature. This form takes the path to a WRF output file as its sole
%   input and returns the absolute temperature in K.
%
%   T = CONVERT_WRF_TEMPERATURE( T, P, PB ) takes the arrays containing the
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

% Add an offset of 300 K to the value of T since it is defined as a
% perturbation to potential temperature.

T = T+300;

% Now convert from potential temperature to absolute temperature
T = T .* ((P+PB) ./ 1e5) .^ 0.2865;
end

