function [ in, quadrants, options ] = subset_WRF_grid_cells( xlon, xlat, domain_lonlim, domain_latlim, options )
%IN = SUBSET_WRF_GRID_CELLS( XLON, XLAT, DOMAIN_LONLIM, DOMAIN_LATLIM)
%   Returns a list of indicies for grid cells that match the given criteria
%   for grid cell location and size. At its most basic, it simply finds all
%   grid cells wholly contained within the domain specified by
%   DOMAIN_LONLIM and DOMAIN_LATLIM. XLON and XLAT are the longitude and
%   latitude of the WRF-Chem grid cells. DOMAIN_LONLIM and DOMAIN_LATLIM
%   must be two element vectors specifying the respective limits of the
%   domain. It will return IN, the linear indicies of the grid cells
%   meeting all criteria.
%
%   In this mode, the function will interactively ask about additional
%   options, such as the ability to subset the data into four quadrants
%   about a point, further limit by distance from that point, and remove
%   edge grid cells by either viewing zenith angle or row.
%
%   IN = SUBSET_WRF_GRID_CELLS( XLON, XLAT, DOMAIN_LONLIM, DOMAIN_LATLIM, [] ) will
%   turn off the interactive setting of options and only subset grid cells by
%   the domain.
%
%   IN = SUBSET_WRF_GRID_CELLS( XLON, XLAT, DOMAIN_LONLIM, DOMAIN_LATLIM, OPTIONS )
%   allows the options to be specified by the structure OPTIONS rather than
%   requiring interactive questions. OPTIONS must contain the fields:
%       quad_bool - true or false, determines whether or not to subset the
%       grid cells into the four quadrants. Will be output as the second output
%       variable QUADRANTS (see below).
%  
%       dist_lim - distance in kilometers from the center lat and center
%       lon to keep grid cells for. The domain limits still need to be
%       specified as these are used to cut down the swaths with simpler
%       operations first to speed up processing. Since 1 km ~ 0.01 degrees,
%       you may compute a domain that contains as few possible grid cells as
%       possible.
%
%       size_lim_type - one of the strings 'vza', 'row', or 'none'. This
%       determines how edge grid cells will be removed. See lim_crit for more
%       information.
%
%       center_lon, center_lat - these must be scalar values between
%       -180,180 and -90,90 respectively. These set the point from which
%       distance and angle (for determination of quadrant) will be
%       calculated.
%
%   [IN, QUADRANTS] = SUBSET_WRF_GRID_CELLS( ___ ) returns a second vector the
%   same size as IN that contains the geometric quadrant each grid cell belongs
%   in relative to the center lon and lat defined by OPTIONS or
%   interactively. This vector, QUADRANTS, will have values 1-4 reflecting
%   the quadrant number: 1 is NE, 2 is NW, 3 is SW, and 4 is SE.
%
%   [IN, QUADRANTS, OPTIONS] = SUBSET_WRF_GRID_CELLS( ___ ) returns the
%   OPTIONS structure generated as well. This can then be fed back into the
%   function in the future. This is useful to both store the exact settings
%   used to subset grid cells or to use the interactive questions once in a
%   long for loop.
%
%   OPTIONS = SUBSET_WRF_GRID_CELLS() When given no input, this function
%   assumes you only want to set the options and get the resulting
%   structure back.
%   
E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


req_options = {'dist_limit','quad_bool','center_lon','center_lat'};
if ~exist('options','var')
    options.quad_bool = ask_multichoice('Divide into quadrants?',{'y','n'});
    options.dist_limit = ask_number('Only use grid cells within x km of Atlanta? 0 if no','default',0,'testfxn',@(x)x>=0,'testmsg','Value must be >= 0');
    options.quad_bool = strcmpi(options.quad_bool,'y');
    if options.quad_bool || options.dist_limit > 0
        options.center_lon = ask_number('Enter the longitude of the point to calculate distance and/or angles from','testfxn',@(x) x>=-180 && x<=180,'testmsg','Value must be between -180 and 180');
        options.center_lat = ask_number('Enter the latitude of the point to calculate distance and/or angles from','testfxn',@(x) x>=-90 && x<=90,'testmsg','Value must be between -90 and 90');
    end
elseif isempty(options)
    % Set default values
    options.quad_bool = false;
    options.dist_limit = 0;
    options.center_lon = 0;
    options.center_lat = 0;
elseif ~isstruct(options) || any(~isfield(options,req_options))
    E.badinput('OPTIONS must be a structure with the fields %s',strjoin(req_options,', '));
end

if ~exist('xlon','var')
    % If no input given, assume we just want to generate the options
    % structure.
    in = options;
    return
end

if ~isnumeric(xlon) || ~isnumeric(xlat) || ~isequal(size(xlon), size(xlat))
    E.badinput('XLON and XLAT must be the same size numeric arrays');
end
if numel(domain_lonlim) ~= 2 || numel(domain_latlim) ~= 2
    E.badinput('DOMAIN_LONLIM and DOMAIN_LATLIM must be be two-element vectors')
end

if nargout < 2 && options.quad_bool
    warning('You have selected to subset grid cells by quadrants but are not receiving the vector with that information which is the second output.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

pp = xlon >= min(domain_lonlim) & xlon <= max(domain_lonlim) & xlat >= min(domain_latlim) & xlat <= max(domain_latlim);
if sum(pp) < 1;
    in = [];
    quadrants = [];
    return
else
    pp = find(pp);
end

lon_in = xlon(pp);
lat_in = xlat(pp);

xx = true(size(lon_in));

if options.dist_limit < 0
    E.badinput('options.dist_limit must be >= 0')
elseif options.dist_limit > 0
    dist_in = zeros(size(lon_in));
    for a=1:numel(lon_in)
        dist_in(a) = m_lldist([options.center_lon, lon_in(a)], [options.center_lat, lat_in(a)]);
    end
    xx(dist_in > options.dist_limit) = false;
end

% Cut down for distance and grid cell size
in = pp(xx);
if isempty(in)
    quadrants = [];
    return
end
lon_in = xlon(in);
lat_in = xlat(in);

if nargout > 1
    if options.quad_bool
        angle_in = atan2d(lat_in - options.center_lat, lon_in - options.center_lon);
        quadrants = size(angle_in);
        quadrants(angle_in >= 0 & angle_in < 90) = 1;
        quadrants(angle_in >= 90) = 2;
        quadrants(angle_in <= -90) = 3;
        quadrants(angle_in > -90 & angle_in < 0) = 4;
    else
        quadrants = ones(size(lon_in));
    end
end

end

