function [ tp_lev , tp_pres] = find_wrf_tropopause( wrf_info, varargin )
%FIND_WRF_TROPOPAUSE Find the model level where the tropopause is
%   The WRF preprocessor determines the tropopause level in the model as
%   being where the average lapse rate over 3 model layers is < 2 K/km.
%   This replicates that calculation.
%
%   [ TP_LEV, TP_PRES ] = FIND_WRF_TROPOPAUSE( WRF_INFO ) returns the model
%   level that the tropopause resides in as TP_LEV and the pressure of that
%   model level and TP_PRES. Both will have size west_east by south_north
%   by time. If it fails to identify the tropopause, then it will return -1
%   for the level and 0 for the pressure.
%
%   WRF_INFO is a structure obtained from ncinfo that points
%   to a netCDF file of WRF output. It must contain information about
%   temperature, elevation, pressure, latitude, and longitude. The first
%   three quantities can be either calculated by calculated_quantities.nco
%   or available in the "raw" form (i.e. potential temperature,
%   geopotential height, and base + perturbation pressure).
%
%   [ __ ] = FIND_WRF_TROPOPAUSE( WRF_INFO, true ) if the tropopause isn't
%   found, assume that it lies above the model domain. This will return the
%   top model layer and pressure for those locations.
%
%   Additional parameters:
%
%       'error_if_missing_units' - by default, this checks the units in the
%       WRF file it is reading from. However, in some cases, units might
%       not have been stored in the output, so setting this parameter to
%       false allowed you to override that and assume the units are
%       correct.
%
%   This function works by starting at the top of each vertical profile and
%   moving downwards, looking for the last group of 3 model layers that
%   have an average lapse rate < 2 K/km.  It will stop at z = 3000 m, to
%   avoid accidentally detecting a temperature inversion.  If it reaches
%   the 3000 m cutoff, it will by default assign a value of -1 as the
%   tropopause level, which should be interpreted as indicating that the
%   algorithm failed to find a tropopause.  By passing "true" as the
%   optional second parameter, it will check to see if the lapse rate was
%   ever < 2 K/km. If not, it will assume that it did not find a tropopause
%   because all layers in that profile are below it, and so the tropopause
%   level will be set to the top-most index. It also does three additional
%   filtering steps:
%
%       1) It reruns the above lapse rate calculation omitting the top
%       three layers, and sees if it gets the same tropopause level. If
%       not, it prefers the lower one. This helps catch unexpectely high
%       tropopause levels if the laps rate near the top dips below 2 K/km.
%
%       2) It looks for big jumps along the zonal direction, assuming that
%       the tropopause pressure will be reasonably invariant along the same
%       latitude (actually the same south_north row in the WRF output). If
%       there is a jump, it finds the region that has a jump of >= 50 hPa
%       compared to the rest of the zonal band and marks that it could not
%       find a valid tropopause there.
%
%       3) Finally, it compares the values of each tropopause pressure to
%       the median in a 50-by-30 chunk and any pressure that exceeds a 70
%       hPa difference compared to the chunk median is marked as "could not
%       find tropopause".

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
p = advInputParser;
p.addOptional('assume_top', false);
p.addParameter('error_if_missing_units', true);
p.parse(varargin{:});
pout = p.Results;

assume_top = pout.assume_top;
error_if_missing_units = pout.error_if_missing_units;

if (~islogical(assume_top) && ~isnumeric(assume_top)) || ~isscalar(assume_top)
    E.badinput('assume_top must be a scalar logical or numeric value, or be left unspecified');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%



% Find the dimensions so we know how big the arrays will be
dims = {wrf_info.Dimensions.Name};

dd = strcmp('bottom_top',dims);
sz_bt = wrf_info.Dimensions(dd).Length;
dd = strcmp('west_east',dims);
sz_we = wrf_info.Dimensions(dd).Length;
dd = strcmp('south_north',dims);
sz_sn = wrf_info.Dimensions(dd).Length;
dd = strcmp('Time',dims);
if sum(dd) == 0
    sz_time = 1;
else
    sz_time = wrf_info.Dimensions(dd).Length;
end

tp_lev = zeros(sz_we, sz_sn, sz_time);
tp_pres = zeros(sz_we, sz_sn, sz_time);
%  The WRF pre-processor defines the tropopause as the first level where the
% average lapse rate over 3 layers is < 2 K/km. So we calculate the lapse
% rate averaged over 3 bins and look for the lowest one that meets the
% criteria.

T = read_wrf_preproc(wrf_info.Filename, 'temperature', 'error_if_missing_units', error_if_missing_units);
% read_wrf_preproc will convert from staggered to centered z levels if
% necessary
z_lev = read_wrf_preproc(wrf_info.Filename, 'z_center', 'error_if_missing_units', error_if_missing_units);
pres = read_wrf_preproc(wrf_info.Filename, 'pressure', 'error_if_missing_units', error_if_missing_units);

wrf_lon = ncread(wrf_info.Filename,'XLONG');
wrf_lat = ncread(wrf_info.Filename,'XLAT');

for x = 1:sz_we
    for y = 1:sz_sn
        for t = 1:sz_time
            lt_2Kkm = false;
            for z = (sz_bt-1):(-1):1
                % Go from the top down. Once we find 1 case where the lapse
                % rate is < 2 K/km, search until we hit one > 2 K/km
                
                % Include the top two lapse rates, but they'll need
                % calculated specially (since there's not 3 layers to
                % average)
                if sz_bt - z == 1 || sz_bt - z == 2
                    end_ind = sz_bt;
                else
                    end_ind = z + 3;
                end
                
                % Calculate the lapse rate in three layer chunks. z is in meters, so
                % convert it to km so the lapse rate is K/km. Also take the negative since
                % lapse rate is defined as -dT/dz.
                lapse = -(T(x,y,end_ind,t) - T(x,y,z,t))/((z_lev(x,y,end_ind,t)-z_lev(x,y,z,t))/1000);
                if ~lt_2Kkm
                    if lapse < 2
                        lt_2Kkm = true;
                    end
                else
                    if lapse > 2
                        tp_lev(x,y,t) = z;
                        tp_pres(x,y,t) = pres(x,y,z,t);
                        break
                    end
                end
                
                % Reject if the pressure is >500 hPa (i.e. below the 500
                % hPa altitude). This is part of the WMO definition of
                % tropopause, which rejects any lapse rates < 2 K/km below
                % this unless it is the only one.  We are likewise going to
                % always reject these because surface temperature
                % inversions will confuse the algorithm.
                if pres(x,y,z,t) > 500
                    % If we never found any point with a lapse rate < 2
                    % K/km at all and the assume_top parameter is set,
                    % assume that we didn't see a tropopause b/c it was
                    % above the top box.  Otherwise, set the level as -1 as
                    % a cue to the user that the conditions were never met.
                    
                    if assume_top && ~lt_2Kkm
                        tp_lev(x,y,t) = sz_bt;
                        tp_pres(x,y,t) = pres(x,y,sz_bt,t);
                    else
                        tp_lev(x,y,t) = nan;
                        tp_pres(x,y,t) = nan;
                    end
                    break
                end
                
            end
            
            % Here we found if searching lapse rate larger than 2 from top
            % down, in some cases the lapse rate in the first three layers
            % are slightly over 2 so the function above will recognize it
            % as tropopause pressure. However, it will cause the sharp
            % change of tropopause pressure in adjacent grid cells. To get
            % rid of this, we omit the top three layers and search lapse
            % rate again. If we find 1 case where the lapse rate is < 2
            % K/km, and also hit one > 2 K/km from the lower layers,
            % tropopause pressure will be replace by lower value. If this
            % results in oddly high (low altitude) tropopause, it should be
            % caught later in the code.
            %
            lt_2Kkm = false;
            for z = (sz_bt-4):(-1):1
                end_ind = z + 3;
                lapse = -(T(x,y,end_ind,t) - T(x,y,z,t))/((z_lev(x,y,end_ind,t)-z_lev(x,y,z,t))/1000);
                
                if ~lt_2Kkm
                    if lapse < 2
                        lt_2Kkm = true;
                    end
                else
                    if lapse > 2
                        tp_lev(x,y,t) = z;
                        tp_pres(x,y,t) = pres(x,y,z,t);
                        break
                    end
                end
                
                if pres(x,y,z,t) > 500
                    break
                end
            end
        end
    end
end

% Considering the weakness of this algorithm, floodfill is intended to find
% the points with absurd tropopause pressure and set the tp_lev to be -1
% and tp_pres to be 0, i.e. not found. This can be interpolated by the
% calling function if desired.
xx_bad_trop = false(size(tp_pres));

% search center points along the altitude, locate the adjacent points
% with sharp changes in tropopause pressure and set the first point as
% center point in the function floodfill
for y = 1:sz_sn
    for t = 1:sz_time
        tp_pres_diff = abs(tp_pres(2:end,y,t)-tp_pres(1:end-1,y,t));
        dp_pres = find(tp_pres_diff >= 50); % 50 hPa jump chosen based on 2012 WRF data?
        for i = 1:numel(dp_pres)
            % With some test, quantile(tp_pres_diff,0.7) is always around
            % 0.5 pa. (Pa or hPa?)
            % Need to keep an eye on this - there might occur a case where
            % the bad tropopause values are at the edge of the domain and I
            % think this might mark the good tropopause as bad because it
            % comes across the jump from bad --> good.
            if ~xx_bad_trop(dp_pres(i),y)
                tolerance_pres = quantile(tp_pres_diff,0.7);
                threshold = @(t) abs(t) < tolerance_pres;
                center_lon = wrf_lon(dp_pres(i),y);
                center_lat = wrf_lat(dp_pres(i),y);
                [in_plume] = floodfill(tp_pres, wrf_lon, wrf_lat, threshold, center_lon, center_lat);
                xx_bad_trop = xx_bad_trop | in_plume;
            end
        end
    end
end


tp_pres(xx_bad_trop) = nan;
tp_lev(xx_bad_trop) = nan;

% second run of filter: devide the map by 50X30 chunks, in each chunk find
% the grid cell that the tropopause pressure is 70hpa larger or lower than
% the median tropopause pressure, set it to be nan; (Why 70 hPa here, why
% 50x30 chunks?)

% This chunk size was selected for WRF domains at 12 km resolution. It may
% need adjusted in the future to be grid spacing-aware.
bulk = [50,30];
s1 = fix(sz_we/bulk(1));
s2 = fix(sz_sn/bulk(2));

for i=1:s1
    for j=1:s2
        if j ==s2
            ybulk = bulk(2)*(j-1)+1:max(bulk(2)*j,sz_sn);
        else
            ybulk = bulk(2)*(j-1)+1:bulk(2)*j;
        end
        if i ==s1
            xbulk = bulk(1)*(i-1)+1:max(bulk(1)*i,sz_we);
        else
            xbulk = bulk(1)*(i-1)+1:bulk(1)*i;
        end
        pres_bulk = tp_pres(xbulk,ybulk);
        median_pres = nanmedian(pres_bulk(:));
        pres_diff = abs(pres_bulk-median_pres);
        indx = pres_diff > 70;
        pres_bulk(indx) = nan;
        tp_pres(xbulk,ybulk) = pres_bulk;
    end
end

tp_pres(isnan(tp_pres)) = 0;
tp_lev(isnan(tp_lev)) = -1;
end
