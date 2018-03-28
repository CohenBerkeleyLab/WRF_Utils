function [ tp_lev , tp_pres] = find_wrf_tropopause( wrf_info, assume_top )
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
%   to a netCDF file of WRF output that contains the calculated quantities
%   TT (actual temperature), z (altitude calculated from geopotential), and
%   pres (pressure). Both of these will be calculated if you use the
%   "calculated_quantities.nco" script (which you can run using
%   (slurm)run_wrf_output.sh - these quantities will not be in regular WRF
%   output).
%
%   [ __ ] = FIND_WRF_TROPOPAUSE( WRF_INFO, true ) if the tropopause isn't
%   found, assume that it lies above the model domain. This will return the
%   top model layer and pressure for those locations.
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
%   level will be set to the top-most index.
%
%   Josh Laughner <joshlaugh5@gmail.com> 20 Jul 2015

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

vars = {wrf_info.Variables.Name};

if ~exist('assume_top','var')
    assume_top = false;
end

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

if ismember('TT',vars)
    T = ncread(wrf_info.Filename, 'TT'); % temperature of each level in K
else
    T = convert_wrf_temperature(wrf_info.Filename);
end

if ismember('z',vars)
    z_lev = ncread(wrf_info.Filename, 'z'); % layer thickness in meters  
else
    z_lev = calculate_wrf_altitude(wrf_info.Filename);
end

if ismember( 'pres',vars)
    pres = ncread(wrf_info.Filename, 'pres'); % model box center pressure in hPa
else
    pres = (ncread(wrf_info.Filename, 'P') + ncread(wrf_info.Filename, 'PB'))/100;

end

% Since T is defined at the layer center and z the edges (staggered
% coordinates) let's convert z to non-staggered coordinates
z_lev = (z_lev(:,:,2:end,:)+z_lev(:,:,1:end-1,:))/2;

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
                if pres(x,y,z,t) > 500;
                    % If we never found any point with a lapse rate < 2
                    % K/km at all and the assume_top parameter is set,
                    % assume that we didn't see a tropopause b/c it was
                    % above the top box.  Otherwise, set the level as -1 as
                    % a cue to the user that the conditions were never met.
                    if assume_top && ~lt_2Kkm
                        tp_lev(x,y,t) = sz_bt;
                        tp_pres(x,y,t) = pres(x,y,sz_bt,t);
                    else
                        tp_lev(x,y,t) = -1;
                        tp_pres(x,y,t) = 0;
                    end
                    break
                end

                    
            end

        end
    end
end


end

