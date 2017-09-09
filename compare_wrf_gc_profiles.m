function [ gc_no2_match, wrf_no2_match, lon, lat, pres, dnums ] = compare_wrf_gc_profiles( gc_no2, gc_pres, gc_dnums, wrf_no2, wrf_pres, wrf_dnums, wrf_lon, wrf_lat, month_avg )
%COMPARE_WRF_GC_PROFILES Compare GEOS-Chem profiles to average WRF-Chem ones
%   [GC_NO2_MATCH, WRF_NO2_MATCH, LON, LAT, DNUMS] =
%   COMPARE_WRF_GC_PROFILES( GC_NO2, GC_PRES, GC_DNUMS, WRF_NO2, WRF_PRES,
%   WRF_DNUMS, WRF_LON, WRF_LAT) will examine the WRF-Chem domain and
%   determine which GEOS-Chem grid cells fall entirely within the WRF
%   domain. It will then average the WRF-Chem profiles to the GEOS-Chem
%   grid cells and return the GEOS-Chem profiles that match the averaged
%   WRF-Chem profiles along with the lon, lat, and dates of the profiles.
%   Average profiles will be returned with the first dimension
%   corresponding to the vertical dimension, the second dimension the
%   spatial dimension, and the third dimension as time. WRF profiles are
%   interpolated to GEOS-Chem pressures before averaging. GC_NO2 and
%   GC_PRES should have the same size, and should be from 2x2.5 degree
%   simulations and their length in the fourth dimension must match the
%   length of GC_DNUMS, a vector of date numbers. Likewise, WRF_NO2 and
%   WRF_PRES must be the same size and their first two dimensions must be
%   the same size as WRF_LON and WRF_LAT and their fourth dimension must
%   match the length of WRF_DNUMS.
%
%   [ ... ] = COMPARE_WRF_GC_PROFILES( ___ , true ) will indicate that the
%   WRF profiles given are monthly averages and the GEOS-Chem profiles
%   should be averaged accordingly.
%
%   The goal of this function is to enable me to compare UT GEOS-Chem and
%   WRF-Chem profiles to see how much of a difference lightning causes in
%   the profile.
%
%   Josh Laughner <joshlaugh5@gmail.com> 17 Jun 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

if ~exist('month_avg','var')
    month_avg = false;
else
    if ~isscalar(month_avg) || ~islogical(month_avg)
        E.badinput('MONTH_AVG (if given) must be a scalar logical')
    elseif month_avg
        E.notimplemented('monthly averages')
    end
end

[gloncorn, glatcorn] = geos_chem_corners;
gloncorn = gloncorn';
glatcorn = glatcorn';

% The GC and WRF NO2 arrays should be 3 or 4 dimensional, no more. Also
% check that the size is what we expect, as far as matching lat/lon
% coordinates.

if ndims(gc_no2) > 5
    E.badinput('GC_NO2 cannot have greater than 4 non-singleton dimensions (lon, lat, alt, time)')
end
if ndims(gc_pres) > 5
    E.badinput('GC_PRES cannot have greater than 4 non-singleton dimensions (lon, lat, alt, time)')
end
if ndims(wrf_no2) > 5
    E.badinput('WRF_NO2 cannot have greater than 4 non-singleton dimensions (lon, lat, alt, time)')
end
if ndims(wrf_pres) > 5
    E.badinput('WRF_PRES cannot have greater than 4 non-singleton dimensions (lon, lat, alt, time)')
end

if size(gc_no2,1) ~= size(gloncorn,1)-1 || size(gc_no2,2) ~= size(gloncorn,2)-1
    E.badinput('GC_NO2 must be from a 2x2.5 deg simulation (144x91 grid cells')
end
if size(gc_pres,1) ~= size(gloncorn,1)-1 || size(gc_pres,2) ~= size(gloncorn,2)-1
    E.badinput('GC_PRES must be from a 2x2.5 deg simulation (144x91 grid cells')
end
if ~isvector(gc_dnums)
    E.badinput('GC_DNUMS must be a vector')
end
if size(gc_no2,4) ~= length(gc_dnums)
    E.badinput('The fourth dimension of GC_NO2 must have the same length as the vector GC_DNUMS')
end
if size(gc_pres,4) ~= length(gc_dnums)
    E.badinput('The fourth dimension of GC_PRES must have the same length as the vector GC_DNUMS')
end
if ~isequal(size(gc_no2), size(gc_pres))
    E.badinput('GC_NO2 and GC_PRES must be the same size')
end

if ndims(wrf_lon) > 3 
    wrf_lon = wrf_lon(:,:,:);
    if any(isnan(wrf_lon(:)))
        E.badinput('WRF_LON should not contain NaNs')
    end
    test_diff = sum(abs(diff(wrf_lon,1,3)));
    if any(test_diff(:) > 0)
        E.badinput('All 2D slices of WRF_LON must be the same if passing a higher-dimensions matrix. Only the first 2D slice is used')
    end
    wrf_lon = wrf_lon(:,:,1);
end
if ndims(wrf_lat) > 3 
    wrf_lat = wrf_lat(:,:,:);
    if any(isnan(wrf_lat(:)))
        E.badinput('WRF_LAT should not contain NaNs')
    end
    test_diff = sum(abs(diff(wrf_lat,1,3)));
    if any(test_diff(:) > 0)
        E.badinput('All 2D slices of WRF_LAT must be the same if passing a higher-dimensions matrix. Only the first 2D slice is used')
    end
    wrf_lat = wrf_lat(:,:,1);
end
if size(wrf_no2,1) ~= size(wrf_lon, 1) || size(wrf_no2,2) ~= size(wrf_no2,2)
    E.badinput('WRF_NO2 must have the same size first two dimensions as WRF_LON and WRF_LAT')
end
if size(wrf_pres,1) ~= size(wrf_lon, 1) || size(wrf_pres,2) ~= size(wrf_no2,2)
    E.badinput('WRF_PRES must have the same size first two dimensions as WRF_LON and WRF_LAT')
end
if ~isvector(wrf_dnums)
    E.badinput('WRF_DNUMS must be a vector')
end
if size(wrf_no2, 4) ~= length(wrf_dnums)
    E.badinput('The fourth dimension of WRF_NO2 must have the same length as the vector WRF_DNUMS')
end
if size(wrf_pres, 4) ~= length(wrf_dnums)
    E.badinput('The fourth dimension of WRF_NO2 must have the same length as the vector WRF_DNUMS')
end
if ~isequal(size(wrf_no2), size(wrf_pres))
    E.badinput('WRF_NO2 and WRF_PRES must be the same size')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% First identify which GC grid cells fall entirely within the WRF domain.
% Assume that the first dimension of the WRF arrays is the west_east one,
% then look for the least western point on the west edge and sim. on the
% east, north, south. This accounts for the use of non-equirectangular map
% projections in WRF by essentially inscribing a box around the coordinates
% guaranteed to be inside the domain.
wrf_lonbdy = [max(min(wrf_lon,[],1),[],2), min(max(wrf_lon,[],1),[],2)];
wrf_latbdy = [max(min(wrf_lat,[],2),[],1), min(max(wrf_lat,[],2),[],1)];

% GEOS-Chem corners are equirectangular so we can just look at one vector
% to find which grid cells are entirely inside the WRF domain.
glonvec = gloncorn(:,1);
glatvec = glatcorn(1,:);

gcxx = glonvec(1:end-1) >= wrf_lonbdy(1) & glonvec(2:end) <= wrf_lonbdy(2);
gcyy = glatvec(1:end-1) >= wrf_latbdy(1) & glatvec(2:end) <= wrf_latbdy(2);

% Now we start by matching days, then we'll loop over days to actually pull
% out profiles

if ~month_avg
    dd = find(ismember(gc_dnums, wrf_dnums));
    dnums = gc_dnums(dd);
end
gc_no2_match = nan(size(gc_no2,3), sum(gcxx)*sum(gcyy), numel(dd));
wrf_no2_match = nan(size(gc_no2,3), sum(gcxx)*sum(gcyy), numel(dd));
pres = nan(size(gc_no2,3), sum(gcxx)*sum(gcyy), numel(dd));

[glon, glat] = geos_chem_centers('2x25');

lon = repmat(glon(gcxx)', 1, sum(gcyy));
lon = lon(:);
lat = repmat(glat(gcyy), sum(gcxx), 1);
lat = lat(:);


gcxxf = find(gcxx);
gcyyf = find(gcyy);

for a=1:numel(dd)
    gctt = gc_dnums == dnums(a);
    gc_no2_tmp = gc_no2(gcxx, gcyy, :, gctt);
    gc_no2_tmp = reshape(gc_no2_tmp, [], size(gc_no2_tmp,3), size(gc_no2_tmp,4));
    gc_no2_tmp = permute(gc_no2_tmp, [2 1 3]);
    gc_no2_match(:,:,a) = gc_no2_tmp;
    
    gc_pres_tmp = gc_pres(gcxx, gcyy, :, gctt);
    gc_pres_tmp = reshape(gc_pres_tmp, [], size(gc_pres_tmp,3), size(gc_pres_tmp,4));
    gc_pres_tmp = permute(gc_pres_tmp, [2 1 3]);
    pres(:,:,a) = gc_pres_tmp;
    
    gc_lonmins = reshape(gloncorn(gcxxf, gcyyf),1,[]);
    gc_lonmaxes = reshape(gloncorn(gcxxf+1, gcyyf+1),1,[]);
    gc_latmins = reshape(glatcorn(gcxxf, gcyyf),1,[]);
    gc_latmaxes = reshape(glatcorn(gcxxf+1, gcyyf+1),1,[]);
    
    wrftt = wrf_dnums == dnums(a);
    
    for b=1:numel(gc_lonmins)
        wrfxx = wrf_lon >= gc_lonmins(b) & wrf_lon < gc_lonmaxes(b);
        wrfyy = wrf_lat >= gc_latmins(b) & wrf_lat < gc_latmaxes(b);
        
        wrf_no2_tmp = permute(wrf_no2(:,:,:,wrftt), [3 1 2]);
        wrf_no2_tmp = wrf_no2_tmp(:, wrfxx & wrfyy);
        
        wrf_pres_tmp = permute(wrf_pres(:,:,:,wrftt), [3 1 2]);
        wrf_pres_tmp = wrf_pres_tmp(:, wrfxx & wrfyy);
        
        wrf_no2_interp = nan(size(gc_no2_match,1), size(wrf_no2_tmp,2));
        for c=1:size(wrf_no2_tmp, 2)
            % Interpolate WRF NO2 profiles to GC pressures in log-log space
            % which I believe accounts for the exponential relationship
            % between mixing ratio concentration and pressure
            
            wrf_no2_interp(:,c) = exp(interp1(log(wrf_pres_tmp(:,c)), log(wrf_no2_tmp(:,c)), log(pres(:,b,a))));
        end
        
        wrf_no2_match(:,b,a) = mean(wrf_no2_interp, 2);
    end
end


end

