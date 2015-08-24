function [ fire_mask, frp, lon, lat ] = read_myd14( hdfi )
%READ_MYD14(HDFI) Read in MODIS thermal anomaly file
%   Function will read in the L2 swath MODIS thermal anomaly file
%   identified by the input hdfi. This can be either a structure returned
%   by hdfinfo or the full path to an MYD14 (or MOD14) file.
%
%   This will return four vectors: the value of the fire mask, fire
%   radiative power, longitude, and latitude of each pixel identified as a
%   fire pixel.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(hdfi)
    if ~isfield(hdfi,'Filename')
        E.badinput('hdfi does not seem to be a structure returned from hdfinfo');
    end
    [~,filename] = fileparts(hdfi.Filename);
    if isempty(regexp(filename,'M[OY]D14','ONCE'))
        E.badinput('hdfi does not appear to point to a MODIS MOD14 or MYD14 file');
    end
elseif ischar(hdfi)
    if ~exist(hdfi,'file')
        E.badinput('%s does not exist',hdfi);
    end
    [~,filename] = fileparts(hdfi);
    if isempty(regexp(filename,'M[OY]D14','ONCE'))
        E.badinput('hdfi does not appear to point to a MODIS MOD14 or MYD14 file');
    end
    hdfi = hdfinfo(hdfi);
end

%%%%%%%%%%%%%%%%%%%
%%%%% READ IN %%%%%
%%%%%%%%%%%%%%%%%%%

% Assuming that the HDF file is structures as it was when I wrote this, so
% that the SDSs are right in the top level of the structure.

fire_mask = hdfread(hdfi.Filename, hdfdsetname(hdfi,'fire mask'));

% Fires are indicated in the fire mask by values of 7, 8, or 9 (7 is the
% lowest confidence, 9 is the highest confidence). If there are no fires,
% the following three quantities cannot be imported without an error, so we
% set them to be empty matrices.

if sum(fire_mask(:)>6) < 1
    frp = [];
    lon = [];
    lat = [];
    fire_mask = [];
else
    frp = hdfread(hdfi.Filename, hdfdsetname(hdfi, 'FP_power'));
    lon = hdfread(hdfi.Filename, hdfdsetname(hdfi, 'FP_longitude'));
    lat = hdfread(hdfi.Filename, hdfdsetname(hdfi, 'FP_latitude'));
    
    % Now handle fire_mask, which comes in as a 2D matrix compared to the other
    % quantities, which are 1-by-n vectors, n is the number of fires. FP_line
    % and FP_sample have the indices of fire_mask that correspond to the other
    % quantities, but the are 0 based (and Matlab is 1 based).
    
    fp_line = hdfread(hdfi.Filename, hdfdsetname(hdfi, 'FP_line'));
    fp_sample = hdfread(hdfi.Filename, hdfdsetname(hdfi, 'FP_sample'));
    
    linind = zeros(size(fp_line));
    for a=1:numel(linind)
        linind(a) = sub2ind(size(fire_mask),fp_line(a)+1,fp_sample(a)+1);
    end
    
    fire_mask = fire_mask(linind);
    
end

end


