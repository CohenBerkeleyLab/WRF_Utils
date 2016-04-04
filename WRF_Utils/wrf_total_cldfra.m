function [ new_info ] = wrf_total_cldfra( ncdf_file, overwrite )
%[ NEW_INFO ] = WRF_TOTAL_CLDFRA( NCDF_FILE, OVERWRITE ) Computes total column cloud fraction
%   [ NEW INFO ] = WRF_TOTAL_CLDFRA( NCDF_FILE ) The first argument can
%   either be a path to the netCDF file to compute total cloud fraction for
%   or a structure returned by ncinfo. Either way, it must contain the
%   variables CLDFRA and zlev. The function will ask how to save the
%   updated netCDF file: should to append to the existing file (overwriting
%   it with the new version) or create a new copy.
%
%   [ NEW_INFO ] = WRF_TOTAL_CLDFRA( NCDF_FILE, OVERWRITE ) The OVERWRITE
%   argument allows you to specify how to save the file; 'y' will append
%   the variable to the original file, 'n' will create a new copy, and 'i'
%   (default) will prompt before appending.
%
%   WRF-Chem outputs cloud fraction on a level-by-level basis. This is
%   great for radiative transfer modeling or other applications where we
%   want specific information on where the clouds are, but is less useful
%   for satellite measurements, as what we need is the overall cloud
%   fraction for the column. Since there is definitely some overlap between
%   levels, this becomes a more complicated issue.
%
%   The simplest approach is to assume maximum overlap; thus the column
%   cloud fraction will be simply the maximum cloud fraction observed at
%   any level. But, this is not terribly realistic. Hogan and Illingworth
%   (Q.J.R. Meteorol. Soc., 2000, 126, 2903-2909) came up with an approach
%   that slightly modifies this by taking maximum overlap and adding some
%   random variation, constrained by an e-folding or correlation length
%   between layers.  Although this is applied to radar measurements, it
%   should be okay for WRF output as well. Li et al. (Q. J. R. Meteorol.
%   Soc., 2005, 131, 1607-1629) expand on this and is a useful read (esp.
%   section 3a) but we will not be using the full complexity of that paper.
%
%   Josh Laughner <joshlaugh5@gmail.com> 16 Mar 2016

E = JLLErrors;
DEBUG_LEVEL = 1;

%%%%%%%%%%%%%%%%%
%%%%% INPUT %%%%%
%%%%%%%%%%%%%%%%%

if ischar(ncdf_file)
    try
        ni = ncinfo(ncdf_file);
    catch err
        if strcmpi(err.identifier,'MATLAB:imagesci:netcdf:unableToOpenFileforRead') || strcmpi(err.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
            E.badinput('NCDF_FILE does not appear to be a valid netCDF file; either it does not exist or is a different format. Check the path and try again or pass in a structure returned from ncinfo instead.');
        else
            rethrow(err)
        end
    end
elseif isstruct(ncdf_file)
    ni = ncdf_file;
    if ~isfield(ni,'Variables')
        E.badinput('The structure NCDF_FILE does not contain the field "Variables" - are you sure it is the output from ncinfo?')
    end
end

% Find which variable is CLDFRA
cc = strcmp('CLDFRA',{ni.Variables.Name});
if sum(cc) ~= 1
    E.badinput('Either cannot find CLDFRA in the variables in NCDF_FILE or this matched multiple variables')
end
zz = strcmp('zlev',{ni.Variables.Name});
if sum(zz) ~= 1
    E.badinput('Either cannot find CLDFRA in the variables in NCDF_FILE or this matched multiple variables. Is this file the result of running slurmrun_wrf_output.sh (see BEHR/WRF_Utils folder)?')
end

if ~exist('overwrite','var')
    overwrite = 'i';
else
    overwrite = lower(overwrite);
    allowed_overwrite = {'i','y','n'};
    if ~ismember(overwrite,allowed_overwrite)
        E.badinput('OVERWRITE must be one of %s (if given)',strjoin(allowed_overwrite,', '));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the correlation length from Table 1 of Hogan and Illingworth.
% We'll use the 1 hr, 720 m one. WRF output is usually every hour or half
% hour, but since it's a model, I'd error on the side of assuming the
% clouds are more correlated. Also at the vertical level where the clouds
% are most prominent, the level spacing is usually ~700 m.
dz0 = 1900; % m

% Load the necessary vars
cldfra = ncread(ni.Filename,'CLDFRA');
zlev = ncread(ni.Filename,'zlev');
cldfra_total = nan(size(cldfra,1), size(cldfra,2), size(cldfra,4));
% Dimensions of both should be lon x lat x alt x time
for t=1:size(cldfra,4)
    if DEBUG_LEVEL > 0; fprintf('Hour %d\n',t); end
    for i=1:size(cldfra,1)
        if DEBUG_LEVEL > 1; fprintf('i=%d ', i); end
        for j=1:size(cldfra,2)
            ck = cldfra(i,j,1,t);
            for k=1:size(cldfra,3)-1
                ck1 = cldfra(i,j,k+1,t);
                alpha = exp(-zlev(i,j,k,t)/dz0);
                O = alpha*min([ck, ck1]) + (1-alpha)*ck*ck1;
                ck = ck + ck1 - O;
            end
            cldfra_total(i,j,t) = ck;
        end
    end
    if DEBUG_LEVEL > 1; fprintf('\n'); end
end

if strcmp(overwrite,'i')
    overwrite = ask_multichoice('Do you wish to append the total cldfrac to the existing file?', {'y','n'});
end

if strcmp(overwrite,'n')
    % Figure out what to name the file so as to not overwrite anything
    [fpath, fname, fext] = fileparts(ni.Filename);
    fname_raw = fname;
    n=1;
    while(exist(fullfile(fpath, [fname, fext]),'file'))
        fname = sprintf('%s-%02d',fname_raw,n);
        n = n+1;
    end
    save_file = fullfile(fpath, [fname, fext]);
    % Make of copy of the existing netCDF file so that we can append to it
    copyfile(ni.Filename, save_file);
else
    save_file = ni.Filename;
end

% Create the new variable in the file and write the data to it.
nccreate(save_file,'CLDFRA_TOT','Dimensions',{'west_east','south_north','hour_index'});
ncwrite(save_file,'CLDFRA_TOT',cldfra_total);
ncwriteatt(save_file,'CLDFRA_TOT','description','Estimated overall column cloud fraction');
ncwriteatt(save_file,'CLDFRA_TOT','units','unitless fraction');
ncwriteatt(save_file,'CLDFRA_TOT','stagger','');
ncwriteatt(save_file,'CLDFRA_TOT','coordinates','XLONG XLAT');

new_info = ncinfo(save_file);
end