function [ noaa_sites ] = load_noaa_isd( data_dir, lonlim, latlim, datelim )
%LOAD_NOAA_ISD Loads all NOAA integrated surface database files
%   LOAD_NOAA_ISD Uses default data directory and loads all sites
%
%   LOAD_NOAA_ISD( DATA_DIR ) Uses DATA_DIR instead of the default
%
%   LOAD_NOAA_ISD( DATA_DIR, LONLIM, LATLIM ) Only load sites in the given
%   longitude/latitude limits.
%
%   LOAD_NOAA_ISD( '', LONLIM, LATLIM ) Only load sites in the given
%   longitude/latitude limits, but use the default data directory

E = JLLErrors;
DEBUG_LEVEL = 1;

if ~exist('data_dir','var') || isempty(data_dir)
    data_dir = '/Volumes/share2/USERS/LaughnerJ/MetObsData/NOAAIntegratedSurfaceDatabase';
elseif ~ischar(data_dir) 
    E.badinput('DATA_DIR must be a string')
elseif ~exist(data_dir, 'dir')
    E.dir_dne(data_dir);
end

if ~exist('lonlim','var')
    lonlim = [-180, 180];
elseif ~isnumeric(lonlim) ||  numel(lonlim) ~= 2 
    E.badinput('LONLIM must be a 2 element numeric vector');
end

if ~exist('latlim','var')
    latlim = [-90, 90];
elseif ~isnumeric(latlim) || numel(latlim) ~= 2
    E.badinput('LATLIM must be a 2 element numeric vector');
end

if ~exist('datelim','var')
    datelim = [];
elseif ~isnumeric(datelim) || numel(datelim) ~= 2 || any(datelim < 0)
    E.badinput('DATELIM must be a 2 element numeric vector with both elements >0');
end

% Find all the files. They have no extension, so we will need to remove the
% . and .. entries
F = dir(data_dir);
F(1:2)=[];

noaa_sites = cell(size(F));
def_site = NOAAISDSite;

for a=1:numel(F)
    if DEBUG_LEVEL > 0; fprintf('Loading %s... ', F(a).name); end
    fname = fullfile(data_dir, F(a).name);
    site = parse_noaa_isd(fname, lonlim, latlim, datelim);
    if strcmp(site.usaf_id, def_site.usaf_id)
        % The site is out of the lon and lat limits or date limits, do not
        % add it
        continue
    end
    if DEBUG_LEVEL > 0; fprintf('Adding site\n'); end
    noaa_sites{a} = site;
end

xx = ~iscellcontents(noaa_sites, @isempty);
noaa_sites = noaa_sites(xx);

end

