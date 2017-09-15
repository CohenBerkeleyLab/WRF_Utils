function [ dnums ] = date_from_wrf_filenames( files )
%DATE_FROM_WRF_FILENAMES Generates an array of datenums corresponding to the given WRF files
%   DNUMS = DATE_FROM_WRF_FILENAMES( FILES ) Returns the array DNUMS with
%   the date numbers corresponding to the name of each of the files in
%   FILES. FILES may be a single filename as a string, a cell array of
%   strings, or a structure returned by DIR(). In the latter two cases,
%   DNUMS will have the same shape as FILES. 
%
%   This function assumes that the date and time is given in the filename
%   as yyyy-mm-dd_HH:MM:SS or yyyy-mm-dd_HH-MM-SS. It also makes the
%   implicit assumption that only one time is represented by the wrfout
%   file, i.e. frames_per_outfile in the &time_control section of the WRF
%   namelist.input file is set to 1. If this is not true, only the time
%   represented in the filename is captured.

E = JLLErrors;

if isstruct(files) && isfield(files,'name')
    files = {files.name};
elseif ischar(files)
    files = {files};
elseif ~iscellstr(files)
    E.badinput('FILES must be a structure with the field "name", a cell array of strings, or a string');
end

dnums = nan(size(files));
for a = 1:numel(files)
    % Match either the format yyyy-mm-dd_HH:MM:SS or yyyy-mm-dd_HH-MM-SS.
    % The :'s are the default format, but they don't play nice on Macs at
    % least, so I usually santize them to -'s.
    dstr = regexp(files{a}, '\d\d\d\d-\d\d-\d\d_\d\d[\-:]\d\d[\-:]\d\d', 'match','once');
    % Homogenize the format
    dstr = strrep(dstr,':','-');
    dnums(a) = datenum(dstr, 'yyyy-mm-dd_HH-MM-SS');
end

end

