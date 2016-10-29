function [ dnums ] = date_from_wrf_filenames( files )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

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

