function [ varargout ] = read_wrf_vars( filepath, filenames, varnames, force_dims, DEBUG_LEVEL )
%READ_WRF_VARS Reads WRF_BEHR .nc files
%   This function will read in arbitrary variables from WRF_BEHR .nc files
%   (the result of processing raw wrfout files using
%   (slurm)run_wrf_output.sh in this folder). Requires 3 arguments:
%
%       FILEPATH - a path to the .nc files.
%
%       FILENAMES - this can either be a structure output from the dir()
%       function (i.e. each file is an entry in the structure, and has its
%       filename in the "name" field) or a cell array of filenames.
%
%       VARNAMES - a cell array of variable names to read in, i.e. the name
%       you would pass as the second argument to ncread.
%
%   There is an optional argument, force_dims, that controls how the output
%   dimensions of the variables are set up. By default, each variable will
%   be concatenated along a new dimension, i.e. if U is a 4-D variable in
%   the files, then the various files will be concatenated along the 5th
%   dimension, whereas if time is a 1-D variable, it will be concatenated
%   along the 2nd dimension. However, if force_dims is set to true, all
%   variables will be concatenated along the same dimension. So in the
%   example above, both U and time will be concatenated along the 5th
%   dimension.
%
%   This will return the variables in the same order they are named.
%
%   Josh Laughner <joshlaugh5@gmail.com> 28 Aug 2015


E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT PARSING %%%%
%%%%%%%%%%%%%%%%%%%%%%%

if ~ischar(filepath)
    E.badinput('filepath must be a string')
elseif ~exist(filepath, 'dir')
    E.badinput('%s is an invalid directory')
end

if isstruct(filenames)
    if isfield(filenames, 'name')
        filenames = {filenames.name};
    else 
        E.badinput('filenames does not contain the field "name"')
    end
end
if iscell(filenames)
    if any(~iscellcontents(filenames,'ischar'))
        E.badinput('All entries in "filenames" must be strings')
    end
else
    E.badinput('filenames must be a struct with field "names" or a cell array of filenames')
end

if ischar(varnames)
    varnames = {varnames};
end

if ~iscell(varnames) || any(~iscellcontents(varnames,'ischar'))
    E.badinput('varnames must be a cell array of strings or a single string')
elseif ~isvector(varnames)
    E.badinput('varnames should be a vector cell array')
end

if ~exist('force_dims','var')
    force_dims = false;
elseif ~isscalar(force_dims) || (~isnumeric(force_dims) && ~islogical(force_dims))
    E.badinput('force_dims, if given, must be a scalar logical or number than can be converted to a logical')
end

if ~exist('DEBUG_LEVEL','var')
    DEBUG_LEVEL = 1;
elseif ~isscalar(DEBUG_LEVEL) || ~isnumeric(DEBUG_LEVEL)
    E.badinput('DEBUG_LEVEL must be a scalar')
end
    
%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN FUNCTION %%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Prepare varargout to hold the variables requested
varargout = cell(size(varnames));
n_days = numel(filenames);

wrfi = ncinfo(fullfile(filepath, filenames{1}));
wrf_vars = {wrfi.Variables.Name};
var_dims = zeros(1,numel(varnames)); % used for concatenation
var_sizes = cell(1,numel(varnames)); % used to check each file's variables

% First figure out how many dimensions each variable has; we'll need this
% whether or not force_dims is set to determine which dimension to
% concatenate along.
for a=1:numel(varnames)
    vv = strcmp(varnames{a}, wrf_vars);
    if sum(vv) < 1
        E.callError('var_not_found','Variable %s cannot be found in the first WRF file',varnames{a})
    elseif sum(vv) > 1
        E.callError('too_many_var','Multiple instances of %s found in the first WRF file',varnames{a})
    end
    
    var_dims(a) = numel(wrfi.Variables(vv).Size)+1;
    var_sizes{a} = wrfi.Variables(vv).Size;
end

max_dims = max(var_dims);
if force_dims
    % Make all concatenate along the same dimension
    var_dims(:) = max_dims;
end

% for a=1:numel(varnames)
%     vv = strcmp(varnames{a}, wrf_vars);
%     if sum(vv) < 1
%         E.callError('var_not_found','Variable %s cannot be found in the first WRF file',varnames{a})
%     elseif sum(vv) > 1
%         E.callError('too_many_var','Multiple instances of %s found in the first WRF file',varnames{a})
%     end
%     
%     if ~force_dims
%         sz = [wrfi.Variables(vv).Size, n_days];
%         varargout{a} = nan(sz);
%     else
%         % This will create a matrix with singleton dimensions in any
%         % dimensions between its last one and the concatenation dimension.
%         sz = wrfi.Variables(vv).Size;
%         sz(max_dims) = n_days;
%         sz(sz==0) = 1;
%         varargout{a} = nan(sz); 
%     end      
% end

% Now go through each file and import each variable into the correct cell
% in varargout, concatenating as we go.

for d=1:n_days
    if DEBUG_LEVEL > 0
        fprintf('Reading day %d of %d: ', d, n_days);
    end
    for a=1:numel(varnames)
        if DEBUG_LEVEL > 0 
            fprintf('%s ', varnames{a});
        end
        var = ncread(fullfile(filepath, filenames{d}),varnames{a});
        
        % Check that the variables have the same dimensions as the first
        % day.
        %if (ndims(var) ~= numel(var_sizes{a}) && var_sizes{a}(end) > 1) || (var_sizes{a}(end) > 1 && any(size(var) ~= var_sizes{a}))
        %    E.callError('import','%s in the file %s is not the same size as in the first .nc file (%s vs. %s in the first)', varnames{a}, filenames{a}, mat2str(size(var)), mat2str(var_sizes{a}));
        %end
        varargout{a} = cat(var_dims(a), varargout{a}, var);
    end
    if DEBUG_LEVEL > 0
        fprintf('\n');
    end
end

end

