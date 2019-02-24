function [ varargout ] = read_wrf_vars( filepath, filenames, varnames, varargin )
%READ_WRF_VARS Reads WRF_BEHR .nc files
%   [ VAR1, VAR2, ...] = READ_WRF_VARS( FILEPATH, FILENAMES, VARNAMES )
%   This function will read in arbitrary variables from WRF_BEHR .nc files
%   (the result of processing raw wrfout files using
%   (slurm)run_wrf_output.sh in this folder). Requires 3 arguments:
%
%       FILEPATH - a path to the netCDF files. Alternatively, pass an empty
%       string to indicate that the full path is contained in each file
%       name.
%
%       FILENAMES - this can either be a structure output from the dir()
%       function (i.e. each file is an entry in the structure, and has its
%       filename in the "name" field) or a cell array of filenames.
%
%       VARNAMES - a cell array of variable names to read in, i.e. the name
%       you would pass as the second argument to read_wrf_preproc. Since
%       this calls read_wrf_preproc internally, you can specify computed
%       quantities like 'temperature' and 'pressure' that read_wrf_preproc
%       understands.
%
%   By default, each variable will be returned as a separate output
%   variable (that is, there will be as many outputs as there were entries
%   in VARNAMES). Each output will be the corresponding variable
%   concatenated along a new dimension. So if P is a 4D variable in each
%   WRF file, it will be returned as a 5D variable, where P(:,:,:,:,1) is
%   from the first file, P(:,:,:,:,2) from the second, and so on.
%
%   [ ___ ] = READ_WRF_VARS( ___, FORCE_DIMS ) controls how the output
%   dimensions of the variables are set up. By default, each variable will
%   be concatenated along a new dimension, i.e. if U is a 4-D variable in
%   the files, then the various files will be concatenated along the 5th
%   dimension, whereas if time is a 1-D variable, it will be concatenated
%   along the 2nd dimension. However, if force_dim can be given one of
%   several values:
%
%       'no' - default behavior 
%
%       'same' - all variables will be concatenated along the same
%       dimension. So in the example above, both U and time will be
%       concatenated along the 5th dimension.
%
%       'squeeze' - all variables will be squeezed before returning,
%       eliminating singleton dimensions. In the example above, if U is 4-D
%       but the 4th dimension is always only 1 long (because there's only
%       one time in each file), that dimension will be removed so that
%       U(:,:,:,1) is from the first file, instead of U(:,:,:,:,1).
%
%   FORCE_DIMS may also be specified as a scalar logical value (true or
%   false) or a scalar number. If it evaluates as false, that is equivalent
%   to a value of 'no', and if it evaluates as true, that is equivalent to
%   'same'. This is preserved for backwards compatibility.
%
%   [ ___ ] = READ_WRF_VARS( ___, FORCE_DIMS, DEBUG_LEVEL ) You may also
%   specify how verbose this function is by the DEBUG_LEVEL argument, which
%   is a scalar number >= 0 or the char array 'visual'. If it is 'visual',
%   a graphical wait bar is used instead.
%
%   [ ___ ] = READ_WRF_VARS( ___, 'as_struct' ) With any of the other
%   syntaxes, including the 'as_struct' flag will return all the variables
%   in a structure, with each variable as a field name, instead of as
%   separate output variables.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT PARSING %%%%
%%%%%%%%%%%%%%%%%%%%%%%

p = advInputParser;
p.addOptional('force_dims', 'no');
p.addOptional('DEBUG_LEVEL', 1);
p.addFlag('as_struct');

p.parse(varargin{:});
pout = p.Results;

force_dims = pout.force_dims;
DEBUG_LEVEL = pout.DEBUG_LEVEL;
as_struct = pout.as_struct;


if ~ischar(filepath)
    E.badinput('filepath must be a string')
elseif ~exist(filepath, 'dir') && ~isempty(filepath)
    E.badinput('%s is an invalid directory', filepath)
end

filenames = files_input(filenames);

if ischar(varnames)
    varnames = {varnames};
end

if ~iscell(varnames) || any(~iscellcontents(varnames,'ischar'))
    E.badinput('varnames must be a cell array of strings or a single string')
elseif ~isvector(varnames)
    E.badinput('varnames should be a vector cell array')
end

if islogical(force_dims) || isnumeric(force_dims)
    if ~isscalar(force_dims)
        E.badinput('If given as a logical or numeric value, FORCE_DIMS must be a scalar value')
    elseif force_dims
        force_dims = 'same';
    else
        force_dims = 'no';
    end
elseif ischar(force_dims)
    allowed_dim_modes = {'no','same','squeeze'};
    if ~ismember(force_dims, allowed_dim_modes)
        E.badinput('If given as a char array, FORCE_DIMS must be one of: %s', strjoin(force_dims, ', '));
    end
end

waitbar_bool = false;
if ischar(DEBUG_LEVEL) && strcmpi('visual',DEBUG_LEVEL)
    if isDisplay
        waitbar_bool = true;
        DEBUG_LEVEL = 0;
    else
        warning('Cannot use waitbar when running in non-GUI mode, defaulting to DEBUG_LEVEL = 1')
        DEBUG_LEVEL = 1;
    end
elseif ~isscalar(DEBUG_LEVEL) || ~isnumeric(DEBUG_LEVEL)
    E.badinput('DEBUG_LEVEL must be a scalar or the string ''visual''')
end
    
%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN FUNCTION %%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Prepare read_data to hold the variables requested
read_data = cell(size(varnames));
n_days = numel(filenames);


wrfi = ncinfo(fullfile(filepath, filenames{1}));
wrf_vars = {wrfi.Variables.Name};
var_dims = zeros(1,numel(varnames)); % used for concatenation
var_sizes = cell(1,numel(varnames)); % used to check each file's variables

% First figure out how many dimensions each variable has; we'll need this
% whether or not force_dims is set to determine which dimension to
% concatenate along.
temperature_strs = {'t','temp','temperature'};
pressure_strs = {'p','pres','pressure'};
elevation_strs = {'z','elev','elevation'};
elev_center_strs = {'z_center'};
ndens_strs = {'ndens','number density'};
no2_ndens_strs = {'no2_ndens'};

for a=1:numel(varnames)
    % Special cases are the precalculated variables available in the subset
    % files. We use aliases for those when loading, so first see if the actual
    % variable we need exists.
    if any(strcmpi(temperature_strs, varnames{a}))
        var_for_size = 'TT';
    elseif any(strcmpi(pressure_strs, varnames{a}))
        var_for_size = 'pres';
    elseif any(strcmpi(elevation_strs, varnames{a}))
        var_for_size = 'z';
    elseif any(strcmpi(elev_center_strs, varnames{a}))
        var_for_size = 'P'; % will remove the stagger, so don't use z
    elseif any(strcmpi(ndens_strs, varnames{a}))
        var_for_size = 'ndens';
    elseif any(strcmpi(no2_ndens_strs, varnames{a}))
        var_for_size = 'no2_ndens';
    else
        var_for_size = varnames{a};
    end

    vv = strcmp(var_for_size, wrf_vars);
    if sum(vv) < 1
        % If we can't find the variable, it might be a preprocessed
        % variable that will be computed by read_wrf_preproc. So here we'll
        % map those preprocessed variables to the underlying variables that
        % we can get the size from. I'm being a bit lazy and not checking
        % that these do exist and that they are not multiply defined. Also,
        % if read_wrf_preproc gets updated, this will need updated too.
        if any(strcmpi(temperature_strs, varnames{a}))
            vv = strcmp('T', wrf_vars);
        elseif any(strcmpi(pressure_strs, varnames{a}))
            vv = strcmp('P', wrf_vars);
        elseif any(strcmpi(elevation_strs, varnames{a}))
            vv = strcmp('PH', wrf_vars);
        elseif any(strcmp(elev_center_strs, varnames{a}))
            vv = strcmp('PH', wrf_vars);
        elseif any(strcmpi(ndens_strs, varnames{a}))
            vv = strcmp('P', wrf_vars);
        elseif any(strcmpi(no2_ndens_strs, varnames{a}))
            vv = strcmp('no2', wrf_vars);
        else
            E.callError('var_not_found','Variable %s cannot be found in the first WRF file',varnames{a})
        end
    elseif sum(vv) > 1
        E.callError('too_many_var','Multiple instances of %s found in the first WRF file',varnames{a})
    end
    
    var_dims(a) = numel(wrfi.Variables(vv).Size)+1;
    var_sizes{a} = wrfi.Variables(vv).Size;
end

max_dims = max(var_dims);
if strcmpi(force_dims, 'same')
    % Make all concatenate along the same dimension
    var_dims(:) = max_dims;
end

% Now go through each file and import each variable into the correct cell
% in read_data, concatenating as we go.

if waitbar_bool
    wb = waitbar(0, 'Reading WRF files');
end

for d=1:n_days
    if waitbar_bool
        waitbar(d/n_days)
    end
    if DEBUG_LEVEL > 0
        fprintf('Reading file %d of %d: ', d, n_days);
    end
    for a=1:numel(varnames)
        if DEBUG_LEVEL > 0 
            fprintf('%s ', varnames{a});
        end
        var = read_wrf_preproc(fullfile(filepath, filenames{d}),varnames{a});
        
        % Check that the variables have the same dimensions as the first
        % day.
        %if (ndims(var) ~= numel(var_sizes{a}) && var_sizes{a}(end) > 1) || (var_sizes{a}(end) > 1 && any(size(var) ~= var_sizes{a}))
        %    E.callError('import','%s in the file %s is not the same size as in the first .nc file (%s vs. %s in the first)', varnames{a}, filenames{a}, mat2str(size(var)), mat2str(var_sizes{a}));
        %end
        read_data{a} = cat(var_dims(a), read_data{a}, var);
    end
    if DEBUG_LEVEL > 0
        fprintf('\n');
    end
end

if waitbar_bool
    close(wb)
end

if strcmpi(force_dims, 'squeeze')
    for i_var = 1:numel(read_data)
        read_data{i_var} = squeeze(read_data{i_var});
    end
end

if as_struct
    varargout{1} = make_struct_from_field_values(varnames, read_data);
else
    varargout = read_data;
end

end

