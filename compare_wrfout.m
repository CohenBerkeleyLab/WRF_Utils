function [  ] = compare_wrfout( old_dir, new_dir, varargin )
%COMPARE_WRFOUT Compare several variables from two sets of WRF-Chem output
%   COMPARE_WRFOUT( OLD_DIR, NEW_DIR ) Reads wrfout_d* files from OLD_DIR
%   and NEW_DIR and plots the average absolute and percent difference in
%   no, no2, o3, U, and V across all files common to both directories. If
%   "plot_slice_gui" exists on your Matlab path, then it will be used to
%   allow you to examine each model level at your leisure. Otherwise,
%   pcolor will be used to plot the top and bottom levels. A list of
%   statistics will also be printed out.
%
%   Additional parameters:
%       'files' - a cell array of strings listing individual files
%       that you wish to load. If any are missing from either OLD_DIR or
%       NEW_DIR, an error is thrown. By default all files in common between
%       the two directories are loaded.
%
%       'vars' - a cell array of strings listing which variables you wish
%       to compare. Default is {'no', 'no2', 'o3', 'U', 'V'}. Variables
%       will automatically be unstagged.
%
%       'filestem' - the beginning of the WRF-Chem output file names, used
%       to filter out unrelated files. Default is 'wrfout_d'.
%
%       'maxvargb' - the maximum size a single variable loaded from two
%       files can have in gigabytes. Because WRF-Chem output can get crazy
%       large, this function estimates how large each variable loaded will
%       be and halts if it exceeds this value to prevent grinding your poor
%       computer to a halt.

E = JLLErrors;

p = inputParser;
p.addParameter('files',[]);
p.addParameter('vars',{'no','no2','o3','U','V'});
p.addParameter('filestem','wrfout_d');
p.addParameter('maxvargb',2);

p.parse(varargin{:});
pout = p.Results;

files = pout.files;
vars = pout.vars;
filestem = strcat(pout.filestem,'*');
maxvarGB = pout.maxvargb;

if ~iscellstr(vars)
    E.badinput('Parameter "vars" must be a cell array of strings')
elseif ~ischar(filestem)
    E.badinput('Parameter "filestem" must be a string')
end

% If a list of files is NOT given, find all wrfout files common to both the
% old and new directories. Make sure that no date/time is repeated.
if isempty(files)
    Fold = dir(fullfile(old_dir,filestem));
    fold_names = {Fold.name};
    Fnew = dir(fullfile(new_dir,filestem));
    fnew_names = {Fnew.name};
    xx = find_common_elements(fold_names, fnew_names, 'nodup');
    files = fnew_names(xx);
elseif ~iscellstr(files)
    E.badinput('The value for parameter "files" must be a cell array of strings')
else
    missing_files_new = {};
    missing_files_old = {};
    for a=1:numel(files)
        if ~exist(fullfile(new_dir,files{a}),'file')
            missing_files_new{end+1} = files{a}; %#ok<AGROW>
        end
        if ~exist(fullfile(old_dir,files{a}),'file')
            missing_files_old{end+1} = files{a}; %#ok<AGROW>
        end
    end
    doerr = false;
    if numel(missing_files_new) > 0
        fprintf('The following files are missing from %s: \n\t%s\n', new_dir, strjoin(missing_files_new, '\n\t'));
        doerr = true;
    end
    if numel(missing_files_old) > 0
        fprintf('The following files are missing from %s: \n\t%s\n', old_dir, strjoin(missing_files_old, '\n\t'));
        doerr = true;
    end
    if doerr
        E.filenotfound('%s', 'Files specified by the parameter "files" were missing. See previous messages for a list.')
    end
end

% Check the array size of the first file, use to issue a warning if too
% much data will be read in
wrf_sz = get_wrf_array_size(fullfile(old_dir, files{1}));
gigabytes = prod(wrf_sz) * 4 * numel(files) * 2 / 1e9; % read in as single precision (4 bytes/value) and each file name is loaded from 2 directories
if gigabytes > maxvarGB
    E.callError('var_too_large','Loading just one variable for all files would require %.2f GB, > %.2f GB (you can change this limit with the "maxvargb" parameter', gigabytes, maxvarGB)
end

xlon = double(ncread(fullfile(old_dir, files{1}),'XLONG'));
xlat = double(ncread(fullfile(old_dir, files{1}),'XLAT'));
coast_dat = load('coast');
if exist('plot_slice_gui','file')
    plot_fxn = @(del,tstr) plot_with_gui(del,tstr);
else
    plot_fxn = @(del,tstr) plot_with_pcolor(del,tstr);
end

% Load each variable; calculate the difference between the old and new WRF
% runs, then plot (either using the plot slice gui or just pcolor)

substruct = struct('mean_diff_by_level',[],'sdev_diff_by_level',[],'mean_perdiff_by_level',[],'sddev_perdiff_by_level',[],...
    'mean_absdiff_by_level',[],'sdev_absdiff_by_level',[],'mean_absperdiff_by_level',[],'sdev_absperdiff_by_level',[]);
stats = make_empty_struct_from_cell(vars,substruct);



for a=1:numel(vars)
    this_var_old = read_wrf_vars(old_dir,files,vars(a),false,'visual');
    this_var_new = read_wrf_vars(new_dir,files,vars(a),false,'visual');
    stagger_id = ncreadatt(fullfile(old_dir,files{1}),vars{a},'stagger');
    
    if ~isequal(size(this_var_new), size(this_var_old))
        E.sizeMismatch('this_var_old',sprintf('this_var_new (%s arrays concatenated from all WRF files)',vars{a}))
    end
    
    del = this_var_new - this_var_old;
    [this_mean, this_std] = diff_stats(del, stagger_id);
    stats.(vars{a}).mean_diff_by_level = this_mean;
    stats.(vars{a}).sdev_diff_by_level = this_std;
    
    [this_mean, this_std] = diff_stats(abs(del), stagger_id);
    stats.(vars{a}).mean_absdiff_by_level = this_mean;
    stats.(vars{a}).sdev_absdiff_by_level = this_std;
    
    plot_fxn(del,sprintf('Absolute difference %s',vars{a}));
    
    del = reldiff(this_var_new, this_var_old)*100;
    [this_mean, this_std] = diff_stats(del, stagger_id);
    stats.(vars{a}).mean_perdiff_by_level = this_mean;
    stats.(vars{a}).sdev_perdiff_by_level = this_std;
    
    [this_mean, this_std] = diff_stats(abs(del), stagger_id);
    stats.(vars{a}).mean_absperdiff_by_level = this_mean;
    stats.(vars{a}).sdev_absperdiff_by_level = this_std;
    
    plot_fxn(del,sprintf('Percent difference %s',vars{a}));
end

% Print the stats
for a=1:numel(vars)
    fprintf('%s\n:',vars{a});
    fns = fieldnames(stats.(vars{a}));
    for b=1:numel(fns)
        fprintf('  %s: %g\n', fns{b}, stats.(vars{a}).(fns{b}));
    end
end

    function plot_with_gui(del, tstr)
        del = double(nanmean(del(:,:,:,:),4)); %average along all time dimensions
        plot_slice_gui(del, xlon, xlat, coast_dat.long, coast_dat.lat, tstr);
    end

    function plot_with_pcolor(del,tstr)
        del = double(nanmean(del(:,:,:,:),4)); %average along all time dimensions
        figure;
        pcolor(xlon, xlat, del(:,:,1));
        line(coast_dat.long, coast_dat.lat,'color','k');
        colorbar;
        title(sprintf('%s - level 1',tstr));
        
        figure;
        pcolor(xlon, xlat, del(:,:,end));
        line(coast_dat.long, coast_dat.lat,'color','k');
        colorbar;
        title(sprintf('%s - level %d',tstr,size(del,3)));
    end

end

function [m,s] = diff_stats(del, unstagger_dim)
if ~isempty(unstagger_dim)
    del = unstagger(del, unstagger_dim);
end
pvec = 1:ndims(del);
pvec(3) = [];
pvec = cat(2,3,pvec);
del = permute(del,pvec);
del = reshape(del,size(del,1),[]);
m = nanmean(del,2);
s = nanstd(del,0,2);
end


