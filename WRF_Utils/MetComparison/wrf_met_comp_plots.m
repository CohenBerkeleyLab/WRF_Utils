function [ varargout ] = wrf_met_comp_plots( plottype, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;
noaa_match_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/WRF/NOAA ISD Obs Comparison/';

switch lower(plottype)
    case 'spatial'
    case 'noaa-timeser-all'
        all_noaa_timeseries();
    case 'noaa-timeser'
        noaa_error_timeseries();
    case 'met-error'
        met_error_plots();
    case 'rmse-trend'
        rmse_trend_from_fig(varargin{1});
% Functions below here aren't really meant to be called externally, but are
% available just the same
    case 'match_wrf_met'
        varargout{1} = match_wrf_met();
end

    function plot_spatial_agreement()
        % This function will make a plot showing the agreement between WRF
        % variables and NOAA observations as a map. It can plot either
        % winds or temperature and as RMSE or correlation
        agreement_type = ask_multichoice('How to measure agreement?', {'RMSE','Correlation'},'list',true);
        quantity = ask_multichoice('Which quantity to plot?', {'Winds', 'Temperature'}, 'list', true);
    end

    function all_noaa_timeseries()
        quantity = ask_multichoice('Which quantity to plot?', {'Winds (U-V)','Winds (vel-dir)', 'Temperature'}, 'list', true);
        
        % Load all matches
        files = dir(fullfile(noaa_match_dir,'*Match*.mat'));
        legstr = cell(size(files));
        Matches = cell(size(files));
        for a=1:numel(files)
            M = load(fullfile(noaa_match_dir, files(a).name));
            fns = fieldnames(M);
            if numel(fns) > 1; warning('Only plotting first variable in %s', files(a).name); end
            Matches{a} = M.(fns{1});
            [~, legstr{a}] = fileparts(files(a).name);
        end
        
        if strcmpi(quantity, 'winds (vel-dir)')
            rmses = cell(numel(Matches),2);
        else
            rmses = cell(numel(Matches),1);
        end
        xdates = cell(numel(Matches,1));
        
        for a=1:numel(Matches)
            [xdates{a}, rmses(a,:), ystr] = noaa_error_timeseries(Matches{a}, quantity);
        end
        
        for b=1:size(rmses,2)
            figure;
            hold on
            for a=1:size(rmses,1)
                plot(xdates{a}, rmses{a,b}, 'o-','linewidth',2)
            end
            datetick('x');
            set(gca,'XTickLabelRotation',45);
            set(gca,'FontSize',16);
            ylabel(ystr{b});
            legend(legstr{:});
        end
        
    end

    function varargout = noaa_error_timeseries(Match, quantity)
        % This function will make a time series plot showing the difference
        % between WRF and observed quantities day by day. In theory, this
        % agreement should get worse over time. If there is an output, it
        % won't plot, but will return the datenums and errors (this way you
        % can use it to create timeseries for multiple cases and plot them
        % in another function).
        
        % Currently unused
        %agreement_type = ask_multichoice('How to measure agreement?', {'RMSE'},'list',true);
        
        if ~exist('quantity','var')
            quantity = ask_multichoice('Which quantity to plot?', {'Winds (U-V)','Winds (vel-dir)', 'Temperature'}, 'list', true);
        end
            
        if ~exist('Match','var')
            % Load the Match object from match_wrf_noaa
            files = dir(fullfile(noaa_match_dir,'*Match*.mat'));
            files = {files.name};
            match_file = ask_multichoice('Choose the match file to plot', files, 'list', true);
            M = load(fullfile(noaa_match_dir,match_file));
            fns = fieldnames(M);
            if numel(fns) > 1; warning('Only plotting first variable in %s', match_file); end
            Match = M.(fns{1}); clear M
        end
        ndays = size(Match.wrf_U,3);
        
        switch lower(quantity)
            case 'winds (u-v)'
                rmse = cell(1,1);
                rmse = wind_rmse(reshape(Match.wrf_U,[],ndays), reshape(Match.wrf_V, [], ndays), reshape(Match.noaa_U, [], ndays), reshape(Match.noaa_V, [], ndays));
                ystr = {'RMSE winds (U+V, m/s)'};
            case 'winds (vel-dir)'
                rmse = cell(1,2);
                [rmse{1}, rmse{2}] = wind_rmse_veldir(reshape(Match.wrf_U,[],ndays), reshape(Match.wrf_V, [], ndays), reshape(Match.noaa_U, [], ndays), reshape(Match.noaa_V, [], ndays));
                ystr = {'RMSE velocity (m/s)', 'RMSE direction (deg)'};
            case 'temperature'
                rmse = cell(1,1);
                rmse = temperature_rmse(reshape(Match.wrf_T,[],ndays), reshape(Match.noaa_T,[],ndays));
                ystr = {'RMSE temperature (K)'};
            otherwise
                E.notimplemented(quantity);
        end
        
        
        
        % Get rid of hour, minute, second in the datenum
        xdates = datenum(datestr(Match.dnums(1,:),'yyyy-mm-dd'));
        
        if nargout == 3
            varargout{1} = xdates;
            varargout{2} = rmse;
            varargout{3} = ystr;
        else
            for a=1:numel(rmse)
                figure;
                plot(xdates, rmse{a}, 'ko-');
                datetick('x','mm/dd');
                xlabel('Date')
                set(gca,'xtickLabelRotation',45);
                ylabel(ystr{a});
            end
        end
    end

    function met_error_plots(Match)
        matchdir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/WRF/NARR Comparison';
        if ~exist('Match','var')
            F=dir(fullfile(matchdir,'*.mat'));
            opts = {F.name};
            opts{end+1} = 'Gen new match file';
            user_ans = ask_multiselect('Which match file to use?', opts);
            if any(strcmpi(user_ans,'gen new match file'))
                fprintf('To generate a new match file, use wrf_met_comp_plots(''match_wrf_met'') and save the structure in\n%s\n',matchdir);
                return
            else
                Match = cell(1,numel(user_ans));
                legstr = cell(1,numel(user_ans));
                for a=1:numel(user_ans)
                    M = load(fullfile(matchdir, user_ans{a}),'Match');
                    Match{a} = M.Match;
                    legstr{a} = user_ans{a};
                end
            end
        end
        % Check if we're loading Match files with vertical information or
        % not
        d3_matches = false(size(Match));
        for a=1:numel(d3_matches)
            d3_matches(a) = ndims(Match{a}.wrf_U)>3;
        end
        if ~all(d3_matches) && ~all(~d3_matches)
            E.callError('inconsistent_matches','Loaded Matches have inconsistent dimensions')
        end
        
        plot_types = {'time series','spatial','spatial slices'};
        plt = ask_multichoice('Which type of plot to make?', plot_types, 'list', true');
        allowed_quantities = {'winds (U+V)','winds (dir+vel)','temperature'};
        quantity = ask_multichoice('Which quantity to plot?', allowed_quantities, 'list', true);
        if all(d3_matches)
            nlevels = size(Match{1}.wrf_U,3);
            levels = ask_number('Enter what model levels to plot the RMSE for (separated by a space)','default',1,'testfxn',@(x) any(x >= 1 & x <= 29 & mod(x,1) == 0),'testmsg',...
                sprintf('All values must be between 1 and %d and must be integers',29));
        else
            levels = 1;
        end
        switch lower(plt)
            case 'time series'
                met_error_timeseries(Match, quantity, legstr, levels);
            case 'spatial'
                if numel(Match) > 1
                    E.badinput('Plotting multiple matched runs is incompatible with a spatial plot')
                elseif iscell(Match)
                    Match = Match{1};
                end
                met_error_spatial(Match, quantity);
            case 'spatial slices'
                met_error_slices(Match, quantity, {F.name});
            otherwise
                E.notimplemented(plt)
        end
    end

    function met_error_slices(Match_in, quantity, matchnames)
        % Call interactively from met_error_plots.
        % First we need to find which Match object has the shortest date
        % vector; we will restrict all others to that.
        
        dvec = Match_in{1}.dvec;
        for a=2:numel(Match_in)
            if numel(Match_in{a}.dvec) < numel(dvec)
                dvec = Match_in{a}.dvec;
            end
        end
        
        [slon, slat] = state_outlines('not', 'ak', 'hi');
        
        % Now get the difference of either wind speed or wind direction
        for a=1:numel(Match_in)
            Match = Match_in{a};
            tt = ismember(Match.dvec, dvec);
            if sum(tt) ~= numel(dvec)
                E.callError('Not all dates in the smallest datevec are present in %s', matchnames{a});
            end
            if strcmpi(quantity, 'winds (dir+vel)')
                wrf_vel = sqrt(Match.wrf_U(:,:,tt) .^ 2 + Match.wrf_V(:,:,tt) .^ 2);
                met_vel = sqrt(Match.met_U(:,:,tt) .^ 2 + Match.met_V(:,:,tt) .^ 2);
                diff_vel = wrf_vel - met_vel;
                
                wrf_dir = atan2d(Match.wrf_V(:,:,tt), Match.wrf_U(:,:,tt));
                met_dir = atan2d(Match.met_V(:,:,tt), Match.met_U(:,:,tt));
                diff_dir = angle_diffd(wrf_dir, met_dir);
            else
                E.notimplemented(quantity)
            end
            
            plot_slice_gui(diff_vel, Match.lon, Match.lat, slon, slat, sprintf('%s - velocity', matchnames{a}));
            plot_slice_gui(diff_dir, Match.lon, Match.lat, slon, slat, sprintf('%s - direction', matchnames{a}));
        end
            
    end

    function met_error_spatial(Match, quantity)
        % Call interactively from met_error_plots to load the match file
        % and get the quantity to plot.
        blank_mat = nan(size(Match.lon));
        if strcmpi(quantity, 'winds (dir+vel)')
            rmse = cell(1,2);
            rmse(:) = {blank_mat};
        else
            rmse = {blank_mat};
        end
        for a=1:size(rmse{1},1)
            for b=1:size(rmse{1},2)
                if ~isempty(regexp(quantity,'winds','once'))
                    wrf_U_slice = squeeze(Match.wrf_U(a,b,:));
                    wrf_V_slice = squeeze(Match.wrf_V(a,b,:));
                    narr_U_slice = squeeze(Match.met_U(a,b,:));
                    narr_V_slice = squeeze(Match.met_V(a,b,:));
                    if strcmpi(quantity, 'winds (U+V)')
                        rmse{1}(a,b) = wind_rmse(wrf_U_slice, wrf_V_slice, narr_U_slice, narr_V_slice);
                        cblabel{1} = 'RMSE wind U+V (m/s)';
                    elseif strcmpi(quantity, 'winds (dir+vel)');
                        [rmse{1}(a,b), rmse{2}(a,b)] = wind_rmse_veldir(wrf_U_slice, wrf_V_slice, narr_U_slice, narr_V_slice);
                        cblabel{1} = 'RMSE wind speed (m/s)';
                        cblabel{2} = 'RMSE wind direction (deg)';
                    else
                        E.notimplemented(quantity)
                    end
                elseif strcmpi(quantity, 'temperature')
                    wrf_T_slice = squeeze(Match.wrf_T(a,b,:));
                    narr_T_slice = squeeze(Match.met_T(a,b,:));
                    rmse{1}(a,b) = temperature_rmse(wrf_T_slice, narr_T_slice);
                    cblabel{1} = 'RMSE temperature (K)';
                else
                    E.notimplemented(quantity)
                end
            end
        end
        
        for f=1:numel(rmse)
            figure;
            pcolor(Match.lon, Match.lat, rmse{f});
            colormap jet
            cb = colorbar;
            cb.Label.String = cblabel{f};
            shading flat
            state_outlines('k','not','ak','hi');
            set(gca,'fontsize',16)
        end
    end

    function met_error_timeseries(Match, quantity, legstr, levels)
        % Call interactively from met_error_plots to load the match file
        % and get the quantity to plot.
        
        if ~exist('levels','var')
            levels = 1;
        end
        
        % Handle 3D (rather than 4D) match structures by creating a
        % singleton dimension 3
        for a=1:numel(Match)
            fns = fieldnames(Match{a});
            for f=1:numel(fns)
                if ndims(Match{a}.(fns{f})) == 3
                    sz = size(Match{a}.(fns{f}));
                    Match{a}.(fns{f}) = reshape(Match{a}.(fns{f}), [sz(1:3),1,sz(3)]);
                end
            end
        end
        
        levels_inds = zeros(numel(Match), numel(levels));
        for a=1:numel(Match)
            if isfield(Match{a}, 'model_levels')
                match_levels = Match{a}.model_levels;
            else
                match_levels = 1:size(Match{a}.wrf_U,3);
            end
            for b=1:numel(levels)
                li = find(levels(b)==match_levels);
                if ~isempty(li)
                    levels_inds(a,b) = li;
                end
            end
        end
        
        if strcmpi(quantity, 'winds (dir+vel)')
            rmse = cell(numel(Match),numel(levels),2);
        else
            rmse = cell(numel(Match),numel(levels),1);
        end
        for a=1:numel(Match)
            ndays = numel(Match{a}.dvec);
            for b=1:numel(levels)
                if levels_inds(a,b) > 0
                    l = levels_inds(a,b);
                else
                    % Level not present in this Match object, skip
                    continue
                end
                switch lower(quantity)
                    case 'winds (u+v)'
                        rmse{a,b} = wind_rmse(reshape(Match{a}.wrf_U(:,:,l,:),[],ndays), reshape(Match{a}.wrf_V(:,:,l,:),[],ndays), reshape(Match{a}.met_U(:,:,l,:),[],ndays), reshape(Match{a}.met_V(:,:,l,:),[],ndays));
                        ystr = {'RMSE in winds (U+V, m/s)'};
                    case 'winds (dir+vel)'
                        [rmse{a,b,1}, rmse{a,b,2}] = wind_rmse_veldir(reshape(Match{a}.wrf_U(:,:,l,:),[],ndays), reshape(Match{a}.wrf_V(:,:,l,:),[],ndays), reshape(Match{a}.met_U(:,:,l,:),[],ndays), reshape(Match{a}.met_V(:,:,l,:),[],ndays));
                        ystr = {'RMSE in wind speed (m/s)', 'RMSE in wind direction (deg)'};
                    case 'temperature'
                        rmse{a,b} = temperature_rmse(reshape(Match{a}.wrf_T(:,:,l,:),[],ndays), reshape(Match{a}.met_T(:,:,l,:),[],ndays));
                        ystr = {'RMSE in temperature (K)'};
                    otherwise
                        E.notimplemented(quantity);
                end
            end
        end
        
        levels_sym = {'o','s','d','v'};
        levels_lstyle = {'-','--',':','-.'};
        starttime_colors = {'b','r',[0 0.5 0],'m'};
        
        for i=1:size(rmse,3)
            l = gobjects(size(rmse,1), size(rmse,2));
            lstr = cell(size(l));
            figure;
            for b=1:size(rmse,2)
                bmod = mod(b-1, numel(levels_sym))+1;
                marker = levels_sym{bmod};
                lstyle = levels_lstyle{bmod};
                for a=1:size(rmse,1)
                    % You'll want to add more colors if you need to plot
                    % more starting times
                    sercol = starttime_colors{a};
                    l(a,b) = line(Match{a}.dvec, rmse{a,b,i}, 'linestyle',lstyle,'marker',marker,'color',sercol,'linewidth',2);
                    lstr{a,b} = sprintf('%s, level %d', legstr{a}, levels(b));
                end
            end
            datetick('x');
            set(gca,'XTickLabelRotation',45);
            set(gca,'FontSize',16);
            ylabel(ystr{i});
            legend(l(:),lstr(:));
        end
    end

    function rmse_trend_from_fig(axes)
        ch = get(axes, 'children');
        for a=1:numel(ch)
            dvec = ch(a).XData;
            x = dvec - floor(dvec(1));
            rmse = ch(a).YData;
            
            figure;
            plot(x, rmse, 'ko');
            % There is no uncertainty in the date, so use y-residual
            % fitting
            plot_fit_line(x, rmse, 'regression', 'y-resid', 'one2one', false);
            set(gca,'fontsize', 16);
            ylabel(axes.YLabel.String);
            title(ch(a).DisplayName);
        end
    end

    function Match = match_wrf_met(wrfdir, metdir, avgday)
        % This function will compare WRF output meteorology to that
        % contained in the met_em files that are derived from NARR
        % directly.
        
        sharedir='/Volumes/share2/USERS/LaughnerJ/WRF/NudgeTest-US/';
        D = dir(sharedir);
        D(1:2) = []; % remove the . and .. entries
        
        if ~exist('wrfdir','var')
            wrfdir = ask_multichoice('Choose the directory with the WRF output.', {D.name}, 'list', true);
            wrfdir = fullfile(sharedir, wrfdir);
        end
        if ~exist('metdir','var')
            metdir = ask_multichoice('Choose the directory with the met_em files or wrfinput_subset files.', {D.name}, 'list', true);
            metdir = fullfile(sharedir, metdir);
        end
        if ~exist('surf_bool','var')
            levels = ask_number('Enter the levels to include in the output, separated by space. Valid levels are 1-29', 'testfxn', @(x) all(x>=1 & x<=29), 'testmsg', 'All level indices must be between 1 and 29');
        end
%         if ~exist('avgday','var')
%             user_ans = ask_multichoice('Average all available hours in a day to a single value?', {'y','n'});
%             avgday = strcmp(user_ans,'y');
%         end
%         if avgday
%             E.notimplemented('averaging to one per day');
%         end
        
        % Met files are usually output every three hours, wrf files every
        % hour to half hour. If we're supposed to be using WRF-BEHR files,
        % the only overlap for each day will be 1800 and 2100 UTC.
        % Otherwise, if we have the wrfout_subset files, they should be
        % available every three hours
        
        wrffiles = dir(fullfile(wrfdir, 'wrfout_subset*'));
        if ~isempty(wrffiles)
            all_hours = true;
        else
            wrffiles = dir(fullfile(wrfdir, 'WRF_BEHR*.nc'));
            all_hours = false;
            if isempty(wrffiles)
                E.filenotfound('wrfout_subset* or WRF_BEHR*.nc')
            end
        end
        
        metfiles = dir(fullfile(metdir, 'wrfinput_subset*'));
        met_varnames.lon = 'XLONG';
        met_varnames.lat = 'XLAT';
        met_varnames.U = 'U';
        met_varnames.V = 'V';
        met_varnames.T = 'TT';
        if isempty(metfiles)
            metfiles = dir(fullfile(metdir, 'met_em*'));
            met_varnames.lon = 'XLONG_M';
            met_varnames.lat = 'XLAT_M';
            met_varnames.U = 'UU';
            met_varnames.V = 'VV';
            met_varnames.T = 'TT';
            if ~isempty(metfiles)
                warning('Could not find wrfinput_subset files, using met_em files. These are not properly interpolated to WRF vertical coordinates, so the comparisons may have greater error than expected.');
            else
                E.filenotfound('wrfinput_subset* or met_em*');
            end
        end
        if ~all_hours
            metfiles = glob({metfiles.name}, '_(18|21)-00-00');
        else
            metfiles = {metfiles.name};
        end
        
        % Go ahead and load lat/lon, can use to get the size of array we
        % need
        xlon = ncread(fullfile(metdir, metfiles{1}),met_varnames.lon);
        xlat = ncread(fullfile(metdir, metfiles{1}),met_varnames.lat);
        
        wrf_size = get_wrf_array_size(fullfile(metdir, metfiles{1}));
        
        blank_mat = nan([wrf_size(1:2), numel(levels), numel(metfiles)]);
        
        wrf_U = blank_mat;
        wrf_V = blank_mat;
        wrf_T = blank_mat;
        met_U = blank_mat;
        met_V = blank_mat;
        met_T = blank_mat;
        dvec = nan(numel(metfiles),1);
        
        % Loop through the met files. For each find the corresponding day's
        % WRF_BEHR file. 
        last_wrf_file = '';
        wrffiles = {wrffiles.name};
        include_inds = true(numel(metfiles),1);
        for a = 1:numel(metfiles)
            dstr = regexp(metfiles{a}, '\d\d\d\d-\d\d-\d\d_\d\d-\d\d-\d\d', 'match', 'once');
            dnum = datenum(dstr, 'yyyy-mm-dd_HH-MM-SS');
            dvec(a) = dnum;
            % Load met data
            [met_U_a, met_V_a, met_T_a, met_cos, met_sin] = read_wrf_vars(metdir, metfiles(a), {met_varnames.U,met_varnames.V,met_varnames.T,'COSALPHA','SINALPHA'});
            [met_U_a, met_V_a] = wrf_winds_transform(met_U_a, met_V_a, met_cos, met_sin);
            
            met_U(:,:,:,a) = met_U_a(:,:,levels);
            met_V(:,:,:,a) = met_V_a(:,:,levels);
            met_T(:,:,:,a) = met_T_a(:,:,levels);
            
            if ~all_hours
                this_wrf_file = glob(wrffiles, datestr(dnum, 'yyyy-mm-dd'));
            else
                this_wrf_file = glob(wrffiles, datestr(dnum, 'yyyy-mm-dd_HH[:-]MM[:-]SS'));
            end
            if isempty(this_wrf_file)
                include_inds(a) = false;
                continue
            end
            this_wrf_file = this_wrf_file{1};
            if ~strcmp(this_wrf_file, last_wrf_file)
                [wrf_U_a, wrf_V_a, wrf_T_a, wrf_cos, wrf_sin, wrf_times] = read_wrf_vars(wrfdir, {this_wrf_file}, {'U','V','TT','COSALPHA','SINALPHA','Times'});
                wrf_cos = wrf_cos(:,:,1);
                wrf_sin = wrf_sin(:,:,1);
                [wrf_U_a, wrf_V_a] = wrf_winds_transform(wrf_U_a, wrf_V_a, wrf_cos, wrf_sin);
                wrf_times = wrf_times';
                last_wrf_file = this_wrf_file;
            end
            
            % Put the proper slices of WRF data in the output
            % arrays. Extraneous if using wrfout_subset files that only
            % include a single time, but if multiple times included in the
            % output, this is necessary.
            xx = strcmp(cellstr(wrf_times), datestr(dnum,'yyyy-mm-dd_HH:MM:SS'));
            if sum(xx) < 1
                E.callError('time_not_found', 'Time %s not found in WRF output file %s', datestr(dnum, 'yyyy-mm-dd_HH:MM:SS'), this_wrf_file);
            end
            
            wrf_U(:,:,:,a) = wrf_U_a(:,:,levels,xx);
            wrf_V(:,:,:,a) = wrf_V_a(:,:,levels,xx);
            wrf_T(:,:,:,a) = wrf_T_a(:,:,levels,xx);
            
        end
        
        % If we want to average down to one value per day, then we need to
        % identify which slices belong to the same day
        if avgday
            E.notimplemented('averaging to one per day');
            %dvec_day = unique(floor(dvec));
            %for d=1:numel(dvec_day)
            %   dd = floor(dvec) == dvec_day(d)
            %   ...
        end
        
        Match.lon = xlon;
        Match.lat = xlat;
        Match.dvec = dvec(include_inds);
        Match.model_levels = levels;
        
        Match.wrf_U = wrf_U(:,:,:,include_inds);
        Match.wrf_V = wrf_V(:,:,:,include_inds);
        Match.wrf_T = wrf_T(:,:,:,include_inds);
        Match.met_U = met_U(:,:,:,include_inds);
        Match.met_V = met_V(:,:,:,include_inds);
        Match.met_T = met_T(:,:,:,include_inds);
        
    end
end

function rmse = wind_rmse(wrf_U, wrf_V, true_U, true_V)
% Computes root mean squared error for the winds (as a sum of the errors in
% the winds). RMSE will be calculated for each column in the inputs.
del_u = wrf_U - true_U;
del_v = wrf_V - true_V;
sq_del = nansum2(del_u .^ 2 + del_v .^ 2, 1);
valid_obs = ~isnan(wrf_U) & ~isnan(wrf_V) & ~isnan(true_U) & ~isnan(true_V);
n_obs = nansum2(valid_obs,1);
rmse = sqrt ( sq_del ./ n_obs );
end

function [rmse_vel, rmse_dir] = wind_rmse_veldir(wrf_U, wrf_V, true_U, true_V)
wrf_vel = sqrt(wrf_U .^ 2 + wrf_V .^ 2);
wrf_dir = atan2d(wrf_V, wrf_U);
true_vel = sqrt(true_U .^ 2 + true_V .^ 2);
true_dir = atan2d(true_V, true_U);

del_vel = wrf_vel - true_vel;
del_dir = angle_diffd(wrf_dir, true_dir);

sq_delvel = nansum2( del_vel .^ 2, 1);
sq_deldir = nansum2( del_dir .^ 2, 1);

valid_obs = ~isnan(wrf_U) & ~isnan(wrf_V) & ~isnan(true_U) & ~isnan(true_V);
n_obs = nansum2(valid_obs,1);

rmse_vel = sqrt( sq_delvel ./ n_obs );
rmse_dir = sqrt( sq_deldir ./ n_obs );
end

function rmse = temperature_rmse(wrf_T, true_T)
% Computes root mean squared error for the temperatures. RMSE will be
% calculated for each column in the inputs.
del_t = wrf_T - true_T;
sq_del = nansum2( del_t .^ 2, 1);
valid_obs = ~isnan(wrf_T) & ~isnan(true_T);
n_obs = nansum2(valid_obs,1);
rmse = sqrt ( sq_del ./ n_obs );
end

