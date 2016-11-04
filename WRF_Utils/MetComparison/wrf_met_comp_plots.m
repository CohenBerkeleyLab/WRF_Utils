function [ varargout ] = wrf_met_comp_plots( plottype )
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
            opts{end+1} = 'All';
            opts{end+1} = 'Gen new match file';
            user_ans = ask_multichoice('Which match file to use?', opts, 'list', true);
            if strcmpi(user_ans,'gen new match file')
                fprintf('Use wrf_met_comp_plots(''match_wrf_met'') and save the structure in\n%s\n',matchdir);
                return
            elseif strcmpi(user_ans,'all')
                Match = cell(size(F));
                for a=1:numel(F)
                    M = load(fullfile(matchdir, F(a).name), 'Match');
                    Match{a} = M.Match;
                end
                legstr = {F.name}';
            else
                M = load(fullfile(matchdir, user_ans),'Match');
                Match = {M.Match};
                legstr = {user_ans};
            end
        end
        
        plot_types = {'time series','spatial'};
        plt = ask_multichoice('Which type of plot to make?', plot_types, 'list', true');
        allowed_quantities = {'winds (U+V)','winds (dir+vel)','temperature'};
        quantity = ask_multichoice('Which quantity to plot?', allowed_quantities, 'list', true);
        switch lower(plt)
            case 'time series'
                met_error_timeseries(Match, quantity, legstr);
            case 'spatial'
                if numel(Match) > 1
                    E.badinput('Plotting multiple matched runs is incompatible with a spatial plot')
                elseif iscell(Match)
                    Match = Match{1};
                end
                met_error_spatial(Match, quantity);
            otherwise
                E.notimplemented(plt)
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

    function met_error_timeseries(Match, quantity, legstr)
        % Call interactively from met_error_plots to load the match file
        % and get the quantity to plot.
        
        if strcmpi(quantity, 'winds (dir+vel)')
            rmse = cell(numel(Match),2);
        else
            rmse = cell(numel(Match),1);
        end
        for a=1:numel(Match)
            ndays = numel(Match{a}.dvec);
            switch lower(quantity)
                case 'winds (u+v)'
                    rmse{a} = wind_rmse(reshape(Match{a}.wrf_U,[],ndays), reshape(Match{a}.wrf_V,[],ndays), reshape(Match{a}.met_U,[],ndays), reshape(Match{a}.met_V,[],ndays));
                    ystr = {'RMSE in winds (U+V, m/s)'};
                case 'winds (dir+vel)'
                    [rmse{a,1}, rmse{a,2}] = wind_rmse_veldir(reshape(Match{a}.wrf_U,[],ndays), reshape(Match{a}.wrf_V,[],ndays), reshape(Match{a}.met_U,[],ndays), reshape(Match{a}.met_V,[],ndays));
                    ystr = {'RMSE in wind speed (m/s)', 'RMSE in wind direction (deg)'};
                case 'temperature'
                    rmse{a} = temperature_rmse(reshape(Match{a}.wrf_T,[],ndays), reshape(Match{a}.met_T,[],ndays));
                    ystr = {'RMSE in temperature (K)'};
                otherwise
                    E.notimplemented(quantity);
            end
        end
        
        for b=1:size(rmse,2)
            figure;
            hold on
            for a=1:size(rmse,1)
                plot(Match{a}.dvec, rmse{a,b}, 'o-');
            end
            datetick('x');
            set(gca,'XTickLabelRotation',45);
            set(gca,'FontSize',16);
            ylabel(ystr{b});
            legend(legstr{:});
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
            metdir = ask_multichoice('Choose the directory with the met_em files.', {D.name}, 'list', true);
            metdir = fullfile(sharedir, metdir);
        end
        if ~exist('avgday','var')
            user_ans = ask_multichoice('Average all available hours in a day?', {'y','n'});
            avgday = strcmp(user_ans,'y');
        end
        if avgday
            E.notimplemented('averaging to one per day');
        end
        
        % Met files are usually output every three hours, wrf files every
        % hour to half hour. Therefore, since we're supposed to be using
        % WRF-BEHR files, the only overlap for each day will be 1800 and
        % 2100 UTC.
        
        metfiles = dir(fullfile(metdir, 'met_em*'));
        metfiles = glob({metfiles.name}, '_(18|21)-00-00');
        
        % Go ahead and load lat/lon, can use to get the size of array we
        % need
        xlon = ncread(fullfile(metdir, metfiles{1}),'XLONG_M');
        xlat = ncread(fullfile(metdir, metfiles{1}),'XLAT_M');
        
        blank_mat = nan([size(xlon), numel(metfiles)]);
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
        wrffiles = dir(fullfile(wrfdir, 'WRF_BEHR*.nc'));
        wrffiles = {wrffiles.name};
        include_inds = true(numel(metfiles),1);
        for a = 1:numel(metfiles)
            dstr = regexp(metfiles{a}, '\d\d\d\d-\d\d-\d\d_\d\d-\d\d-\d\d', 'match', 'once');
            dnum = datenum(dstr, 'yyyy-mm-dd_HH-MM-SS');
            dvec(a) = dnum;
            % Load met data
            [met_U_a, met_V_a, met_T_a, met_cos, met_sin] = read_wrf_vars(metdir, metfiles(a), {'UU','VV','TT','COSALPHA','SINALPHA'});
            [met_U_a, met_V_a] = wrf_winds_transform(met_U_a, met_V_a, met_cos, met_sin);
            % Get the surface winds/temperature
            met_U(:,:,a) = met_U_a(:,:,1);
            met_V(:,:,a) = met_V_a(:,:,1);
            met_T(:,:,a) = met_T_a(:,:,1);
            
            this_wrf_file = glob(wrffiles, datestr(dnum, 'yyyy-mm-dd'));
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
            % arrays
            xx = strcmp(cellstr(wrf_times), datestr(dnum,'yyyy-mm-dd_HH:MM:SS'));
            wrf_U(:,:,a) = wrf_U_a(:,:,1,xx);
            wrf_V(:,:,a) = wrf_V_a(:,:,1,xx);
            wrf_T(:,:,a) = wrf_T_a(:,:,1,xx);
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
        Match.wrf_U = wrf_U(:,:,include_inds);
        Match.wrf_V = wrf_V(:,:,include_inds);
        Match.wrf_T = wrf_T(:,:,include_inds);
        Match.met_U = met_U(:,:,include_inds);
        Match.met_V = met_V(:,:,include_inds);
        Match.met_T = met_T(:,:,include_inds);
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

