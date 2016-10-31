function [ varargout ] = wrf_noaa_comp_plots( plottype )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;

switch lower(plottype)
    case 'spatial'
    case 'timeser'
        error_timeseries();
end

    function plot_spatial_agreement()
        % This function will make a plot showing the agreement between WRF
        % variables and NOAA observations as a map. It can plot either
        % winds or temperature and as RMSE or correlation
        agreement_type = ask_multichoice('How to measure agreement?', {'RMSE','Correlation'},'list',true);
        quantity = ask_multichoice('Which quantity to plot?', {'Winds', 'Temperature'}, 'list', true);
    end

    function varargout = error_timeseries()
        % This function will make a time series plot showing the difference
        % between WRF and observed quantities day by day. In theory, this
        % agreement should get worse over time. If there is an output, it
        % won't plot, but will return the datenums and errors (this way you
        % can use it to create timeseries for multiple cases and plot them
        % in another function).
        
        % Currently unused
        %agreement_type = ask_multichoice('How to measure agreement?', {'RMSE'},'list',true);
        quantity = ask_multichoice('Which quantity to plot?', {'Winds', 'Temperature'}, 'list', true);
        
        % Load the Match object from match_wrf_noaa
        match_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/WRF/NOAA ISD Obs Comparison/';
        files = dir('*Match*.mat');
        files = {files.name};
        match_file = ask_multichoice('Choose the match file to plot', files, 'list', true);
        M = load(fullfile(match_dir,match_file));
        fns = fieldnames(M);
        if numel(fns) > 1; warning('Only plotting first variable in %s', match_file); end
        Match = M.(fns{1}); clear M
        
        ndays = size(Match.wrf_U,3);
        
        switch lower(quantity)
            case 'winds'
                del_u = Match.wrf_U - Match.noaa_U;
                del_v = Match.wrf_V - Match.noaa_V;
                sq_del = nansum2(reshape(del_u .^ 2,[],ndays) + reshape(del_v .^ 2, [], ndays), 1);
                valid_obs = ~isnan(Match.wrf_U) & ~isnan(Match.wrf_V) & ~isnan(Match.noaa_U) & ~isnan(Match.noaa_V);
                ystr = 'RMSE winds (U+V, m/s)';
            case 'temperature'
                del_t = Match.wrf_T - Match.noaa_T;
                sq_del = nansum2( reshape(del_t .^ 2, [], ndays), 1);
                valid_obs = ~isnan(Match.wrf_T) & ~isnan(Match.noaa_T);    
                ystr = 'RMSE temperature (K)';
            otherwise
                E.notimplemented(quantity);
        end
        
        n_obs = nansum2(reshape(valid_obs,[],ndays),1);
        rmse = sqrt ( sq_del ./ n_obs );
        
        % Get rid of hour, minute, second in the datenum
        xdates = datenum(datestr(Match.dnums(1,:),'yyyy-mm-dd'));
        
        if nargout == 2
            varargout{1} = xdates;
            varargout{2} = rmse;
        else
            figure;
            plot(xdates, rmse, 'ko-');
            datetick('x','mm/dd');
            xlabel('Date')
            set(gca,'xtickLabelRotation',45);
            ylabel(ystr);
        end
    end

end

