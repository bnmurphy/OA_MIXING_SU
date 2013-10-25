%Make a bunch of MFR plots

for iday = 1:ndays+1
    if iday <= ndays
        y = squeeze(mfr(:,:,iday));
        ii = strmatch(scenario.date(iday), FAME_DAYS,'exact');
        y = [FAME_TDC(1:24,ii)'; y];
    else
        y = nanmean(mfr,3);
        y = [nanmean(FAME_TDC(1:24,:),2)'; y];
    end
    
    if ndays ~= 1 || iday == 1
        
        t = [1:1:24];
        trend_names = {'FAME-08'};
        for a = accom
            trend_names = [trend_names {['alpha = ', num2str(a)]}];
        end
        psize = plot_timeseries(y,t,'Mass Fraction Remaining (MFR)',xtitle,...
            [0, 1], trend_names, ['Time Series at Finokalia: ' cell2mat(scenario.tag(ll))], 1);
        
        if iday <= ndays
            plotname = ['MFR_', cell2mat(scenario.date(iday))];
        else
            plotname = ['MFR_AVG'];
        end
        save_fig
    end

    
end

        