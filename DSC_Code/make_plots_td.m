
%Make plot of varying alpha
y = [];
for tag = tags
    itag = strmatch(tag,tags,'exact');
    for dhvap = dhvaps
        idh = find(dhvaps == dhvap,1);
        for day = days
            iday = strmatch(day,days,'exact');
            
            %Start with FAME Observation
            trend_names = {'FAME-08'};
            if strcmp(days(iday),{'avg'})
                y(1,1:24) = nanmean(FAME_TDC,2);
            else
                dd = strmatch(days(iday),FAME_DAYS,'exact');
                y(1,1:24) = FAME_TDC(1:24,dd)';
            end
            
            %Cycle through all alphas
            iy = 1;
            for a = accom
                trend_names = [trend_names {['alpha = ', num2str(a)]}];
                iy = iy + 1;
                iaccom = find(accom==a,1); 
                y(iy,1:24) = squeeze(mfr(itag,idh,iaccom,1:24,iday));
            end
            t = [1:1:24];

            if strcmp(days(iday),{'avg'})
                xtitle = ['Local Time (Hour)'];
            else
                xtitle = ['HOUR OF DAY', cell2mat(days(iday))];
            end
            
            psize = plot_timeseries(y,t,[{'Mass Fraction'} {'Remaining (MFR)'}],xtitle,...
                [0, 1], trend_names, ['Time Series at Finokalia: ' cell2mat(ctags(itag))], 1);
            
            scen2 = [cell2mat(age),'/', cell2mat(tag)];
            plotname = ['MFR_AVG_dHvap',num2str(dhvap)];
            save_fig_td
        end
    end
end
            
%Make plot of varying Delta Hvap
y = [];
for tag = tags
    itag = strmatch(tag,tags,'exact');
    for a = accom
        iaccom = find(accom==a,1);
        for day = days
            iday = strmatch(day,days,'exact');
            
            %Start with FAME Observation
            trend_names = {'FAME-08'};
            if strcmp(days(iday),{'avg'})
                y(1,1:24) = nanmean(FAME_TDC,2);
            else
                dd = strmatch(days(iday),FAME_DAYS,'exact');
                y(1,1:24) = FAME_TDC(1:24,dd)';
            end
            
            %Cycle through all alphas
            iy = 1;
            for dhvap = dhvaps
                idh = find(dhvaps == dhvap,1);
                trend_names = [trend_names {['\DeltaH_{vap} = ', num2str(dhvap), ' kJ mol^{-1}']}];
                iy = iy + 1;
                y(iy,1:24) = squeeze(mfr(itag,idh,iaccom,1:24,iday));
            end
            t = [1:1:24];

            if strcmp(days(iday),{'avg'})
                xtitle = ['Local Time (Hour)'];
            else
                xtitle = ['HOUR OF DAY', cell2mat(days(iday))];
            end
            psize = plot_timeseries(y,t,[{'Mass Fraction'} {'Remaining (MFR)'}],xtitle,...
                [0, 1], trend_names, ['Time Series at Finokalia: ' cell2mat(ctags(itag))], 1);
            
            scen2 = [cell2mat(age),'/', cell2mat(tag)];
            plotname = ['MFR_AVG_alpha=',num2str(a)];
            save_fig_td
            
        end
    end
end
