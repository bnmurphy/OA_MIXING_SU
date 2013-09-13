%Write data out

if size(scenario.date) > 1, date = 'avg'; else date = cell2mat(scenario.date); end
scen2 = [cell2mat(scenario.age),'/', cell2mat(scenario.tag),'/dHvap=' num2str(dhvap)];

if ~exist(['../Figs/', scen2],'dir'); mkdir(['../Figs/', scen2]); end;

for iaccom = 1:length(accom)
    fname = ['../Figs/', scen2, '/', 'td.data.alpha_', num2str(accom(iaccom))];
    fid = fopen(fname,'w');
    
    for iday = 1:ndays+1
        for ihr = 1:24
            if iday <= ndays
                y = squeeze(mfr(iaccom,:,iday));
                sday = cell2mat(scenario.date(iday));
            else
                y = nanmean(mfr(iaccom,1:24,:),3);
                sday = 'avg';
            end
            fprintf(fid, '\n%s  %2.0f  %6.5f',sday,ihr,y(ihr));
        end
    end
    fclose(fid);
end
