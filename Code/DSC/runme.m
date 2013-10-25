% Wrapper for running Ilona's thermodenuder dynamic model with the results
% from the trajectory model.


%Define Scenario and Run Days
scenario.tags = {'2bin_ksivoc_1stgen'};
dhvaps = [50, 75, 100];


for dhvap = dhvaps
%     if dhvap == 50
%         scenario.tags = {'FUNCNEIL_nohet_nowet',...
%             'FRAG_0.8_nohet_nowet','FRAG_VARY_nohet_nowet',...
%             'FRAG_0.3_nohet_nowet','FRAG_0.5_bsoa_nohet_nowet'};
%     else
%         scenario.tags = {'FUNC100_nohet_nowet','2bin_nohet_nowet',...
%             'FUNCNEIL_nohet_nowet',...
%             'FRAG_0.8_nohet_nowet','FRAG_VARY_nohet_nowet',...
%             'FRAG_0.3_nohet_nowet','FRAG_0.5_bsoa_nohet_nowet'};
%     end
for tag = scenario.tags
    
    dhvap
    tag
    
    scenario.tag = tag;
    scenario.age = {'FULL_HYSP'};
    scenario.loc = {'fino'};
    scenario.date = {'130','133','136','137','138',...
        '147','148','149','150'};
    
    nvol = 9;
    nsize = 5;
    ndays = size(scenario.date,2);
    
    set_meas
    
    spec_list = [];
    for ivol = 1:nvol
        for isize = 1:nsize
            spec_list = [spec_list, {['TD_0',num2str(ivol),num2str(isize)]}];
        end
        spec_list = [spec_list, {['TD_0',num2str(ivol),'V']}];
    end
    
    %Get All Size and Volatility Data
    get_size_dist
    
    %Set Thermodenuder Inputs
    vol_dist = zeros(nvol,6,ndays);
    size_dist = zeros(nsize,6,ndays);
    for dd = 1:ndays
        for ivol = 1:nvol
            for isize = 1:nsize
                ii = strmatch({['TD_0',num2str(ivol),num2str(isize)]},spec_list);
                
                vol_dist(ivol,1:6,dd) = vol_dist(ivol,1:6,dd) + traj_td(ii,:,dd);
                size_dist(isize,1:6,dd) = size_dist(isize,1:6,dd) + traj_td(ii,:,dd);
            end
            %Calculate Aerosol Mass Fraction
            ii = strmatch({['TD_0',num2str(ivol),'V']},spec_list);
            %         vol_dist(ivol,1:24,dd) = vol_dist(ivol,1:24,dd) + traj_td(ii,:,dd);
        end
        
        %Convert Volatility Distribution to mass fraction
        for ihr = 1:6
            vol_dist(:,ihr,dd) = vol_dist(:,ihr,dd) ./ sum(vol_dist(:,ihr,dd));
        end
    end
    
    %Average the Volatility and Size Distributions Across Hours And Days for
    %now
    % size_dist_use = nanmean(nanmean(size_dist,3),2)'
    % vol_dist_use = nanmean(nanmean(vol_dist,3),2)'
    % size_dist_use = size_dist_use .* 1e-9;  %ug/m3 -> kg/m3
    % mfr_dist = model_size_dist(vol_dist_use, size_dist_use)
    
    
    %Run Model
    accom = [0.01, 0.1, 1];
    for iaccom = 1:length(accom);
        for iday = 1:ndays
            for ihr = 1:6
                vol_dist_use = vol_dist(:,ihr,iday)';
                size_dist_use = size_dist(:,ihr,iday)';
                T_TD = FAME_TDT(ihr,iday);
                iday
                ihr
                accom(iaccom)
                
                if ~isnan(T_TD)
                    size_dist_use = size_dist_use .* 1e-9;  %ug/m3 -> kg/m3
                    mfr_dist = model_size_dist(vol_dist_use, size_dist_use, accom(iaccom), T_TD, dhvap)
                else
                    mfr_dist = NaN;
                end
                hr_map = [1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 1 1];
                indx = find(hr_map == ihr);
                for kk = indx
                    mfr(iaccom,kk,iday) = mfr_dist;
                end
            end
        end
    end
    
    write_data
    %Make Plots
%     if dom_vars.date == 1
%         xtitle = ['HOUR OF DAY', cell2mat(scenario.date)];
%     else
%         xtitle = ['HOUR OF DAY'];
%     end
    
    %Timeseries of MFR at TD Temp
%     make_plots
%     close all

end
end





