function plot_2D(x, y, lmeas, xstr, ystr, tstr, xlims, ylims, ...
    plotDir, plotName, legstr, lylog, lsave)
%This function makes a 2-D plot.
%
%It is not dependent on any outside information except VERY relevant
%information for the plot. It needs to know what the data is and what
%labels to put on. It does the rest.
%
%Ben Murphy Jan 2013


%Retrieve Font and Axis Properties
set_plot_props

%Create Figure
full_page = [2.5 1 26 19];
half_page = [2.5 1 14.25 10];
quarter_page = [2.5 1 9 8.5];
tall_page = [2.5 1 7 13.5];
long_page = [1 1 26.5 5];
figure('units','normalized',...
             'position',[0.2 0.3 0.6 0.35],...
             'visible','on',...
             'paperorientation','portrait',...
             'papertype','a4',...
             'paperunit','centimeters',...
             'paperposition',quarter_page);

%Draw Plot
line_color = {'k',[1 0 0],[0 0 1],[0 145/255 0],[153 0 153]./255};

for itrend = 1:length(x)
    x1 = cell2mat(x(itrend));
    y1 = cell2mat(y(itrend));
    lc = cell2mat(line_color(itrend));
    lm = cell2mat(lmeas(itrend));
    
    h(itrend) = plot(x1,y1, 'LineWidth',3,'color',lc);
    hold on
    if lm,
        set(h(itrend),'linestyle','none','marker','o','markersize',2,...
            'markerfacecolor',lc);
    end
end

%Set Axis Limits, Tick Marks, etc
ax = gca;
set(ax,'xlim',xlims,'ylim',ylims,'tickdir','out', ...
    'LineWidth',axis_lnwd,'box','on',...
    'fontname',sfont,'fontsize',axis_fsz,'fontweight',axis_style);

if lylog, set(ax,'yscale','log'); end

%Label the Plot and Create Legend
xlabel(xstr,'fontname',sfont,'fontsize',axis_lbl_fsz,'fontweight',axis_style)
ylabel(ystr,'fontname',sfont,'fontsize',axis_lbl_fsz,'fontweight',axis_style)
title(tstr,'fontname',sfont,'fontsize',title_fsz,'fontweight',title_style,...
    'interpreter','none','horizontalalignment','center')

if ~isempty(legstr); legend(h,legstr,'Location','NorthWest'); end


%Save the Figure
% Make sure there is a place for figures and Save as PDF and Tiff
if lsave
    if ~exist(plotDir,'dir')
        mkdir(plotDir);
    end
    saveas(gcf, [plotDir plotName '.pdf'],'pdf');
    print (gcf, '-dtiff', '-r600', [plotDir plotName]);
end



end