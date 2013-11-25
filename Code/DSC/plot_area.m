function plot_area(x, y, xstr, ystr, tstr, xlims, ylims, ylog,...
    plotDir, plotName, legstr, lsave)
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
quarter_page = [2.5 1 9 8.5];  %Good for a square
tall_page = [2.5 1 7 13.5];
long_page = [1 1 11 5];
figure('units','normalized',...
             'position',[0.2 0.3 0.7 0.3],...
             'visible','on',...
             'paperorientation','portrait',...
             'papertype','a4',...
             'paperunit','centimeters',...
             'paperposition',long_page);

%Draw Plot
line_color = {'k',[1 0 0],[0 0 1],[0 145/255 0],[153 0 153]./255};

h = area(x, y);


%Set Axis Limits, Tick Marks, etc
ax = gca;
set(ax,'xlim',xlims,'ylim',ylims,'tickdir','out', ...
    'LineWidth',axis_lnwd,'box','on',...
    'fontname',sfont,'fontsize',6,'fontweight',axis_style);

if ylog, set(ax,'yscale','log'); end

%Label the Plot and Create Legend
xlabel(xstr,'fontname',sfont,'fontsize',axis_lbl_fsz,'fontweight',axis_style)
ylabel(ystr,'fontname',sfont,'fontsize',axis_lbl_fsz,'fontweight',axis_style)
title(tstr,'fontname',sfont,'fontsize',title_fsz,'fontweight',title_style,...
    'interpreter','none','horizontalalignment','center')

if ~isempty(legstr); legend(h,legstr,'Location','NorthWest','fontsize',6); end


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