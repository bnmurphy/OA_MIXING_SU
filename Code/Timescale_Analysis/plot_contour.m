function plot_contour(x, y, z, zlines, lzlog, xstr, ystr, tstr, tinterpreter,...
    xlims, ylims, xticks, yticks, xticklbl, yticklbl,...
    plotDir, plotName, legstr, lxlog, lylog, lsave)
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
             'position',[0.2 0.3 0.7 0.6],...
             'visible','on',...
             'paperorientation','portrait',...
             'papertype','a4',...
             'paperunit','centimeters',...
             'paperposition',half_page);

%Draw Plot
[cs, h] = contour(x, y, z, zlines,'linewidth',3);


%Set Axis Limits, Tick Marks, etc
ax = gca;
set(ax,'xlim',xlims,'ylim',ylims,'xtick',xticks,'ytick',yticks,...
    'tickdir','out', ...
    'LineWidth',axis_lnwd,'box','on',...
    'fontname',sfont,'fontsize',axis_fsz,'fontweight',axis_style,...
    'position',[0.13 0.16, 0.775, 0.70]);

if lxlog, set(ax,'xscale','log'); end
if lylog, set(ax,'yscale','log'); end

%Format Contour Labels
hlbl = clabel(cs,h,'labelspacing',144,'fontsize',axis_fsz,'fontname',sfont,...
    'fontweight','bold','fontangle','oblique','edgecolor','k','rotation',0,...
    'backgroundcolor',[1 1 1]);
if lzlog  %Make contour labels exponential
    for ilbl = 1:length(hlbl)
        temp_str = get(hlbl(ilbl),'string');
        set(hlbl(ilbl),'string', ['10^{' temp_str '}']);
    end
end

%Format Tick Labels
[hx, hy] = format_ticks(gca, xticklbl, yticklbl, xticks, yticks,[],[],[]);

%Format Axes Labels
hxlbl = xlabel(xstr,'fontname',sfont,'fontsize',axis_lbl_fsz,'fontweight',axis_style,...
    'units','normalized');
pos = get(hxlbl,'position');
set(hxlbl,'position',[pos(1) -0.09 0]);

hylbl = ylabel(ystr,'fontname',sfont,'fontsize',axis_lbl_fsz,'fontweight',axis_style,...
    'units','normalized');
pos = get(hylbl,'position');
set(hylbl,'position',[-0.09 pos(2) 0]);

%Format Title
title(tstr,'fontname',sfont,'fontsize',title_fsz,'fontweight',title_style,...
    'interpreter',tinterpreter,'horizontalalignment','center','units','normalized',...
    'position',[0.5 1.04 0.0]);


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