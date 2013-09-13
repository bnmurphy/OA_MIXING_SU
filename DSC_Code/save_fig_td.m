%Save_Fig gets the dimensions of the current figure and saves it
%
%BNM - 10/1/10



% Make PDF Output Pretty
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'PaperPosition',[0.1 0.1 pos(3) pos(4)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', psize);
 

% Make sure there is a place for figures and save as a pdf
if ~exist(['../Figs/', scen2],'dir'); mkdir(['../Figs/', scen2]); end;

fname = ['../Figs/', scen2 '/' plotname];
print (gcf, '-dtiff', '-r300', fname);
% print (gcf, '-dpdf', '-r600', fname);
saveas(gcf,[fname '.pdf'],'pdf');
