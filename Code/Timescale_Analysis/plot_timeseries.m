function psize = plot_timeseries(y, t, ystr, tstr, ylimits, legstr, titlestr, newfig)

if newfig == 1, figure, end
set(gcf,'Position',[250 200 1000 400])
set(gca,'Position',[0.1 0.2 0.85 0.65])

y3 = y(1,:) + 0.1;
y4 = y(1,:) - 0.1;
a1 = area([1:24],[y4',y3'-y4']);
set(a1(2),'facecolor',[160,160,160]./255)
set(a1(1),'facecolor',[1,1,1],'edgecolor','none')
set(a1,'edgecolor','none')
hold on
set(gca,'layer','top')

h = plot(t,y, 'LineWidth',4);

set(h(1), 'color','k','linestyle','none','marker','o','markersize',4,...
    'markerfacecolor','k');
if size(h,1) > 1, set(h(2), 'color',[1 0 0]); end
if size(h,1) > 2, set(h(3), 'color',[0 0 1]); end
if size(h,1) > 3, set(h(4), 'color',[0 145/255 0]); end
if size(h,1) > 4, set(h(5), 'color',[153/255 0 153/255]); end

ylim(ylimits);
xlim([t(1),t(end)]);
sfont = 'Arial';

title(titlestr, 'FontName',sfont,'FontSize',20,'FontWeight','bold')
ylabel(ystr, 'FontName',sfont,'FontSize',16,'FontWeight','bold')
xlabel(tstr, 'FontName',sfont,'FontSize',16,'FontWeight','bold')
legend(h,legstr,'location','Southeast')
set(gca, 'FontName',sfont,'FontSize',12,'linewidth',2)
box on

psize = [11 4.5];