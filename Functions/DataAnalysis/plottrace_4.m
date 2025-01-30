function plottraces_4(data,xvals,xstring,ystring,ylimits,smoothoption,drugspike,tracestats)
figure, hold on;
for i=1:size(data,1)
    plotdata=data(i,:);
    if smoothoption
        realframes=find(~isnan(data(i,:)));
        plotdata(realframes)=smooth(data(i,realframes));
    end
    line(xvals,plotdata,'color','r');
%     line(xvals+0.2,plotdata,'color','b','linestyle','--','linewidth',3);
%     hold on;
%     scatter(tracestats(i,1)/5,plotdata(tracestats(i,1)),100,'markerfacecolor','w','markeredgecolor','b','linewidth',4);
end

line([0 0],ylimits,'color','k','linestyle','--','linewidth',2);

ylim(ylimits);
xlabel(xstring); xlabel('Time (hrs)');
ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'D:\Downloads\FigAllLines.jpg');