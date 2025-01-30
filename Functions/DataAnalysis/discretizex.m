function [mu,std]=discretizex(graphoption,stdorsem,x,y,stepsize,minbin,maxbin,colorchar)
if stepsize==2
    modval=mod(minbin,2);
    if mod(maxbin,2)~=modval
        maxbin=maxbin-1;
    end
end
bins=minbin:stepsize:maxbin;
numbins=length(bins);
mu=ones(numbins,1)*NaN; std=ones(numbins,1)*NaN;
for i=1:numbins
    b=bins(i);
    binsample=y(x>=b-stepsize/2 & x<b+stepsize/2);
    if isempty(binsample)
        continue;
    end
    %[mu(i),std(i)]=normfit(binsample);
    mu(i)=nanmean(binsample);
%     fprintf('%0.2f\n',mu(i));
    switch stdorsem
        case 'std'
            std(i)=nanstd(binsample);
        case 'sem'
            std(i)=nanstd(binsample)/sqrt(numel(binsample));
    end
end

if strcmp(colorchar,'r')
    colorcode=[1 0 0];
elseif strcmp(colorchar,'g')
    colorcode=[0 1 0];
elseif strcmp(colorchar,'b')
    colorcode=[0 0 1];
elseif strcmp(colorchar,'k')
    colorcode=[0 0 0];
elseif strcmp(colorchar,'c')
    colorcode=[0 1 1];
end

switch graphoption
    case 'errorbars'
        line(bins',mu,'color',colorcode,'linewidth',1);
        hold on;
        h=errorbar(bins,mu,std,'.','linewidth',1,'color',colorcode);
        %errorbar(bins,mu,std,'o','linewidth',2,'color','k','markerfacecolor','w');
        h.CapSize=10;
    case 'shading'
        binfill=[bins fliplr(bins)];
        fadedcolor=colorcode+.7*(1-colorcode);
        line(bins',mu,'color',colorcode,'linewidth',1);
        hold on;
        auto='off';
        switch auto
            case 'on'
                fill(binfill',[mu'+std' fliplr(mu'-std')],fadedcolor,'edgecolor',fadedcolor,'FaceAlpha', 0.4);
            case 'off' %make in illustrator
                line(bins',mu+std);
                line(bins',mu-std);
                line([bins(1) bins(1)],[mu(1)-std(1) mu(1)+std(1)]);
                line([bins(end) bins(end)],[mu(end)-std(end) mu(end)+std(end)]);
        end
end