function discretizexpercy(x,y,stepsize,minbin,maxbin,colorchar)
if stepsize==2
    modval=mod(minbin,2);
    if mod(maxbin,2)~=modval
        maxbin=maxbin-1;
    end
end
bins=minbin:stepsize:maxbin;
numbins=length(bins);
perc=ones(numbins,1)*NaN;
for i=1:numbins
    b=bins(i);
    binsample=y(x>=b-stepsize/2 & x<b+stepsize/2);
    if isempty(binsample)
        continue;
    end
    %[mu(i),std(i)]=normfit(binsample);
    %perc(i)=100*sum(binsample>thresh)/numel(binsample);
    perc(i)=100*sum(binsample)/numel(binsample);
end

if strcmp(colorchar,'r')
    colorcode=[1 0 0];
elseif strcmp(colorchar,'g')
    colorcode=[0 1 0];
elseif strcmp(colorchar,'b')
    colorcode=[0 0 1];
elseif strcmp(colorchar,'k')
    colorcode=[0 0 0];
end

line(bins',perc,'color',colorcode,'linewidth',2);
