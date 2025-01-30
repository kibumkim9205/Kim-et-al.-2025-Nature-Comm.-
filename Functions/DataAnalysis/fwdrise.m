function avgrise=fwdrise(series,step)
numsteps=length(step);
steplength=length(series)-max(step);
rise=zeros(numsteps,steplength);
for cc=1:numsteps
    i=step(cc);
    temprise=series(1+i:end)-series(1:end-i);
    rise(cc,:)=temprise(1:steplength);
end
avgrise=mean(rise,1);
avgrise=[avgrise,ones(1,max(step))*avgrise(end)];