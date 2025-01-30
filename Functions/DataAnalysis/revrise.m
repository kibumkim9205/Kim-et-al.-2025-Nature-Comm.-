function avgrise=revrise(series,step)
numsteps=length(step);
steplength=length(series)-max(step);
rise=zeros(numsteps,steplength);
for cc=1:numsteps
    i=step(cc);
    temprise=series(1:end-i)-series(1+i:end);
    rise(cc,:)=temprise(end-steplength+1:end);
end
avgrise=mean(rise,1);
avgrise=[ones(1,max(step))*avgrise(1),avgrise];