function avgtheta=revrise(series,step)
numsteps=length(step);
steplength=length(series)-max(step);
theta=zeros(numsteps,steplength);
for cc=1:numsteps
    i=step(cc);
    temp=atan2(series(1:end-i)-series(1+i:end),i);
    theta(cc,:)=temp(end-steplength+1:end);
end
avgtheta=mean(theta,1);
avgtheta=[ones(1,max(step))*avgtheta(1),avgtheta];