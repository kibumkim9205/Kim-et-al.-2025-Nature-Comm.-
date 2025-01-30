function ksyn=calcksyn(traces,tracestats,Sdelay)
maxframe=size(traces,2);
numtraces=size(traces,1);
tracestep=diff(traces,1,2);
Sksyn=nan(numtraces,1);
for i=1:numtraces
    if tracestats(i,1)+Sdelay<maxframe
        tracestepsmooth=conv(tracestep(i,:),normpdf(-3:3,0,2),'same');
        Sksyn(i)=tracestepsmooth(tracestats(i,1)+Sdelay);
    end
end

keyboard;
histogram(Sksyn,30);
end