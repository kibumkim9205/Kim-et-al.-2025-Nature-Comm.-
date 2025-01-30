function normsignals=normalizetraces_5(signals,tracestats,motherstats,refpoint)
numtraces=size(signals,1);
%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numtraces
    %realframes=find(~isnan(signals(i,:)));
    %signals(i,realframes)=smooth(signals(i,realframes));
    signals(i,:)=smoothignorenans(signals(i,:),3);
end
%%% get maxval to normalize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxval=NaN(numtraces,1);
switch refpoint
    case 'mother'
        for i=1:numtraces
            maxval(i)=max(signals(i,motherstats(i,1):motherstats(i,2)));
        end
    case 'daughter'
        for i=1:numtraces
            maxval(i)=max(signals(i,tracestats(i,1):tracestats(i,2)));
        end
end
%%% normalize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normsignals=ones(size(signals))*NaN;
for i=1:numtraces
    normsignals(i,:)=signals(i,:)/maxval(i);
end