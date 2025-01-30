function [risetimeabs,badtraces]=getGemininfeatures_1(apctraces,samplestats,gemtraces)
[samplesize,tracelength]=size(apctraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
risetimeabs=ones(samplesize,1)*NaN;
risetimerel=risetimeabs;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore=sigstore;
prebuffer=6;
for i=1:samplesize
    signal=apctraces(i,:);
    gemsig=gemtraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    numframes=lastframe(i)-firstframe(i)+1;
    signal=signal(firstframe(i):lastframe(i));
    %signal=smooth(signal(firstframe(i):lastframe(i)))'; %incoming signal always raw (unsmoothened)
    sigstore(i,firstframe(i):lastframe(i))=signal;
    altstore(i,firstframe(i):lastframe(i))=gemsig(firstframe(i):lastframe(i));
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% general requirements
    postbuffer=10;
    gate=zeros(1,numframes);
    for j=1:numframes-postbuffer
        presentheight=signal(j)<0.5;
        maxfutureheight=max(signal(j+1:j+postbuffer))<signal(j);
        gate(j)=presentheight & maxfutureheight;
    end
    
    %sig_revslope=10*getslope_reverse(signal,1:10);
    %sig_fwdslope=10*getslope_forward_avg(signal,1:5);
    sig_fwdslope=10*getslope_forward_avg(signal,1:postbuffer);
    sig_time=(1:length(signal))/tracelength;
    %filterscore=2*sig_fwdslope-4*signal+sig_time+6;
    filterscore=signal;
    filterscore=filterscore.*gate;
    %altstore(i,firstframe(i):lastframe(i))=sig_fwdslope;
    %tempsearch=find(sig_fwdslope>0.05 & abs(signal)<0.03,1,'last');
    %tempsearch=find(filterscore>0.05,1,'first');
    filtermax=max(filterscore);
    tempsearch=find(filterscore==filtermax,1,'first');
    %if isempty(tempsearch) || signal(end)<0.05 %0.05
    if ~isempty(tempsearch) && filtermax>0
        risetimerel(i)=tempsearch;
        risetimeabs(i)=tempsearch+firstframe(i)-1; %return absolute POI rather than relative to mitosis
    %elseif signal(end)<0.1
    elseif max(signal(prebuffer:end))<0.1
        risetimeabs(i)=NaN;
    else
        badtraces(i)=1;
    end
end
keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
traceIDs=1:96;
ylims=[0 1.2];
POIpanel(traceIDs,firstframe,lastframe,sigstore,risetimerel,badtraces,ylims,altstore);
%}