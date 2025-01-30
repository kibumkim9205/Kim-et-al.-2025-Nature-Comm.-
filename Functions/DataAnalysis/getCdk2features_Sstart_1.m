function [risetime,badtraces]=getCdk2features_2(sampletraces,samplestats)
%slopewindow=min([minlength 10]); %default 40
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
%minval=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
risetimedb=risetime;
%riseslope=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore1=ones(samplesize,tracelength)*NaN;
altvar1=ones(samplesize,1)*NaN;
prebuffer=10;
Sstartthresh=0.85;
for i=1:samplesize
    signal_total=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    %signal=signal(firstframe(i):lastframe(i));
    signal=smooth(signal_total(firstframe(i):lastframe(i)))'; %incoming signal always raw (unsmoothened)
    numframes=lastframe(i)-firstframe(i)+1;
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %badtraces(i)=min(signal(1:slopewindow))>1; %will remove permanently high signals
    %badtraces(i)=mean(signal(1:slopewindow))>1;
    if mean(signal(1:5))>1
        badtraces(i)=1;
        continue;
    end
    %%%%%% build filter
    gate=zeros(1,numframes);
    for j=prebuffer+2:numframes
        pastheight=max(signal(prebuffer:j-1))<Sstartthresh;
        presentheight=signal(j)>=Sstartthresh;
        %futureheight=max(signal(j:j+gatewindow))>signal(j)+0.05; %0.15
        %gate(j)=pastheight && presentheight && futureheight;
        gate(j)=pastheight && presentheight;
    end
    %%% find cdk2start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if max(gate)==0
        risetime(i)=NaN;
        continue;
    else
        tempsearch=find(gate,1,'first');
        risetimedb(i)=tempsearch;
        risetime(i)=tempsearch+firstframe(i)-1;
    end
end
%keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(1:length(frames),sigstore(i,frames));
    ylim([0.3 1.8]);
    hold on;
    %plot(1:length(frames),altstore(i,frames),'r');
    if badtraces(i)==1
        framemid=round(median(frames));
        plot(framemid-frames(1)+1,sigstore(i,framemid),'rx','markersize',20);
        continue;
    elseif isnan(risetime(i))
        %framemid=round(median(frames));
        %plot(framemid-frames(1)+1,sigstore(i,framemid),'bx','markersize',20);
        continue;
    end
    plot(risetimedb(i),sigstore(i,frames(risetimedb(i))),'go','markerfacecolor','g','markersize',6);
end
%}