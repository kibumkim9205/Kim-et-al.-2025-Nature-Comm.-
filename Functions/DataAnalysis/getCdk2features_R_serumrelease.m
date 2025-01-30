function [risetime,badtraces]=getCdk2features_R_serumrelease(sampletraces,samplestats)
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
for i=1:samplesize
    signal_total=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    signal=smooth(signal_total(firstframe(i):lastframe(i)))'; %incoming signal always raw (unsmoothened)
    numframes=lastframe(i)-firstframe(i)+1;
    sigstore(i,firstframe(i):lastframe(i))=signal;

    sig_fwdslope=getslope_forward_avg(signal,4:6);
    sig_revslope=getslope_reverse_avg(signal,4:6);
    %filter=[diff(signal)>0.5 0];
    %filter=sig_fwdslope>0.005 & sig_revslope>-0.002;
    %filter=signal<0.45 & sig_fwdslope>0.005;
    sig_time=(1:length(signal))/tracelength;
    filter=2*sig_fwdslope-4*signal+sig_time+6;
    filter(1:25)=0;
    
%     if max(filter)==1
%         risetime(i)=find(filter,1,'first');
%     else
%         risetime(i)=NaN;
%         continue;
%     end
    filtermax=max(filter);
    risetime(i)=find(filter==filtermax,1,'first');
    if filtermax<=0 || risetime(i)>numframes-5
        risetime(i)=NaN;
        continue;
    end
    
    
    %altvar1(i)=2*sig_fwdslope(risetime(i))-4*signal(risetime(i));
    risetimedb(i)=risetime(i);
    risetime(i)=risetime(i)+firstframe(i); %return absolute POI rather than relative to mitosis
    %altstore1(i,firstframe(i):lastframe(i))=filter;
    
%     if risetime(i)-firstframe(i)<length(signal)-20;
%         if signal_total(risetime(i)+20)-signal_total(risetime(i))<0.14;
%             risetime(i)=NaN;
%         end
%     end

end
%keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%altstore1(altstore1==0)=NaN;
%altstore_norm=normalizeMyTracesGeminin_alt2(altstore1,0.01);
for i=1:96
%for i=1:samplesize
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(1:length(frames),sigstore(i,frames));
    axis([1 length(frames) 0.2 1]);

    if badtraces(i)==1 || isnan(risetimedb(i))
        continue;
    end
    hold on;
    plot(risetimedb(i),sigstore(i,frames(risetimedb(i))),'go','markerfacecolor','g','markersize',6);
    %plot(1:length(frames),altstore1(i,frames)-1,'r');
    %plot(risetimedb(i),altstore_norm(i,frames(risetimedb(i)))+0.5,'go','markerfacecolor','g','markersize',6);
end
%}