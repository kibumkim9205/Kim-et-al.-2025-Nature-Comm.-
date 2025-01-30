function bg=ThreshImage(ImageOld)

TempSeries=ImageOld(:);
% nbin=numel(tbin);
% tmax=max(tbin);
% tmin=min(tbin);
% tstep=(tmax-tmin)/(nbin-1);

%%% get background value: search intensity distribution for most convex point %%%%%%%%%%%%%%%%%%%%%%%
tmax=max(TempSeries);
tmin=min(TempSeries);
nbin=200;
tstep=(tmax-tmin)/nbin;
tmin=tmin+tstep/2;tmax=tmax-tstep/2;
[n,xout]=ksdensity(TempSeries,tmin:tstep:tmax);
gp=max([2,ceil(nbin/50)]);  %gp=4
ng=getcurvature(n,gp);      %returns the angle change of intensity histogram over one bin at each point

negpeak=regionprops(bwlabel(ng<0),'PixelIdxList','Area'); %returns indices and area (# of indices) for each contiguous set of bins with negative (convex) curvature
negth=prctile(ng,25)-1.5*iqr(ng);                       %threshold for low outliers
Ibg00=zeros(1,size(negpeak,1));AC=Ibg00;DD=Ibg00;
for cc=1:size(negpeak,1)
    AC(cc)=negpeak(cc).Area;
    tempeakvalues=ng(negpeak(cc).PixelIdxList);         %curvature for each index of negpeak group
    [DD(cc),I]=min(tempeakvalues);                      %most negative curvature
    Ibg00(cc)=negpeak(cc).PixelIdxList(I);              %index for this most negative curve for this negpeak group
end
neglo=((DD<=negth)|(AC>=gp));       %each peak legitimate if: curvature is low outlier or area > 4 bins
Ibg00=Ibg00(neglo);DD=DD(neglo);    %remove illegitimate peaks
if isempty(Ibg00)
    [DD,Ibg00]=min(ng);             %if set of peaks is empty, set to index with lowest global curvature
end
Ibg00=[Ibg00,size(ng,2)-gp];        %pad row vector with last intensity bin (for threshold calc contingency)
[dummy,I_DD]=min(DD);               %I_DD = negpeak group # with the lowest global curvature
bg0=xout(Ibg00(I_DD));              %intensity bin of peak with lowest global curvature
bg=bg0;

% TempSeries2=TempSeries((TempSeries>(bg0-5*gp*tstep))&(TempSeries<(bg0+5*gp*tstep)));  %narrow search for peak
% [n_bg,xout_bg]=ksdensity(TempSeries2);
% [dummy,Ibg]=max(n_bg);              %index of peak
% xfit=xout_bg(Ibg-1:Ibg+1)';yfit=n_bg(Ibg-1:Ibg+1)';
% det_n=det([xfit.^2,xfit,[1;1;1]]);
% a=det([yfit,xfit,[1;1;1]])/det_n;
% b=det([xfit.^2,yfit,[1;1;1]])/det_n;
% bg=-b/a/2;                          %refined index of peak --> background intensity
