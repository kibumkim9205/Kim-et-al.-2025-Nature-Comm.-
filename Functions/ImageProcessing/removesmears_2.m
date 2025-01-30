function [correctedimage,stainflag]=removesmears_2(rawimage,foregroundthresh,areathresh,bgprctile)
smear_mask=rawimage>foregroundthresh;
smear_mask=bwareaopen(smear_mask,areathresh);
if max(smear_mask(:))>0
    stainflag=1;
else
    stainflag=0;
end
smear_mask=imdilate(smear_mask,strel('disk',100,0));
foreground=rawimage(~smear_mask);
background=prctile(foreground(:),bgprctile);
correctedimage=rawimage;
correctedimage(smear_mask)=background;
end