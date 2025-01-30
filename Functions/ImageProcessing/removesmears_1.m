function correctedimage=removesmears_1(rawimage,foregroundthresh,areathresh,bgprctile)
smear_mask=rawimage>foregroundthresh;
smear_mask=bwareaopen(smear_mask,areathresh);
smear_mask=imdilate(smear_mask,strel('disk',100,0));
foreground=rawimage(~smear_mask);
background=prctile(foreground(:),bgprctile);
correctedimage=rawimage;
correctedimage(smear_mask)=background;
end