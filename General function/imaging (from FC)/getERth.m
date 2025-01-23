function [th0_1a,th1_1a,th0_2a,th1_2a]=getERth(tempMIX0,sz1b,sz1s)

%%
tempMIX1=bgsub(tempMIX0,sz1b,0.05);
tempMIX1a=imadjust(mat2gray(tempMIX1));
tempMIX2=bgsub(tempMIX0,sz1s,0.05);
tempMIX2a=imadjust(mat2gray(tempMIX2));

%%
tempseries_2a=tempMIX2a((tempMIX2a>0)&(tempMIX2a<1));
th0_2a=graythresh(tempseries_2a);
[~,th1_2a]=ThreshImage(tempseries_2a);
th0_2a=1/mean([1/th1_2a,1/th0_2a]);

%%
tempseries_1a=tempMIX1a((tempMIX1a>0)&(tempMIX1a<1));
th0_1a=graythresh(tempseries_1a);
[~,th1_1a]=ThreshImage(tempseries_1a);
th0_1a=1/mean([1/th1_1a,1/th0_1a]);