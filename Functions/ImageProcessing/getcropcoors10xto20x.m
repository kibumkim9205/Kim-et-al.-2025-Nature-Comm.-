function [cropcoors]=getcropcoors10xto20x(dims10x,dims20x,jitx,jity,subsite)
%%% Use: only 1 10x image followed by 1 20x image %%%%%%%%%%%%%%%%%%%%%%%%%
x20x=dims20x(1); y20x=dims20x(2);
x10x=dims10x(1); y10x=dims10x(2);
%%% Establish limits and jitter relative to the correct quadrant %%%%%%%%%%
switch subsite
    case 1
        minx=1; maxx=x20x;
        miny=1; maxy=y20x;
    case 2
        minx=x20x+1; maxx=x10x;
        miny=1; maxy=y20x;
        jitx=jitx-x20x;
    case 3
        minx=1; maxx=x20x;
        miny=y20x+1; maxy=y10x;
        jity=jity-y20x;
    case 4
        minx=x20x+1; maxx=x10x;
        miny=y20x+1; maxy=y10x;
        jitx=jitx-x20x;
        jity=jity-y20x;
end
numshifts=numel(jitx);
cropcoors=ones(numshifts+1,4)*NaN;
%%% get edge cropping coordinates of 10x image %%%%%%%%%%%%%%%%%%%%%%%%%%%%
topidx10x=miny+jity;
topidx10x=max([1;topidx10x]);
maxtopshift=topidx10x-miny;
bottomidx10x=maxy+jity;
bottomidx10x=min([y10x;bottomidx10x]);
maxbottomshift=bottomidx10x-maxy;
leftidx10x=minx+jitx;
leftidx10x=max([1;leftidx10x]);
maxleftshift=leftidx10x-minx;
rightidx10x=maxx+jitx;
rightidx10x=min([x10x;rightidx10x]);
maxrightshift=rightidx10x-maxx;
cropcoors(1,:)=[miny+maxtopshift maxy+maxbottomshift minx+maxleftshift maxx+maxrightshift];

%%% get edge cropping coordinates of 20x image %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cropcoors(2,:)=[1 y20x 1 x20x];
if topidx10x==1, cropcoors(2,1)=1-jity; end
if bottomidx10x==y10x, cropcoors(2,2)=y20x-jity; end
if leftidx10x==1, cropcoors(2,3)=1-jitx; end
if rightidx10x==x10x, cropcoors(2,4)=x20x-jitx; end