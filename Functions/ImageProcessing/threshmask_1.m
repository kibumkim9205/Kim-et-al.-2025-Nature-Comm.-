function mask=threshmask_1(varargin)

image = varargin{1};
blurradius = varargin{2};
blur=imfilter(image,fspecial('disk',blurradius),'symmetric'); %10x:3 20x:6
logblur=log(blur);
normlog=mat2gray(logblur);
if numel(varargin)>2
    minval=min(logblur(:));
    maxval=max(logblur(:));
    rangeval=maxval-minval;
    thresh=(log(varargin{3})-minval)/rangeval;
else
    thresh=graythresh(normlog);
end
mask=im2bw(normlog,thresh);
mask=imfill(mask,'holes');
end