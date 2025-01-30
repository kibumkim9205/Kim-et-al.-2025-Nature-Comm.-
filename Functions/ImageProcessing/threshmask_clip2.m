function mask=threshmask_clip2(image,blurradius,threshmin)
if blurradius>0
    blur=imfilter(image,fspecial('disk',blurradius),'symmetric'); %10x:3 20x:6
else
    blur=image;
end
blur(blur<1)=1;
normlog=mat2gray(log(blur));
% normlog(normlog<threshmin)=0;
% normlog2=mat2gray(blur);
thresh=graythresh(normlog);
mask=im2bw(normlog,thresh);
mask=imfill(mask,'holes');
end