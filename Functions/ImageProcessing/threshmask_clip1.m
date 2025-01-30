function mask=threshmask_clip1(image,blurradius)
if blurradius>0
    blur=imfilter(image,fspecial('disk',blurradius),'symmetric'); %10x:3 20x:6
else
    blur=image;
end
blur(blur<1)=1;
normlog=mat2gray(log(blur));
%normlog=mat2gray(blur);
[thresh thresh2]=graythresh(normlog);
mask=im2bw(normlog,thresh);
mask=imfill(mask,'holes');
end