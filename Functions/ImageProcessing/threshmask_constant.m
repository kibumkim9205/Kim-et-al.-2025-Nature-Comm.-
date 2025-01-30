function mask=threshmask(image,blurradius,thresh)
if blurradius>0
    blur=imfilter(image,fspecial('disk',blurradius),'symmetric'); %10x:3 20x:6
else
    blur=image;
end
%normlog=mat2gray(log(blur));
%thresh=graythresh(normlog);
mask=imbinarize(blur,thresh);
mask=imfill(mask,'holes');
end