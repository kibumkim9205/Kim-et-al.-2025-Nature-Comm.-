function img=sub2img(dims,coorvec)
boundarymarkersidx=sub2ind(dims,coorvec(:,2),coorvec(:,1));
img=zeros(dims);
img(boundarymarkersidx)=1;
end