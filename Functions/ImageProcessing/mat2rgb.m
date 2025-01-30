function RGB=mat2rgb(mat,cmap)
mat=mat*63; mat=mat+1;
[height,width]=size(mat);
template=zeros(height,width);
r=template; g=template; b=template;
r(:)=cmap(round(mat(:)),1);
g(:)=cmap(round(mat(:)),2);
b(:)=cmap(round(mat(:)),3);
RGB=template;
RGB(:,:,1)=r;
RGB(:,:,2)=g;
RGB(:,:,3)=b;
end