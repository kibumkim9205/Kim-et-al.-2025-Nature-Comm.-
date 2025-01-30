function bs=bgsubmasked_local(mask,nuc_info,raw,winrad,compression)
mask=imresize(mask,compression,'nearest');
raw=imresize(raw,compression,'nearest');

numcells=length(nuc_info);
[height,width]=size(mask);
raw(mask>0)=NaN;
bs=raw;
for i=1:numcells
    center=round(nuc_info(i).Centroid*compression);
    minrow=center(2)-winrad; minrow=max([minrow 1]);
    maxrow=center(2)+winrad; maxrow=min([maxrow height]);
    mincol=center(1)-winrad; mincol=max([mincol 1]);
    maxcol=center(1)+winrad; maxcol=min([maxcol width]);
    nucwindow=raw(minrow:maxrow,mincol:maxcol);
    background=nanmedian(nucwindow(:));
    bs(nuc_info(i).PixelIdxList)=background;
    %imshow(imadjust(mat2gray(YFP_raw(minrow:maxrow,mincol:maxcol))));
end
