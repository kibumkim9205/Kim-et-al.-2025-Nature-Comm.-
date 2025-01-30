function [bs,bg]=bgsubring(nuc_info,ring_info,raw,bgprc)

numcells=length(nuc_info);
bg=raw;
for i=1:numcells
    %background=nanmedian(raw(ring_info(i).PixelIdxList));
    background=prctile(raw(ring_info(i).PixelIdxList),bgprc);
    bg(nuc_info(i).PixelIdxList)=background;
    %imshow(imadjust(mat2gray(YFP_raw(minrow:maxrow,mincol:maxcol))));
end
bs=raw-bg;