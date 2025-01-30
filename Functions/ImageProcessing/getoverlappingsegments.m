function overlapsegs=getoverlappingsegments(mask,segs)
seglabel=bwlabel(segs);
overlaplabel=seglabel.*mask;
uniquesegs=unique(overlaplabel); uniquesegs(1)=[];
overlapsegs=ismember(seglabel,uniquesegs);
end