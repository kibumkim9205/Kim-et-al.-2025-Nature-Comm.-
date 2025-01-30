function ringvals=sigringcalc(cellid,nuc_label,innerrad,outerrad,numgood,img,ringthresh)
if isempty(cellid)
   cellid=1:numgood;     
end
ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,img);
ring_info=regionprops(ring_label,'PixelIdxList');
ringvals=ones(numgood,1)*NaN;
for n=1:numgood
    cc=cellid(n);
    if cc>numel(ring_info)
        break;
    end
    ring2all=img(ring_info(cc).PixelIdxList);
    %ring2all(ring2all>prctile(ring2all,98))=[];
    ring2foreground=ring2all(ring2all>ringthresh);
    %if numel(ring2foreground)<100
    %     ring2foreground=ring2all;
    %end
    if numel(ring2all)>100
         ringvals(n)=nanmedian(ring2foreground);
    end
end
end