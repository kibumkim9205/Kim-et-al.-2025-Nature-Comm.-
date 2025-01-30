function ringvals=sigringcalc(nuc_label,innerrad,outerrad,numgood,img,ringthresh,prctilethresh)
ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,img);
ring_info=regionprops(ring_label,'PixelIdxList');
ringvals=ones(numgood,1)*NaN;
for cc=1:numgood
    if cc>numel(ring_info)
        break;
    end
    ring2all=img(ring_info(cc).PixelIdxList);
    %ring2all(ring2all>prctile(ring2all,98))=[];
    ring2foreground=ring2all(ring2all>ringthresh);
    if numel(ring2foreground)<100
         ring2foreground=ring2all;
    end
    if numel(ring2all)>100
         %ringvals(cc)=nanmedian(ring2foreground);
         ringvals(cc)=prctile(ring2foreground,prctilethresh);
    end
end
end