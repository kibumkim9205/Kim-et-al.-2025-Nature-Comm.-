function ringvals=sigringcalc_cdk4(cellid,nuc_label,innerrad,outerrad,numgood,img,Minthresh,MAXthresh)
if isempty(cellid)
   cellid=1:numgood;     
end
ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,img);
%{
ring_mask=logical(ring_label);
ring_intensity=img(ring_mask);
figure; histogram(ring_intensity);

testIm=img; testIm(testIm>5000)=0;
extractmask_ring=bwmorph(ring_mask,'remove');
[height width]=size(extractmask_ring);
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(testIm));
tempframe(:,:,2)=extractmask_ring;
figure,imshow(tempframe);
%}

ring_info=regionprops(ring_label,'PixelIdxList');
ringvals=ones(numgood,1)*NaN;
for n=1:numgood
    cc=cellid(n);
    if cc>numel(ring_info)
        break;
    end
    ring2all=img(ring_info(cc).PixelIdxList);    
    ring2foreground=ring2all(ring2all>Minthresh);
    ring2all(ring2all>MAXthresh)=[];
    %if numel(ring2foreground)<100
    %     ring2foreground=ring2all;
    %end
    if numel(ring2all)>100
         ringvals(n)=nanmedian(ring2foreground);
    end
end
end