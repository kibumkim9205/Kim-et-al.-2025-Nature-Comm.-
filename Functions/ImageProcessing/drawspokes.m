function mask=drawspokes(mask,coorset1,coorset2)
[height,width]=size(mask);
numcoors=size(coorset1,1);
%%%%%% generate lines (spokes) connecting coorset1 and coorset2 %%%%%%%%%%%
spokes=zeros(height,width);
for pc=1:numcoors
    [bx,by]=bridge(coorset1(pc,:),coorset2(pc,:));
    for boi=1:length(bx)
        spokes(by(boi),bx(boi))=1;
    end
end
spokes=logical(spokes);
%%%%%% define bins by splitting with spokes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thickspokes=imdilate(spokes,strel('disk',1));
mask(thickspokes)=0;