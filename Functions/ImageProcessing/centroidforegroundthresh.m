function newmask=centroidforegroundthresh(nuc_mask,centers,raw,nucr,localoffsetratio,foregroundoffsetratio,resizeratio)
nuc_mask=imresize(nuc_mask,resizeratio,'nearest');
raw=imresize(raw,resizeratio,'nearest');

localoffset=round(resizeratio*localoffsetratio*nucr);
foregroundoffset=round(resizeratio*foregroundoffsetratio*nucr);
%regionlength=2*localoffset+1;
[height,width]=size(nuc_mask);
bgraw=raw;
bgraw(nuc_mask>0)=NaN;

centers=centers*resizeratio;
numcells=size(centers,1);
x=centers(:,1); y=centers(:,2);

bgprctile=50;
[~,localminrow,localmaxrow,localmincol,localmaxcol,pixelindices]=calclocal(bgraw,x,y,localoffset,bgprctile);
% %%% calculate local boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% onesmatrix=ones(numcells,1);
% localminrow=round(y-localoffset); localminrow=max(localminrow,onesmatrix);
% localmaxrow=round(y+localoffset); localmaxrow=min(localmaxrow,onesmatrix*height);
% localmincol=round(x-localoffset); localmincol=max(localmincol,onesmatrix);
% localmaxcol=round(x+localoffset); localmaxcol=min(localmaxcol,onesmatrix*width);
% %%% calculate local background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% localheightmatrix=localmaxrow-localminrow;
% localwidthmatrix=localmaxcol-localmincol;
% templatecols=ones(regionlength,1)*[0:height:(regionlength-1)*height];
% templaterows=[1:regionlength]'*ones(1,regionlength);
% templateindices=templatecols+templaterows;
% bgraw(1)=NaN;
% linearindices=ones(regionlength*regionlength,1);
% linearindices=linearindices*ones(1,numcells);
% firstindex=sub2ind([height width],localminrow,localmincol);
% for i=1:numcells
%     indices=templateindices(1:localheightmatrix(i),1:localwidthmatrix(i))+firstindex(i)-1;
%     linearindices(1:numel(indices),i)=indices(:);
% end
% localbgraw=bgraw(linearindices);
% localbgcalc=prctile(localbgraw,prcthresh)'; %ignores NaNs
%%% calculate foreground sampling boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fgprctile=5;
[foregroundvals,~,~,~,~]=calclocal(raw,x,y,foregroundoffset,fgprctile);
% foregroundminrow=round(y-foregroundoffset); foregroundminrow=max(foregroundminrow,onesmatrix);
% foregroundmaxrow=round(y+foregroundoffset); foregroundmaxrow=min(foregroundmaxrow,onesmatrix*height);
% foregroundmincol=round(x-foregroundoffset); foregroundmincol=max(foregroundmincol,onesmatrix);
% foregroundmaxcol=round(x+foregroundoffset); foregroundmaxcol=min(foregroundmaxcol,onesmatrix*width);
%%% for each region, calculate mask pixels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
onesmatrix=ones(numcells,1);
relx=round(x-localmincol+1); relx=max(relx,onesmatrix);
rely=round(y-localminrow+1); rely=max(rely,onesmatrix);
maskindices=[];
for i=1:numcells
    %relx=round(x(i)-localmincol(i)+1);
    %rely=round(y(i)-localminrow(i)+1);
    localraw=raw(localminrow(i):localmaxrow(i),localmincol(i):localmaxcol(i));
    localmask=localraw>foregroundvals(i);
    localmask=imfill(localmask,'holes');
    centermask=sub2img([localmaxrow(i)-localminrow(i)+1 localmaxcol(i)-localmincol(i)+1],[relx(i) rely(i)]);
    localmask=getoverlappingsegments(centermask,localmask);
    fgindices=find(localmask==1);
    localindices=pixelindices(fgindices,i);
    maskindices=[maskindices;localindices];
end
newmask=zeros(height,width);
newmask(maskindices)=1;
end

function [localval,minrow,maxrow,mincol,maxcol,pixelindices]=calclocal(raw,x,y,offset,prcthresh)
[height,width]=size(raw);
numcells=numel(x);
regionlength=2*offset+1;
%%% calculate local boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
onesmatrix=ones(numcells,1);
minrow=round(y-offset); minrow=max(minrow,onesmatrix);
maxrow=round(y+offset); maxrow=min(maxrow,onesmatrix*height);
mincol=round(x-offset); mincol=max(mincol,onesmatrix);
maxcol=round(x+offset); maxcol=min(maxcol,onesmatrix*width);
%%% calculate local background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% heightmatrix=maxrow-minrow;
% widthmatrix=maxcol-mincol;
heightmatrix=maxrow-minrow+1;
widthmatrix=maxcol-mincol+1;
templatecols=ones(regionlength,1)*[0:height:(regionlength-1)*height];
templaterows=[1:regionlength]'*ones(1,regionlength);
templateindices=templatecols+templaterows;
raw(1)=NaN;
pixelindices=ones(regionlength*regionlength,1);
pixelindices=pixelindices*ones(1,numcells);
firstindex=sub2ind([height width],minrow,mincol);
for i=1:numcells
    indices=templateindices(1:heightmatrix(i),1:widthmatrix(i))+firstindex(i)-1;
    pixelindices(1:numel(indices),i)=indices(:);
end
localraw=raw(pixelindices);
localval=prctile(localraw,prcthresh)'; %ignores NaNs
end



%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
[th,tw]=size(localraw);
rgb=zeros(th,tw,3);
rgb(:,:,1)=imadjust(mat2gray(localraw));
% a=zeros(height,width);
% a(pixelindices(:,51))=1;
% rgb(:,:,2)=bwperim(a);
rgb(:,:,2)=bwperim(localmask);
imshow(rgb);
%}
