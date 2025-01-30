function [FRET_ratio_image] = FRETimage(f,sig1,sig2,jitters,height,width,numlivecells,cellid,nuc_info,nuc_mask,thresh)
FRET_ratio=sig2./sig1;
%hist(FRET_ratio,100); %FRET range
if isempty(thresh)
    thresh=[0.8 2.7];
end

FRET_ratio_image=zeros(height,width);
for n=1:numlivecells
    cc=cellid(n);
    if cc>numel(nuc_info)
        break;
    end
    if FRET_ratio(n)<=thresh(1)
        FRET_ratio_image(nuc_info(cc).PixelIdxList)=thresh(1)+0.01;
    elseif FRET_ratio(n)>=thresh(2)
        FRET_ratio_image(nuc_info(cc).PixelIdxList)=thresh(2);
    else
        FRET_ratio_image(nuc_info(cc).PixelIdxList)=FRET_ratio(n);
    end
end
%%% generate FRET ratio image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FRET_ratio_image(~nuc_mask)=nan;
FRET_ratio_image(isinf(FRET_ratio_image))=nan;
%FRET_ratio_image(FRET_ratio_image<0.5)=0.51;

%%% jitter correct FRET ratio image
tempjity=round(abs(jitters(f,2)));
pady=zeros(tempjity,width);
if jitters(f,2)>0
    FRET_ratio_image=[pady;FRET_ratio_image]; %add top
    FRET_ratio_image(end-tempjity+1:end,:)=[]; %remove bottom
else
    FRET_ratio_image(1:tempjity,:)=[]; %remove top
    FRET_ratio_image=[FRET_ratio_image;pady]; %add bottom
end
tempjitx=round(abs(jitters(f,1)));
padx=zeros(height,tempjitx);
if jitters(f,1)>0
    FRET_ratio_image=[padx FRET_ratio_image]; %add left
    FRET_ratio_image(:,end-tempjitx+1:end)=[]; %remove right
else
    FRET_ratio_image(:,1:tempjitx)=[]; %remove left
    FRET_ratio_image=[FRET_ratio_image padx]; %add right
end

%imwrite(uint16(FRET_ratio),[maskdir,'\','FRETratio_',num2str(f),'.tif']);
end