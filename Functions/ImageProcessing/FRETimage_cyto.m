function [FRET_ratio_image] = FRETimage_cyto(f,imRatio,jitters)
%hist(FRET_ratio,100); %FRET range
[height width]=size(imRatio);
FRET_ratio_image=imRatio;
%%% generate FRET ratio image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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