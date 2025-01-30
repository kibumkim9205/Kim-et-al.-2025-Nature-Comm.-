function bridgeindices = bridgevec(pos1,pos2)
% format for pos and bridgeindices is [x,y]
diffx=pos2(:,1)-pos1(:,1);
diffy=pos2(:,2)-pos1(:,2);
longerside=max([abs(diffx) abs(diffy)],[],2);
stepx=diffx./longerside;
stepy=diffy./longerside;
% bridgex=zeros(longerside+1,1);bridgey=zeros(longerside+1,1);
% for bs=0:longerside
%     bridgex(1+bs)=pos1(1)+round(bs*stepx);
%     bridgey(1+bs)=pos1(2)+round(bs*stepy);
% end
numpos=size(pos1,1);
bridgeindices=[];
for i=1:numpos
    if longerside(i)>0
        tempbridgex=pos1(i,1)+ round((0:longerside(i))*stepx(i));
        tempbridgey=pos1(i,2)+ round((0:longerside(i))*stepy(i));
        bridgeindices=[bridgeindices;[tempbridgex' tempbridgey']];
    end
end

%%% debug: stop in for-loop to confirm proper bridging %%%%%%%%%%%%%%%%%%%%
%{
posimg=zeros(850);
posimg(pos1(i,2),pos1(i,1))=1;
posimg(pos2(i,2),pos2(i,1))=1;
figure;imshow(posimg);
bridgeimg=zeros(850);
bridgeindices2=sub2ind([850 850],bridgeindices(:,2),bridgeindices(:,1));
bridgeimg(bridgeindices2)=1;
figure;imshow(bridgeimg);
combinedimg=posimg | bridgeimg;
figure;imshow(combinedimg);
%}