function [intersectcoors,intersectidx]=getanglecoors(centroid,boundaryset,angles)
numboundarypix=size(boundaryset,1);
centroidx=ones(numboundarypix,1)*centroid(1);
centroidy=ones(numboundarypix,1)*centroid(2);
centroidset=[centroidx centroidy];
centroiddiff=boundaryset-centroidset;
centroidangles=atan2(centroiddiff(:,2),centroiddiff(:,1));
intersectidx=ones(numel(angles),1)*NaN;
for ac=1:numel(angles)
    angleset=ones(numboundarypix,1)*angles(ac);
    %[~,intersectidx(ac)]=min(abs(centroidangles-angleset));
    [~,intersectidx(ac)]=min(mod(centroidangles-angleset,2*pi));
end
intersectcoors=boundaryset(intersectidx,:);
end