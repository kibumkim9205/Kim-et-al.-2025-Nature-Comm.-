function makeheatmaps_2(datasorted,POI,time,heatmapmax,heatmapmin,markoption,framesperhr)
numtraces=size(datasorted,1);
datasorted(datasorted>heatmapmax)=heatmapmax;
datasorted(datasorted<heatmapmin)=heatmapmin;
datasorted(isnan(datasorted))=heatmapmin; %make invalid data black
traceid=1:numtraces;
figure;
imagesc(time,traceid,datasorted);
cmap=colormap(parula);
%cmap(1,:)=[0 0 0];
colormap(cmap);
%%% overlay markers %%%%%%%%%%%%%%%%%%%%%%%%%
if markoption
    hold on;
    scatter((POI(~isnan(POI)))/framesperhr,traceid(~isnan(POI)),8,'r','*');
end