function array=shiftarray(array,padsize)
reljitx=padsize(2);
if reljitx>0
    array=padarray(array,[0 abs(reljitx)],'replicate','pre');
    array=array(:,1:end-abs(reljitx));
else
    array=padarray(array,[0 abs(reljitx)],'replicate','post');
    array=array(:,abs(reljitx)+1:end);
end
reljity=padsize(1);
if reljity>0
    array=padarray(array,[abs(reljity) 0],'replicate','pre');
    array=array(1:end-abs(reljity),:);
else
    array=padarray(array,[abs(reljity) 0],'replicate','post');
    array=array(abs(reljity)+1:end,:);
end
end