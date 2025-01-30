function boundaryset=getboundarycoors(mask)
maskboundary=bwboundaries(mask,'noholes');
boundaryset=maskboundary{1};
boundaryset=[boundaryset(end:-1:1,2) boundaryset(end:-1:1,1)]; %reverse order
end