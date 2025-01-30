function RGB=seteachchannel(RGB,mask,setvalue)
tempmat=RGB(:,:,1); tempmat(mask)=setvalue; RGB(:,:,1)=tempmat;
tempmat=RGB(:,:,2); tempmat(mask)=setvalue; RGB(:,:,2)=tempmat;
tempmat=RGB(:,:,3); tempmat(mask)=setvalue; RGB(:,:,3)=tempmat;
end