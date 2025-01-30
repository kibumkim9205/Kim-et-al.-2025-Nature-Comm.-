function normhist(vals,bin)
bstep=bin(2)-bin(1);
binedges=bin(1)-0.5*bstep:bstep:bin(end)+0.5*bstep;
pdfvals=histcounts(vals,binedges);
pdfvals=100*pdfvals/sum(pdfvals);
bar(bin,pdfvals);
end