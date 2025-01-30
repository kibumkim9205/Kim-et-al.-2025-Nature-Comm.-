function [pdfvals,cdfvals]=pdfcdf(vals,bin,varargin)
pdfvals=histc(vals,bin);
pdfvals=100*pdfvals/sum(pdfvals);
cdfvals=cumsum(pdfvals);
cdfvals=100*cdfvals/max(cdfvals);
cdfvals(cdfvals==max(cdfvals))=99;
if ~isempty(varargin)
    if strcmp(varargin{1},'plot')
        figure;
        line(bin,pdfvals);
        hold on;
        line(bin,cdfvals,'color','r');
    end
end
end