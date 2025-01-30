function [filenames,filedates]=getFilenamesandtimes(path)
% returns a list of the filenames in a specified directory

d0=dir(path);
files=d0([d0.isdir]==0);
filenames={files.name};
filenames=filenames(:);
filedates={files.date};
filedates=filedates(:);
[~,nondbidx]=setdiff(filenames,{'Thumbs.db','thumbs.db'});
filenames=filenames(nondbidx,:);
filedates=filedates(nondbidx,:);

