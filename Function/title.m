function [hh,hhsub] = title(varargin)
%TITLE  Graph title.
%   TITLE('txt') adds the specified title to the axes or chart returned by
%   the gca command. Reissuing the title command causes the new title to 
%   replace the old title.
%
%   TITLE('txt','subtxt') adds a subtitle in addition to the title.
%
%   TITLE(...,'Property1',PropertyValue1,'Property2',PropertyValue2,...)
%   sets the values of the specified properties of the title, and subtitle
%   if specified.
%
%   TITLE(target,...) adds the title to the specified target object.
%
%   H = TITLE(...) returns the handle to the text object used as the title.
%
%   See also XLABEL, YLABEL, ZLABEL, TEXT, SUBTITLE.

%   Copyright 1984-2020 The MathWorks, Inc.

% if the input has a title property which is a text object, use it to set
% the title on.

if nargout>0
    [isaxarr,hh]=matlab.graphics.chart.internal.objArrayDispatch(@title,varargin{:});
else
    isaxarr=matlab.graphics.chart.internal.objArrayDispatch(@title,varargin{:});
end
if isaxarr
    % Warn when an array of targets and a second output argument
    if nargout>1
        hhsub=[];
        warning(message('MATLAB:title:MultiOutputMultiTarget'))
    end
    return
end

[ax,args,nargs] = labelcheck('Title',varargin);

if nargs == 0 
  error(message('MATLAB:title:InvalidNumberOfInputs'))
end

if isempty(ax)
    ax = gca;
    % Chart subclass support
    % Invoke title method with same number of outputs to defer output arg
    % error handling to the method.
    if isa(ax,'matlab.graphics.chart.Chart')
        if(nargout == 1)
            hh = title(ax,args{:});
        else
            title(ax,args{:});
        end
        return
    end
end

dosubtitle=false;
if (nargs > 1 && (rem(nargs-1,2) ~= 0))
    subtitlestr=string(args{2});
    args(2)=[];
    dosubtitle=true;
end

titlestr = args{1};
if isempty(titlestr), titlestr=''; end
pvpairs = args(2:end);

% get-set does not support strings as of now
pvpairs = matlab.graphics.internal.convertStringToCharArgs(pvpairs);

%---Check for bypass option
if isappdata(ax,'MWBYPASS_title')       
   h = mwbypass(ax,'MWBYPASS_title',titlestr,pvpairs{:});

%---Standard behavior      
else
   matlab.graphics.internal.markFigure(ax);
   h = get(ax,'Title');
   set(h, 'String', titlestr, pvpairs{:});
   
    if dosubtitle && isprop(ax,'Subtitle')
        hSub = get(ax,'Subtitle');
        set(hSub, 'String', subtitlestr, pvpairs{:});
    end
end

if nargout > 0
  hh = h;
end
if nargout > 1 
    if exist('hSub','var')
        hhsub=hSub;
    elseif isprop(ax,'Subtitle')
        hhsub=get(ax,'Subtitle');
    else
        hhsub=[];
    end
end
