function im2 = imclearborder(varargin)
%IMCLEARBORDER Suppress light structures connected to image border.
%   IM2 = IMCLEARBORDER(IM) suppresses structures that are lighter than
%   their surroundings and that are connected to the image border.  IM can
%   be an intensity or binary image.  The output image, IM2, is intensity or
%   binary, respectively.  The default connectivity is 8 for two dimensions,
%   26 for three dimensions, and CONNDEF(NDIMS(BW),'maximal') for higher
%   dimensions.
%
%   For intensity images, IMCLEARBORDER tends to reduce the overall
%   intensity level in addition to suppressing border structures.
%
%   IM2 = IMCLEARBORDER(IM,CONN) specifies the desired connectivity.
%   CONN may have the following scalar values:  
%
%       4     two-dimensional four-connected neighborhood
%       8     two-dimensional eight-connected neighborhood
%       6     three-dimensional six-connected neighborhood
%       18    three-dimensional 18-connected neighborhood
%       26    three-dimensional 26-connected neighborhood
%
%   Connectivity may be defined in a more general way for any dimension by
%   using for CONN a 3-by-3-by- ... -by-3 matrix of 0s and 1s.  The 1-valued
%   elements define neighborhood locations relative to the center element of
%   CONN.  CONN must be symmetric about its center element.
%
%   Class Support
%   -------------
%   IM can be a numeric or logical array of any dimension, and it must be
%   nonsparse and real.  IM2 has the same class as IM.
%
%   Note
%   ----
%   A pixel on the edge of the input image might not be considered to be
%   a "border" pixel if a nondefault connectivity is specified.  For
%   example, if CONN = [0 0 0; 1 1 1; 0 0 0], elements on the first and
%   last row are not considered to be border pixels because, according to
%   that connectivity definition, they are not connected to the region
%   outside of the image.
%
%   References: 
%   -----------
%   P. Soille, Morphological Image Analysis: Principles and Applications,
%   Springer, 1999, pp. 164-165.
%
%   Example 1
%   ---------
%   % This example shows how to clear the border of an intensity image:
%
%       I = imread('rice.png');
%       I2 = imclearborder(I);
%       montage({I,I2})
%
%   Example 2
%   ---------
%   % This example shows how to clear the border of a binary image:
%
%       BW = imbinarize(imread('rice.png'));
%       BW2 = imclearborder(BW);
%       montage({I,I2})
%
%   Example 3
%   ---------
%   % This example shows how to retain some borders while clearing:
% 
%       I = imread('rice.png');
%       % Specify border, in order from: top, right, bottom, left.
%       bordersToKeep = [0 1 0 1]; % keep right and left
%       % Top and left
%       Ipad = padarray(I, [bordersToKeep(1), bordersToKeep(4)], 'pre');
%       % Bottom and right
%       Ipad = padarray(Ipad, [bordersToKeep(3), bordersToKeep(2)], 'post');
%       I2 = imclearborder(Ipad);
%       % Unpad
%       topLeft = [bordersToKeep(1) bordersToKeep(4)]+1;
%       botRight = topLeft + [size(I,1), size(I,2)] - 1;
%       I2 = I2(topLeft(1):botRight(1), topLeft(2):botRight(2));
%       montage({I,I2})
%
%   See also IMRECONSTRUCT.

%   Copyright 1993-2019 The MathWorks, Inc.

[im,conn] = parse_inputs(varargin{:});
conn = conn2array(conn);

marker = im;

% Now figure out which elements of the marker image are connected to the
% outside, according to the connectivity definition.
im2 = true(size(marker));
im2 = padarray(im2, ones(1,ndims(im2)), 0, 'both');
im2 = imerode(im2,conn);
idx = cell(1,ndims(im2));
for k = 1:ndims(im2)
    idx{k} = 2:(size(im2,k) - 1);
end
im2 = im2(idx{:});

% Set all elements of the marker image that are not connected to the
% outside to the lowest possible value.
if islogical(marker)
    marker(im2) = false;
else
    marker(im2) = -Inf;
end

im2 = imreconstruct(marker, im, conn);
if islogical(im2)
    im2 = im & ~im2;
else
    im2 = im - im2;
end

%%%
%%% parse_inputs
%%%
function [im,conn] = parse_inputs(varargin)


narginchk(1,2);

im = varargin{1};
validateattributes(im, {'numeric' 'logical'}, {'nonsparse' 'real'}, ...
           mfilename, 'IM', 1);

if nargin < 2
    conn = conndef(ndims(im),'maximal');
else
    conn = varargin{2};
    iptcheckconn(conn,mfilename,'CONN',2);
end

% Skip NaN check here; it will be done by imreconstruct if input
% is double.
