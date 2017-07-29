function [e, thresh]=canny_edge(Im,th,sigma)
%canny_edge:Find edges in intensity image.
%canny_edge takes an intensity image I as its input, and returns a binary image
% BW of the same size as I, with 1's where the function finds edges in I
% and 0's elsewhere.
%
% The Canny method finds edges by looking for local maxima of the
% gradient of I. The gradient is calculated using the derivative of a
% Gaussian filter. The method uses two thresholds, to detect strong
% and weak edges, and includes the weak edges in the output only if
% they are connected to strong edges. This method is therefore less
% likely than the others to be "fooled" by noise, and more likely to
% detect true weak edges.
%
%
% Canny Method
% ----------------------------
% BW = canny_edge(I) specifies the Canny method.
%
% BW = canny_edge(I,THRESH) specifies sensitivity thresholds for the
% Canny method. THRESH is a two-element vector in which the first element
% is the low threshold, and the second element is the high threshold. If
% you specify a scalar for THRESH, this value is used for the high
% threshold and 0.4*THRESH is used for the low threshold. If you do not
% specify THRESH, or if THRESH is empty ([]), EDGE chooses low and high
% values automatically.
%
% BW = canny_edge(I,THRESH,SIGMA) specifies the Canny method, using
% SIGMA as the standard deviation of the Gaussian filter. The default
% SIGMA is sqrt(2); the size of the filter is chosen automatically, based
% on SIGMA.
%
% [BW,thresh] = canny_edge(I,...) returns the threshold values as a
% two-element vector.
%
% Class Support
% -------------
% I can be of class uint8, uint16, or double. BW is of class uint8.
%
% Example
% -------
% Find the edges of the rice.tif image using the Prewitt and Canny
% methods:
%
% I = imread('rice.png');
% BW1 = edge(I,'prewitt');
% BW2 = canny_edge(I);
% figure, imshowpair(BW1,BW2,'montage');
%
% Clay M. Thompson 10-8-92
% Revised by Chris Griffin, 1996,1997
% Revised by Pengfei Zhang, 2016
% Copyright 1993-1998 The MathWorks, Inc. All Rights Reserved.

assert(ismatrix(Im),'The 1st input is assumed to be gray scale image (matrix)!');
if isa(Im, 'uint8') || isa(Im, 'uint16'),
    Im = im2double(Im);
end;
[m,n]=size(Im);
e = false( m, n);

PercentNotEdges=0.7;%用于计算边缘门限
ThresholdRatio=0.4;%设置两个门限的比例

if nargin<3 || isempty(sigma),
    sigma = sqrt(2);
end;

% bigeps=0.001;  %设定高斯函数消失门限
% pw = 1:30;
% n = find(exp(-pw.^2/(2*ssq))>bigeps,1,'last');
% if isempty(n)
%     n = 1;  %当使用者键入很小的sigma时
% end;

% gau2 = gau'*gau;     % 二维高斯滤波器
% dgau2dx = -x(ones(2*a+1,1), :).*gau2/ssq;  %高斯滤波器x方向的一阶导数
% negVals = dgau2dx < 0;
% posVals = dgau2dx > 0;
% dgau2dx(posVals) = dgau2dx(posVals)/sum(dgau2dx(posVals));
% dgau2dx(negVals) = dgau2dx(negVals)/abs(sum(dgau2dx(negVals)));
% % gradient of gaussian image
% dx = imfilter(Im, dgau2dx, 'conv','replicate');
% dy = imfilter(Im, dgau2dx', 'conv','replicate');

[dx, dy] = gaussGradient(Im, sigma);
% Calculate Magnitude of Gradient
magGrad = hypot(dx, dy);

% Normalize for threshold selection
magmax = max(magGrad(:));
if magmax > 0,
    magGrad = magGrad / magmax;
end;

if nargin<2 || isempty(th),
    counts=imhist(magGrad, 64);
    highThresh = find(cumsum(counts) > PercentNotEdges*m*n, 1,'first') / 64;
    lowThresh = ThresholdRatio*highThresh;
elseif length(th)==1,
    if th>=1,
        error('Threshold Must Be Less Than One!');
    end;
    highThresh = th;
    lowThresh = ThresholdRatio*th;
elseif length(th)==2,
    lowThresh = th(1);
    highThresh = th(2);
    if (lowThresh >= highThresh) || (highThresh >= 1),
        error('Threshold Out Of Range!');
    end;
end;

thresh=[lowThresh,highThresh];

 % The next step is to do the non-maximum supression.
% We will accrue indices which specify ON pixels in strong edgemap
% The array e will become the weak edge map.
idxStrong = [];
for dir = 1:4,
    idxLocalMax =  cannyFindLocalMaxima(dir,dx,dy,magGrad);
    idxWeak = idxLocalMax(magGrad(idxLocalMax) > lowThresh);
    e(idxWeak)=1;
    idxStrong = [idxStrong; idxWeak(magGrad(idxWeak) > highThresh)];
end
rstrong = rem(idxStrong-1, m)+1;
cstrong = floor((idxStrong-1)/m)+1;
e = bwselect(e, cstrong, rstrong, 8);
% e = bwmorph(e, 'thin', inf);

if nargout==0,
    imshow(e);
end;
return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : cannyFindLocalMaxima
%
function idxLocalMax = cannyFindLocalMaxima(direction,gx,gy,mag)
%
% This sub-function helps with the non-maximum suppression in the Canny
% edge detector.  The input parameters are:
%
%   direction - the index of which direction the gradient is pointing,
%               read from the diagram below. direction is 1, 2, 3, or 4.
%   gx        - input image filtered by derivative of gaussian along x
%   gy        - input image filtered by derivative of gaussian along y
%   mag       - the gradient magnitude image
%
%    there are 4 cases: boundary perpendicular to vector [gx, gy]
%
%             3              2
%       O----0----0
%     4 |                           | 1
%        |                           |
%        O            X            O
%        |                           |
%    (1)|                          |(4)
%       O----O----O
%             (2)           (3)
%
%     The X marks the pixel in question, and each of the quadrants for the
%     gradient vector fall into two cases, divided by the 45 degree line.
%     In one case the gradient vector is more horizontal, and in the other
%     it is more vertical.  There are eight divisions, but for the
%     non-maximum suppression we are only worried about 4 of them since we
%     use symmetric points about the center pixel.

[m,n] = size(mag);

% Find the indices of all points whose gradient (specified by the vector
% (gx,gy)) is going in the direction we're looking at.

switch direction
    case 1
        idx = find((gy<=0 & gx>-gy)  | (gy>=0 & gx<-gy));
    case 2
        idx = find((gx>0 & -gy>=gx)  | (gx<0 & -gy<=gx));
    case 3
        idx = find((gx<=0 & gx>gy) | (gx>=0 & gx<gy));
    case 4
        idx = find((gy<0 & gx<=gy) | (gy>0 & gx>=gy));
end

% Exclude the exterior pixels
if ~isempty(idx),
    v = mod(idx,m);
    extIdx = (v==1 | v==0 | idx<=m | (idx>(n-1)*m));
    idx(extIdx) = [];
end;

ixv = gx(idx);
iyv = gy(idx);
gradmag = mag(idx);

% Do the linear interpolations for the interior pixels
switch direction
    case 1
        d = abs(iyv./ixv);
        gradmag1 = mag(idx+m).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx-m).*(1-d) + mag(idx-m+1).*d;
    case 2
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx-m+1).*d;
    case 3
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx+m+1).*d;
    case 4
        d = abs(iyv./ixv);
        gradmag1 = mag(idx-m).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+m).*(1-d) + mag(idx+m+1).*d;
end
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2);
return;
