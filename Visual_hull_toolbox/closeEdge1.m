function bw = closeEdge1(bw, winsize, flag)
% CloseEdge: close gaps in a binary edge map image
%
% Usage: bw2 = CloseEdge(bw, gapsize)
%
% Arguments:    bw - Binary edge image
%          winsize - The window size for function to fill big gaps.
%          flag - 1: Delete short piece of edges and fill gaps all at once (default).
%                     0: Keep short piece of edges and fill gaps iteratively with growing window size.
%
% Returns:     bw2 - The binary edge image with gaps filled.
%
%
% Strategy: extend endings in rough direction of edge if the stretching can
% reach another edge in the window.
%
% See also: canny_edge.

% Pengfei Zhang, May 2016

if nargin<3 || isempty(flag),
    flag = 1;
else
    flag = ~~flag;
end;

if nargin<2 || isempty(winsize),
    winsize = 10;
else
    winsize = max(round(winsize),4);
end;

assert(islogical(bw),'The first input is assumed to be binary image (logical)!');
[rows, cols] = size(bw);
bw = bwmorph(bw, 'thin', inf);
% Set up a look up table to find endings.
lut = makelut(@ending, 3);
ends = bwlookup(bw, lut);
[re,ce] = find(ends);

% find distance of xx less than threshold
xx = [re,ce]';
n = length(re);
a2 =  repmat(sum(xx.*xx,1),[n,1]);
d = sqrt(a2'+a2 - 2*(xx'*xx));
jj = repmat(1:n,n,1);
d(jj'<=jj)=Inf;
[d,id1] = min(d,[],1);
id2 = find(d<winsize/2);
id1 = id1(id2);
% link specified ends
for i=1:length(id1),
    ii = id1(i); jj = id2(i);
    ri = re(ii); ci = ce(ii); rj = re(jj); cj = ce(jj);
    l = cross([ri,ci,1], [rj,cj,1]);
    l = l/norm(l(1:2));
    if abs(l(1))>=abs(l(2)),
        if cj>=ci,
            y = ci:cj;
        else
            y = cj:ci;
        end;
        x = round(-(l(2)*y+l(3))/l(1));
    else
        if rj>=ri,
            x = ri:rj;
        else
            x = rj:ri;
        end;
        y = round(-(l(1)*x+l(3))/l(2));
    end;
    bw((y-1)*rows+x) = 1;
end;

bw = bwmorph(bw, 'thin', inf);  % thin again
ends = bwlookup(bw, lut);
[re,ce] = find(ends);
% generate connector
blob = true(3);
% fill small gaps by place square blob at every endpoint
for n = 1:length(re),
    r = re(n); c = ce(n);
    if  r>=2 && r<=rows-1 && c>=2 && c<=cols-1,
        bw(r-1:r+1, c-1:c+1) = bw(r-1:r+1, c-1:c+1) | blob;
    end;
end;
bw = bwmorph(bw, 'thin', inf);  % thin again

if flag,
    % delete small piece of edges
    bw = bwareaopen(bw,winsize);
else
    % keep small pieces and fill gaps iteratively
    winsize = 4:winsize;
end;

for win = winsize,
    ends = bwlookup(bw, lut);
    [re,ce] = find(ends);
    %   imshow(bw); hold on; plot(ce,re,'g+');
    % Place a circular blob at every endpoint and isolated pixel
    for n = 1:length(re),
        r = re(n); c = ce(n);
        lu = max([r-win, c-win], [1,1]);
        rb = min([r+win, c+win], [rows,cols]);
        BW = fillgap(bw(lu(1):rb(1), lu(2):rb(2)), r-lu(1)+1, c-lu(2)+1);
        if ~isempty(BW),
            bw(lu(1):rb(1), lu(2):rb(2)) = BW;
        end;
    end;
    bw = bwmorph(bw, 'thin', inf);
end;
return;


%----------------------------------------------------------------------
% Function to test whether the centre pixel within a 3x3 neighbourhood is an ending.
% The centre pixel must be set and the number of transitions/crossings between 0
% and 1 as one traverses the perimeter of the 3x3 region must be 2.
%
% Pixels in the 3x3 region are numbered as follows
%
%       1 4 7
%       2 5 8
%       3 6 9
function b = ending(x)
a = [x(1) x(2) x(3) x(6) x(9) x(8) x(7) x(4)]';
b = [x(2) x(3) x(6) x(9) x(8) x(7) x(4) x(1)]';
crossings = sum(abs(a-b));

b = x(5) && crossings == 2;
return;


%----------------------------------------------------------------------
% Function to extend ending in the edge direction.
function BW = fillgap(bw,r,c)
[rn, cn] = size(bw);
L = labelmatrix(bwconncomp(bw));
BW = bw;
[x0,y0] = find(L==L(r,c));
% checkerboard distance == 1
id = find(max(abs(x0-r),abs(y0-c))==1,1);
r1 = x0(id)-r;
c1 = y0(id)-c;
% linear regression of line close to points (x,y)
xx = [x0,y0,ones(size(x0))];
[~,~,l] = svd(xx'*xx);
l = l(:,3);
l = l/norm(l(1:2));
% make sure the line pass through [r,c]
l(3) = -[r,c]*l(1:2);
if abs(l(1))>=abs(l(2)),
    y = [1:c-1,c+1:cn];     % exclude point [r, c]
    x = -(l(2)*y+l(3))/l(1);
    x = [floor(x), ceil(x)];
    y = [y,y];
    % along the inverse direction of ending edge
    id = x>=1 & x<=rn & (x-r)*r1+(y-c)*c1<=0;
    x = x(id);
    y = y(id);
    id = find(bw((y-1)*rn+x));
    n = length(id);
    if n==0,
        BW = [];
        return;
    else
        x1 = x(id);
        y1 = y(id);
        if n>1,
            [~,id] = min((x1-r).^2+(y1-c).^2);
            % x1 = x1(id);
            y1 = y1(id);
        end;
    end;
    id = find((y-y1).*(y-c)<=0);
else
    x = [1:r-1,r+1:rn];     % exclude point [r, c]
    y = -(l(1)*x+l(3))/l(2);
    y = [floor(y),ceil(y)];
    x = [x,x];
    % along the inverse direction of ending edge
    id = y>=1 & y<=cn & (x-r)*r1+(y-c)*c1<=0;
    x = x(id);
    y = y(id);
    id = find(bw((y-1)*rn+x));
    n = length(id);
    if n==0,
        BW = [];
        return;
    else
        x1 = x(id);
        y1 = y(id);
        if n>1,
            [~,id] = min((x1-r).^2+(y1-c).^2);
            x1 = x1(id);
            % y1 = y1(id);
        end;
    end;
    id = find((x-x1).*(x-r)<=0);
end;

x1 = x(id);
y1 = y(id);
BW((y1-1)*rn+x1) = 1;

return;