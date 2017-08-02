% DIJKSTRA Calculate Minimum Costs and Paths using Dijkstra's Algorithm
%
% Reference: Joseph Kirk's dijkstra.m
%
% Inputs:
%         A:  is a NxN adjacency matrix, where A(I,J) is nonzero (weight or distance)
%               from point I to point J. Works for both symmetric and asymmetric matrix.
%        SID: (optional) 1xL vector of starting points
%               if unspecified, the algorithm will calculate the minimal path from
%               all N points to the finish point(s) (automatically sets SID = 1:N)
%        FID: (optional) 1xM vector of finish points
%               if unspecified, the algorithm will calculate the minimal path from
%               the starting point(s) to all N points (automatically sets FID = 1:N)
%
% Outputs:
%     Dijkstra(int v0  costs:   is an LxM matrix of minimum cost values for the minimal paths
%      paths:   is an LxM cell array containing the shortest path arrays
%
% Usage:
%     [costs,paths] = dijkstra(A)
%         -or-
%     [costs,paths] = dijkstra(A,SID,FID)
%
% Example:
%     % Calculate the (all pairs) shortest distances and paths
%     n = 7; A = zeros(n); xy = 10*rand(n,2)
%     tri = delaunay(xy(:,1),xy(:,2));
%     I = tri(:); J = tri(:,[2 3 1]); J = J(:);
%     IJ = I + n*(J-1); A(IJ) = 1;
%     a = 1:n; b = a(ones(n,1),:);
%     C = reshape(sqrt(sum((xy(b,:) - xy(b',:)).^2,2)),n,n);
%     AC = A.*C;
%     [costs,paths] = dijkstra(AC)
%
% Example:
%     % Calculate the shortest distance and path between two points
%     n = 1000; A = zeros(n); xy = 10*rand(n,2);
%     tri = delaunay(xy(:,1),xy(:,2));
%     I = tri(:); J = tri(:,[2 3 1]); J = J(:);
%     IJ = I + n*(J-1); A(IJ) = 1;
%     a = 1:n; b = a(ones(n,1),:);
%     C = reshape(sqrt(sum((xy(b,:) - xy(b',:)).^2,2)),n,n);
%     A(C>0.75) = 0;
%     AC = A.*C;
%     [costs,paths] = dijkstra(AC,1,n);
%     gplot(A,xy,'k.:'); hold on;
%     plot(xy(paths,1),xy(paths,2),'ro-','LineWidth',2); hold off
%     title(sprintf('Distance from %d to %d = %1.3f',1,n,costs))
%
% Web Resources:
%   <a href="http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm">Dijkstra's Algorithm</a>
%   <a href="http://en.wikipedia.org/wiki/Graph_%28mathematics%29">Graphs</a>
%   <a href="http://en.wikipedia.org/wiki/Adjacency_matrix">Adjacency Matrix</a>
%
% See also: gplot, gplotd, gplotdc.

function [costs,paths] = dijkstra(A,SID,FID)

narginchk(1,3);
[n,nc] = size(A);
assert(n==nc,'The first input is assumed to be an adjacency matrix!');
if nargin < 3,
    FID = 1:n;
    if nargin < 2,
        SID = 1:n;
    end;
end;

A(A==inf) = 0;
L = length(SID);
M = length(FID);
% Reverse the algorithm if it will be more efficient
if L > M,
    A = A';
    tmp = L;
    L = M;
    M = tmp;
    tmp = SID;
    SID = FID;
    FID = tmp;
    flag = true;
else
    flag = false;
end;

% Initialize output variables
costs = zeros(L,M);
paths = cell(L,M);

% Find the minimum costs and paths using Dijkstra's Algorithm
for k = 1:L,
    i = SID(k);
    isok = false(1,n);
    isok(i) = true;
    dist = inf(1,n);
    dist(i) = 0;
    road = cell(1,n);
    road{i} = i;
    % Execute Dijkstra's Algorithm for this point
    while any(~isok(FID)),
        % Calculate the costs to the neighbor nodes and record paths
        for j = find(A(i,:)),
            if ~isok(j),
                c = dist(i) + A(i,j);
                if dist(j) > c,
                    dist(j) = c;
                    if flag,
                        road{j} = [j, road{i}];
                    else
                        road{j} = [road{i}, j];
                    end;
                end;
            end;
        end;
        % solve the minimum distance from source
        ind = find(~isok);
        [c,id] = min(dist(ind));
        if isinf(c),
            break;
        else
            i = ind(id);
            isok(i) = true;
        end;
    end;
    % Store costs and paths
    costs(k,:) = dist(FID);
    paths(k,:) = road(FID);
end;

% Reformat outputs if algorithm was reversed
if flag,
    costs = costs';
    paths = paths';
end;

% Pass the path as an array if only one source/destination were given
if L == 1 && M == 1,
    paths = paths{1};
end;

return;



%% Test
n = 7; A = zeros(n); xy = 10*rand(n,2)
tri = delaunay(xy(:,1),xy(:,2));
I = tri(:); J = tri(:,[2 3 1]); J = J(:);
IJ = I + n*(J-1); A(IJ) = 1;
a = 1:n; b = a(ones(n,1),:);
C = reshape(sqrt(sum((xy(b,:) - xy(b',:)).^2,2)),n,n);
AC = A.*C;
[costs,paths] = dijkstra(AC,[1 3 4],[2 3 5 7]);
[costs1,paths1] = dijkstra1(A,xy,[1 3 4],[2 3 5 7]);
