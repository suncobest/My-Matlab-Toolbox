% DIJKSTRA Calculate Minimum Costs and Paths using Dijkstra's Algorithm
%
% Reference: Joseph Kirk's dijkstra.m
%
% Inputs:
%         WA:  is a NxN adjacency matrix, where A(I,J) is nonzero (weight or distance)
%               from point I to point J. Works for both symmetric and asymmetric matrix.
%        SID: (optional) 1xL vector of starting points
%               if unspecified, the algorithm will calculate the minimal path from
%               all N points to the finish point(s) (automatically sets SID = 1:N)
%        FID: (optional) 1xM vector of finish points
%               if unspecified, the algorithm will calculate the minimal path from
%               the starting point(s) to all N points (automatically sets FID = 1:N)
%
% Outputs:
%      costs:   is an LxM matrix of minimum cost values for the minimal paths
%      paths:   is an LxM cell array containing the shortest path arrays
%
% Usage:
%     [costs,paths] = dijkstra(WA)
%         -or-
%     [costs,paths] = dijkstra(WA,SID,FID)
%
% Example:
%     % Calculate the (all pairs) shortest distances and paths
%     n = 7; A = zeros(n); xy = 10*rand(n,2)
%     tri = delaunay(xy(:,1),xy(:,2));
%     I = tri(:); J = tri(:,[2 3 1]); J = J(:);
%     IJ = I + n*(J-1); A(IJ) = 1;
%     a = 1:n; b = a(ones(n,1),:);
%     WA = round(reshape(sqrt(sum((xy(b,:) - xy(b',:)).^2,2)),n,n));
%     WA = WA.*A;
%     [costs,paths] = dijkstra(WA)
%
% Example:
%     % Calculate the shortest distance and path between two points
%     n = 1000; A = zeros(n); xy = 10*rand(n,2);
%     tri = delaunay(xy(:,1),xy(:,2));
%     I = tri(:); J = tri(:,[2 3 1]); J = J(:);
%     D = sqrt(sum((xy(I,:)-xy(J,:)).^2,2));
%     I(D > 0.75,:) = []; J(D > 0.75,:) = [];
%     IJ = I + n*(J-1); A(IJ) = 1;
%     a = 1:n; b = a(ones(n,1),:);
%     WA = round(reshape(sqrt(sum((xy(b,:) - xy(b',:)).^2,2)),n,n));
%     WA = WA.*A;
%     [cost,path] = dijkstra(WA,1,n);
%     gplot(A,xy,'k.:'); hold on;
%     plot(xy(path,1),xy(path,2),'ro-','LineWidth',2); hold off
%     title(sprintf('Distance from %d to %d = %1.3f',1,n,cost))
%
% Web Resources:
%   <a href="http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm">Dijkstra's Algorithm</a>
%   <a href="http://en.wikipedia.org/wiki/Graph_%28mathematics%29">Graphs</a>
%   <a href="http://en.wikipedia.org/wiki/Adjacency_matrix">Adjacency Matrix</a>
%
% See also: gplot, gplotd, gplotdc.

function [costs,paths] = dijkstra(WA,SID,FID)

narginchk(1,3);
[n,nc] = size(WA);
assert(n==nc,'The first input is assumed to be an adjacency matrix!');
if nargin < 3,
    FID = 1:n;
    if nargin < 2,
        SID = 1:n;
    end;
end;

WA(WA==0) = inf;
L = length(SID);
M = length(FID);
% Reverse the algorithm if it will be more efficient
if L > M,
    WA = WA';
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
for i = SID,
    iTable = NaN(n,1);
    minCost = Inf(n,1);
    isSettled = false(n,1);
    minCost(I) = 0;
    iTable(I) = 0;
    isSettled(I) = true;
    paths = {I};

    % Execute Dijkstra's Algorithm for this point
    while any(~isSettled(FID)),
        % Update the table
        jTable = iTable;
        iTable(I) = NaN;
        nodeIndex = find(E(:,1) == I);

        % Calculate the costs to the neighbor nodes and record paths
        for kk = 1:length(nodeIndex),
            J = E(nodeIndex(kk),2);
            if ~isSettled(J),
                c = cost(I,J);
                empty = isnan(jTable(J));
                if empty || (jTable(J) > (jTable(I) + c)),
                    iTable(J) = jTable(I) + c;
                    if flag,
                        path{J} = [J path{I}];
                    else
                        path{J} = [path{I} J];
                    end;
                else
                    iTable(J) = jTable(J);
                end;
            end;
        end;

        % Find values in the table
        K = find(~isnan(iTable));
        if isempty(K),
            break
        else
            % Settle the minimum value in the table
            [~,N] = min(iTable(K));
            I = K(N);
            minCost(I) = iTable(I);
            isSettled(I) = true;
        end;
    end;

    % Store costs and paths
    costs(k,:) = minCost(FID);
    paths(k,:) = path(FID);
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



%% test
