function q = convex_hull(p)
% reurn the subset of points that lie on the convex hull for 2 or 3 dimensions
% see also convhull, convhulln, convexHull, delaunayTriangulation

m = size(p,1);
if m== 2,
    k = convhull( p(1,:), p(2,:) );
    q = p(:, k(1:end-1) );
elseif m == 3,
    k = convhull( p' );
    q = p(:, unique(k) );
else
    error('Input matrix must be 2 or 3 rows!');
end;

return;

%% Test against delaunayTriangulation convexHull
p = 10*randn(3,10000);
tic; q1 = convex_hull(p); t1 = toc

tic;
DT = delaunayTriangulation(p');
k = convexHull(DT);
q2 = p(:,unique(k));
t2 = toc

diff = norm(q1-q2)