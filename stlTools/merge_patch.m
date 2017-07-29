function [v, f] = merge_patch( v1,f1,v2,f2 )
%MERGE_PATCH merge two patches into one piece.
%   v1: vertices of input patch 1;
%   f1: faces of input patch 1;
%   v2: vertices of input patch 2;
%   f2: faces of input patch 2;
%   v: vertices of output patch;
%   f: faces of output patch;

assert(ismatrix(v1) && ismatrix(v2) && ismatrix(f1) && ismatrix(f2), 'Input must be matrix!');
[m,nv] = size(v1);
n = size(v2,2);
assert(nv==n, 'Vertice of the two patches must have same dimension!');
assert(nv==2 || nv==3, 'Input vertice must be dimension of 2 and 3!');
nf = size(f1,2);
n = size(f2,2);
assert(nf==n, 'The two patches must have the same number of vertices on each face!');
v = [v1;v2];
f2 = f2+m;
f = [f1; f2];
[v, ~, indexn] =  unique(v, 'rows');
f = indexn(f);
return;


%%  test
[verts1,faces1] = stlRead('front_tower.stl');
[verts2,faces2] = stlRead('front_tower_top.stl');
[verts, faces] = merge_patch(verts1,faces1,verts2,faces2);
stlPlot(verts,faces, 'Tower');

