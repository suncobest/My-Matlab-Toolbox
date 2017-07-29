function test = is3D(X)

[m,Np] = size(X);
if m~=3,
    test = 0;
    return;
end;

% Check for planarity of the structure:
X_mean = mean(X,2);
Y = X - (X_mean*ones(1,Np));
YY = Y*Y';

[~,S,~] = svd(YY);
r = S(3,3)/S(2,2);
test = (r > 1e-3);

