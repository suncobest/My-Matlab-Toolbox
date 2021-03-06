function [xd,dxddk,dxddx] = apply_distortion(x,k)

assert(ismatrix(x) && size(x,1)==2,'Unexpected dimension of the 1st argument!');

% Complete the distortion vector if you are using the simple distortion model:
m = length(k);
if m<5,
    flag = 1;
    k = [k(:); zeros(5-m,1)];
else
    flag = 0;
end;

[~,n] = size(x);
% Add distortion:
r2 = x(1,:).^2 + x(2,:).^2;
r4 = r2.^2;
r6 = r2.^3;

% Radial distortion:
cdist = 1 + k(1) * r2 + k(2) * r4 + k(5) * r6;
xd1 = x .* (ones(2,1)*cdist);

% tangential distortion:
a1 = 2.*x(1,:).*x(2,:);
a2 = r2 + 2*x(1,:).^2;
a3 = r2 + 2*x(2,:).^2;

delta_x = [k(3)*a1 + k(4)*a2 ;
   k(3) * a3 + k(4)*a1];

xd = xd1 + delta_x;

if nargout > 1,
    dcdistdk = [r2', r4', zeros(n,2), r6'];
    dxd1dk = zeros(2*n,5);
  	dxd1dk(1:2:end,:) = (x(1,:)'*ones(1,5)) .* dcdistdk;
  	dxd1dk(2:2:end,:) = (x(2,:)'*ones(1,5)) .* dcdistdk;
  	ddelta_xdk = zeros(2*n,5);
  	ddelta_xdk(1:2:end,3) = a1';
  	ddelta_xdk(1:2:end,4) = a2';
  	ddelta_xdk(2:2:end,3) = a3';
  	ddelta_xdk(2:2:end,4) = a1';
    dxddk = dxd1dk + ddelta_xdk ;
    if flag,
        dxddk = dxddk(:,1:m);
    end;
    if nargout > 2,
        d1 = 1 + k(1).*a2' + k(2).*r2'.*(a2+2*x(1,:).^2)' + k(5).*r4'.*(a2+4*x(1,:).^2)' + 2*k(3).*x(2,:)' + 6*k(4).*x(1,:)';
        d2 = 1 + k(1).*a3' + k(2).*r2'.*(a3+2*x(2,:).^2)' + k(5).*r4'.*(a3+4*x(2,:).^2)' + 6*k(3).*x(2,:)' + 2*k(4).*x(1,:)';
        d3 = (k(1) + 2*k(2).*r2' + 3*k(5).*r4').*a1' + 2*k(3).*x(1,:)' + 2*k(4).*x(2,:)';
        i = [1:2:2*n, 1:2:2*n, 2:2:2*n, 2:2:2*n];
        j = [1:2:2*n, 2:2:2*n, 1:2:2*n, 2:2:2*n];
        s = [d1', d3', d3', d2'];
        dxddx = sparse(i, j, s, 2*n, 2*n);
    end
end;

return;


%% Test of the Jacobians:
n = 100;
x = 10*randn(2,n);
k = 0.5*randn(5,1);
[xd,dxddk,dxddx] = apply_distortion(x,k);

% Test on k: OK!!
dk = 0.001 * norm(k)*randn(5,1);
k2 = k + dk;
[x2] = apply_distortion(x,k2);
x_pred = xd + reshape(dxddk * dk,2,n);
norm(x2-xd)/norm(x2-x_pred)

% Test on x:
dx = 0.000001 *randn(2,n);
x2 = x + dx;
xd2 = apply_distortion(x2,k);
x_pred = xd + reshape(dxddx*dx(:),2,n);
norm(xd2-xd)/norm(xd2-x_pred)
