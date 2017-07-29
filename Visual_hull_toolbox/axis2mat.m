function	out=axis2mat(in,outdims)

% AXIS2MAT	Transform rotation vectors into rotation matrices.
% in : rotation vectors; dimension of (3,n);
% out: rotation matrices; dimension of (3,3,n) or (3,3*n)
% outdims: dimentions of output; default and 3---(3,3,n); 2---(3,3*n)

[m,n] = size(in);
if m~=3   %% in is a rotation vector list
    error('Omega (rotation) have to be 3 dimensional! (3 rows)')
end

theta = sqrt(sum(in.^2));   % 沿z方向的矢量
logid = theta<eps;
theta(logid)=inf;  % 将趋于零的除数赋值为inf
omega= in./repmat(theta,3,1);  
theta(logid)=0;     % 若theta接近于0，则对应的omega和theta都取0（初值为0）

costheta = reshape(cos(theta),1,1,[]);
sintheta = reshape(sin(theta),1,1,[]);
omega = reshape(omega,3,1,[]);  

alpha = omega(1,1,:);
beta = omega(2,1,:);
gama = omega(3,1,:);

A = repmat(omega,1,3);
B = A.*permute(A,[2 1 3]) .* repmat(1-costheta,3,3);
A = zeros(1,1,n);
C = [A, -gama, beta;
     gama, A, -alpha;
     -beta, alpha, A] .* repmat(sintheta,3,3);
A = repmat(eye(3),1,1,n) .* repmat(costheta,3,3);

if nargin==1||outdims==3
    out = A+B+C;
elseif outdims==2
    out = reshape(A+B+C,3,[]);
else
    error('Unexpected dimension!');
end

return;

   