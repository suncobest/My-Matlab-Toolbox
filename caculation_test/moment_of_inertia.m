% Xnew = R*X
% Xnew'*Inew*Xnew = X'*I*X

syms Ix Iy Ixy theta real
I = [Ix,-Ixy;-Ixy,Iy]
R=[cos(theta),sin(theta);-sin(theta),cos(theta)]
Inew = simplify(R*I*R')
