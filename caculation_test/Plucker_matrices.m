%% Plucker matrices

symvect=@(name,length) sym([name '%d'],[length,1]);   % name is a string
symmat=@(name,m,n) sym([name '%d%d'],[m,n]);
symreal=@(x) sym(x,'real');

plucker_m2v = @(m) [m(1,2),m(1,3),m(1,4),m(2,3),m(4,2),m(3,4)];
plucker_v2m = @(v) [0,v(1),v(2),v(3);-v(1),0,v(4),-v(5);-v(2),-v(4),0,v(6);-v(3),v(5),-v(6),0];

% A=sym('A%d',[4 1]);
% A=sym(A,'real');
% B=sym('B%d',[4 1]);
% B=sym(B,'real');
% a=sym('a%d',[4 1]);
% a=sym(a,'real');
% b=sym('b%d',[4 1]);
% b=sym(b,'real');

A=symreal(symvect('A',4));
B=symreal(symvect('B',4));
a=symreal(symvect('a',4));
b=symreal(symvect('b',4));


W=[A';B'];
L=A*B'-B*A';
l=a*b'-b*a';

Lv=plucker_m2v(L);

ABab=det(cat(2,A,B,a,b));
biprod=@(x,y) x(1,2)*y(3,4)+y(1,2)*x(3,4)+x(1,3)*y(4,2)+y(1,3)*x(4,2)+x(1,4)*y(2,3)+y(1,4)*x(2,3);
Ll=biprod(L,l);
simplify(ABab-Ll)
