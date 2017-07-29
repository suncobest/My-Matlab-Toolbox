% twist test

% syms w1 w2 w3 n1 n2 n3 x1 x2 x3 v1 v2 v3 theta real
% a=sym(zeros(4));
% b=a;
% a(1:3,1:3)=[0 -w3 w2
%             w3 0 -w1
%            -w2 w1 0]
% b(1:3,4)=[x1;x2;x3]
% c=a+b
% e=subs(c,[w1 w2 w3 x1 x2 x3],[n1 n2 n3 v1 v2 v3])
% simplify(c*e-e*c)

wt=randn(3,1)
theta=norm(wt);
v=randn(3,1)
kesi=zeros(4);
kesi(1:3,1:3)=skew3(wt);
kesi(1:3,4)=v*theta
R=rodrigues(wt)
Rt=expm(kesi)
t=(eye(3)-expm(skew3(wt)))*(cross(wt/theta,v))+v'*wt/theta*wt
T=(1-cos(theta))*cross(wt/theta,v)+(theta-sin(theta))*v'*wt/theta^2*wt+sin(theta)*v


% [v d]=eig(a)
% simplify(v\a*v)
% d=simplify(subs(d,w1^2+w2^2+w3^2,theta^2));
% d=simplify(subs(d,(-theta^2)^(1/2),i*theta))
% 
% ert=simplify(expm(c))
% ert=simplify(subs(ert,w1^2+w2^2+w3^2,theta^2));
% ert=simplify(subs(ert,(-theta^2)^(1/2),i*theta))
