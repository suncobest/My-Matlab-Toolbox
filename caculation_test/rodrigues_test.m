clear all, clc
syms theta w1 w2 w3 real

omega =[w1 w2 w3]';
n=norm([w1 w2 w3]);

alpha = cos(theta);
beta = sin(theta);
gamma = 1-cos(theta);
omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
A = omega*omega';

R = eye(3)*alpha + omegav*beta + A*gamma

% R=subs(R,[w1 w2 w3 theta],[w1/n w2/n w3/n n]);
% RRT=simplify(R*R')

dRdin = [0 0 0;
    0 0 1;
    0 -1 0;
    0 0 -1;
    0 0 0;
    1 0 0;
    0 1 0;
    -1 0 0;
    0 0 0]


%m1 = [alpha;beta;gamma;omegav;A];
dRdm1 = zeros(9,21);
dRdm1(:,2) = omegav(:);
dRdm1([1 5 9],1) = ones(3,1);
dRdm1(:,4:12) = beta*eye(9);
dRdm1(:,3) = A(:);
dRdm1(:,13:21) = gamma*eye(9);

%m2 = [omega;theta]
dm1dm2 = zeros(21,4);
dm1dm2(1,4) = -sin(theta);
dm1dm2(2,4) = cos(theta);
dm1dm2(3,4) = sin(theta);
dm1dm2(4:12,1:3) = [0 0 0 0 0 1 0 -1 0;
    0 0 -1 0 0 0 1 0 0;
    0 1 0 -1 0 0 0 0 0]';


dm1dm2(13:21,1) = [2*w1;w2;w3;w2;0;0;w3;0;0];
dm1dm2(13: 21,2) = [0;w1;0;w1;2*w2;w3;0;w3;0];
dm1dm2(13:21,3) = [0;0;w1;0;0;w2;w1;w2;2*w3];

%m3 = [in,theta]
dm2dm3 = [eye(3)/theta -omega/theta; zeros(1,3) 1];


dm3din = [eye(3);[w1 w2 w3]];

dout = dRdm1 * dm1dm2 * dm2dm3 * dm3din
r1=reshape(dout(:,1),3,3);
r2=reshape(dout(:,2),3,3);
r3=reshape(dout(:,3),3,3);

% r1=subs(r1,[w1 w2 w3 theta],[w1/n w2/n w3/n n]);
% r2=subs(r2,[w1 w2 w3 theta],[w1/n w2/n w3/n n]);
% r3=subs(r3,[w1 w2 w3 theta],[w1/n w2/n w3/n n]);
% differ=simplify(R-r1*r2*r3)

% subs(r1,[w1 w2 w3],theta*[1 0 0])
% subs(r2,[w1 w2 w3],theta*[1 0 0])
% subs(r3,[w1 w2 w3],theta*[1 0 0])

r1=subs(r1,[w1 w2 w3],[w1/n w2/n w3/n]);
r2=subs(r2,[w1 w2 w3],[w1/n w2/n w3/n]);
r3=subs(r3,[w1 w2 w3],[w1/n w2/n w3/n]);

a1=subs(r1,[w1 w2 w3],[1 0 0])
b1=subs(r2,[w1 w2 w3],[1 0 0])
c1=subs(r3,[w1 w2 w3],[1 0 0])

a2=subs(r1,[w1 w2 w3],[0 1 0])
b2=subs(r2,[w1 w2 w3],[0 1 0])
c2=subs(r3,[w1 w2 w3],[0 1 0])

a3=subs(r1,[w1 w2 w3],[0 0 1])
b3=subs(r2,[w1 w2 w3],[0 0 1])
c3=subs(r3,[w1 w2 w3],[0 0 1])





