% H is the homography between m and M
% m is the measured 2D points
% M is the model 2D points

m=sym('m%d',[2 1]); 
m=sym(m,'real');
m=[m;1];
M=sym('M%d',[2 1]);
M=sym(M,'real');
M=[M;1];

% The homography
H=sym('h%d%d',[3 3]);
H=sym(H,'real');

%%% m=H*M (ideal formula)

% Let us reproject now model points on the image plane
mp = H*M;
mpn = mp/mp(3);

% calculate the geometric distance
dm=m-mpn;

% Evaluate the symbolic expression of the Jacobian w.r.t. the estimated parameters
param=reshape(H',9,1);
J=jacobian(dm(1:2),param);



% m1=m(1);
% m2=m(2);
% M1=M(1);
% M2=M(2);
% h11=H(1,1);
% h12=H(1,2);
% h13=H(1,3);
% h21=H(2,1);
% h22=H(2,2);
% h23=H(2,3);
% h31=H(3,1);
% h32=H(3,2);
% h33=H(3,3);


%% another configuration
m=sym('m%d',[3 1]); 
m=sym(m,'real');
mn=m/norm(m);
M=sym('M%d',[3 1]);
M=sym(M,'real');
Mn=M/norm(M);
% The homography
H=sym('h%d%d',[3 3]);
H=sym(H,'real');

%%% m=H*M (ideal formula)

% Let us reproject now model points on the image plane
mp = H*M;
mpn = mp/norm(mp);

% calculate the geometric distance
dm=mn-mpn;

% Evaluate the symbolic expression of the Jacobian w.r.t. the estimated parameters
param=reshape(H',9,1);
J=jacobian(dm,param);



syms lmp mp1 mp2 mp3 real
J=subs(J,norm(mp),lmp);
J=simplify(J);
J=subs(J,mp(1),mp1);
J=subs(J,mp(2),mp2);
J=subs(J,mp(3),mp3);
J=simplify(J);