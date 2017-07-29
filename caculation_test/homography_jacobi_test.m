%%% test of jacobi
% mp=H*M (ideal formula)
% dm=m-mp;
% param=reshape(H',9,1);
% J=jacobian(dm(1:2),param);
% 则m对param的导数为-J


H = 10*randn(3);
M = 5*randn(2,1);
M = [M;1];
mp = H*M;
m = mp(1:2)/mp(3);
dH = randn(3)/1000;
mnew = (H+dH)*M;
mnew = mnew(1:2)/mnew(3);
dm = mnew-m;


hr1M = mp(1);  % hr1M = h13 + M1*h11 + M2*h12;
hr2M = mp(2);  % hr2M = h23 + M1*h21 + M2*h22;
hr3M = mp(3);  % hr3M = h33 + M1*h31 + M2*h32;

J(1, 1:3) = -M'/hr3M;
J(1, 4:6) = 0;
J(1, 7:9) = M'*hr1M/hr3M^2;
J(2, 1:3) = 0;
J(2, 4:6) = -M'/hr3M;
J(2, 7:9) = M'*hr2M/hr3M^2;

m2 = m-J*reshape(dH',9,1);
dm2 = mnew-m2;
tol = norm(dm)/norm(dm2);  % 一阶残差/二阶残差

%%% 测试何时失效
% 测试结论：当H将M投影到无穷远时，即m(3)=0时，J的某些项趋于无穷大

% for i=1:10000
% homography_jacobi_test
% if tol<5
% i
% M
% H
% dH
% J
% m
% mnew
% break
% end
% end


%% 统计
N_data =10000;
H = 10*randn(3,3);
M = 5*randn(2,N_data);
M = [M;ones(1,N_data)];
mp = H*M;
m = mp(1:2,:)./(ones(2,1)*mp(3,:));
dH = randn(3)/1000;
mnew = (H+dH)*M;
mnew = mnew(1:2,:)./(ones(2,1)*mnew(3,:));
dm = mnew-m;   %f(x+dx)-f(x)


hr1M = mp(1,:);  % hr1M = h13 + M1*h11 + M2*h12;
hr2M = mp(2,:);  % hr2M = h23 + M1*h21 + M2*h22;
hr3M = mp(3,:);  % hr3M = h33 + M1*h31 + M2*h32;

J= zeros(2*N_data,9);
J(1:2:2*N_data, 1:3) = -(M./hr3M(ones(3,1),:))';
J(1:2:2*N_data, 7:9) = (M.*(ones(3,1)*(hr1M./hr3M.^2)))';
J(2:2:2*N_data, 4:6) = -(M./hr3M(ones(3,1),:))';
J(2:2:2*N_data, 7:9) = (M.*(ones(3,1)*(hr2M./hr3M.^2)))';

m2 = m(:)-J*reshape(dH',9,1);    % f(x)+f'(x)dx; f'(x)=-J
dm2 = mnew-reshape(m2,2,N_data);   % f(x+dx)-f(x)-f'(x)dx
tol = sqrt((dm(1,:).^2+dm(2,:).^2)./(dm2(1,:).^2+dm2(2,:).^2));  % 一阶残差/二阶残差
min_ratio = min(tol)
max_ratio = max(tol)
id_min = find(tol==min_ratio)
mp_id = mp(:,id_min)
lM_id = H(3,:)*M(:,id_min)


%%  单位归一方式: mn=m/norm(m)
N_data =100000;
H = 10*randn(3,3);
while cond(H)<50
    H = 10*randn(3,3);
end
M = 5*randn(3,N_data);

% LM = sqrt(sum(M.^2));
% id = find(LM>1e-5);
% N_data = length(id);
% mp = H*M(:,id);

mp = H*M;
lmp = sqrt(sum(mp.^2));
m = mp./(ones(3,1) * lmp);
dH = randn(3)/1000;
mnew = (H+dH)*M;
mnew = mnew./(ones(3,1)*sqrt(sum(mnew.^2)));
dm = mnew-m;   %f(x+dx)-f(x)


hr1M = mp(1,:);  % hr1M = h13 + M1*h11 + M2*h12;
hr2M = mp(2,:);  % hr2M = h23 + M1*h21 + M2*h22;
hr3M = mp(3,:);  % hr3M = h33 + M1*h31 + M2*h32;

% J =
% [ -(M1*(lmp^2 - mp1^2))/lmp^3, -(M2*(lmp^2 - mp1^2))/lmp^3, -(M3*(lmp^2 - mp1^2))/lmp^3,          (M1*mp1*mp2)/lmp^3,          (M2*mp1*mp2)/lmp^3,          (M3*mp1*mp2)/lmp^3,          (M1*mp1*mp3)/lmp^3,          (M2*mp1*mp3)/lmp^3,          (M3*mp1*mp3)/lmp^3]
% [          (M1*mp1*mp2)/lmp^3,          (M2*mp1*mp2)/lmp^3,          (M3*mp1*mp2)/lmp^3, -(M1*(lmp^2 - mp2^2))/lmp^3, -(M2*(lmp^2 - mp2^2))/lmp^3, -(M3*(lmp^2 - mp2^2))/lmp^3,          (M1*mp2*mp3)/lmp^3,          (M2*mp2*mp3)/lmp^3,          (M3*mp2*mp3)/lmp^3]
% [          (M1*mp1*mp3)/lmp^3,          (M2*mp1*mp3)/lmp^3,          (M3*mp1*mp3)/lmp^3,          (M1*mp2*mp3)/lmp^3,          (M2*mp2*mp3)/lmp^3,          (M3*mp2*mp3)/lmp^3, -(M1*(lmp^2 - mp3^2))/lmp^3, -(M2*(lmp^2 - mp3^2))/lmp^3, -(M3*(lmp^2 - mp3^2))/lmp^3]


J= zeros(3*N_data,9);
J(1:3:3*N_data, 1:3) = -(M.*(ones(3,1)*((lmp.^2-hr1M.^2)./lmp.^3)))';
J(1:3:3*N_data, 4:6) = (M.*(ones(3,1)*(hr1M.*hr2M./lmp.^3)))';
J(1:3:3*N_data, 7:9) = (M.*(ones(3,1)*(hr1M.*hr3M./lmp.^3)))';
J(2:3:3*N_data, 1:3) = (M.*(ones(3,1)*(hr1M.*hr2M./lmp.^3)))';
J(2:3:3*N_data, 4:6) = -(M.*(ones(3,1)*((lmp.^2-hr2M.^2)./lmp.^3)))';
J(2:3:3*N_data, 7:9) = (M.*(ones(3,1)*(hr2M.*hr3M./lmp.^3)))';
J(3:3:3*N_data, 1:3) = (M.*(ones(3,1)*(hr1M.*hr3M./lmp.^3)))';
J(3:3:3*N_data, 4:6) = (M.*(ones(3,1)*(hr2M.*hr3M./lmp.^3)))';
J(3:3:3*N_data, 7:9) = -(M.*(ones(3,1)*((lmp.^2-hr3M.^2)./lmp.^3)))';

m2 = m(:)-J*reshape(dH',9,1);    % f(x)+f'(x)dx; f'(x)=-J
dm2 = mnew-reshape(m2,3,N_data);   % f(x+dx)-f(x)-f'(x)dx
tol = sqrt(sum(dm.^2)./sum(dm2.^2));  % 一阶残差/二阶残差
cond_H = cond(H)
min_ratio = min(tol)
max_ratio = max(tol)