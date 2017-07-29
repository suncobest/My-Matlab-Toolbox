% u_t = a*u_xx, 0<x<1,t>0,   a=1
% u(x,0)=sin(pi*x), 0<x<1;
% u(0,t)=u(1,t)=0.

lamda=0.1;  %0.1
dx=0.01;
dt=lamda*dx^2;
x = 0:dx:1;
t = 0:dt:0.5;
n=length(x);
m=length(t);

a = 1;

%% forward --- O(dt+dx^2)                      稳定条件： lamda<=0.5
% U(i+1,j)=U(i,j)+a*lamda*(U(i,j+1)-2*U(i,j)+U(i,j-1))

%  this method is faster but need more memory
tic;
U=zeros(m,n);
U(1,:)=sin(pi*x);
for i=1:m-1,
    for j=2:n-1,
        U(i+1,j)=a*lamda*(U(i,j-1)+U(i,j+1))+(1-2*a*lamda)*U(i,j);
    end;
end;
FTCS_t1 = toc
FTCS_err1 = U(m,:) - exp(-pi^2*t(m))*U(1,:);
FTCS_err1=FTCS_err1*FTCS_err1'

%% forward
tic;
U = sin(pi*x);
Ut = U;
temp = zeros(1,n);
for i=1:m-1,
    for j=2:n-1,
        temp(j)=a*lamda*(Ut(j-1)+Ut(j+1))+(1-2*a*lamda)*Ut(j);        
    end;
    Ut = temp;
end;
FTCS_t2=toc
FTCS_err2 = Ut - exp(-pi^2*t(m))*U;
FTCS_err2=FTCS_err2*FTCS_err2'

%% forward
coeff = [a*lamda,1-2*a*lamda,a*lamda];
tic;
U = sin(pi*x);
Ut = U;
for i=1:m-1,
    Ut(2:n-1) = conv(Ut(2:n-1),coeff,'same');
end;
FTCS_t3=toc
FTCS_err3 = Ut - exp(-pi^2*t(m))*U;
FTCS_err3=FTCS_err3*FTCS_err3'

%% backward --- O(dt+dx^2)
% U(i,j)=U(i-1,j)+a*lamda*(U(i,j+1)-2*U(i,j)+U(i,j-1));
% A*U(i+1,:)' = U(i,:)'; 

tic;
U = sin(pi*x);
Ut = U';
A = diag([1,(1+2*a*lamda)*ones(1,n-2),1])-diag([0,a*lamda*ones(1,n-2)],1)-diag([a*lamda*ones(1,n-2),0],-1);
for i=1:m-1,
    Ut = A\Ut;  
end;
BTCS_t1=toc
BTCS_err1 = Ut' - exp(-pi^2*t(m))*U;
BTCS_err1=BTCS_err1*BTCS_err1'


%% backward --- tridiag_lu

tic;
U = sin(pi*x);
Ut = U;

Vb = [1,(1+2*a*lamda)*ones(1,n-2),1];
Va = [-a*lamda*ones(1,n-2),0];
Vc = [0,-a*lamda*ones(1,n-2)];
Ve = zeros(1,n);
Ve(1)=Vb(1);
for i=2:n,
    Ve(i)=Vb(i)-Va(i-1)*Vc(i-1)/Ve(i-1);
end;
Vd=Va./Ve(1:n-1);

for i=1:m-1,
    for j = 2:n,
        Ut(j)=Ut(j)-Vd(j-1)*Ut(j-1);
    end;
    Ut(n)=Ut(n)/Ve(n);
    for j = n-1:-1:1,
        Ut(j)=(Ut(j)-Vc(j)*Ut(j+1))/Ve(j);
    end;    
end;
BTCS_t2=toc
BTCS_err2 = Ut - exp(-pi^2*t(m))*U;
BTCS_err2=BTCS_err2*BTCS_err2'


%% Crank-Nicolson --- O(dt^2+dx^2)
% (1+a*lamda)*U(i+1,j)-a*lamda/2*(U(i+1,j+1)+U(i+1,j-1))=(1-a*lamda)*U(i,j)+a*lamda/2*(U(i,j+1)+U(i,j-1))
% A*U(i+1,:)' = B*U(i,:)';

tic;
U = sin(pi*x);
Ut = U';

A = diag([1,(1+a*lamda)*ones(1,n-2),1])-diag([0,a*lamda/2*ones(1,n-2)],1)-diag([a*lamda/2*ones(1,n-2),0],-1);
B = diag([1,(1-a*lamda)*ones(1,n-2),1])+diag([0,a*lamda/2*ones(1,n-2)],1)+diag([a*lamda/2*ones(1,n-2),0],-1);
for i=1:m-1,
    Ut = A\B*Ut;
end;
Crank_t1=toc
Crank_err1 = Ut' - exp(-pi^2*t(m))*U;
Crank_err1=Crank_err1*Crank_err1'

%% Crank-Nicolson --- tridiag_lu
% (1+a*lamda)*U(i+1,j)-a*lamda/2*(U(i+1,j+1)+U(i+1,j-1))=(1-a*lamda)*U(i,j)+a*lamda/2*(U(i,j+1)+U(i,j-1))

tic;
U = sin(pi*x);
Ut = U;
temp = zeros(1,n);

BVb = [1,(1-a*lamda)*ones(1,n-2),1];
BVa = [a*lamda/2*ones(1,n-2),0];
BVc = [0,a*lamda/2*ones(1,n-2)];

AVb = [1,(1+a*lamda)*ones(1,n-2),1];
AVa = [-a*lamda/2*ones(1,n-2),0];
AVc = [0,-a*lamda/2*ones(1,n-2)];

AVe = zeros(1,n);
AVe(1)=AVb(1);
for i=2:n,
    AVe(i)=AVb(i)-AVa(i-1)*AVc(i-1)/AVe(i-1);
end;
AVd=AVa./AVe(1:n-1);

for i=1:m-1,
    for j = 2:n-1,
        temp(j)= BVa(j-1)*Ut(j-1)+BVb(j)*Ut(j)+BVc(j)*Ut(j+1); 
    end;
    Ut=temp;
    for j = 2:n,
        Ut(j)=Ut(j)-AVd(j-1)*Ut(j-1);
    end;
    Ut(n)=Ut(n)/AVe(n);
    for j = n-1:-1:1,
        Ut(j)=(Ut(j)-AVc(j)*Ut(j+1))/AVe(j);
    end;    
end;

Crank_t2=toc
Crank_err2 = Ut - exp(-pi^2*t(m))*U;
Crank_err2=Crank_err2*Crank_err2'


%% Du Fort-Frankel	--- O(dt^2+dx^2+dt^2/dx^2)
% U(i+1,j) = U(i-1,j) + 2*a*lamda*(U(i,j+1)-(U(i+1,j)+U(i-1,j))+U(i,j-1))
% i.e.
% (1+2*a*lamda)*U(i+1,j)=(1-2*a*lamda)*U(i-1,j)+2*a*lamda*(U(i,j+1)+U(i,j-1));

tic;
U=zeros(m,n);
U(1,:)=sin(pi*x);
for j=2:n-1,
    U(2,j)=a*lamda*(U(1,j-1)+U(1,j+1))+(1-2*a*lamda)*U(1,j);
end;
for i=2:m-1,
    for j=2:n-1,
       U(i+1,j) = ((1-2*a*lamda)*U(i-1,j)+2*a*lamda*(U(i,j+1)+U(i,j-1)))/(1+2*a*lamda);
    end;
end;
Du_t1 = toc
Du_err1 = U(m,:) - exp(-pi^2*t(m))*U(1,:);
Du_err1=Du_err1*Du_err1'


%%  Du Fort-Frankel

tic;
U=sin(pi*x);
Ut=zeros(2,n);
Ut(1,:)=U;
temp = zeros(1,n);
for j=2:n-1,
    Ut(2,j)=a*lamda*(Ut(1,j-1)+Ut(1,j+1))+(1-2*a*lamda)*Ut(1,j);
end;
for i=2:m-1,
    for j=2:n-1,
       temp(j) = ((1-2*a*lamda)*Ut(1,j)+2*a*lamda*(Ut(2,j+1)+Ut(2,j-1)))/(1+2*a*lamda);
    end;
    Ut=Ut([2 1],:);
    Ut(2,:)=temp;
end;
Du_t2 = toc
Du_err2 = Ut(2,:) - exp(-pi^2*t(m))*U;
Du_err2=Du_err2*Du_err2'


%% frog  --- equivalent to Du Fort-Frankel
% forward: even --- mod(i+j,2)=0; backward: odd --- mod(i+j,2)=1.

tic;
U = sin(pi*x);
Ut=zeros(2,n);
Ut(1,:)=U;
temp = zeros(1,n);
for j=2:2:n-1,
    Ut(2,j)=a*lamda*(Ut(1,j-1)+Ut(1,j+1))+(1-2*a*lamda)*Ut(1,j);    
end;
for i=2:2:m-1,
    for j=3:2:n-1,
        Ut(2,j)=(Ut(1,j)+a*lamda*(Ut(2,j+1)+Ut(2,j-1)))/(1+2*a*lamda);
        temp(j)=2*Ut(2,j)-Ut(1,j);
    end;
    Ut=Ut([2 1],:);
    Ut(2,:)=temp;
    for j=2:2:n-1,
        Ut(2,j)=(Ut(1,j)+a*lamda*(Ut(2,j+1)+Ut(2,j-1)))/(1+2*a*lamda);
        temp(j)=2*Ut(2,j)-Ut(1,j);
    end;
    Ut=Ut([2 1],:);
    Ut(2,:)=temp;
end;

if mod(m,2)==1,
    Ut = Ut(1,:);
else
    for j=3:2:n-1,
        Ut(2,j)=(Ut(1,j)+a*lamda*(Ut(2,j+1)+Ut(2,j-1)))/(1+2*a*lamda);
    end;
    Ut = Ut(2,:);
end;

frog_t1=toc
frog_err1 = Ut - exp(-pi^2*t(m))*U;
frog_err1=frog_err1*frog_err1'


%% frog  --- equivalent to Du Fort-Frankel
% forward: odd --- mod(i+j,2)=1; backward: even --- mod(i+j,2)=0.

tic;
U = sin(pi*x);
Ut=zeros(2,n);
Ut(1,:)=U;
temp = zeros(1,n);
for j=3:2:n-1,
    Ut(2,j)=a*lamda*(Ut(1,j-1)+Ut(1,j+1))+(1-2*a*lamda)*Ut(1,j);    
end;
for i=2:2:m-1,
    for j=2:2:n-1,
        Ut(2,j)=(Ut(1,j)+a*lamda*(Ut(2,j+1)+Ut(2,j-1)))/(1+2*a*lamda);
        temp(j)=2*Ut(2,j)-Ut(1,j);
    end;
    Ut=Ut([2 1],:);
    Ut(2,:)=temp;
    for j=3:2:n-1,
        Ut(2,j)=(Ut(1,j)+a*lamda*(Ut(2,j+1)+Ut(2,j-1)))/(1+2*a*lamda);
        temp(j)=2*Ut(2,j)-Ut(1,j);
    end;
    Ut=Ut([2 1],:);
    Ut(2,:)=temp;
end;

if mod(m,2)==1,
    Ut = Ut(1,:);
else
    for j=2:2:n-1,
        Ut(2,j)=(Ut(1,j)+a*lamda*(Ut(2,j+1)+Ut(2,j-1)))/(1+2*a*lamda);
    end;
    Ut = Ut(2,:);
end;
frog_t2=toc
frog_err2 = Ut - exp(-pi^2*t(m))*U;
frog_err2=frog_err2*frog_err2'
