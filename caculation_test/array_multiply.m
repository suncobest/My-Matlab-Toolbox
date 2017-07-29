clear; echo off;

% [i,j,k]=ind2sub([3 3 3],1:27);
% name=['a'*ones(27,1),num2str(i'),num2str(j'),num2str(k')];
% 
% for n=1:27
%     a(n)=sym(name(n,:),'real');
% end
% 
% a=a';
% 
% a=reshape(a,3,3,3)
% a1=reshape(a,[],3)

% 如果要生成一维或二维符号数组，可以用sym('a%d%d',[M N])命令，此命令不能生成3维以上的符号数组

A=sym('a%d%d',[3 5]);
B=sym('b%d%d',[5 4]);

[p,n] = size(A); [n2,q] = size(B);

if n2 ~= n,
   error(' A and B must have equal inner dimensions');
end;

AB=A*B;

dABdA = sym(zeros(numel(AB),numel(A)));
dABdB = sym(zeros(numel(AB),numel(B)));

for i = 1:numel(AB)
    for j = 1:numel(A)
        dABdA(i,j) = diff(AB(i), A(j));
    end
end

for i = 1:numel(AB)
    for j = 1:numel(B)
        dABdB(i,j) = diff(AB(i), B(j));
    end
end

% 观察结果后，采用赋值算法

