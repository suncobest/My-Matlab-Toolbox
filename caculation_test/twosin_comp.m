function [rt,f,g] = twosin_comp(A,u0,v0,m,n)
%TWOSIN compares for loops and vectorization.
%   The comparision is based upon the function 'f(x,y)=Asin(u0x+v0y)'. 
%   'm' and 'n' stand for the grid-point numbers of axes x and y respectively.
%   'rt' is the runing time ratio of loops to vetorization.
%   'f' and 'g' are the result matrixes of the loop method and
%   vectorization method respectively.

%   the loop method
tic   % start timing
for r=1:m
    u0x=u0*(r-1);
    for c=1:n
        v0y=v0*(c-1);
        f(r,c)=A*sin(u0x+v0y);
    end
end
t1=toc;   % end timing

% the vectorization method
tic   % start timing
r=0:m-1;
c=0:n-1;
[C,R]=meshgrid(c,r);
g=A*sin(u0*R+v0*C);
t2=toc;   % end timing

% compute the ratio of the two method running time.
rt=t1/(t2+eps);   % use eps in case t2 is close to 0.

end

