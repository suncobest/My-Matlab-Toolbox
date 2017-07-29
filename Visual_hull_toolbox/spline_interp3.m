function [yy, coefs] = spline_interp3(x,y,xx,kind)
% SPLINE_INTERP3 - cubic spline interpolation
%   COEFS = SPLINE_INTERP3(X,Y) provides the piecewise polynomial form of the
%   cubic spline curve to interpolate the values Y at the break sites X.
%   If Y is a vector, then Y(j) is taken as the value to be matched at X(j),
%   hence Y must be of the same length as X.
%   If Y is a matrix, then Y(:,j) is taken as the value to be matched at X(j),
%   hence the columns of Y must equal length(X).
%   COEFS((j-1)*ndim+1:j*ndim,:) is the j-th polynomial of the spline curve.
%   coefficient [a_i,b_i,c_i,d_i]:   s_i = a_i*t_i^3+b_i*t_i^2+c_i*t_i+d_i, t_i = (x-x_i)/h(i);
%   See the exception in case 2.
%   Let n = length(x); If n==3, then the curve is degenerated into quadratic;
%   if n==2, then the curve is degenerated into linear spline;
%
%   KIND describe which kind of spline you want, There are four cases: 
%   1, Natural; 2, Clamped; 3, Not-a-knot. 4, Quadratic bound.  Kind is a
%   character: 'N' stands for case 1; 'C' stands for case 2; 'K' stands for
%   case 3; 'Q' stands for case 4.
%   Ordinarily, the not-a-knot end conditions are used. As a result,
%   the default value of argument KIND is 'K'--case 3. 
%
%   Case1, Natural spline. There is no force at the two ends. So the second
%   derivative at the 1st and last knots are 0, then b1=0,
%   b_n-2/h_n-2+2b_n-1(h_n-1+h_n-2)/h_n-1^2=3(delta_n-1-delta_n-2).
%   Y must be of the same length as X in this case.
%
%   Case2, Clamped spline. Given the tangents of the 1st and last knots, so
%   s1_x(0)=c1/h1=A, s_n-1_x(1)=(3a_n-1+2b_n-1+c_n-1)/h_n-1=B.
%   In this case, tangent elements must be added  before and after all knots in Y.
%   If Y have the same length as X, then the routine add zero tangent
%   automatically at the two ends.
%
%   Case3, Not-a-knot spline. Assuming the 1st and 2nd curve pieces are the
%   same curve, and so are the last and 2nd last ones. Then the 2nd and 2nd
%   last breaks are not knots anymore, the whole curve have continuous 3rd
%   derivative there. Y must have the same length as X in this case.
%
%   Case4, Quadratic bound spline. Assuming the 1st and last curve pieces are
%   degenerated into quadratic curve, then a_1=a_n-1=0. Y must have the
%   same length as X in this case. 
%
%   For 3-order Bezier spline curve, each piece is generated by four points
%   (P_i, A_i, B_i, P_i+1). P_i and P_i+1 are knots, while A_i and B_i are
%   auxiliary points, which can be caculated under C2 continuity assumption.
%   On the boundary, it is assumed that A_1=P_1, B_n-1=P_n. This curve
%   is a special case of clamped spline (Case 2), with boundary tangents A=B=0.
%
%   YY = SPLINE_INTERP3(X,Y,XX) provide cubic spline interpolation YY at the
%   values of XX.
%
%   See also spline_interp2, spline_interp1.

% by zpf, form BIT, 2015-7-2

assert(isvector(x),'The first argument must be a vector!');
x = x(:)';
n = length(x);
assert(n>=2,'It takes at least two points to make a spline curve!');
h = diff(x);
assert(all(h>0) || all(h<0),'The first argument must be monotonic!');
[ndim,ny]=size(y);
if ny==1,
    y=y'; ny=ndim; ndim=1;
end;
if h(1)<0,
    n1 = n:-1:1;
    x = x(n1);
    h = -h(n1(2:end));
    y = y(:,ny:-1:1);
end;
h2 = h.^2;

% switch to turn off piecewise interpolation (only output coefficients)
SW = 0;
if nargin==2,
    SW = 1;
    kind = 'k';
elseif nargin==3,
    if ischar(xx),
        SW = 1;
        kind = xx;
    else
        kind = 'k';
    end;
else
    assert(ischar(kind),'The 4th argument must be char!');
end;
kind = lower(kind);
kind = kind(1);

nd1 = ones(ndim,1);
nd0 = zeros(ndim,1);
switch kind,
    case 'n',   % Case 1, Natural spline.
        assert(ny==n,'Number of breaks do not match!');
        dy = diff(y,[],2);
        if n==2,
            a = zeros(ndim, 1);
            b = a;
        else
            delta = dy ./ h(nd1,:);
            ddelta = diff(delta,[],2);
            
            % tridiag lu factorization
            va = [0,1./h(1:n-2)];
            vb = [1,2*(h(1:n-2)+h(2:n-1))./h2(2:n-1)];
            vc = [0,h(2:n-2)./h2(3:n-1),0];
            A = [va', vb', vc'];
            b = 3*[nd0, ddelta]';
            [~,~,~,b] = tridiaglu(A,b);
            b = b';
            a = zeros(ndim, n-1);
            A = h2(1:n-2)./h2(2:n-1);
            a(:, 1:n-2) = (b(:, 2:n-1).*A(nd1,:)-b(:,1:n-2))/3;
            a(:, n-1) = -b(:, n-1)/3;
        end;
        d = y(:,1:n-1);
        c = dy-a-b;
        
    case 'c',   % Case 2, Clamped spline.
        if ny == n,
            y = [nd0 , y, nd0];
            ny = n+2;
        end;
        assert(ny==n+2,'Number of breaks do not match!');
        dy = diff(y(:,2:n+1),[],2);
        if n==2,
            c = y(:,1).*h;
            b = 3*dy-2*c-y(:,end).*h;
        else
            delta = dy ./ h(nd1,:);
            ddelta = diff(delta,[],2);
            
            % tridiag lu factorization
            va = [0, 1./h(1:n-3), 2/h(n-2)];
            vb = [2/h(1), 2*(h(1:n-3)+h(2:n-2))./h2(2:n-2), (4*h(n-2)+3*h(n-1))/h2(n-1)];
            vc = [h(1:n-2)./h2(2:n-1), 0];
            A = [va', vb', vc'];
            b = 3*[delta(:,1)-y(:,1), ddelta(:,1:n-3), 3*delta(:,n-1)-2*delta(:,n-2)-y(:,end)]';
            [~,~,~,b] = tridiaglu(A, b);
            b = b';
            c = zeros(ndim, n-1);
            A = h2(1:n-2)./h2(2:n-1);
            c(:,1:n-2) = dy(:,1:n-2)-(b(:, 2:n-1).*A(nd1, :)+2*b(:,1:n-2))/3;
            A = h(n-2)/h(n-1);
            c(:, n-1) = (b(:, n-2)/3+dy(:, n-2))/A+2*A*b(:, n-1)/3;
        end;
        d = y(:,2:n);
        a = dy-b-c;
        
    case 'k',   % Case 3, Not-a-knot spline.
        assert(ny==n,'Number of breaks do not match!');
        dy = diff(y,[],2);
        delta = dy ./ h(nd1,:);
        if n==2,
            a = zeros(ndim, 1);
            b = a;
        elseif n==3,
            a = zeros(ndim, 2);
            c = (delta(:,2)-delta(:,1))/(h(1)+h(2));
            b(:,1) = c*h2(1);
            b(:,2) = c*h2(2);
        else           
            ddelta = diff(delta,[],2);           
            % band-diagonal matrix lu factorization
            va = [0, 1./h(1:n-3), 1/h(n-2)-h2(n-1)/h(n-2)^3];
            vb = [h2(2)/(h2(1)*h(1))-1/h(1), 2*(h(1:n-3)+h(2:n-2))./h2(2:n-2), 1/h(n-2)+3/h(n-1)+2*h(n-2)/h2(n-1)];
            vc = [-(2*h(1)+3*h(2))/h2(2)-1/h(1), h(2:n-2)./h2(3:n-1), 0];
            A = [va', vb', vc'];
            b = 3*[-ddelta(:,1), ddelta]';
            [~,~,~,~,b] = bandlu(A, 1, 1, b);       %  tridiag lu factorization error if h(1)==h(2), then vb(1)=0
            b = b';
            a = zeros(ndim, n-1);
            A = h2(1:n-2)./h2(2:n-1);
            a(:, 1:n-2) = (b(:, 2:n-1).*A(nd1, :)-b(:,1:n-2))/3;
            A = (h(n-1)/h(n-2))^3;
            a(:, n-1) = a(:, n-2)*A;
        end;
        d = y(:,1:n-1);
        c = dy-a-b;
        
    case 'q',   % Case 4, Quadratic bound spline.
        assert(ny==n,'Number of breaks do not match!');
        dy = diff(y,[],2);
        if n==2,
            a = zeros(ndim, 1);
            b = a;
        else
            delta = dy ./ h(nd1,:);
            ddelta = diff(delta,[],2);
            
            %  tridiag lu factorization
            va = [0,1./h(1:n-2)];
            vb = [1, 2*(h(1:n-3)+h(2:n-2))./h2(2:n-2), (3*h(n-1)+2*h(n-2))/h2(n-1)];
            vc = [-h2(1)/h2(2), h(2:n-2)./h2(3:n-1), 0];
            A = [va', vb', vc'];
            b = 3*[nd0, ddelta]';
            [~,~,~,b] = tridiaglu(A, b);
            b = b';
            a = zeros(ndim, n-1);
            A = h2(2:n-2)./h2(3:n-1);
            a(:, 2:n-2) = (b(:, 3:n-1).*A(nd1, :)-b(:, 2:n-2))/3;
        end;
        d = y(:,1:n-1);
        c = dy-a-b;
        
    otherwise,
        error('Unexpected input!');    
end; 
    
% coefficients of splines
coefs =[a(:),b(:),c(:),d(:)];

if SW,
    yy = coefs;
else
    assert(isvector(xx),'The third argument must be a vector!');
    m = size(xx,2);
    if m==1,
        xx = xx';
    end;
    Np = length(xx);
    [~,index] = histc(xx,[-inf,x(2:n-1),inf]);    % histogram filter, if [n,bin] = histc(x,edges), n(k) = sum(bin==k).
    % now go to local coordinates
    t = (xx-x(index))./h(index);

    if ndim>1,  % replicate t and index in case Points are vector
        t = reshape(t(nd1,:),Np*ndim,1);
        index = ndim*index; temp = (-ndim:-1)';
        index = reshape(1+index(nd1,:)+temp(:,ones(1,Np)), Np*ndim,1);
    end;

    v = coefs(index,1);
    for i = 2:4,
        v = v.*t(:)+coefs(index,i);
    end;

    yy = reshape(v,ndim,Np);
    if  ndim==1 && m==1,
        yy = yy';
    end;

end;

return;



%% test
n1 = 6;
n2 = 200;
interval  = [-12 12 -1.5 1.5];
x1=linspace(-10, 10, n1);
xx1=linspace(-11, 11, n2);
y1=sin(x1);
z1=cos(x1.^3);
yz = [y1; z1];
kinds = 'NCKQ';
tname = {'Natural spline', 'Clamped spline', 'Not-a-knot spline', 'Quadratic bound spline'};

figure(1);
for ii=1:4,
    yy1 = spline_interp3(x1,yz,xx1,kinds(ii));
    subplot(2,2,ii);
    plot(x1,y1,'ro', x1,z1,'yo', xx1,yy1(1,:),'b-', xx1,yy1(2,:),'g-');
    title(tname{ii});
    axis(interval);
end;

% [yy1,coef] = spline_interp3(x1,[1,y1,1;-0.5,z1,-1],xx1,'c'); 
[yy1,coef] = spline_interp3(x1,yz,xx1,'c');      % clamp the curves with 0 tangent

%%  see draw_jointing_bezier: spline_interp3- bezier curve gif 