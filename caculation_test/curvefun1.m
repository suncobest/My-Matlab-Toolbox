function f = curvefun1(x,tdata)
assert(isvector(x) && length(x)==2 && isvector(tdata), 'Unexpected input!');
f=x(1)*x(2)*(0.064+0.008*exp(-tdata/0.008/x(2))-0.072*exp(-tdata/0.072/x(2))); %其中x1=E1,x2=theta1
return;

%% test
tdata=[0,0.0214,0.0549,0.0763,0.0983,0.1203,0.1422,0.1636,0.1856,0.2075,0.2294,0.2509,0.2728,0.3];
cdata=[ 0,0.0007,0.0017,0.0018,0.0021,0.0022,0.0023,0.0024,0.0025,0.0030,0.0032,0.0036,0.00399,0.004];
x0=[0.2,0.05];
x=lsqcurvefit('curvefun1',x0,tdata,cdata)
f=curvefun1(x,tdata);