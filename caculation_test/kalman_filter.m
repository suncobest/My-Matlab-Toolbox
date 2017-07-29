int i ;
i=1
f=[1,1;0,1];
p0=[10,0;0,10];
w=1;
v=1;
m=10;
h=[1 0];
g=[1/2;1];
while i<=m
    m1=f*p0*transpose(f)+g*w*transpose(g);
    k1=m1*transpose(h)/(h*m1*transpose(h)+v)
    p1=m1-k1*h*m1;
    p0=p1;
    i=i+1
end