syms pusi theta fi

a=fi
b=theta
c=pusi

rx=[1,0,0;0,cos(a),-sin(a);0,sin(a),cos(a);]
ry=[cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b);]
rz=[cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1;]

% rx=[1,0,0;0,cos(-a),-sin(-a);0,sin(-a),cos(-a);]
% ry=[cos(-b),0,sin(-b);0,1,0;-sin(-b),0,cos(-b);]
% rz=[cos(-c),-sin(-c),0;sin(-c),cos(-c),0;0,0,1;]
% 
% xyz=rx*ry*rz

zyx=rz*ry*rx
zyx.'