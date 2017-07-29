syms V1 V2 V3 V43 V44

V43=0
V44=pi/2

b=V43 % attack angle
c=V44 % stroke angle

% rx=[1,0,0;0,cos(a),-sin(a);0,sin(a),cos(a);]
% ry=[cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b);]
% rz=[cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1;]

% rx=[1,0,0;0,cos(-a),-sin(-a);0,sin(-a),cos(-a);]
ry=[cos(-b),0,sin(-b);0,1,0;-sin(-b),0,cos(-b);]
rz=[cos(-c),-sin(-c),0;sin(-c),cos(-c),0;0,0,1;]

T=ry*rz*[0,-1,0;1,0,0;0,0,1;]
r=T*[V1;V2;V3;]