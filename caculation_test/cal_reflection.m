syms cx cy cz tx ty tz theta real
rx = [1 0 0 0;0 cz/sqrt(cy^2+cz^2) -cy/sqrt(cy^2+cz^2) 0;0 cy/sqrt(cy^2+cz^2) cz/sqrt(cy^2+cz^2) 0;0 0 0 1]
ry=[sqrt(cy^2+cz^2) 0 -cx 0; 0 1 0 0; cx 0 sqrt(cy^2+cz^2) 0;0 0 0 1]
rtheta = [cos(theta) -sin(theta) 0 0;sin(theta) cos(theta) 0 0;0 0 1 0;0 0 0 1]
T = [1 0 0 -tx;0 1 0 -ty;0 0 1 -tz;0 0 0 1]
Tre = [1 0 0 tx;0 1 0 ty;0 0 1 tz;0 0 0 1]

rotation = simplify(rx'*ry'*rtheta*ry*rx);
rotation = simplify(subs(rotation,cx^2,1-cy^2-cz^2));
rotation = simplify(subs(rotation,cy^2+cz^2,1-cx^2))

rotation = simplify(Tre*rotation*T)

reflect=simplify(rx'*ry'*[1 0 0 0;0 1 0 0;0 0 -1 0;0 0 0 1]*ry*rx);
reflect = simplify(subs(reflect,cx^2,1-cy^2-cz^2));
reflect = simplify(subs(reflect,cy^2+cz^2,1-cx^2))

reflect = simplify(Tre*reflect*T)