function ht = f(xp1,x)
Rx = [1 0 0; 0 cos(x(4)) -sin(x(4)); 0 sin(x(4)) cos(x(4))];
Ry = [cos(x(5)) 0 sin(x(5)); 0 1 0; -sin(x(5)) 0 cos(x(5))];
Rz = [cos(x(6)) -sin(x(6)) 0; sin(x(6)) cos(x(6)) 0; 0 0 1];
Rot = Rz*Ry*Rx;

ht = [Rot'*(xp1(1:3)-x(1:3)); xp1(4:6)-x(4:6)];
end