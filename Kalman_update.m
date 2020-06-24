function [X, PX] = Kalman_update(X, PX, Xold, PXold, z, Pz)

%% Process zest = h(X)
zest = h(X, Xold);

% Calcul de l'innovation
Innovation = z-zest;
Jh_Xold = Jh_xold(X, Xold);
Jh_X = Jh_x(X);

% Covariance de l'innovation
S = Jh_X*PX*Jh_X'+Jh_Xold*PXold*Jh_Xold'+Pz;
% Gain de Kalman
K = PX*Jh_X'/S;
X = X + K*Innovation;
PX = PX - K*Jh_X*PX;

end






function ht = h(xp1,x)
Rx = [1 0 0; 0 cos(x(4)) -sin(x(4)); 0 sin(x(4)) cos(x(4))];
Ry = [cos(x(5)) 0 sin(x(5)); 0 1 0; -sin(x(5)) 0 cos(x(5))];
Rz = [cos(x(6)) -sin(x(6)) 0; sin(x(6)) cos(x(6)) 0; 0 0 1];
Rot = Rz*Ry*Rx;

ht = [Rot'*(xp1(1:3)-x(1:3)); xp1(4:6)-x(4:6)];
end



function H = Jh_x(x)
% jacobienne de h par rapport à x
H = [                                             cos(conj(x(5)))*cos(conj(x(6))),                                             cos(conj(x(5)))*sin(conj(x(6))),              -sin(conj(x(5))), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    cos(conj(x(6)))*sin(conj(x(4)))*sin(conj(x(5))) - cos(conj(x(4)))*sin(conj(x(6))), cos(conj(x(4)))*cos(conj(x(6))) + sin(conj(x(4)))*sin(conj(x(5)))*sin(conj(x(6))), cos(conj(x(5)))*sin(conj(x(4))), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    sin(conj(x(4)))*sin(conj(x(6))) + cos(conj(x(4)))*cos(conj(x(6)))*sin(conj(x(5))), cos(conj(x(4)))*sin(conj(x(5)))*sin(conj(x(6))) - cos(conj(x(6)))*sin(conj(x(4))), cos(conj(x(4)))*cos(conj(x(5))), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0,                                                                       0,                           0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0,                                                                       0,                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0,                                                                       0,                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
end

function H = Jh_xold(xp1,x)
% jacobienne de h par rapport à xp1
H = [                                              -cos(conj(x(5)))*cos(conj(x(6))),                                              -cos(conj(x(5)))*sin(conj(x(6))),                sin(conj(x(5))),                                                                                                                                                                                                                       0,                                           cos(conj(x(5)))*(x(3) - xp1(3)) + cos(conj(x(6)))*sin(conj(x(5)))*(x(1) - xp1(1)) + sin(conj(x(5)))*sin(conj(x(6)))*(x(2) - xp1(2)),                                                                                               cos(conj(x(5)))*sin(conj(x(6)))*(x(1) - xp1(1)) - cos(conj(x(5)))*cos(conj(x(6)))*(x(2) - xp1(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
    cos(conj(x(4)))*sin(conj(x(6))) - cos(conj(x(6)))*sin(conj(x(4)))*sin(conj(x(5))), - cos(conj(x(4)))*cos(conj(x(6))) - sin(conj(x(4)))*sin(conj(x(5)))*sin(conj(x(6))), -cos(conj(x(5)))*sin(conj(x(4))), (cos(conj(x(6)))*sin(conj(x(4))) - cos(conj(x(4)))*sin(conj(x(5)))*sin(conj(x(6))))*(x(2) - xp1(2)) - (sin(conj(x(4)))*sin(conj(x(6))) + cos(conj(x(4)))*cos(conj(x(6)))*sin(conj(x(5))))*(x(1) - xp1(1)) - cos(conj(x(4)))*cos(conj(x(5)))*(x(3) - xp1(3)), sin(conj(x(4)))*sin(conj(x(5)))*(x(3) - xp1(3)) - cos(conj(x(5)))*cos(conj(x(6)))*sin(conj(x(4)))*(x(1) - xp1(1)) - cos(conj(x(5)))*sin(conj(x(4)))*sin(conj(x(6)))*(x(2) - xp1(2)),   (cos(conj(x(4)))*cos(conj(x(6))) + sin(conj(x(4)))*sin(conj(x(5)))*sin(conj(x(6))))*(x(1) - xp1(1)) + (cos(conj(x(4)))*sin(conj(x(6))) - cos(conj(x(6)))*sin(conj(x(4)))*sin(conj(x(5))))*(x(2) - xp1(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
    - sin(conj(x(4)))*sin(conj(x(6))) - cos(conj(x(4)))*cos(conj(x(6)))*sin(conj(x(5))),   cos(conj(x(6)))*sin(conj(x(4))) - cos(conj(x(4)))*sin(conj(x(5)))*sin(conj(x(6))), -cos(conj(x(4)))*cos(conj(x(5))), (cos(conj(x(4)))*cos(conj(x(6))) + sin(conj(x(4)))*sin(conj(x(5)))*sin(conj(x(6))))*(x(2) - xp1(2)) - (cos(conj(x(4)))*sin(conj(x(6))) - cos(conj(x(6)))*sin(conj(x(4)))*sin(conj(x(5))))*(x(1) - xp1(1)) + cos(conj(x(5)))*sin(conj(x(4)))*(x(3) - xp1(3)), cos(conj(x(4)))*sin(conj(x(5)))*(x(3) - xp1(3)) - cos(conj(x(4)))*cos(conj(x(5)))*cos(conj(x(6)))*(x(1) - xp1(1)) - cos(conj(x(4)))*cos(conj(x(5)))*sin(conj(x(6)))*(x(2) - xp1(2)), - (cos(conj(x(6)))*sin(conj(x(4))) - cos(conj(x(4)))*sin(conj(x(5)))*sin(conj(x(6))))*(x(1) - xp1(1)) - (sin(conj(x(4)))*sin(conj(x(6))) + cos(conj(x(4)))*cos(conj(x(6)))*sin(conj(x(5))))*(x(2) - xp1(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0,                                                                         0,                            0,                                                                                                                                                                                                                      -1,                                                                                                                                                       0,                                                                                                                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0,                                                                         0,                            0,                                                                                                                                                                                                                       0,                                                                                                                                                      -1,                                                                                                                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0,                                                                         0,                            0,                                                                                                                                                                                                                       0,                                                                                                                                                       0,                                                                                                                                                                              -1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
end
