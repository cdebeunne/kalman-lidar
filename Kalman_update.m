function [X, PX] = Kalman_update(X, PX, Xold, PXold, z, Pz)

%% Process zest = h(X)
zest = h(X, Xold);

% Calcul de l'innovation
Innovation = z-zest;
Jh_Xold = Jh_xold(X, Xold);
Jh_X = Jh_x(X, Xold);

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

Rx = [1 0 0; 0 cos(xp1(4)) -sin(xp1(4)); 0 sin(xp1(4)) cos(xp1(4))];
Ry = [cos(xp1(5)) 0 sin(xp1(5)); 0 1 0; -sin(xp1(5)) 0 cos(xp1(5))];
Rz = [cos(xp1(6)) -sin(xp1(6)) 0; sin(xp1(6)) cos(xp1(6)) 0; 0 0 1];
Rot1 = Rz*Ry*Rx;

Rot2 = Rot1*Rot';
sy = sqrt(Rot2(1,1)^2+Rot2(2,1)^2);
if sy<1e-6
    yaw = atan2(Rot2(3,2),Rot2(3,3));
    pitch = atan2(-Rot2(3,1), sy);
    roll = atan2(Rot2(2,1), Rot2(1,1));
else
    yaw = atan2(-Rot2(2,3), Rot2(2,2));
    pitch = atan2(-Rot2(3,1), sy);
    roll = 0;
end
ang = [yaw;pitch;roll];
ht = [Rot*(xp1(1:3)-x(1:3)); ang];
end



function H = Jh_xold(xp1,x)
% jacobienne de h par rapport à x
dt = 0.0001;
delta_x = [0.001; zeros(14,1)];
H1 = (h(xp1,x)-h(xp1,x+delta_x))/dt;
delta_y = [0; 0.001; zeros(13,1)];
H2 = (h(xp1,x)-h(xp1,x+delta_y))/dt;
delta_z = [0;0;0.001; zeros(12,1)];
H3 = (h(xp1,x)-h(xp1,x+delta_z))/dt;
delta_yaw = [zeros(3,1); 0.005; zeros(11,1)];
H4 = (h(xp1,x)-h(xp1,x+delta_yaw))/dt;
delta_pitch = [zeros(4,1); 0.005; zeros(10,1)];
H5 = (h(xp1,x)-h(xp1,x+delta_pitch))/dt;
delta_roll = [zeros(5,1); 0.005; zeros(9,1)];
H6 = (h(xp1,x)-h(xp1,x+delta_roll))/dt;
H = [H1,H2,H3,H4,H5,H6,zeros(6,9)];
end

function H = Jh_x(xp1,x)
% jacobienne de h par rapport à xp1
dt = 0.0001;
delta_x = [0.001; zeros(14,1)];
H1 = (h(xp1,x)-h(xp1+delta_x,x))/dt;
delta_y = [0; 0.001; zeros(13,1)];
H2 = (h(xp1,x)-h(xp1+delta_y,x))/dt;
delta_z = [0;0;0.001; zeros(12,1)];
H3 = (h(xp1,x)-h(xp1+delta_z,x))/dt;
delta_yaw = [zeros(3,1); 0.005; zeros(11,1)];
H4 = (h(xp1,x)-h(xp1+delta_yaw,x))/dt;
delta_pitch = [zeros(4,1); 0.005; zeros(10,1)];
H5 = (h(xp1,x)-h(xp1+delta_pitch,x))/dt;
delta_roll = [zeros(5,1); 0.005; zeros(9,1)];
H6 = (h(xp1,x)-h(xp1+delta_roll,x))/dt;
H = [H1,H2,H3,H4,H5,H6,zeros(6,9)];
end
