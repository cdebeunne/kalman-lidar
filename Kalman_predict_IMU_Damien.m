function [X, PX] = Kalman_predict_IMU_Damien(X, PX, imuACC, imuGYR, Q, dt)


RImu_lidar = [1 0 0; 0 -1 0; 0 0 -1];
u = [RImu_lidar*imuACC'; RImu_lidar*imuGYR']; %[ax ay az gx gy gz] dans le repère LIDAR

% X = [x y z theta phi psi, vx vy vz,  bgx bgy bgz bax, bay baz]

Rx = [1 0 0; 0 cos(X(4)) -sin(X(4)); 0 sin(X(4)) cos(X(4))];
Ry = [cos(X(5)) 0 sin(X(5)); 0 1 0; -sin(X(5)) 0 cos(X(5))];
Rz = [cos(X(6)) -sin(X(6)) 0; sin(X(6)) cos(X(6)) 0; 0 0 1];
Rot = Rz*Ry*Rx; % LIDAR 2 World

pos = 1:3; ang = 4:6; vit = 7:9; bg = 10:12; ba = 13:15;
ua = 1:3; ug = 4:6;


% IMU in world coordinates
uw = u;
uw(ug) = (u(ug)-X(bg));
uw(ua) = Rot*(u(ua)-X(ba));



X(pos) = X(pos) + X(vit)*dt + uw(ua)*dt^2/2;
X(ang) = X(ang) + uw(ug)*dt;
X(vit) = X(vit) + uw(ua)*dt;
X(bg) = X(bg);
X(ba) = X(ba);

Jf_X = eye(15);
Jf_X(pos,vit) = (eye(3)*dt);
Jf_X(pos,ang(1)) = Rz*Ry*[1 0 0; 0 -sin(X(4)) -cos(X(4)); 0 cos(X(4)) -sin(X(4))]*(u(ua)-X(ba))*dt^2/2;
Jf_X(pos,ang(2)) = Rz*[-sin(X(5)) 0 cos(X(5)); 0 1 0; -cos(X(5)) 0 -sin(X(5))]*Rx*(u(ua)-X(ba))*dt^2/2;
Jf_X(pos,ang(3)) = [-sin(X(6)) -cos(X(6)) 0; cos(X(6)) -sin(X(6)) 0; 0 0 1]*Ry*Rx*(u(ua)-X(ba))*dt^2/2;

Jf_X(pos,ba) = -Rot.*(eye(3)*dt^2/2);
Jf_X(ang,bg) = -(eye(3)*dt);
Jf_X(vit,ba) = -Rot.*(eye(3)*dt);


Jf_U = eye(15, 6);
Jf_U(pos, ua) = Rot.*(eye(3)*dt^2);
Jf_U(ang, ug) = eye(3)*dt;
Jf_U(vit, ua) = eye(3)*dt;

PX = Jf_X*PX*Jf_X' + Jf_U*Q*Jf_U';

end

