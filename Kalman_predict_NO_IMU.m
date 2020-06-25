function [X, PX] = Kalman_predict_NO_IMU(X, PX, IMU, dt)
RImu_lidar = [1 0 0; 0 -1 0; 0 0 -1];
Rx = [1 0 0; 0 cos(X(4)) -sin(X(4)); 0 sin(X(4)) cos(X(4))];
Ry = [cos(X(5)) 0 sin(X(5)); 0 1 0; -sin(X(5)) 0 cos(X(5))];
Rz = [cos(X(6)) -sin(X(6)) 0; sin(X(6)) cos(X(6)) 0; 0 0 1];
Rot = Rz*Ry*Rx;
X(1:3) = X(1:3) + Rot*X(7:9)*dt;
X(4:6) = X(4:6) + RImu_lidar*IMU*dt;

% to do : calculer la covariance et la jacobienne
end
