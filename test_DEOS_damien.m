% load FILT_IMU_DAT.mat;
% load FILT_VEL_SCAN_DAT.mat;

imuTime = filtered_ROSTime;
imuACC = filtered_ACC;
imuGYR = filtered_GYR;
time = filtered_scantime;

% point cloud analysis parameters
detector_params.c_edge = 0.2;
detector_params.c_plane = 0.05;
detector_params.distThresholdEdge = 0.6;
detector_params.minClusterSizeEdge = 5;
detector_params.barycenterThresholdEdge = 1.5;
detector_params.distThresholdPlane = 1;
detector_params.minClusterSizePlane = 30;
detector_params.barycenterThresholdPlane = 10;


%% KALMAN intialisation
X = zeros(15,1);
X(15) = 9.81; % gravitation biais
PX = [zeros(9,15);
    zeros(6,9),0.001*eye(6,6)];
Xold = X;
PXold = PX;
Q = [1000*eye(3,3) zeros(3,3);
    zeros(3,3), 5e-4*eye(3,3)];
fail = 0;
err = [];
posList = [0;0;0];


%% Time index for data reading in order
idx.imu = 2;
idx.lidar = 2; % 2 parcequ'on peut rien faire du premier tout seul
while idx.imu < length(imuTime) && idx.lidar < length(time)
    
    % If next measure is IMU
    if time(idx.lidar) > imuTime(idx.imu)
        %% Predict pose with IMU Kalman until IMU time
        dt = imuTime(idx.imu)-imuTime(idx.imu-1);
        if dt > 0.1
            dt = 0.01;
        end
        [X, PX] = Kalman_predict_IMU(X, PX, imuACC(idx.imu,:), imuGYR(idx.imu,:), Q, dt);
        
        %% Set pointer to next IMU data
        idx.imu = idx.imu + 1;
        
        % If next measure is LIDAR
    else
        %% Predict pose with state until LIDAR time
        dt = time(idx.lidar)-imuTime(idx.imu);
        [X, PX] = Kalman_predict_NO_IMU(X, PX, dt); % Ou on suppose que la
        %pose est Ok à ce stade, sinon pour faire ca propre il faut les
        %vitesses dans le vecteur d'état pour s'auto predire
        
        
        %% Process LIDAR scan : get transform + uncertainty of the transform
        [z, Pz, fail] = get_relative_transform_LIDAR(filtered_traj{idx.lidar-1}, filtered_traj{idx.lidar}, detector_params);
        if fail
            continue % un des scan est vide ou opti a foiré, on prend la prochaine mesure
        end
        
        %% Update Kalman
        Pz = [0.01*eye(3,3), zeros(3,3);
            zeros(3,3), 0.05*eye(3,3)];
        [X, PX] = Kalman_update(X, PX, Xold, PXold, z, Pz);
        Xold = X;
        PXold = PX;
        disp(X');
        
        %% Set pointer to next LIDAR data
        idx.lidar = idx.lidar + 1;
    end
    
    %% Save pose for display
    posList = [posList, X(1:3)];
end

% load GNSS_DAT.mat

pos = groundtruth(filtered_posEnu, -3.2*pi/10);

% display the results

figure(1);
plot(posList(1,:), posList(2,:));
hold on;
plot(pos(:,1),pos(:,2));
axis equal;
legend('Edge & plane odometry', 'Groundtruth');
title('Position comparison');

% save results

save('results.mat', 'tpp', 'posList');