%load KITTI_VEL_SCAN2.mat

% point cloud analysis parameters
c_edge = 0.2;
c_plane = 0.05;
distThresholdEdge = 0.6;
minClusterSizeEdge = 5;
barycenterThresholdEdge = 1.5;
distThresholdPlane = 1;
minClusterSizePlane = 30;
barycenterThresholdPlane = 10;

% estimator parameters

% estimation of x, y and psi
x = [0,0,0,0,0,0];
xWorld = [0;0;0];
posList = [xWorld];
theta = 0;
phi = 0;
psi = 0;
tpp = [theta, phi, psi];
R = eye(3,3);

for k=1:600
    disp(k);
    
    % in case of empty clouds
    if size(filtered_traj{k},1) == 0 || size(filtered_traj{k+1},1)==0
        theta = theta + x(4);
        phi = phi + x(5);
        psi = psi + x(6);
        R = eul2rotm([theta, phi, psi], 'XYZ');
        dxWorld = R*x(1:3)';
        xWorld = xWorld + dxWorld;
        posList = [posList, xWorld];
        tpp = [tpp; [theta, phi, psi]];
        continue
    end
    
    % first filtering
    filteredCloud_1 = cloudFilter(filtered_traj{k},"VLP16");
    [edgeIdx_1, planeIdx_1, labelCloud_1, smoothnessCloud_1] =...
        edgePlaneDetector(filteredCloud_1.Location, c_edge, c_plane);
    
    filteredCloud_2 = cloudFilter(filtered_traj{k+1},"VLP16");
    [edgeIdx_2, planeIdx_2, labelCloud_2, smoothnessCloud_2] =...
        edgePlaneDetector(filteredCloud_2.Location, c_edge, c_plane);
    
    size1 = size(filteredCloud_1.Location, 1);
    size2 = size(filteredCloud_1.Location, 2);
    
    
    
    %----------------------------------------------------------------------
    % evaluate the corespondence
    %----------------------------------------------------------------------
    
    
    
    
    % creating the edgeClouds
    
    edgeCloud_1 = select(filteredCloud_1, ~edgeIdx_1, 'OutputSize', 'full');
    edgeCloud_2 = select(filteredCloud_2, ~edgeIdx_2, 'OutputSize', 'full');
    
    % clustering the edge clouds
    
    [edgePoints_1, barycenterEdge_1, directionsEdge_1, eigenEdge_1]...
        = clusteringEdge(edgeCloud_1, distThresholdEdge, minClusterSizeEdge);
    [edgePoints_2, barycenterEdge_2, directionsEdge_2, eigenEdge_2]...
        = clusteringEdge(edgeCloud_2, distThresholdEdge, minClusterSizeEdge);
    
    % match the subclouds
    
    corespondencesEdge = matchingEdge(edgePoints_1, edgePoints_2,...
        barycenterEdge_1, barycenterEdge_2, eigenEdge_1, eigenEdge_2, barycenterThresholdEdge);
    
    
    % creating the planeCloud
    
    planeCloud_1 = select(filteredCloud_1, ~planeIdx_1, 'OutputSize', 'full');
    planeCloud_2 = select(filteredCloud_2, ~planeIdx_2, 'OutputSize', 'full');
    
    % clustering the plane clouds
    
    [planePoints_1, barycenterPlane_1, normalsPlane_1,...
        normalsStd_1, normalsList_1, labelsPlane_1, validLabels_1]...
        = clusteringPlane(planeCloud_1, distThresholdPlane, minClusterSizePlane);
    [planePoints_2, barycenterPlane_2, normalsPlane_2,...
        normalsStd_2, normalsList_2, labelsPlane_2, validLabels_2]...
        = clusteringPlane(planeCloud_2, distThresholdPlane, minClusterSizePlane);
    
    
    % match the plane clouds
    
    corespondencesPlane = matchingPlane(planePoints_1, planePoints_2,...
        normalsPlane_1, normalsPlane_2, barycenterPlane_1, barycenterPlane_2, barycenterThresholdPlane);
    
    
    
    %--------------------------------------------------------------------------
    % finding the correct rigid transform with Levenberg and Marquardt algorithm
    %--------------------------------------------------------------------------
    
    % finding dx, dy and dpsi with the edges
    
    x0 = [0, 0, 0];
    f = @(x)costEdge(corespondencesEdge, barycenterEdge_1, barycenterEdge_2, x);
    
    % remove outliers
    firstEval = f(x0);
    inliers = ~isoutlier(firstEval);
    inliers = logical(inliers(:,1).*inliers(:,2));
    corespondencesEdge = corespondencesEdge(inliers,:);

    
    y0 = [0,0,0,0,0,0];
    f = @(x)costPlane(corespondencesPlane, normalsPlane_1, normalsPlane_2, barycenterPlane_1, barycenterPlane_2, x);
    
    %remove outliers
    firstEval = f(y0);
    inliers = ~isoutlier(firstEval);
    inliers = logical(inliers(:,1).*inliers(:,2));
    corespondencesPlane = corespondencesPlane(inliers,:);
    
    % optimisation
    
    x0 = [0,0,0,0,0,0];
    lb = [-1.5, -0.02, -0.02, -0.01, -0.01, -pi/6];
    ub = [1.5, 0.02, 0.02, 0.01, 0.01, pi/6];
    f = @(x)globalCost_orth(corespondencesEdge, corespondencesPlane,...
    edgePoints_1, directionsEdge_1, barycenterEdge_2, directionsEdge_2,...
    normalsPlane_1, normalsPlane_2, x);
    
    try
        options = optimoptions('lsqnonlin','FunctionTolerance', 0.001, 'MaxFunctionEvaluations', 1000);
        [x, ~] = lsqnonlin(f,x,lb,ub,options);
    catch
        warning('optimisation failure')
    end    
    
    %----------------------------------------------------------------------
    % adding the new pose in world coordinates
    %----------------------------------------------------------------------
    
    % x,y and psi
    theta = theta + x(4);
    phi = phi + x(5);
    psi = psi + x(6);
    R = eul2rotm([theta, phi, psi], 'XYZ');
    dxWorld = R*x(1:3)';
    xWorld = xWorld + dxWorld;
    disp(xWorld);
    posList = [posList, xWorld];
    tpp = [tpp; [theta, phi, psi]];
end

% load KITTI_OSTX.mat
pos = groundtruth(filtered_posEnu, -3.2*pi/10);

% display the results

figure(1);
plot(posList(1,:), posList(2,:));
hold on;
axis equal;
plot(pos(:,1),pos(:,2));
legend('Edge & plane odometry', 'Groundtruth');
title('Position comparison');

figure(2)
plot(tpp(:,1));
hold on;
plot(att(:,1));
legend('Edge & plane odometry', 'Groundtruth');
title('Attitude comparison');

% save results

save('results.mat', 'tpp', 'posList');