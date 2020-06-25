function [z, Pz, fail] = get_relative_transform_LIDAR(filtered_traj1, filtered_traj2, detector_params)

% if size(filtered_traj1,1) == 0 || size(filtered_traj2,1)==0
%     fail = 1;
%     return;
% end
    

% first filtering
filteredCloud_1 = cloudFilter(filtered_traj1,"VLP16");
[edgeIdx_1, planeIdx_1, ~, ~] =...
    edgePlaneDetector(filteredCloud_1.Location, detector_params.c_edge, detector_params.c_plane);

filteredCloud_2 = cloudFilter(filtered_traj2,"VLP16");
[edgeIdx_2, planeIdx_2, ~, ~] =...
    edgePlaneDetector(filteredCloud_2.Location, detector_params.c_edge, detector_params.c_plane);



%----------------------------------------------------------------------
% evaluate the corespondence
%----------------------------------------------------------------------

% creating the edgeClouds

edgeCloud_1 = select(filteredCloud_1, ~edgeIdx_1, 'OutputSize', 'full');
edgeCloud_2 = select(filteredCloud_2, ~edgeIdx_2, 'OutputSize', 'full');

% clustering the edge clouds

[edgePoints_1, barycenterEdge_1, directionsEdge_1, eigenEdge_1, ~, ~]...
    = clusteringEdge(edgeCloud_1, detector_params.distThresholdEdge, detector_params.minClusterSizeEdge);
[edgePoints_2, barycenterEdge_2, directionsEdge_2, eigenEdge_2, ~, ~]...
    = clusteringEdge(edgeCloud_2, detector_params.distThresholdEdge, detector_params.minClusterSizeEdge);

% match the subclouds

corespondencesEdge = matchingEdge(edgePoints_1, edgePoints_2,...
    barycenterEdge_1, barycenterEdge_2, eigenEdge_1, eigenEdge_2, detector_params.barycenterThresholdEdge);


% creating the planeCloud

planeCloud_1 = select(filteredCloud_1, ~planeIdx_1, 'OutputSize', 'full');
planeCloud_2 = select(filteredCloud_2, ~planeIdx_2, 'OutputSize', 'full');

% clustering the plane clouds

[planePoints_1, barycenterPlane_1, normalsPlane_1]...
    = clusteringPlane(planeCloud_1, detector_params.distThresholdPlane, detector_params.minClusterSizePlane);
[planePoints_2, barycenterPlane_2, normalsPlane_2]...
    = clusteringPlane(planeCloud_2, detector_params.distThresholdPlane, detector_params.minClusterSizePlane);


% match the plane clouds

corespondencesPlane = matchingPlane(planePoints_1, planePoints_2,...
    normalsPlane_1, normalsPlane_2, barycenterPlane_1, barycenterPlane_2, detector_params.barycenterThresholdPlane);


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
try
    corespondencesEdge = corespondencesEdge(inliers,:);
catch
    warning('not enough matches');
end

y0 = [0,0,0,0,0,0];
f = @(x)costPlane(corespondencesPlane, normalsPlane_1, normalsPlane_2, barycenterPlane_1, barycenterPlane_2, x);

%remove outliers
firstEval = f(y0);
inliers = ~isoutlier(firstEval);
inliers = logical(inliers(:,1).*inliers(:,2));
corespondencesPlane = corespondencesPlane(inliers,:);

% optimisation
x0 = [0,0,0,0,0,0];
lb = [-1.5, -0.05, -0.02, -0.01, -0.01, -pi/6];
ub = [1.5, 0.05, 0.02, 0.01, 0.01, pi/6];
f = @(x)globalCost_orth(corespondencesEdge, corespondencesPlane,...
    edgePoints_1, directionsEdge_1,...
    barycenterEdge_2, directionsEdge_2,...
    normalsPlane_1, normalsPlane_2, x);
try
    options = optimoptions('lsqnonlin','FunctionTolerance', 0.001, 'MaxFunctionEvaluations', 1000);
    [x,resnorm,~,~,~,~,jacobian] = lsqnonlin(f,x0,lb,ub,options);
    z = x';
    Pz = resnorm\(jacobian'*jacobian);
    fail = 0;
catch
    warning('optimisation failure')
    z = zeros(6,1);
    Pz = 10000*eye(6,6);
    fail = 1;
end
end