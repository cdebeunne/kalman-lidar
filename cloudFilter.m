function filteredCloud = cloudFilter(pc, type)
% Format the cloud and remove the ground 

% reshape point cloud

if type == "VLP16"
    filteredCloud = pointcloudMatrixVLP16(pc);
elseif type == "HDL64"
    filteredCloud = pointcloudMatrixHDL64(pc);
end

end

