
function [indices]=block_indices(V,bsize)
%V is Nx3, a point cloud, each row are the xyz coordinates of a point, 
% each coordinate x,y,z is an integer
% this assumes point cloud is morton ordered


    V_coarse=floor(V/bsize)*bsize;

    variation=sum(abs(V_coarse(2:end,:)-V_coarse(1:end-1,:)),2);

    variation=[1;variation];

    indices=find(variation);

  


end

