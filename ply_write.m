%MSR file
function ply_write(filename,V,C,F)
% Write the mesh to a ply file
% Usage: ply_write(filename,V,F,C)
% Input:
%   filename: name of ply file in current directory
%   V: Nx3 matrix of 3D coordinates of total N vertices
%   C: Nx3 matrix of R,G,B colors on the N vertices
%   F: Mx3 matrix of vertex indices of M triangles
%

N = size(V,1);
M = size(F,1);

fid = fopen(filename,'w','n','UTF-8');
fprintf(fid,'ply\n');
fprintf(fid,'format ascii 1.0\n');
fprintf(fid,'element vertex %d\n',N);
fprintf(fid,'property float x\n');
fprintf(fid,'property float y\n');
fprintf(fid,'property float z\n');
fprintf(fid,'property uchar red\n');
fprintf(fid,'property uchar green\n');
fprintf(fid,'property uchar blue\n');
if M > 0
    fprintf(fid,'element face %d\n',M);
    fprintf(fid,'property list uchar int vertex_index\n');
end
fprintf(fid,'end_header\n');
fprintf(fid,'%g %g %g %u %u %u\n',[V,double(C)]');
%fprintf(fid,'%f %f %f\n',V');
if M > 0
    fprintf(fid,'3 %u %u %u\n',F');
end
fclose(fid);