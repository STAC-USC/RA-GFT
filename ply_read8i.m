function [V,C,J] = ply_read8i(filename)
% Read vertices and their colors from a ply file
% Usage: [V,C] = ply_write(filename)
% Input:
%   filename: name of ply file in current directory
% Output:
%   V: Nx3 matrix of x,y,z coordinates of total N vertices
%   C: Nx3 matrix of R,G,B colors on the N vertices
%   J: voxel depth

fid = fopen(filename,'r','n','UTF-8');
fscanf(fid,'ply\n');
fscanf(fid,'format ascii 1.0\n');
%N = fscanf(fid,'element vertex %d\n',1);
fscanf(fid,'comment Version 2, Copyright 2017, 8i Labs, Inc.\n');
fscanf(fid,'comment frame_to_world_scale %g\n',1);
fscanf(fid,'comment frame_to_world_translation %g %g %g\n',3);
w=fscanf(fid,'comment width %d\n',1);
N=fscanf(fid,'element vertex %d\n',1);
%
fscanf(fid,'property float x\n');
fscanf(fid,'property float y\n');
fscanf(fid,'property float z\n');
fscanf(fid,'property uchar red\n');
fscanf(fid,'property uchar green\n');
fscanf(fid,'property uchar blue\n');
s = fgetl(fid);
% if strcmp(s(1:7),'element')
%     M = sscanf(s,'element face %d\n',1);
%     fscanf(fid,'property list uchar int vertex_index\n');
% else % strcmp(s(1:10),'end_header')
%     M = 0;
% end
A = fscanf(fid,'%g %g %g %u %u %u\n',[6,N]);
V = A(1:3,:)';
C = A(4:6,:)';
J=log2(w+1);
% if M > 0
%     F = fscanf(fid,'3 %u %u %u\n',[3,M])';
% else
%     F = [];
% end
fclose(fid);