% Authors - Eduardo Pavez <eduardo.pavez.carvelli@gmail.com, pavezcar@usc.edu>
% Copyright Eduardo Pavez,  University of Southern California
%Los Angeles, USA, 05/30/2020
% E. Pavez, B. Girault, A. Ortega, and P. A. Chou. 
%"Region adaptive graph Fourier transform for 3D point clouds". 
%IEEE International Conference on Image Processing (ICIP), 2020
%https://arxiv.org/abs/2003.01866
%% script to test RA-GFT and its inverse
%encoder RA-GFT
clear;
filename = 'longdress_vox10_1051.ply';
[V,Crgb,J] = ply_read8i(filename);
N = size(V,1);
C = RGBtoYUV(Crgb); %transform to YUV

%%
bsize=[ 2 2 2 2 2 2  2 2 2 2];
param.V=V;
param.J=J;
param.bsize = bsize;
param.isMultiLevel=1;
tic;

step = 64;
%C = ones(N,1);
[Coeff, Gfreq, weights]  = Region_Adaptive_GFT( C, param );
toc;
Y = Coeff(:,1);
Coeff_quant = round(Coeff/step)*step;
%%
tic;
[ start_indices, end_indices, V_MR, Crec ] = iRegion_Adaptive_GFT( Coeff_quant, param );
toc;

Crgb_rec = double(YUVtoRGB(Crec));

psnr_Y = -10*log10(norm(Y - Coeff_quant(:,1))^2/(N*255^2))

% ply_write('PC_original.ply',V,Crgb,[]);
% ply_write('PC_coded.ply',V,Crgb_rec,[]);

