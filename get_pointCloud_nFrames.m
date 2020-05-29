% Authors - Eduardo Pavez <eduardo.pavez.carvelli@gmail.com, pavezcar@usc.edu>
% Copyright Eduardo Pavez,  University of Southern California
%Los Angeles, USA, 05/30/2020
% E. Pavez, B. Girault, A. Ortega, and P. A. Chou. 
%"Region adaptive graph Fourier transform for 3D point clouds". 
%IEEE International Conference on Image Processing (ICIP), 2020
%https://arxiv.org/abs/2003.01866
function [ nFrames ] = get_pointCloud_nFrames(dataset, sequence )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
switch dataset
    
    case '8iVFBv2'
        
        switch sequence
            
            case 'redandblack'
                startFrame= 1450;
                endFrame = 1749;
            case 'soldier'
                startFrame= 536;
                endFrame = 835;
            case 'longdress'
                startFrame= 1051;
                endFrame = 1350;
            case  'loot'
                startFrame= 1000;
                endFrame = 1299;
            otherwise
                warning('the provided sequence %s does not belong in dataset %s', sequence, dataset);
        end
          
    case 'MVUB'
        
        switch sequence
            case 'andrew9'
                startFrame= 0;
                endFrame = 317;
                
            case 'david9'
                startFrame= 0;
                endFrame = 215;
            case 'phil9'
                startFrame= 0;
                endFrame = 244;
            case 'ricardo9'
                startFrame= 0;
                endFrame = 215;
            case 'sarah9'
                startFrame= 0;
                endFrame = 206;
            otherwise
                warning('the provided sequence %s does not belong in dataset %s', sequence, dataset);
                
        end
        
        
        J = 9;
       
        
        
    otherwise
        warning('%s is not proper dataset', dataset);
        
end
nFrames = endFrame - startFrame +1;
end

