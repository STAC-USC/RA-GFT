% Authors - Eduardo Pavez <eduardo.pavez.carvelli@gmail.com, pavezcar@usc.edu>
% Copyright Eduardo Pavez,  University of Southern California
%Los Angeles, USA, 05/30/2020
% E. Pavez, B. Girault, A. Ortega, and P. A. Chou. 
%"Region adaptive graph Fourier transform for 3D point clouds". 
%IEEE International Conference on Image Processing (ICIP), 2020
%https://arxiv.org/abs/2003.01866
function [ V,C,J ] = get_pointCloud(dataset, sequence,  frame )
%
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
          getframe = startFrame -1 + frame;
        if(getframe>endFrame)
            warning('the frame number %d does not exist', frame);
        end
        filename = sprintf('8iVFBv2/%s/Ply/%s_vox10_%04d.ply',sequence, sequence,getframe);
        [V,C,J] = ply_read8i(filename);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        getframe = startFrame -1 + frame;
        if(getframe>endFrame)
            warning('the frame number %d does not exist', frame);
        end
        filename = sprintf('MVUB/%s/ply/frame%04d.ply', sequence,getframe);
        [V,C] = ply_readMVUB(filename);
        
        
    otherwise
        warning('%s is not proper dataset', dataset);
        
end

end

