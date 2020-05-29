% Authors - Eduardo Pavez <eduardo.pavez.carvelli@gmail.com, pavezcar@usc.edu>
% Copyright Eduardo Pavez,  University of Southern California
%Los Angeles, USA, 05/30/2020
% E. Pavez, B. Girault, A. Ortega, and P. A. Chou. 
%"Region adaptive graph Fourier transform for 3D point clouds". 
%IEEE International Conference on Image Processing (ICIP), 2020
%https://arxiv.org/abs/2003.01866
function [ Ahat, freqs, weights ] = Region_Adaptive_GFT( A, params )
%This function implements the Region adaptive graph fourier transform
%(RA-GFT) for point cloud attributes of voxelized point clouds
%A: attribute
%params.V pointcloud coordinates
%params.bsize = block size
%params.J, depth of octree
V           = params.V;
b           = params.bsize; %b = [b_1, b_2,...,b_L] if multilevel, or b scalar. b_1*b_2*...*b_L = 2^J, 
J           = params.J;
isMultiLevel = params.isMultiLevel;
N = size(V,1); %number of points

%% Check consistency of block sizes, resolution levels, and octree depth
if(length(b)==1)
    
    if(isMultiLevel)%basically all levels have the same block size
        
        base_bsize = log2(b);
        if(floor(base_bsize)~= base_bsize)%make sure that block size is power of 2
            error('block size bsize should be a power of 2');
        end
        L = J/base_bsize;
        
        if( L ~= floor(L))%make sure number of levels is an integer
            error('block size do not match number of levels');
        end
        bsize = ones(L,1)*b; %block size at each level is the same
        
    else
        base_bsize = log2(b);
        if(floor(base_bsize)~= base_bsize)%make sure that block size is power of 2
            error('block size bsize should be a power of 2');
        end
        L=1;
        bsize = b;
        
    end
else
    bsize =b;
    L = length(bsize);
    
    %check all entries of bsize are powers of 2
    base_bsize = log2(b);
    if(sum(base_bsize ~= floor(base_bsize)))
        error('entries of block size should be a power of 2');
    end
    %check if block sizes are consistent with octree depth
    if(sum(base_bsize)>J)
        error('block sizes do not match octree depth J');
    end
    
end
%%
Ahat = [];
Vcurr = V;
Acurr = A;
Qin = ones(N,1);
Gfreq_curr = zeros(N,1);
freqs =[];
weights=[];
for level=L:-1:1
    
    %%% block level processing
    %get block indices
    start_indices = block_indices(Vcurr,bsize(level)); %start index of blocks
    Nlevel = size(Vcurr,1);                   %number of points at curr level
    end_indices = [start_indices(2:end)-1;Nlevel];
    %get blocks with more than 1 point
    ni = end_indices - start_indices +1;
    %unchanged =  find(ni==1);%indices of blocks with single point
    to_change = find(ni ~=1); %indices of blocks that have more than 1 point
    
    Acurr_hat = Acurr;
    Qout=Qin;
    Gfreq_curr = zeros(size(Qin));
    %
    for currblock =1:length(to_change)
        
        first_point = start_indices(to_change(currblock));
        last_point  = end_indices(to_change(currblock));
        Vblock = Vcurr(first_point:last_point,:);
        Qin_block = Qin(first_point:last_point);
        Ablock =Acurr(first_point:last_point,:);
        
        [Ahatblock, Gfreq_block,weights_block] = block_coeffs(Vblock,Ablock,Qin_block,bsize(level));
        
        Acurr_hat(first_point:last_point,:) = Ahatblock;
        Qout(first_point:last_point) = weights_block;
        Gfreq_curr(first_point:last_point) = Gfreq_block;
    end
    %%% prepare for next level
    Vcurr = floor( Vcurr(start_indices,:)/bsize(level));
    Acurr = Acurr_hat(start_indices,:);
    Qin = Qout(start_indices);
    
    %%% prepare outputs, high pass coeffs, correspondig weights and graph
    %%% frequencies
    Coeff_high = Acurr_hat; %high pass coeffs
    Qout_high = Qout;       %weights of high pass coeffs
    Gfreq_high = Gfreq_curr;%graph frequency of high pass coefficients
       
     Coeff_high(start_indices,:)=[];
     Qout_high(start_indices) =[];
     Gfreq_high(start_indices)=[];
 
    Ahat = [Coeff_high ; Ahat];
    freqs = [Gfreq_high; freqs];
    weights=[Qout_high;weights];
    
    if(level==1)
        %since this is the last level, put in output DC coeffs, ans the
        %corresponding frequencies and weights Q.
        Gfreq_low = Gfreq_curr(start_indices);
        Ahat = [Acurr ; Ahat];
        freqs = [Gfreq_low ; freqs];
        weights=[Qin ; weights];
    end
    
end



end
%% Block level transform
function [Ahat, Gfreq,weights] = block_coeffs(Vblock,A,Q,bsize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%[W,~] = compute_graph_MSR(Vblock);
[W,~] = compute_graph_gaussian(Vblock);
if(sum(isnan(W),'all'))
          disp(['nana']);  
end
if (bsize == 2)
    %do standard RA-GFT with a connected graph
    [Ahat, Gfreq, weights] = RAGFT_connected_graph(W,A,Q);
else
    %check of graph is connected
    [p, ~, r, ~] = dmperm( W + eye(size(W)));
    numConnComp = size( r, 2 ) - 1;
    if (numConnComp==1)%graph is connected
        %do standard RA-GFT with a connected graph
        [Ahat, Gfreq, weights] = RAGFT_connected_graph(W,A,Q);
    else
        %if graph is disconnected,
        [Ahat, Gfreq, weights] = RAGFT_disconnected_graph(W,A,Q,Vblock,numConnComp,p,r);
    end
end
end
function [Coeff, Gfreq, weights] = RAGFT_connected_graph(W,A,Q)
[ GFT, Gfreq ] = compute_GFT( W, Q );
weights = repmat(sum(Q),size(A,1),1);
Coeff = GFT*A;
end
function [Coeff, Gfreq, weights] = RAGFT_disconnected_graph(Wcurr,A,Qcurr,Vblock,numDCs,p,r)
%Wcurr = W;
%Qcurr = Q;
%first level
U=[];
isDC=[];
Gfreq_level = [];
weights_level=[];
Vblock_new = zeros(numDCs,3);
for comp=1:numDCs
    %compute GFT
    idx=p(r(comp):r(comp+1)-1);
    [ GFT, Gfreq_tmp ] = compute_GFT( Wcurr(idx,idx), Qcurr(idx) );
    Utmp=zeros(size(Wcurr,1),length(idx));
    Utmp(idx,:)=GFT';
    U=[U,Utmp];
    isDCtmp=zeros(length(idx),1);
    isDCtmp(1)=1;
    isDC=[isDC;isDCtmp];
    Gfreq_level = [Gfreq_level; Gfreq_tmp];
    weights_level = [weights_level ; ones(length(idx),1)*sum(Qcurr(idx))];
    %compute average of points in block
    Vblock_new(comp,:) =  sum(diag(Qcurr(idx))*Vblock(idx,:),1)/sum(Qcurr(idx));
    %    
end
Ahat_1 = U'*A;
isDC_index = find(isDC);
notDC_index = find(~isDC);
Ahat_low = Ahat_1(isDC_index,:);%low pass coeffs for further processing
Ahat_high=Ahat_1(notDC_index,:);%high pass coeffs
%level 2
Qnew = weights_level(isDC_index);
Wnew = complete_graph(Vblock_new);

[GFT_new, Gfreq_new ] = compute_GFT( Wnew, Qnew );
Coeff = [GFT_new*Ahat_low;Ahat_high];
Gfreq = [Gfreq_new ;Gfreq_level(notDC_index) ];
weights = [ones(length(Qnew),1)*sum(Qnew) ; weights_level(notDC_index)];
end
