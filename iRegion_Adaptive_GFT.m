% Authors - Eduardo Pavez <eduardo.pavez.carvelli@gmail.com, pavezcar@usc.edu>
% Copyright Eduardo Pavez,  University of Southern California
%Los Angeles, USA, 05/30/2020
% E. Pavez, B. Girault, A. Ortega, and P. A. Chou. 
%"Region adaptive graph Fourier transform for 3D point clouds". 
%IEEE International Conference on Image Processing (ICIP), 2020
%https://arxiv.org/abs/2003.01866
function [ starti, endi, V_MR, Coeff_aux ] = iRegion_Adaptive_GFT( Coeff, params )

V           = params.V;
b           = params.bsize;
J           = params.J;
isMultiLevel = params.isMultiLevel;
N = size(V,1);
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

starti = cell(L,1);  
Q      = cell(L,1);
endi =  cell(L,1);
V_MR = cell(L,1);
Vcurr = V;
V_MR{L} = Vcurr;
Q{L}=ones(N,1);
for level = L : -1 :1
    
    start_indices = block_indices(Vcurr,bsize(level)); %start index of blocks
    Nlevel = size(Vcurr,1);                   %number of points at curr level
    end_indices = [start_indices(2:end)-1;Nlevel];
    Vcurr = floor( Vcurr(start_indices,:)/bsize(level));
    starti{level} = start_indices;
    endi{level}   = end_indices;
    if(level >1)
        %compute pointcloud at lower resolution
        V_MR{level -1} = Vcurr;
        %collect Q matrices
        Q{level-1} = compute_Q_lower_resolution(Q{level}, start_indices, end_indices);
    end
    
end



Coeff_aux = Coeff;


for level =1:L
    start_indices = starti{level};
    end_indices   =endi{level};
    level_indices = start_indices(1):end_indices(end);
    Nlevel = length(level_indices);
    
    %get blocks with more than 1 point
    ni = end_indices - start_indices +1;
    %unchanged =  find(ni==1);%indices of blocks with single point
    to_change = find(ni ~=1); %indices of blocks that have more than 1 point
    
    coeff_level = Coeff_aux(level_indices,:);
    coeff_level_permuted = invert_level_permutation(start_indices,  coeff_level);
    
    % iterate over all blocks of level 
    Vcurr = V_MR{level};
    Qin = Q{level};
    for currblock =1:length(to_change)
        
        first_point = start_indices(to_change(currblock));
        last_point  = end_indices(to_change(currblock));
        Vblock = Vcurr(first_point:last_point,:);
        Qin_block = Qin(first_point:last_point);
        Ahatblock =coeff_level_permuted(first_point:last_point,:);
        
        Ablock = invert_block_coeffs(Vblock,Ahatblock,Qin_block,bsize(level));
        coeff_level_permuted(first_point:last_point,:)=Ablock;
%         Acurr_hat(first_point:last_point,:) = Ahatblock;
%         Qout(first_point:last_point) = weights_block;
%         Gfreq_curr(first_point:last_point) = Gfreq_block;
    end
    Coeff_aux(level_indices,:) = coeff_level_permuted;
end



end
function [coeff_level_permuted] = invert_level_permutation(start_indices,  coeff_level)
%
coeff_level_permuted = zeros(size(coeff_level));

ndc = length(start_indices);
level_indices = 1:size(coeff_level,1);
level_indices_high = level_indices;
level_indices_high(start_indices)=[];
coeff_level_permuted(start_indices,:) = coeff_level(1:ndc,:);
coeff_level_permuted(level_indices_high,:) = coeff_level(ndc+1:end,:);

end
function Qout = compute_Q_lower_resolution(Qin, start_indices, end_indices)

Nl=length(start_indices);
Qout = zeros(Nl,1);

for i=1:Nl
    starti = start_indices(i);
    endi   = end_indices(i);
    Qout(i) = sum(Qin(starti:endi));
end
end
function [A] = invert_block_coeffs(Vblock,Ahat,Q,bsize)

[W,~] = compute_graph_MSR(Vblock);
%[W,~] = compute_graph_gaussian(Vblock);

if (bsize == 2)
    %do standard RA-GFT with a connected graph
    A = inverse_RAGFT_connected_graph(W,Ahat,Q);
else
    %check of graph is connected
    [p, ~, r, ~] = dmperm( W + eye(size(W)));
    numConnComp = size( r, 2 ) - 1;
    if (numConnComp==1)%graph is connected
        %do standard RA-GFT with a connected graph
        A = inverse_RAGFT_connected_graph(W,Ahat,Q);
    else
        %if graph is disconnected,
        A = inverse_RAGFT_disconnected_graph(W,Ahat,Q,Vblock,numConnComp,p,r);
    end
end

end
function [A] = inverse_RAGFT_connected_graph(W,Ahat,Q)
[ GFT, ~ ] = compute_GFT( W, Q );
%weights = repmat(sum(Q),size(Ahat,1),1);
A = GFT'*Ahat;
end
function [A] = inverse_RAGFT_disconnected_graph(Wcurr,Ahat,Qcurr,Vblock,numDCs,p,r)

%first level, compute the transform for all connected comopnents
U=[];
isDC=[];
Gfreq_level = [];
weights_level=[];
Vblock_new = zeros(numDCs,3);
for comp=1:numDCs
    %compute GFT
    idx=p(r(comp):r(comp+1)-1);
    [ GFT, ~ ] = compute_GFT( Wcurr(idx,idx), Qcurr(idx) );
    Utmp=zeros(size(Wcurr,1),length(idx));
    Utmp(idx,:)=GFT';
    U=[U,Utmp];
    isDCtmp=zeros(length(idx),1);
    isDCtmp(1)=1;
    isDC=[isDC;isDCtmp];
    %Gfreq_level = [Gfreq_level; Gfreq_tmp];
    weights_level = [weights_level ; ones(length(idx),1)*sum(Qcurr(idx))];
    %compute average of points in block
    Vblock_new(comp,:) =  sum(diag(Qcurr(idx))*Vblock(idx,:),1)/sum(Qcurr(idx));
    %
    
end

% Ahat_1 = U'*A;
 %isDC_index = find(isDC);
% notDC_index = find(~isDC);
% Ahat_low = Ahat_1(isDC_index,:);%low pass coeffs for further processing
% Ahat_high=Ahat_1(notDC_index,:);%high pass coeffs
%level 2
%compute transform (complete graph) for DC coeffs
Qnew = weights_level(isDC==1);
Wnew = complete_graph(Vblock_new);

[ GFT_new, ~ ] = compute_GFT( Wnew, Qnew );

%invert level 2 of the transform
A1 = GFT_new'*Ahat(1:numDCs,:);

%reorganize coeffs (scramble DCs and AC coeffs)

%DC_indices = find(isDC);
%AC_indices = find(~isDC);

A1_permutted = zeros(size(A1));
A1_permutted(isDC==1,:) = A1;
A1_permutted(isDC==0,:) = Ahat(numDCs+1:end,:);


%invert level  1 of transform
 A = U * A1_permutted;


% Coeff = [GFT_new*Ahat_low;Ahat_high];
% Gfreq = [Gfreq_new ;Gfreq_level(notDC_index) ];
% weights = [ones(length(Qnew),1)*sum(Qnew) ; weights_level(notDC_index)];

end
% function W = complete_graph(V)
% % V: nx3. n points, vi is the i-th row of V
% %computes complete graph with edge weights 1/sqrt(distance(vi,vj))
% N = size(V,1);
% %
% %compute EDM
% squared_norms = sum(V.^2,2);
% D = sqrt(repmat(squared_norms,1,N) + repmat(squared_norms',N,1) - 2*(V*V'));
% % D = squareform(pdist(coords, 'euclidean')); % pairwise distances, n-by-n
% % matrix% only use pdist if have the statistics/ML toolbox
% iD = D.^(-1);
% 
% iD(find(D==0))   =0;
% W=iD' + iD;
% 
% %idx = find(iD~=0);
% 
% %[I, J] = ind2sub( size(D), idx );
% 
% %edge = [I, J];
% end