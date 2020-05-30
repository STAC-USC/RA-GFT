function [ GFT, Gfreq ] = compute_GFT( Adj, Q )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Qm = diag(Q.^(-1/2)); %assume Q is a vector
L = full(Qm*w2l(Adj)*Qm);

[GFT,D] = eig(L);
[Gfreq, idxSorted ] = sort( diag(D), 'ascend' );
GFT = GFT(:, idxSorted);
GFT(:,1)=abs(GFT(:,1));
GFT = GFT';
Gfreq(1) = abs(Gfreq(1));
end

