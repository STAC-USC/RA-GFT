function L=w2l(W)
% W2L weight matrix to Laplacian matrix
% 
% L=w2l(W)
% 
% by: KS Lu
% 20170712
%
if any(W<0)
    %error('W is not a valid weight matrix');
end
L=diag(sum(W))-W+diag(diag(W));
end