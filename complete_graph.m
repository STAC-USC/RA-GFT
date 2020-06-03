function W = complete_graph(V)
% V: nx3. n points, vi is the i-th row of V
%computes complete graph with edge weights 1/sqrt(distance(vi,vj))
N = size(V,1);
%
%compute EDM
squared_norms = sum(V.^2,2);
D = sqrt(repmat(squared_norms,1,N) + repmat(squared_norms',N,1) - 2*(V*V'));
% D = squareform(pdist(coords, 'euclidean')); % pairwise distances, n-by-n
% matrix% only use pdist if have the statistics/ML toolbox
iD = D.^(-1);

iD(find(D==0))   =0;
W=iD' + iD;

% A = spones(W);
% 
% nnz_distances = 1./nonzeros(W);
% sigma = sqrt(mean(nnz_distances.^2)/3);
% 
% W = A.*exp(-0.5*D.^2/sigma^2);
% Deginv = diag(1./sum(W));
% W = Deginv*W*Deginv;
end