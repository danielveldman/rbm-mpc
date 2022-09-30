function norms = stateL2norm(X, dx)
norms = zeros(1,size(X,2));
for ii = 1:size(X,2)
    norms(ii) = norm(X(:,ii))*sqrt(dx);
end