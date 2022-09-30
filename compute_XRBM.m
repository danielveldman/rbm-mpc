function X = compute_XRBM(ARBM, X0, B, U, tgrid, Mass, batches)

% solves the forward dynamics
% Mass*\dot{x}(t) = Ax(t) + Bu(t),      x(0) = X0,
% on the time grid tgrid by the Crank-Nicholson scheme 
% using the random batch method (RBM). 

X = zeros(length(X0), length(tgrid));
X(:,1) = X0;
dt = diff(tgrid);
ind = 1:length(X0);
for ii = 2:length(tgrid)
    indii = circshift(ind, -batches(ii-1));
    X(indii,ii) = (Mass-dt(ii-1)/2*ARBM)\(Mass*X(indii,ii-1)+dt(ii-1)/2*ARBM*X(indii,ii-1) + dt(ii-1)*B(indii,:)*U(:,ii-1));
end