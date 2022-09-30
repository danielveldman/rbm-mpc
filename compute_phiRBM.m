function phi = compute_phiRBM(A, Q, X, xd, tgrid, Mass, batches)

% solves the backward dynamics
% Mass.'*\dot{phi}(t) = A.'phi(t) + Q*(x(t) - xd(t)),      phi(T) = 0,
% with a scheme that leads to discretely consistent gradients when the
% RBM Crank-Nicholson scheme is used for the forward dynamics.

% see:
% Apel and Flaig, Crank--Nicolson Schemes for Optimal Control Problems
% with Evolution Equations, SIAM journal on Numerical Analysis, 50(3), 2012. 

dt = diff(tgrid); ndt = length(dt);
phi = zeros(size(X,1), length(dt)); N = size(X,1);
dx = X - xd(tgrid.'); I = speye(N);

ind = 1:N;
indii = circshift(ind, -batches(end));
phi(indii,end) = (I-dt(end)/2*A.')\((dt(end)*Q*dx(indii,end)/2));
for ii = ndt-1:-1:1
    indii = circshift(ind, -batches(ii));
    phi(indii,ii) = (Mass-dt(ii)/2*A.')\(Mass*phi(indii,ii+1)+dt(ii)/2*A.'*phi(indii,ii+1) + dt(ii)*Q*dx(indii,ii+1));
end