function [uMPC, XMPC, duration] = compute_MPCcontrolRBM(A, ARBM, X0, B, u0, Q, R, xd, tgrid, Mass, T, tau, Nx, Nh)

tstart = tic;
nt = length(tgrid);
XMPC = zeros(length(X0), length(tgrid)); XMPC(:,1) = X0;
uMPC = zeros(1, length(tgrid)-1);

ind1 = 1; ii = 0;
while ind1 < nt
    % compute optimal control on a finite horizon
    ind2 = find(tgrid <= ii*tau + T  , 1, 'last');
    batches = choose_batches(ind2-ind1, Nx, Nh);
    uii = compute_controlRBM(ARBM, XMPC(:,ind1), B, u0(:,1:ind2-ind1), Q, R, xd, tgrid(ind1:ind2), Mass, batches);
    
    % compute the resulting state trajectory and assign applied control
    ind3 = find(tgrid <= ii*tau + tau, 1, 'last');
    tgridii = tgrid(ind1:ind3); uii = uii(:,1:ind3-ind1);
    XMPC(:,ind1:ind3) = compute_X(A, XMPC(:,ind1), B, uii, tgridii, Mass);
    uMPC(:,ind1:ind3-1) = uii;
    
    % update for next iteration
    ind1 = ind3; ii = ii+1;
end
duration = toc(tstart);