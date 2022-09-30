clear variables
close all
clc

%% Grids
kappa = 0.001; % diffusivity
W = 0.1;   % width of the actuator

% spatial grid
L = 1; % length
Nx = 101; % total number of spatial nodes
xgrid = linspace(0,L,Nx);
dx = L/(Nx-1); % step in position

% temporal grid
T_max = 200; % max simulation time
N_tot = 200 + 1;
tgrid = linspace(0,T_max,N_tot);
dt = tgrid(2) - tgrid(1); % time step
tgrid2 = tgrid(1:end-1) + diff(tgrid)/2;

%% Set Matrices
A = sparse(Nx-1,Nx-1);
Ae = kappa/(dx^2)*[-1, 1; 1 -1];
for ii = 1:Nx-2
    ind = [ii,ii+1];
    A(ind,ind) = A(ind,ind) + Ae;
end
ARBM = (Nx-1)/(Nx-2)*A;
ind = [Nx-1,1];
A(ind,ind) = A(ind,ind) + Ae;

Mass = speye(Nx-1, Nx-1);
X0 = 10*exp(-2*(1 + sin(2*pi*xgrid/L))); X0 = X0(1:end-1);
B = (xgrid <= W/2 | 1-W/2 <= xgrid).';
B = B(1:end-1) / W;
Q = 100*speye(Nx-1,Nx-1)*dx;
R = 1;
xd =@(t) 0;

Nh = 1;   % number of time steps in \Delta t
tau = 10; % control horizon in MPC
Nk = 20;  % number of considered realizations

%% compute reference solution
u0 = zeros(1, N_tot-1);
[~, ~, muinfty] = reference(A,B,Q,R,X0,tgrid);
[uref, Xref, Jref, duration_ref] = compute_control(A, X0, B, u0, Q, R, xd, tgrid, Mass);
Xref_norm = stateL2norm(Xref, dx);
figure(1); plot(tgrid2, uref, 'k')
figure(2); plot(tgrid, Xref_norm, 'k')

% %% compute RBM solution
% for kk = 1:Nk
%     batches = choose_batches(N_tot, Nx, Nh);
%     [uRBM, XRBM, JRBM, duration_RBM(kk)] = compute_controlRBM(ARBM, X0, B, u0, Q, R, xd, tgrid, Mass, batches);
%     XRBM_norm = stateL2norm(XRBM, dx);
%     eRBM_norm = stateL2norm(XRBM - Xref, dx);
%     eRBM(kk) = max(eRBM_norm);
%     euRBM(kk) = max(abs(uRBM-uref));
%     figure(1); hold on; plot(tgrid2, uRBM, 'g')
%     figure(2); hold on; plot(tgrid, XRBM_norm, 'g')
%     figure(3); hold on; plot(tgrid2, abs(uRBM-uref), 'g')
%     figure(4); hold on; plot(tgrid, eRBM_norm, 'g')
% end

T_list = [10:2:30];
for ll = 1:length(T_list)
    T = T_list(ll);   % prediction horizon in MPC
    %% compute MPC solution
    [uMPC, XMPC, duration_MPC(ll)] = compute_MPCcontrol(A, X0, B, u0, Q, R, xd, tgrid, Mass, T, tau);
    XMPC_norm = stateL2norm(XMPC, dx);
    eMPC_norm = stateL2norm(XMPC - Xref, dx);
    eMPC(ll) = max(eMPC_norm);
    euMPC(ll) = norm(uMPC - uref)*sqrt(dt); %max(abs(uMPC - uref));
    figure(1); hold on; plot(tgrid2, uMPC, 'r')
    figure(2); hold on; plot(tgrid, XMPC_norm, 'r')
    figure(3); hold on; plot(tgrid2, abs(uMPC-uref), 'r')
    figure(4); hold on; plot(tgrid, eMPC_norm, 'r')
    
    %% compute RBM-MPC solution
    for kk = 1:Nk
        [uRBMMPC, XRBMMPC, duration_RBMMPC(ll,kk)] = compute_MPCcontrolRBM(A, ARBM, X0, B, u0, Q, R, xd, tgrid, Mass, T, tau, Nx, Nh);
        XRBMMPC_norm  = stateL2norm(XRBMMPC, dx);
        eRBMMPC_norm  = stateL2norm(XRBMMPC - Xref, dx);
        eRBMMPC2_norm = stateL2norm(XRBMMPC - XMPC, dx);
        eRBMMPC(ll,kk) = max(eRBMMPC_norm);
        eRBMMPC2(ll,kk) = max(eRBMMPC2_norm);
        euRBMMPC(ll,kk) = norm(uRBMMPC - uref)*sqrt(dt);% max(abs(uRBMMPC - uref));
        euRBMMPC2(ll,kk) = norm(uRBMMPC - uMPC)*sqrt(dt);% max(abs(uRBMMPC - uMPC));
        figure(1); hold on; plot(tgrid2, uRBMMPC, 'b')
        figure(2); hold on; plot(tgrid, XRBMMPC_norm, 'b')
        figure(3); hold on; plot(tgrid2, abs(uRBMMPC-uref), 'b')
        figure(4); hold on; plot(tgrid, eRBMMPC_norm, 'b')
        figure(5); hold on; plot(tgrid, eRBMMPC2_norm, 'b')
    end
end

% eRBM_mean     = mean(eRBM,2);     eRBM_std     = std(eRBM,[],2);
eRBMMPC_mean  = mean(eRBMMPC,2);  eRBMMPC_std  = std(eRBMMPC,[],2);
eRBMMPC2_mean = mean(eRBMMPC2,2); eRBMMPC2_std = std(eRBMMPC2,[],2);

% euRBM_mean     = mean(euRBM,2);     euRBM_std     = std(euRBM,[],2);
euRBMMPC_mean  = mean(euRBMMPC,2);  euRBMMPC_std  = std(euRBMMPC,[],2);
euRBMMPC2_mean = mean(euRBMMPC2,2); euRBMMPC2_std = std(euRBMMPC2,[],2);

% timeRBM_mean     = mean(duration_RBM,2);     timeRBM_std     = std(duration_RBM,[],2);
timeRBMMPC_mean  = mean(duration_RBMMPC,2);  timeRBMMPC_std  = std(duration_RBMMPC,[],2);

fig = figure(6);
plot(T_list, eMPC, 'DisplayName', 'MPC')
hold on
plot(T_list, eMPC(1)*exp(2*muinfty*(T_list-10)))
plot(T_list, euMPC, 'DisplayName', 'uMPC')
plot(T_list, euMPC(1)*exp(2*muinfty*(T_list-10)))
% errorbar(T_list, eRBM_mean*ones(size(T_list)), 2*eRBM_std*ones(size(T_list)), 'DisplayName', 'RBM')
% errorbar(T_list, euRBM_mean*ones(size(T_list)), 2*euRBM_std*ones(size(T_list)), 'DisplayName', 'uRBM')
errorbar(T_list, eRBMMPC_mean, 2*eRBMMPC_std, 'DisplayName', 'RBMMPC')
errorbar(T_list, euRBMMPC_mean, 2*euRBMMPC_std, 'DisplayName', 'uRBMMPC')
errorbar(T_list, eRBMMPC2_mean, 2*eRBMMPC2_std, 'DisplayName', 'RBMMPC2')
errorbar(T_list, euRBMMPC2_mean, 2*euRBMMPC2_std, 'DisplayName', 'uRBMMPC2')
xlabel 'T'
ylabel 'error'
savefig(fig,'C3_Fig6')

fig = figure(7);
plot(T_list, duration_MPC, 'DisplayName', 'MPC')
hold on
% errorbar(T_list, timeRBM_mean*ones(size(T_list)), 2*timeRBM_std*ones(size(T_list)), 'DisplayName', 'RBM')
errorbar(T_list, timeRBMMPC_mean, 2*timeRBMMPC_std, 'DisplayName', 'RBMMPC')
xlabel 'T'
ylabel 'duration'
savefig(fig,'C3_Fig7')