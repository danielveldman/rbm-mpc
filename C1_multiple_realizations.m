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
T = 15;   % prediction horizon in MPC
tau = 10; % control horizon in MPC
Nk = 20;  % number of considered realizations

%% compute reference solution
u0 = zeros(1, N_tot-1);
[uref, Xref, Jref, duration_ref] = compute_control(A, X0, B, u0, Q, R, xd, tgrid, Mass);
Xref_norm = stateL2norm(Xref, dx);
figure(1); plot(tgrid2, uref, 'k')
figure(2); plot(tgrid, Xref_norm, 'k')

%% compute RBM solution
for kk = 1:Nk
    batches = choose_batches(N_tot, Nx, Nh);
    [uRBM, XRBM, JRBM, duration_RBM(kk)] = compute_controlRBM(ARBM, X0, B, u0, Q, R, xd, tgrid, Mass, batches);
    XRBM_norm = stateL2norm(XRBM, dx);
    eRBM_norm = stateL2norm(XRBM - Xref, dx);
    eRBM(kk) = max(eRBM_norm);
    euRBM(kk) = max(abs(uRBM - uref));
    figure(1); hold on; plot(tgrid2, uRBM, 'g')
    figure(2); hold on; plot(tgrid, XRBM_norm, 'g')
    figure(3); hold on; plot(tgrid2, abs(uRBM-uref), 'g')
    figure(4); hold on; plot(tgrid, eRBM_norm, 'g')
end

%% compute MPC solution
[uMPC, XMPC, duration_MPC] = compute_MPCcontrol(A, X0, B, u0, Q, R, xd, tgrid, Mass, T, tau);
XMPC_norm = stateL2norm(XMPC, dx);
eMPC_norm = stateL2norm(XMPC - Xref, dx);
eMPC = max(eMPC_norm);
euMPC = max(abs(uMPC - uref));
figure(1); hold on; plot(tgrid2, uMPC, 'r')
figure(2); hold on; plot(tgrid, XMPC_norm, 'r')
figure(3); hold on; plot(tgrid2, abs(uMPC-uref), 'r')
figure(4); hold on; plot(tgrid, eMPC_norm, 'r')

%% compute RBM-MPC solution
for kk = 1:Nk
    [uRBMMPC, XRBMMPC, duration_RBMMPC(kk)] = compute_MPCcontrolRBM(A, ARBM, X0, B, u0, Q, R, xd, tgrid, Mass, T, tau, Nx, Nh);
    XRBMMPC_norm  = stateL2norm(XRBMMPC, dx);
    eRBMMPC_norm  = stateL2norm(XRBMMPC - Xref, dx);
    eRBMMPC2_norm = stateL2norm(XRBMMPC - XMPC, dx);
    eRBMMPC(kk) = max(eRBMMPC_norm);
    eRBMMPC2(kk) = max(eRBMMPC2_norm);
    euRBMMPC(kk) = max(abs(uRBMMPC - uref));
    euRBMMPC2(kk) = max(abs(uRBMMPC - uMPC));
    figure(1); hold on; plot(tgrid2, uRBMMPC, 'b')
    figure(2); hold on; plot(tgrid, XRBMMPC_norm, 'b')
    figure(3); hold on; plot(tgrid2, abs(uRBMMPC-uref), 'b')
    figure(4); hold on; plot(tgrid, eRBMMPC_norm, 'b')
    figure(5); hold on; plot(tgrid, eRBMMPC2_norm, 'b')
end

eRBM_mean     = mean(eRBM);     eRBM_std     = std(eRBM);
eRBMMPC_mean  = mean(eRBMMPC);  eRBMMPC_std  = std(eRBMMPC);
eRBMMPC2_mean = mean(eRBMMPC2); eRBMMPC2_std = std(eRBMMPC2);

euRBM_mean     = mean(euRBM);     euRBM_std     = std(euRBM);
euRBMMPC_mean  = mean(euRBMMPC);  euRBMMPC_std  = std(euRBMMPC);
euRBMMPC2_mean = mean(euRBMMPC2); euRBMMPC2_std = std(euRBMMPC2);

timeRBM_mean     = mean(duration_RBM);     timeRBM_std     = std(duration_RBM);
timeRBMMPC_mean  = mean(duration_RBMMPC);  timeRBMMPC_std  = std(duration_RBMMPC);