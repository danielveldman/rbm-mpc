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

%% compute reference solution
u0 = zeros(1, N_tot-1);
[uref, Xref, Jref, duration_ref] = compute_control(A, X0, B, u0, Q, R, xd, tgrid, Mass);
[Xref2, uref2, muinfty] = reference(A,B,Q,R,X0,tgrid);
Xref_norm = stateL2norm(Xref, dx);
Xref2_norm = stateL2norm(Xref2, dx);
figure(1); plot(tgrid2, uref, tgrid2, uref2)
figure(2); plot(tgrid, Xref_norm, tgrid, Xref2_norm)

%% compute RBM solution
batches = choose_batches(N_tot, Nx, Nh);
[uRBM, XRBM, JRBM, duration_RBM] = compute_controlRBM(ARBM, X0, B, u0, Q, R, xd, tgrid, Mass, batches);
XRBM_norm = stateL2norm(XRBM, dx);
figure(1); hold on; plot(tgrid2, uRBM)
figure(2); hold on; plot(tgrid, XRBM_norm)

%% compute MPC solution
[uMPC, XMPC, duration_MPC] = compute_MPCcontrol(A, X0, B, u0, Q, R, xd, tgrid, Mass, T, tau);
XMPC_norm = stateL2norm(XMPC, dx);
figure(1); hold on; plot(tgrid2, uMPC)
figure(2); hold on; plot(tgrid, XMPC_norm)

%% compute RBM-MPC solution
[uRBMMPC, XRBMMPC, duration_RBMMPC] = compute_MPCcontrolRBM(A, ARBM, X0, B, u0, Q, R, xd, tgrid, Mass, T, tau, Nx, Nh);
XRBMMPC_norm = stateL2norm(XRBMMPC, dx);
figure(1); hold on; plot(tgrid2, uRBMMPC); xlabel 'time'; ylabel 'control'
figure(2); hold on; plot(tgrid, XRBMMPC_norm); xlabel 'time'; ylabel 'state'