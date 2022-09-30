function [Xref, uref, muinfty] = reference(A,B,Q,R,X0,tgrid)

[~,K] = icare(full(A),B,full(Q),R);
% K = R\(B.'*Pinfty);
muinfty = max((eig(A-B*K)));

Nx = length(A);
E = speye(Nx, Nx);
Xref = zeros(Nx,length(tgrid));
Xref(:,1) = X0;
dt = tgrid(2) - tgrid(1);
for tt = 2:length(tgrid)
    Xref(:,tt) = (E - dt/2*(A-B*K))\((E+dt/2*(A-B*K))*Xref(:,tt-1));
end

u = -K*Xref;
uref = zeros(size(B,2), length(tgrid)-1);
for tt = 1:length(tgrid)-1 % get u on the intermediate points
    uref(tt) = (u(tt)+u(tt+1))/2;
end

