function batches = choose_batches(NT, Nx, Nh)
batches = zeros(1, NT);
Nb = ceil(NT/Nh);
rand_ind = randi([1 Nx-1],1,Nb);

for kk = 0:Nb-2
    batches(kk*Nh+1:(kk+1)*Nh) = rand_ind(kk+1);
end
batches((Nb-1)*Nh+1:end) = rand_ind(end);