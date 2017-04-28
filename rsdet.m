function N_sons = rsdet(weights)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Determine the offsping indices for the resample process 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(weights);

%   Preallocate the memory for offsprings indices
N_sons = zeros(1,N);

%   Construct the cumulative distribution
dist = cumsum(weights); % cumsum(1:5) = [1 3 6 10 15]

%  Construct the offspring indices, following the equation in the tutorial
%  from Arulampalam et al., PP180 Algorithm 2: u_j = u_1 + 1/Ns(j-1)

%   Generate indicator array
aux  = rand;
u = aux:1:(N -1 + aux);
u = u./N;

j = 1;

for i = 1:N
    while (u(1,i) > dist(1,j))
        j = j + 1;
    end
    N_sons(1,j) = N_sons(1,j) + 1; 
end

























