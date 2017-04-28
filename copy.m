function [theta,x] = copy(thetas,xs,N_sons)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Resample the particles with replacement by using offspring indices from
%   previous function rsdet.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%   Calculate the integer factor M
M           = size(xs,1);

%   Pre-allocate the memory for resampled parameter and state
theta      = zeros(size(thetas));
x          = zeros(size(xs));



%   Mark the index for non-zero weights
indexNonzero = find(N_sons > 0);



% for i = 1:length(indexNonzero)
%     for ind = index:indexNonzero(i)
%         theta_test(:,ind) = thetas(:,indexNonzero(i));
%     end
%     index = indexNonzero(i)+1;
% end

%   Initialize index
index       = 1;
%   Resampling process for parameter space, for which the 
%   dimension of data is 2-D, i.e. 
%   parameter(parameter dimension, number of particles)
for i = 1:length(indexNonzero) 
    temp = thetas(:,indexNonzero(i));
    theta(:,index:index+N_sons(indexNonzero(i))-1) = temp(:,ones(1,N_sons(indexNonzero(i))));
    index = index + N_sons(indexNonzero(i));
end


%   Resampling process for state space, which 
%   is more complicate than the resampling for parameter,
%   since the integer factor M is considered for state space,
%   therefore, the dimension for state data is 4-D, i.e.
%   state(integer factor, particle number, state dimension, 
%           sample points for one time interval).

for j = 1:M
    index = 1;
    for i = 1:length(indexNonzero)
        temp = xs(j,indexNonzero(i),:,:);
        x(j,index:index+N_sons(indexNonzero(i))-1,:,:) = temp(:,ones(1,N_sons(indexNonzero(i))),:,:);
        index = index + N_sons(indexNonzero(i));
    end   
end





























