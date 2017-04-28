function d = distanceFun(x1,x2,x3,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The distance function for measuring distance the 
%   pseudo observations and true observations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Pre-allocate the memory for storing distances
%   the dimension of distance matrix is M * N, in which
%   M is the integer factor and N is the number of particles

d = zeros(size(x1,1),size(x1,2));

for i = 1:size(x1,2)    %   Loop for number of particles N    
    for j = 1:size(x1,1)    %   Loop for integer factor M
            %   the dimension of true state is 2-D
            %   y(state dimension, data points for each time interval).
            %   the dimension of synthetic data is 4-D
            %   x(integer factor, number of particles, state dimension, data points for each time interval)
            %   In calculation, the first summation is for summing
            %   all states in the system, the second summation is 
            %   for data points each state in each time interval.
            d(j,i) = sqrt(sum( sum(((reshape(x1(j,i,:,:),1,size(y,2)) - y(1,:))/max(reshape(x1(j,i,:,:),1,size(y,2)))).^2)...
                                     +sum(((reshape(x2(j,i,:,:),1,size(y,2))- y(2,:))/max(reshape(x2(j,i,:,:),1,size(y,2)))).^2)...
                                     +sum((((reshape(x3(j,i,:,:),1,size(y,2)) - y(3,:)))/max(reshape(x3(j,i,:,:),1,size(y,2)))).^2)));
    end
end








