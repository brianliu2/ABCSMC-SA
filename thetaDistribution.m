%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   This is the sub-function for calculating the distribution for the
%   transformed theta by considering their distribution that is assumed to
%   be satisfied
%
%   Input:
%           theta:      The transformed theta
%           type:       The type of distribution is assumed to satisfy by
%                       theta
%
%   Output:
%           dist:       The distribution of theta
%   Author: Xin Liu
%
%   Date:21/11/12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = thetaDistribution(theta,type)

global pmax;
global pmin;
global pmean;
global pstd;

paraDimSA = size(theta,2);
dist = zeros(size(theta,1),paraDimSA);

switch lower(type)  % Convert string to lowercase
    case {'unif'}
        for d = 1:paraDimSA
            dist(:,d) = (theta(:,d) .* (pmax(d) - pmin(d))) + pmin(d);
        end
    case {'normal'}
        for d = 1:paraDimSA
            dist(:,d) = norminv(theta(:,d),pmean(d),pstd(d));
        end
    case {'lognormal'}
        for d = 1:paraDimSA
            dist(:,d) = norminv(theta(:,d),log(pmean(d)),pstd(d));
        end
    otherwise
        disp('Undefined PDF');
end




