%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   This is the sub-function for calculating the frequencies for the
%   complementary elements, because the frequency for the underlying
%   parameter has been set as the OmegaMax
%
%   Input:
%           paraDim:                  The dimension of unknown parameters
%           OmegaMaxComp:    The maximum frequency can be set for the
%                                           complementary elements
%           underlyingIndex:       Index of underlying unknown parameter
%
%   Output:
%           OmegaComp:       A vector for representing the frequencies of
%                                       the complementary elements.
%
%   Author: Xin Liu
%
%   Date:21/11/12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OmegaComp = setFreq(paraDim, OmegaMaxComp, underlyingIndex)

OmegaComp = zeros(1,paraDim);

if paraDim==1
    
    OmegaComp = 1;
    
elseif OmegaMaxComp == 1
    
    OmegaComp = ones(1, paraDim);
    
else
    
    if(OmegaMaxComp < paraDim)   
        
    % the maximum allowable frequency (OMciMAX) is less than
    % the dimension of unknown parameters, then the
    % INFD is chose as the smaller one (i.e. the OMciMAX)
    % if the OMciMAX is larger than the dimension, then
    % INFD is used as the number of dimension, this is
    % because the INFD is used to decide the step of
    % frequency sequence, and the thumb-up rule is the
    % step should be as large as possible.
    
        INFD = OmegaMaxComp;
        
    else
        
        INFD = paraDim;
        
    end
    
    %  ISTEP is the step for the
    %  frequency selection    
    ISTEP = round((OmegaMaxComp-1)/(INFD-1)); 
   
    OTMP = 1:ISTEP:INFD*ISTEP;
    
    fl_INFD = floor(INFD);  %   find the small integer
    
    for i=1:paraDim
        j = mod(i-1,fl_INFD)+1;
        OmegaComp(i) = OTMP(j);
    end
end

%   frequency for the underlying parameter is set
%   as Omega_max, thus, it is eliminate from the 
%   frequencies for complementary elements
OmegaComp(underlyingIndex)=[];






















