function [epsilon_hist,paraSpace] = ABCSMC_Fcn(paraSamples,Options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Description:    Implementation is for running ABC-SMC algorithm
%                           as a function which needs User to define four
%                           input arguments
%
%   Input:
%                        paraSamples:     Initial samples for unknown parameters
%                        M:                       Integer factor
%                        alphaSMC:          Factor for weighting scheme
%                        N:                        Number of samples used
%                        delta:                   Discount factor for the transition
%                                                    kernel.
%                       options:                Options defined in the main
%                                                    function for running
%                                                    the ABC-SMC
%
%   Output:         
%                       paraSpace:  3-D matrix M*N*T contains estimations
%                                           at every iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global stateDim;
global paraTrue;
global stateTrue;
global sensLab;
global runCount;
%   Pass environmental settings 
tspan = Options.tspan;
alphaSMC = Options.alphaSMC;
%   Dimension of state
initTrue = stateTrue;
%   number of particles N and integer factor M
N = Options.particleNum;
M = Options.M;

%   Counter for resampling
counter = 0;          

%   Threshold for the resampling
essthresh = N/2;      

%   Pre-allocate the memory for weights 
weights = zeros(1,N);

%   Alive indicator vector for samples of current population, following the
%   I() function in pseudo-code for algorithm
npa = zeros(1,N);

%   Alive indicator vector for samples of previous population
pa  = zeros(1,N);       

%   Initial tolerance
epsilon = Options.epsilonInit;
epsilontarget = Options.epsilontarget;
tic

%   Parametrizations for the transition kernel, obtained by the
%   user-specific discount factor delta
delta = Options.delta;
a = (3*delta - 1)/(2*delta);
covParaCoeff = 1 - a^2;
%%  Initialization

%    Define the prior distributions as the accepted population at 1st
%    iteration
paraSpace(:,:,1) = paraSamples;


%   Accuracy options for ODE solver
options   = odeset('RelTol', 1e-02);

%   Define  the prior importance weights in 1st iteration to 1/N
weights = ones(1,N)./N;

timeLength = Options.timeLength;

epsilon_hist(1) = epsilon;
%%  Main loop

loop = 1;

%   Generate true dataset
sol = ode45(@(t,y) hs_ode(t,y,paraTrue),[0 timeLength],initTrue,options);
xtrue = deval(sol,tspan);
xtrue = xtrue + sqrt(0.001) * randn(size(xtrue));

%   For the first time instance
for i = 1:M        
        parfor j = 1:N
            %   Sovle the ODEs of the repressilator system 
            solSMC = ode45(@(t,y) hs_ode(t,y,...
                             paraSamples(:,j)),[0 timeLength],initTrue,options);

            %   Use the last column in the solution as the underlying
            %   state, if the solver is underflow.
            if(solSMC.x(end) >= timeLength)
                temp = deval(solSMC,tspan);
            else
                temp = deval(solSMC,tspan);
            end

            data_ParaLoop(i,j,:,:) = temp;

            fprintf('%d particles in %d M of %d loop and recent epsilon is %d \n',j,i,loop,epsilon);
       end        
end

while (epsilon > epsilontarget)
    
     % Compute the current distance between the synthetic dataset and the
     % real dataset
    distance = distanceFun( data_ParaLoop(:,:,1,:), data_ParaLoop(:,:,2,:), data_ParaLoop(:,:,3,:), xtrue);   
    
    % Calculate the next tolerance by following the condition presented in
    % the pseudo-code as PA() < alpha * PA()
    epsilon_old = epsilon;
    reflevel = alphaSMC * tpa(epsilon_old, N, distance);
    epsilon  = fzero(@(epsilon) tpa(epsilon, N, distance) - reflevel, 0, epsilon_old);
    epsilon_hist(loop+1) = epsilon;
    %   Determine if the current tolerance is less than the ideal tolerance,
    %   terminate the program if so
    if (epsilon < epsilontarget)
        epsilon = epsilontarget;
        epsilon_hist(loop+1) = epsilon;
    end    
  
    %   Compute alive indicator vector considering the previous tolerance
    npa_old = sum((distance < epsilon_old),1);
    
    %   Compute alive indicator vector considering the current tolerance
    npa = sum((distance < epsilon),1);
    
    %   Find the indices for the alive or inactive samples 
    actIndex = find(npa_old > 0);
    unactIndex = find(npa_old == 0);
    
    %   Weight calculation according to the equations stated in the
    %   pseudo-code of ABC-SMC
    weights(1,actIndex) = weights(1,actIndex) .* npa(1,actIndex)./npa_old(1,actIndex);
    
    %   Assign the zero weights to inactive samples
    weights(1,unactIndex) = zeros(1,length(unactIndex));
    
    %   Normalized the weights
    weights = weights./sum(weights);     
    
    %   Resample if neccessary
    if (sum(weights.^2) * essthresh > 1)
        N_sons = rsdet(weights);
         [paraSamples,data_ParaLoop] = copy(paraSamples,data_ParaLoop,N_sons);       
        counter = counter + 1;
        counter
    end
        
      
    %   Move the samples according to the transition kernel which is
    %   presented as Eqn.8 in the manuscript
    for c = 1:size(paraSamples,1)
        if runCount == 1
            paraSamples(c,:) = a * paraSamples(c,:)...
                                    + ones(1,N)*(1-a) * mean((paraSpace(c,:,loop)'))...
                                    + sqrt(covParaCoeff*var((paraSpace(c,:,loop)')))...
                                    * randn(1,N);
        elseif runCount == 2
         if strcmpi(sensLab{c},'stiff')
            paraSamples(c,:) = a * paraSamples(c,:)...
                                    + ones(1,N)*(1-a) * mean((paraSpace(c,:,loop)'))...
                                    + sqrt(covParaCoeff*var((paraSpace(c,:,loop)')))...
                                    * randn(1,N);
         end
        end
    end
    
    %   Generate the synthetic dataset by using updated samples for
    %   parameters
    for i = 1:M        
        parfor j = 1:N
            
            %   Sovle the ODEs of the repressilator system 
            solSMC = ode45(@(t,y) hs_ode(t,y,...
                             paraSamples(:,j)),[0 timeLength],initTrue,options);
                         
            %   Use the last column in the solution as the underlying
            %   state, if the solver is underflow.
            if(solSMC.x(end) >= timeLength)
                temp = deval(solSMC,tspan);
            else
                temp = deval(solSMC,tspan);
            end
            
            data_ParaLoop(i,j,:,:) = temp;
            
            fprintf('%d particles in %d M of %d loop and recent epsilon is %d \n',j,i,loop,epsilon);
        end        
    end
    
    loop
    loop = loop + 1;   
    paraSpace(:,:,loop) = paraSamples;
end




