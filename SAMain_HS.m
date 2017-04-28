%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   This is the main file for the global sensitivity analysis for computing
%   the first order and total effect indices for the Heat Shock responce
%   model, while the used method is Extended Fourier Amplitude Sensitivity
%   Test (eFAST).
%
%   Author: Xin Liu
%
%   Date:21/11/12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

global sumSi;
%%  Initial settings for the eFAST

reqN = nR * nS * paraDimSA;   % the required samples for the SA

timeLength = timeLengthSA;
%   The time settings for ODE solver
tspan = linspace(0,timeLengthSA,dataPointsSA);
%   The time points for executing the SA
time_point = tspan;    
%   Initialization of the states
xInit = [0;0;0]; 

%   Compute the highest frequence of the underlying unknown parameter
OmegaMax = floor((nS- 1) / (2 * M) );

%   Assign the minimum sample size to be used in the FAST from considering
%   the highest frequence
minSampSize = 2 * M * OmegaMax + 1;

%   Pre-allocate the memory for the data

data(minSampSize, length(time_point), stateDim, paraDimSA, nR) = 0;

opts= odeset('RelTol',1e-6) ;

% %   Define the number of processor used for parallel computing
% Num_processor = 4;
% 
% %   Activate the parallel session
% matlabpool('open', Num_processor);
%%  Main Loop for the SA
for i = 1:paraDimSA
    %   1. Select the frequencies for the complementary elements, as the
    %   frequency of the underling parameter has been assigned as the
    %   omega_max.
    OmegaComp = setFreq(paraDimSA, OmegaMax*M/2, i);
    
    %   Loop over the search curves nR
    for R = 1:nR
        
        %   2. Assign the frequency vector for all parameter, each search
        %   curve has the different frequency vector, therefore loop over
        %   all search curves
        compInd = 1;    %   Index in frequencies for the complementary elements
        omega = zeros(1,paraDimSA);
        for paraInd = 1:paraDimSA
            if (paraInd == i)
                omega(i) = OmegaMax;
            else
                omega(paraInd) = OmegaComp(compInd);
                compInd = compInd + 1;
            end
        end   
        
        %   3. Compute the parameters value in frequency domain, such
        %   transformation has been done by the Eqn.20 in the literautre
        %   theta = 0.5 + 1/pi * arcsin(sin(omega*s + xi))
        
        %   Calculation of xi, ranging in [0 2*pi]
        randPhaShft = 2 * pi * rand(1, paraDimSA);
        %   Calculation of S_{i}, in literature it is computed as
        %   S_{k} = pi/Ns * (2*k - nS -1) where k = 1...nS
        sRange = pi * (2 * (1:minSampSize) - minSampSize - 1)/minSampSize;
        
        %   Omega vector is obtained by the step.2, then we are going to
        %   calculate the angle here. Before the calculation, we need to
        %   extense the random phase shift vector xi to the matrix so that
        %   satisfy the dimension requirement.
        randPhaShftMat = randPhaShft(ones(minSampSize,1),1:paraDimSA);
        angle = omega' * sRange + randPhaShftMat';
        
        theta(:,:,i,R) = 0.5 + asin(sin(angle'))/pi;
        %   Transform the theta according to the corresponding PDF of which
        %   is assume to be followed
        currentTheta = thetaDistribution(theta(:,:,i,R),'unif')';
        
        %   4. Do the model evaluation by using the theta obtained at the
        %   step.3
        for sampleNumSA = 1:minSampSize
            %   Track the progress
            disp( [num2str(i) 'th parameter and ' num2str(sampleNumSA) ' samples in '...
                      num2str(R) 'th search curve'] );
               
            %   Solve the ODEs in order to generate the simulations
            sol = ode45(@(t,y) hs_ode(t,y, currentTheta(:,sampleNumSA)),[0 timeLength], xInit);
            y = deval(sol,tspan);
            %   Only save the outputs of the system at the time point of
            %   interest
            data(sampleNumSA,:,:,i,R) = y';
        end    % End of the samples in each search curve    
    end  % End of the search curves
end % End of the parameters


%   Close the memory pool for parallel computing
% matlabpool('close')
% Number_Workers = matlabpool('size');


% save('HS_eFAST_withAllTimeSeries.mat','-v7.3');


% Calculate the sensitivity index
Si = efast(data,OmegaMax,M,time_point,1:stateDim);

%   Normalize SA index
for d = 1:size(Si,3)
    for c = 1:size(Si,2)
        Si(:,c,d) = Si(:,c,d)/sum(Si(:,c,d));
    end
end

% Set the SA index at time zero to zero
Si(:,1,1) = [0;0;0;0;0;0;0];
Si(:,1,2) = [0;0;0;0;0;0;0];
Si(:,1,3) = [0;0;0;0;0;0;0];

%   Calculate the overall SA index
sumSi = sum(Si,2)/size(Si,2);

%   Plot the pie chart shown in Supplementary
figure
pie(sumSi(:,1,1))

%   Plot the complete SA in time series
figure,
plot(Si(1,2:end,1),'r--','LineWidth',3);
hold on;
plot(Si(2,2:end,1),'b--','LineWidth',3);
hold on;
plot(Si(3,2:end,1),'m--','LineWidth',3);
hold on;
plot(Si(4,2:end,1),'g--','LineWidth',3);
hold on;
plot(Si(5,2:end,1),'c--','LineWidth',3);
hold on;
plot(Si(6,2:end,1),'c--','LineWidth',3);
hold on;
plot(Si(7,2:end,1),'k--','LineWidth',3);
grid;
xlabel('Time','FontSize',20);
ylabel('SA','FontSize',20);
ylim([-0.1 1]);
title('Sensitivities of parameters for 1^{st} state S_{t}','FontSize',20);
set(gca,'FontSize',20);
legend('k_{d}', '\alpha_{d}', '\alpha_{0}', '\alpha_{s}', 'k_{s}','k_{u}','dummy');

figure,
plot(Si(1,2:end,2),'r--','LineWidth',3);
hold on;
plot(Si(2,2:end,2),'b--','LineWidth',3);
hold on;
plot(Si(3,2:end,2),'m--','LineWidth',3);
hold on;
plot(Si(4,2:end,2),'g--','LineWidth',3);
hold on;
plot(Si(5,2:end,2),'c--','LineWidth',3);
hold on;
plot(Si(6,2:end,2),'c--','LineWidth',3);
hold on;
plot(Si(7,2:end,2),'k--','LineWidth',3);
grid;
xlabel('Time','FontSize',20);
ylabel('SA','FontSize',20);
ylim([-0.1 1]);
title('Sensitivities of parameters for 2^{nd} state D_{t}','FontSize',20);
set(gca,'FontSize',20);
legend('k_{d}', '\alpha_{d}', '\alpha_{0}', '\alpha_{s}', 'k_{s}','k_{u}','dummy');

figure,
plot(Si(1,2:end,3),'r--','LineWidth',3);
hold on;
plot(Si(2,2:end,3),'b--','LineWidth',3);
hold on;
plot(Si(3,2:end,3),'m--','LineWidth',3);
hold on;
plot(Si(4,2:end,3),'g--','LineWidth',3);
hold on;
plot(Si(5,2:end,3),'c--','LineWidth',3);
hold on;
plot(Si(6,2:end,3),'c--','LineWidth',3);
hold on;
plot(Si(7,2:end,3),'k--','LineWidth',3);
grid;
xlabel('Time','FontSize',20);
ylabel('SA','FontSize',20);
ylim([-0.1 1]);
title('Sensitivities of parameters for 3^{rd} state U_{f}','FontSize',20);
set(gca,'FontSize',20);
legend('k_{d}', '\alpha_{d}', '\alpha_{0}', '\alpha_{s}', 'k_{s}','k_{u}','dummy');













