%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Run ABC-SMC on the repressilator model as a function.
%       User needs to define the prior distribution for generating samples
%       for parameters, integer factor M, number of samples used N,
%       weighting factor alphaSMC and discount factor delta.
%
%
%       Author: Xin Liu
%       Data:    15/03/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
echo off
clc
global paraDim;
global runCount;
global sensLab;

%   Define the number of processor used for parallel computing
Num_processor = 4;

%   Activate the parallel session
matlabpool('open', Num_processor);

%%  First ABC run
%      Initializations of prior distributions
paraSMC_Init = mu(:,ones(1,Options1.particleNum)) + ...
                           sqrt(sigma(:,ones(1,Options1.particleNum))).*...
                           randn(paraDim,Options1.particleNum);

%   Initialize the counter of ABC being performed
runCount = 1;

%   First run of ABC-SMC
[epsilon_hist,paraSpace_1] = ABCSMC_Fcn(paraSMC_Init,Options1);


%%  Second ABC run

%     Turn the flag variable to indicate the ABC is going to run second time
runCount = runCount+1;

%   Reset the initial population according to stiff/sloppy
paraSMC_Init = zeros(paraDim,Options2.particleNum);
for i = 1:paraDim
    if strcmp(sensLab{i},'stiff')
        paraSMC_Init(i,:) = mu(i)*ones(1,Options2.particleNum) + sqrt(sigma(i))*randn(1,Options2.particleNum);
    else
        paraSMC_Init(i,:) = mean(paraSpace_1(i,:,end))*ones(1,Options2.particleNum);
    end
end

%   First run of ABC-SMC
[epsilon_hist_2,paraSpace_2] = ABCSMC_Fcn(paraSMC_Init,Options2);


%%  Save the inference data, closing memory pool and plot figures

save('ABCSMC+SA_HS_Auto.mat');
% % %   Close the memory pool for parallel computing
matlabpool('close')
Number_Workers = matlabpool('size');

%%      Ploting the figure
figure,
subplot(2,2,1)
hist(paraSpace_2(1,:,end),40)
grid;
xlabel('\alpha_{0}, true value is 1','FontSize',20);
ylabel('PDF','FontSize',20);
title('Estimation of \alpha_{0}','FontSize',20);
set(gca,'FontSize',20);


subplot(2,2,2)
hist(paraSpace_2(2,:,end),40)
grid;
xlabel('n, true value is 2','FontSize',20);
ylabel('PDF','FontSize',20);
title('Estimation of n','FontSize',20);
set(gca,'FontSize',20);

subplot(2,2,3)
hist(paraSpace_1(3,:,end),40)
grid;
xlabel('\beta, true value is 5','FontSize',20);
ylabel('PDF','FontSize',20);
title('Estimation of \beta','FontSize',20);
set(gca,'FontSize',20);

subplot(2,2,4)
hist(paraSpace_1(4,:,end),40)
grid;
xlabel('\alpha, true value is 1000','FontSize',20);
ylabel('PDF','FontSize',20);
title('Estimation of \alpha','FontSize',20);
set(gca,'FontSize',20);










