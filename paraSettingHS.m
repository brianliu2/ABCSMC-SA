%% Parameter initializations for SA
global pmax;
global pmin;
global stateDim;
global paraTrue;
global stateTrue;
global paraDim;
%   Dimension definition
nS = 2049;          % the number of the samples for each search curve
nR = 5;                 % the number of the search curves

paraDimSA = 6+1;  % the dimension of unknown parameters, the real unknown 
                                 % parameterss plus one dummy parameter

stateDim = 3;     % the dimension of states in this system                           
                            
M = 4;                  % Maximum number of fourier coefficients, which is always 
                            % chose as 4 or 6

%   data regime setting for SA 
timeLengthSA = 200;
deltaTSA = 0.2;
dataPointsSA = timeLengthSA/deltaTSA;

 %  Boundary set for the SA (scale up process)                           
pmin = [1, %   kd
0.001, %   ad
0.001, %   a0
1, %   as
0.001,  %   ks
0.001,  %   ku
1]';%   dummy parameter

pmax = [5, %   kd
0.1, %   ad
0.1, %   a0
5, %   as
0.1,    %   ks
0.1,    %   ku
6]';%   dummy parameter

% Parameter base lines
kd = 3; 
ad = 0.015;
a0 = 0.03;
as = 3;
ks = 0.05;
ku = 0.0254;
dummy = 1;


%%  Parameter initialization for ABC SMC

%     Setting for the two ABC runs
Options1.M                               = 10;
Options2.M                               = 25;
Options1.alphaSMC                 = 0.99;
Options2.alphaSMC                 = 0.98;
Options1.particleNum               = 1000;
Options2.particleNum               = 2500;
Options1.delta                          = 0.98;
Options2.delta                          = 0.99;

%     Setting for generating data in ABC
paraTrue = [3;0.015;0.03;3;0.05;0.0254];
stateTrue = [0;0;0];
Options1.timeLength = 80;
Options2.timeLength = Options1.timeLength;
deltaT = 0.2;
Options1.samplePoints = Options1.timeLength/deltaT;
Options2.samplePoints = Options1.samplePoints;
Options1.tspan = linspace(0,Options1.timeLength,Options1.samplePoints);
Options2.tspan = Options1.tspan;

%     Tolerance for running ABC 
Options1.epsilonInit = 400;
Options2.epsilonInit = 400;

%   Final target epsilon
Options1.epsilontarget = 200;   
Options2.epsilontarget = 80;

%   Set prior to ABCSMC
priorABCSMC;
paraDim = 6;



















