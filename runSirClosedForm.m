%% runHybridModel.m
%
% Script M-file for generating explisit closed-form solutions for SIR on a
% rooted tree. By editing the parameters in this script, you can get an
% output of the exact solution on a chain as sChain, iChain, and rChain.
%
% 2021-08-16 based on various earlier versions


%% Main parameters

% Number of nodes
numNodes = 10;

% Infection rate
lambda = 1;

% Recovery rate
gamma = 0.1;


%% Time output parameters

% Maximum time for recording data
maxTime = 15;

% Resolution of time in output
timeResolution = 0.01;

% Construct time vector
maxTime = timeResolution*ceil(maxTime/timeResolution);
t = (0:timeResolution:maxTime);
numTimes = numel(t);


%% Explicit solutions
% Calculated for S and I and then R is computed from these

% Initialisation for solutions
sChain = zeros(numNodes+1,numTimes);
iChain = zeros(numNodes+1,numTimes);

% Initialisation for terms inside sum
sumsForSChain = zeros(numNodes,numTimes);
sumsForIChain = zeros(numNodes,numTimes);

% Compute summand (note that n refers to n in the paper and n+2 is used
% because array number starts at 1)
for n = 0:(numNodes-1)
    sumsForSChain(n+2,:) = ((lambda+gamma)*t).^(n)/factorial(n);
    sumsForIChain(n+2,:) = (lambda*t).^(n)/factorial(n);
end

% Compute sums (note that the k+1 term corresponds to a sum up to k-1)
sumsForSChain = cumsum(sumsForSChain,1);
sumsForIChain = cumsum(sumsForIChain,1);
    
% Compute solutions
for k = 0:numNodes
    sChain(k+1,:) = 1 - lambda^k/(lambda+gamma)^k ...
        + lambda^k/(lambda+gamma)^k ...
        * exp(-(lambda+gamma)*t) ...
        .* sumsForSChain(k+1,:);
    iChain(k+1,:) = exp(-gamma*t) ...
        - exp(-(lambda+gamma)*t) ...
        .* sumsForIChain(k+1,:);
end

rChain = 1 - sChain - iChain;