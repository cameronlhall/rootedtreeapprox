%% runSeirGillespieModel.m
%
% Script M-file for running Gillespie simulations of SEIR models in case
% where parameters are constant (same at every node/edge). By editing the
% parameters in this script, you can get an output of the average of a
% given number of Gillespie simulations as probS, probE, probI, probR.
%
% 2021-08-16 based on various earlier versions


%% Graph construction
% Number of nodes
numNodes = 100;

% Flag to indicate whether to generate a new adjacency matrix even if Adj
% is already present in the workspace. If this is false, then an existing
% or imported matrix called Adj will be used as the adjacency matrix.
% (Useful for doing multiple runs with the same randomly generated graph.)
replaceAdj = false;

% Type of graph to generate (if generating a new adjacency matrix)
% Note that a new adjacency matrix will always be generated if Adj is
% currently empty. See generateAdj.m for details of graphParams
graphType = 'ErdosRenyi';
graphParams = {0.05};

% Generate adjacency matrix
if replaceAdj || exist('Adj','var')==0
    [Adj,edgeArray] = generateAdj(numNodes,graphType,graphParams);
else
    % Otherwise, ignore the number of nodes in numNodes and regenerate
    % edgeArray (in case it has not yet been generated)
    numNodes = size(Adj,1);
    [edgeRows, edgeCols] = find(Adj);
    edgeArray = [edgeRows edgeCols];
end


%% Contagion dynamic parameters

% Infection rate
lambda = 1;

% Probability of going to E state on infection
phi = 0.8;

% Rate of E to I transition
mu = 1.2;

% Rate of E to R transition
nu = 0.05;

% Recovery rate
gamma = 0.1;



% Parameters
params = [lambda phi mu nu gamma];


%% Time output parameters

% Maximum time for recording data
maxTime = 12;

% Resolution of time in output
timeResolution = 0.01;

% Construct time vector
maxTime = timeResolution*ceil(maxTime/timeResolution);
t = (0:timeResolution:maxTime);
numTimes = numel(t);


%% Initial conditions

% Standard initial conditions where first node is exposed/infected and all
% others are susceptible.
s0 = ones(numNodes,1);
s0(1) = 0;
e0 = zeros(numNodes,1);
e0(1) = 0;
i0 = zeros(numNodes,1);
i0(1) = 1-e0(1);

initConds = [s0 e0 i0];


%% Gillespie model parameters

% Number of runs of Gillespie model
numRuns = 10^5;


%% Run Gillespie algorithm

[probS,probE,probI,probR,numRuns] = ...
    seirGillespie(...
    Adj, ...                Adjacency matrix (numNodes by numNodes)
    params, ...             Model parameters [lambda, gamma, alpha]
    initConds, ...          Initial conditions for all nodes (numNodes by 2)
    t, ...                  Vector of times (1 by numTimes)
    numRuns ...             Number of runs of Gillespie model
    );
