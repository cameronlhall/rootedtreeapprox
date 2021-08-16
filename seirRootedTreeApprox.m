%% seirRootedTreeApprox.m
% 
% Code to solve differential equations of SEIR rooted tree approximation in
% the case where model parameters are the same at every node/edge.
%
% Input arguments:
%
% edgeArray
% Each row of edgeArray corresponds to an edge of the graph where the node
% indicated by the first column can spread disease to the node indicated by
% the second column. For an undirected graph, there are two rows of
% edgeArray for every edge. edgeArray can be generated from the adjacency
% matrix using find.
%
% params
% Vector of [lambda, phi, mu, nu, gamma]
%
% initConds
% Each row of initConds gives the probability (between 0 and 1) of a node
% being susceptible (first column) or infected (second column) at time 0.
%
% t
% Vector containing list of times for plotting
%
% Authors: C L Hall, B A Siebert
% Date: 2021-08-16
% Note: Minor modifications of 2021-07-10 version (seirHybridModel.m),
% based earlier on sirHybridModel.m

%%
function [...
    sSol, ...               Solution for susceptibles (numNodes by numTimes)
    eSol, ...               Solution for exposeds (numNodes by numTimes)
    iSol, ...               Solution for infectious (numNodes by numTimes)
    rSol ...                Solution for recovereds (numNodes by numTimes)
    ] = seirRootedTreeApprox(...
    edgeArray, ...          List of edges (numDirectedEdges by 2)
    params, ...             Model parameters [lambda, phi, mu, nu, gamma]
    initConds, ...          Initial conditions for all nodes (numNodes by 3)
    t ...                   Vector of times (1 by numTimes)
    )


%% Extracting additional parameters

% Number of nodes in system
numNodes = size(initConds,1);

% Initial conditions
sInit = initConds(:,1);
eInit = initConds(:,2);
y0 = [sInit; eInit; initConds(:,3)];

% Maximum time
maxTime = t(end);

%% Solve differential equations

% Solve sirSpecDE using ODE45
sol = ode45(@(tForFn,y) ...
    seirSpecDE(tForFn,y,params,edgeArray(:,1),edgeArray(:,2),sInit,eInit), ...
    [0 maxTime], y0);

% Evaluate solution at specified times
ySol = deval(sol,t);

% Convert solution to solutions for susceptibles, infecteds, recovereds
sSol = ySol(1:numNodes,:);
eSol = ySol((numNodes+1):2*numNodes,:);
iSol = ySol((2*numNodes+1):3*numNodes,:);
rSol = 1 - sSol - eSol - iSol;

end

%% Nested function to specify differential equation.
function dydt = seirSpecDE(~,y,params,mainNode,neighbourNode,sInit,eInit)

% Note:
% (mainNode,neighbourNode) corresponds to the (i,j) coefficient of one of
% the 1s in the adjacency matrix. They are used in treeModelRates in order
% to avoid doing too many unnecessary calculations of things that will be
% multiplied by zero.

%% Preliminaries

% Parameter extraction
numNodes = numel(y)/3;
lambda = params(1);
phi = params(2);
mu = params(3);
nu = params(4);
gamma = params(5);

% Initialisation
dydt = zeros(size(y));
treeModelRates = zeros(numNodes);
s = y(1:numNodes);
e = y((numNodes+1):2*numNodes);
i = y((2*numNodes+1):3*numNodes);

%% Calculate rates

% In order to calculate the rate of infection at each node, we construct
% matrices (initialised above) where the (i,j) element corresponds to the
% rate of infection at i due to a neighbour at j. These matrices are only
% nonzero where there is an edge between i and j, and so the calculations
% are only done where [i,j] = [mainNode,neighbourNode], based on the
% nonzero elements of the adjacency matrix.

% adjIndex is the index (in single number form) of the nonzero elements of
% the adjancency matrix.
adjIndex = (neighbourNode-1)*numNodes+mainNode;

% Infection rate according to the model that is exact for trees.
treeModelRates(adjIndex) = ...
    max(...
	(lambda+gamma)*s(mainNode) ...
	- lambda*(1-phi)*s(neighbourNode) ...
	- lambda*mu/(mu+nu)*(phi*s(neighbourNode) + e(neighbourNode)) ...
	- gamma*sInit(mainNode) ...
	- lambda*nu/(mu+nu)*(phi*sInit(neighbourNode) + eInit(neighbourNode)) ...
	,0);

% Infection rate adding up over nodes
totalInfRate = sum(treeModelRates,2);

% Rate of change of susceptible probability is -totalInfRate.
% Rate of change of exposed probability is phi*totalInfRate - mu*e - nu*e.
% Rate of change of infected probability is +totalInfRate - recoveryRate.
dydt(1:numNodes) = -totalInfRate;
dydt((numNodes+1):2*numNodes) = phi*totalInfRate - (mu+nu)*e;
dydt((2*numNodes+1):3*numNodes) = (1-phi)*totalInfRate + mu*e - gamma*i;

end
