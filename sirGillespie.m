%% sirGillespie
%
% Code to run Gillespie algorithm simulations of the SIR model in the case
% where model parameters (lambda and gamma) are the same at every
% node/edge.
%
% Input arguments:
%
% Adj
% Adjacency matrix for underlying graph
%
% params
% Vector containing lambda and gamma
%
% initCondProbs
% Each row of initConds gives the probability (between 0 and 1) of a node
% being susceptible (first column) or infected (second column) at time 0.
%
% t
% Vector containing list of times for data storage
%
% numruns
% Number of Gillespie simulations to run and average
%
% Authors: C L Hall, B A Siebert
% Date: 2021-08-16
% History:
% 2021-08-16 Redo with Bram's faster method for storing results each round
% 2021-07-10 Error checking added
% 2020-04-09a Investigating a faster alternative to randsample
% 2020-04-09 Original version
% Note: Speed might be improved by only keeping averages

%%
function [...
    probS, ...              Empirical prob of susceptible (numNodes by numTimes)
    probI, ...              Empirical prob of infecteds (numNodes by numTimes)
    probR, ...              Empirical prob of recovereds (numNodes by numTimes)
    numRuns ...             Number of runs of Gillespie model after correction for rounding
    ] = sirGillespie(...
    Adj, ...                Adjacency matrix (numNodes by numNodes)
    params, ...             Model parameters [lambda, gamma]
    initCondProbs, ...      Initial conditions for all nodes (numNodes by 2)
    t, ...                  Vector of times (1 by numTimes)
    numRuns ...             Number of runs of Gillespie model
    )

%% Extract key parameters

numTimes = size(t,2);
maxTime = t(numTimes);
timeResolution = t(2);
numNodes = size(Adj,1);
lambda = params(1);
gamma = params(2);

%% Convert initial conditions

% Check that initial conditions are plausible (all nonnegative and add to no more than 1)
totalProbability = sum(initCondProbs,2);
if any(totalProbability > 1)
	error('Initial state probabilities cannot add to more than one.')
end
if any(initCondProbs(:) < 0)
	error('Initial state probabilities must be nonnegative')
end

% In which cases are the initial conditions deterministic?
deterministicInitConds = (initCondProbs == 0) | (initCondProbs == 1);

% Check if everything in deterministicInitConds is true (i.e. all initial
% conditions are deterministic not probabilistic)
if ~any(~deterministicInitConds(:))
    
    % initState interpretation:
    % 0 = suceptible
    % 1 = infected
    % 2 = recovered

    % Default to recovered
    initState = 2*ones(numNodes,1);
    
    % Indicate where inital state is susceptible
    initState(initCondProbs(:,1)==1) = 0;
    
    % Inidcate where initial state is infected
    initState(initCondProbs(:,2)==1) = 1;
    
    % No need to generate new initial conditions each round
    recalculateICs = false;
    
else
    
    % Initial conditions need to be randomly resampled each round
    recalculateICs = true;
    
    % Cumulative probabilities for initial condition states (used in
    % recalculate ICs part later)
    cumInitProbs = cumsum(initCondProbs,2);

end

% For doing many many runs without running out of storage.
maxRunsWithoutAverage = 10^3;
if numRuns < maxRunsWithoutAverage
    numRunsOuter = 1;
    numRunsInner = numRuns;
else
    numRunsOuter = ceil(numRuns/maxRunsWithoutAverage);
    numRunsInner = maxRunsWithoutAverage;
end
numRuns = numRunsOuter * numRunsInner;

% State storage for outer loop
probSOuter = zeros(numNodes,numTimes,numRunsOuter);
probIOuter = zeros(numNodes,numTimes,numRunsOuter);
probROuter = zeros(numNodes,numTimes,numRunsOuter);


%% Loop through model runs

for kRunsOuter = 1:numRunsOuter
    
    % Clear stored states for inner loop
    storedStatesInner = zeros(numNodes,numTimes,numRunsInner);
    
    for kRunsInner = 1:numRunsInner
        
        % If necessary, generate initial conditions randomly (treating all
        % nodes as independent)
        if recalculateICs
            
            % initState interpretation: 
            % 0 = suceptible 
            % 1 = infected 
            % 2 = recovered
            
            % Generate random numbers for each node.
            nodeInitStateRand = rand(numNodes,1);
            
            % Default to recovered
            initState = 2*ones(numNodes,1);
            
            % Indicate where initial state is susceptible
            initState(nodeInitStateRand<cumInitProbs(:,1)) = 0;
            
            % Indicate where initial state is infected
            initState(nodeInitStateRand<cumInitProbs(:,2)) = 1;
            
        end
        
        % Initialisations for time and state
        currTime = 0;
        currState = initState;
        tPosition = 1;
        storedStatesInner(:,tPosition,kRunsInner) = currState;

        % Extra details for infected nodes (useful for making chosing a
        % recovering node efficient)
        infNodes = (initState == 1);
        numInfNodes = sum(infNodes);
        infNodesList = zeros(numNodes,1);
        infNodesList(1:numInfNodes) = find(infNodes);        
        
        
        % Initialise matrix of possible links. This should be the adjacency
        % matrix, but any links into infected or recovered nodes should be
        % removed
        
        % List of all of the nodes that are initially infected or recovered
        initAffNodeLocs = [infNodesList(1:numInfNodes); find(initState==2)];
        
        % Construct possLinks as adjacency matrix
        possLinks = Adj;
        
        % Remove any edges going into infected/recovered nodes
        possLinks(initAffNodeLocs,:) = 0;
        
        % Remove any edges going out of recovered nodes
        possLinks(:,find(initState==2)) = 0; %#ok
        
        % Find susceptible nodes that are vulnerable to infection (with
        % weights corresponding to how likely they are to be infected).
        vulNodes = possLinks*infNodes;
        
        % Run Gillespie algorithm
        while currTime < maxTime
            
            %% Time step and event choice
            
            % Old state for storage
            oldState = currState;
            
            % Calculate rates at which events can take place
            rateInf = lambda*sum(vulNodes);
            rateRec = gamma*sum(infNodes);
            rateEvent = rateInf + rateRec;
            
            % Calculate time step
            timeStep = exprnd(1/rateEvent);
            
            % Decide whether event will be a new infection
            if rand < rateInf/rateEvent
                
                % Newly infected node chosen using vulNodes weightings
                newlyInfNode = randsample(numNodes,1,true,vulNodes);
                
                % Update currState, possLinks, infNodes, vulNodes
                % Note, possLinks changes to remove any incoming links to
                % an infected node since it cannot be infected again.
                currState(newlyInfNode) = 1;
                possLinks(newlyInfNode,:) = 0;
                infNodes(newlyInfNode) = true;
                infNodesList(numInfNodes+1) = newlyInfNode;
                numInfNodes = numInfNodes+1;
                vulNodes = possLinks*infNodes;
                
                
            % Otherwise, choose which infected individual recovers
            else
                
                % Newly infected node chosen using uniformly from the
                % infNodesList
                chosenNodeInList = randi(numInfNodes);
                newlyRecNode = infNodesList(chosenNodeInList);
                
                % Update currState, possLinks, infNodes, vulNodes
                % Note possLinks changes to remove any outgoing links from
                % a recovered node since it cannot infect anyone again.
                currState(newlyRecNode) = 2;
                possLinks(:,newlyRecNode) = 0;
                infNodes(newlyRecNode) = false;
                infNodesList(chosenNodeInList) = infNodesList(numInfNodes);
                infNodesList(numInfNodes) = 0;
                numInfNodes = numInfNodes-1;
                vulNodes = possLinks*infNodes;
                
            end
            
            %% Time advance and storage
            
            % Advance time
            currTime = currTime + timeStep;
            
            % Avoid going over the maximum time and fill all the way to the
            % end if maxTime is reached
            if currTime > maxTime
                currTime = maxTime;
                newTPosition = numTimes;
            else
                newTPosition = floor(currTime/timeResolution);
            end
            
            % If appropriate, store results
            if newTPosition > tPosition
                oldStateRepeated = oldState*ones(1,newTPosition-tPosition);
                storedStatesInner(:,(tPosition+1):newTPosition,kRunsInner) ...
                    = oldStateRepeated;
                tPosition = newTPosition;
            end
            
            % Check if current state is steady
            if sum(infNodes) < 1
                currTime = maxTime;
            end
            
        end
        
        % Fill in remaining details if steady state reached
        if tPosition < numTimes
            storedStatesInner(:,tPosition+1:numTimes,kRunsInner) = ...
                currState*ones(1,numTimes-tPosition);
        end
        
    end
    
    % Take average of stored states in inner loop and store them in outer
    % loop
    probSOuter(:,:,kRunsOuter) = sum(storedStatesInner==0,3)/numRunsInner;
    probIOuter(:,:,kRunsOuter) = sum(storedStatesInner==1,3)/numRunsInner;
    probROuter(:,:,kRunsOuter) = sum(storedStatesInner==2,3)/numRunsInner;
    
end

% Take average of results from outer loop
probS = mean(probSOuter,3);
probI = mean(probIOuter,3);
probR = mean(probROuter,3);


