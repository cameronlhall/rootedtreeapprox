%% seirGillespie
%
% Code to run Gillespie algorithm simulations of the SEIR model in the case
% where model parameters are the same at every node/edge.
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
% 2021-08-16 Further polishing
% 2021-07-10 Adapted from SIR code (with Bram's speed improvement) 
% Note: Speed might be improved further by only keeping averages

%%
function [...
    probS, ...              Empirical prob of susceptible (numNodes by numTimes)
    probE, ...              Empirical prob of exposeds (numNodes by numTimes)
    probI, ...              Empirical prob of infectious (numNodes by numTimes)
    probR, ...              Empirical prob of recovereds (numNodes by numTimes)
    numRuns ...             Number of runs of Gillespie model after correction for rounding
    ] = seirGillespie(...
    Adj, ...                Adjacency matrix (numNodes by numNodes)
    params, ...             Model parameters [lambda, phi, mu, nu, gamma, alpha]
    initCondProbs, ...          Initial conditions for all nodes (numNodes by 3) giving probability of initially being in S, E or I state.
    t, ...                  Vector of times (1 by numTimes)
    numRuns ...             Number of runs of Gillespie model
    )

%% Extract key parameters

numTimes = size(t,2);
maxTime = t(numTimes);
timeResolution = t(2);
numNodes = size(Adj,1);
lambda = params(1);
phi = params(2);
mu = params(3);
nu = params(4);
gamma = params(5);

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
    % 1 = exposed
    % 2 = infectious
    % 3 = recovered

    % Default to recovered
    initState = 3*ones(numNodes,1);
    
    % Indicate where inital state is susceptible
    initState(initCondProbs(:,1)==1) = 0;
    
    % Inidcate where initial state is exposed
    initState(initCondProbs(:,2)==1) = 1;

    % Inidcate where initial state is infectious
    initState(initCondProbs(:,3)==1) = 2;

    recalculateICs = false;
    
else
    
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
probEOuter = zeros(numNodes,numTimes,numRunsOuter);
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
            % 1 = exposed
            % 2 = infectious
            % 3 = recovered
            
            % Generate random numbers for each node
            nodeInitStateRand = rand(numNodes,1);
            
            % Default to recovered
            initState = 3*ones(numNodes,1);
            
            % Indicate where initial state is infectious
            initState(nodeInitStateRand<cumInitProbs(:,3)) = 2;
            
            % Indicate where initial state is exposed
            initState(nodeInitStateRand<cumInitProbs(:,2)) = 1;
            
            % Indicate where initial state is susceptible
            initState(nodeInitStateRand<cumInitProbs(:,1)) = 0;
        end

        
        % Initialisations for time and state
        currTime = 0;
        currState = initState;
        tPosition = 1;
        storedStatesInner(:,tPosition,kRunsInner) = currState;

        % Extra details for exposed nodes (useful for making chosing a
        % node efficient)
        expNodes = (initState == 1);
        numExpNodes = sum(expNodes);
        expNodesList = zeros(numNodes,1);
        expNodesList(1:numExpNodes) = find(expNodes);        

        % Extra details for infectious nodes (useful for making chosing a
        % node efficient)
        infNodes = (initState == 2);
        numInfNodes = sum(infNodes);
        infNodesList = zeros(numNodes,1);
        infNodesList(1:numInfNodes) = find(infNodes);        
        
        
        % Initialise matrix of possible links. This should be the adjacency
        % matrix, but any links into exposed, infected or recovered nodes should be
        % removed
        
        % List of all of the nodes that are initially infected or recovered
        initAffNodeLocs = [expNodesList(1:numExpNodes); ...
			infNodesList(1:numInfNodes); ...
			find(initState==3)];
        
        % Construct possLinks as adjacency matrix
        possLinks = Adj;
        
        % Remove any edges going into infected/recovered nodes
        possLinks(initAffNodeLocs,:) = 0;
        
        % Remove any edges going out of recovered nodes
        possLinks(:,find(initState==3)) = 0; %#ok
        
        % Find susceptible nodes that are vulnerable to infection (with
        % weights corresponding to how likely they are to be infected).
        vulNodes = possLinks*infNodes;
        
        % Run Gillespie algorithm
        while currTime < maxTime
            
            %% Time step and event choice
            
            % Old state for storage
            oldState = currState;
            
            % Calculate rates at which events can take place
			allRates = [ ...
				phi*lambda*sum(vulNodes); ... Sus to Exp
				(1-phi)*lambda*sum(vulNodes); ... Sus to Inf
				mu*numExpNodes; ... Exp to Inf
				nu*numExpNodes; ... Exp to Rec
				gamma*numInfNodes ... Inf to Rec
			];
            

            % Calculate time step
            timeStep = exprnd(1/sum(allRates));

            % Cumulative sum of all rates
            cumRates = cumsum(allRates);
                        
            % Event probabilities
            eventProbs = cumRates/cumRates(end);
            
            % Choose a type of event
            eventType = sum(rand<eventProbs);
            % Note: 
            % 5 corresponds to Sus to Exp, 
            % 4 corresponds to Sus to Inf,
            % 3 corresponds to Exp to Inf,
            % 2 corresponds to Exp to Rec,
            % 1 corresponds to Inf to Rec,
            
            switch eventType
                
                case 5
                % Sus to Exp

                % Newly infected node chosen using vulNodes weightings
                newlyExpNode = randsample(numNodes,1,true,vulNodes);
                
                % Update currState, possLinks, expNodes, vulNodes
                % Note, possLinks changes to remove any incoming links to
                % an exposed node since it cannot be infected again.
                currState(newlyExpNode) = 1;
                
                possLinks(newlyExpNode,:) = 0;
                
                expNodes(newlyExpNode) = true;
                expNodesList(numExpNodes+1) = newlyExpNode;
                numExpNodes = numExpNodes+1;
                
                vulNodes = possLinks*infNodes;
                
                
                case 4
                % Sus to Inf
                
                % Newly infected node chosen using vulNodes weightings
                newlyInfNode = randsample(numNodes,1,true,vulNodes);
                
                % Update currState, possLinks, infNodes, vulNodes
                % Note, possLinks changes to remove any incoming links to
                % an infectious node since it cannot be infected again.
                currState(newlyInfNode) = 2;
                
                possLinks(newlyInfNode,:) = 0;
                
                infNodes(newlyInfNode) = true;
                infNodesList(numInfNodes+1) = newlyInfNode;
                numInfNodes = numInfNodes+1;
                
                vulNodes = possLinks*infNodes;
                
                case 3
                % Exp to Inf
                
                % Newly infected node chosen using uniformly from the
                % expNodesList
                chosenNodeInList = randi(numExpNodes);
                newlyInfNode = expNodesList(chosenNodeInList);
                
                % Update currState, expNodes, infNodes, vulNodes
                currState(newlyInfNode) = 2;
                
                expNodes(newlyInfNode) = false;
                expNodesList(chosenNodeInList) = expNodesList(numExpNodes);
                expNodesList(numExpNodes) = 0;
                numExpNodes = numExpNodes-1;

                infNodes(newlyInfNode) = true;
                infNodesList(numInfNodes+1) = newlyInfNode;
                numInfNodes = numInfNodes+1;
                
                vulNodes = possLinks*infNodes;                

                
                case 2
                % Exp to Rec
                
                % Newly recovered node chosen using uniformly from the
                % expNodesList. Note possLinks changes to remove any
                % outgoing links from a recovered node since it cannot
                % infect anyone again.
                chosenNodeInList = randi(numExpNodes);
                newlyRecNode = expNodesList(chosenNodeInList);
                
                % Update currState, possLinks, expNodes (no list of
                % recovered nodes to update).
                currState(newlyRecNode) = 3;

                possLinks(:,newlyRecNode) = 0;
                
                expNodes(newlyRecNode) = false;
                expNodesList(chosenNodeInList) = expNodesList(numExpNodes);
                expNodesList(numExpNodes) = 0;
                numExpNodes = numExpNodes-1;
                
                case 1
                % Inf to Rec
                
                % Newly recovered node chosen using uniformly from the
                % infNodesList. Note possLinks changes to remove any
                % outgoing links from a recovered node since it cannot
                % infect anyone again.
                chosenNodeInList = randi(numInfNodes);
                newlyRecNode = infNodesList(chosenNodeInList);
                
                % Update currState, possLinks, infNodes, vulNodes.
                currState(newlyRecNode) = 3;

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
            if (sum(infNodes)+sum(expNodes)) < 1
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
    probEOuter(:,:,kRunsOuter) = sum(storedStatesInner==1,3)/numRunsInner;
    probIOuter(:,:,kRunsOuter) = sum(storedStatesInner==2,3)/numRunsInner;
    probROuter(:,:,kRunsOuter) = sum(storedStatesInner==3,3)/numRunsInner;
    
end

% Take average of results from outer loop

probS = mean(probSOuter,3);
probE = mean(probEOuter,3);
probI = mean(probIOuter,3);
probR = mean(probROuter,3);


