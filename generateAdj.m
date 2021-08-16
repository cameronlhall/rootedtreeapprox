%% generateAdj.m
% Generate adjacency matrix and list of edges
%
% Current options for graphType
% 'complete' = Complete graph
% 'ErdosRenyi' = ER random graph
% 'line' = Chain of nodes
% 'ring' = Ring of nodes where each node is connected to its N nearest
% neighbours on each side
% 'tree' = Random tree generated from a random Prufer sequence
% 'tree with extras' = Random tree with additional edges
% 'lattice' = square lattice
% 'latent space' = Hoff-Raftery spatial graph 
% 'latent gridlike space' = Hoff-Raftery spatial graph where nodes are
% initialised at random points in the cells of a grid

% Corresponding contents of graphParams
% NOTE: [xSize,ySize] only has purpose for displaying graphs except for
% latent space models.
% 'complete' -> {}
% 'ErdosRenyi' -> {pEdge,[xSize,ySize]}
% 'line' -> {}
% 'ring' -> {numNeigh, radius}
% 'tree' -> {[xSize,ySize]}
% 'tree with extras' -> {numExtraEdges, [xSize,ySize]}
% 'lattice' -> {[xSize,ySize]}
% 'latent space' -> {alpha,scaling,[xSize,ySize]}
% 'latent gridlike space' -> {alpha,scaling,[xSize,ySize]}

% Interpretation of these parameters
% numNeigh = number of nearest neighbours connected on either side
% numExtraEdges = number of additional edges added to tree
% pEdge = probability of edge existing between any two nodes
% alpha = log odds of nodes being connected if no distance between them
% scaling = rate of decrease in log odds with Euclidean distance
% xSize = size of lattice (or graph visualisation) in x direction
% ySize = size of lattice (or graph visualisation) in y direction
%
% Original version: 2021-01-20
% Current version: 2021-08-16

function [...
    Adj, ...                Adjacency matrix (numNodes by numNodes - sparse when appropriate)
    edgeArray, ...          List of edges (numDirectedEdges by 2)
    nodePos ...             Node positions (numNodes by 2)
    ] = generateAdj(...
    numNodes, ...           Number of Nodes
    graphType, ...          String to indicate graph type
    graphParams ...         Cell array of any relevant parameters for graph
    )



%% Generate Adj
switch graphType
    
    case 'complete'
        
        % Adjacency matrix with all but self-links
        Adj = ones(numNodes) - eye(numNodes);
        
        % Node locations evenly spaced around a circle
        evenlySpaced = linspace(0,1,numNodes+1)';
        evenlySpaced(end) = [];
        nodePos = [cos(2*pi*evenlySpaced),sin(2*pi*evenlySpaced)];
        
    case 'ErdosRenyi'
        
        % Extract p
        pEdge = graphParams{1};
        
        % For backwards compatibility, no requirement to give [xSize,ySize]
        % for plotting.
        if numel(graphParams) == 1
            nodePosFlag = false;
            nodePos = [];
        else
            nodePosFlag = true;
            xSize = graphParams{2}(1);
            ySize = graphParams{2}(2);
        end
        
%         if pEdge < 2*log(numNodes)/numNodes
%             disp('Warning: Graph unlikely to be connected')
%         end

        % Generate upper triangular matrix of random numbers with zeros on diagonal
        adjRand = triu(rand(numNodes),1);
        
        % Use this to create adjacency matrix. (Note the need to be careful because
        % of the zero entries in adjRand.
        Adj = adjRand>(1-pEdge);
        Adj = Adj + Adj';
        Adj = sparse(Adj);
        
        % Random node locations (could probably calculate something better)
        % nodePos = rand(numNodes,2);
        % nodePos = [rand(numNodes,1)*13,rand(numNodes,1)*7];
        
        % Find node positions from Matlab's plot function and rescale them
        if nodePosFlag
            figure
            h = plot(graph(Adj),'Layout','force');
            nodeXPos = h.XData';
            nodeYPos = h.YData';
            
            % Rescale appropriately
            nodeXPos = nodeXPos - min(nodeXPos);
            nodeYPos = nodeYPos - min(nodeYPos);
            nodeXPos = nodeXPos*(xSize-1)/max(nodeXPos)+1;
            nodeYPos = nodeYPos*(ySize-1)/max(nodeYPos)+1;
            
            nodePos = [nodeXPos nodeYPos];
        end
        
    case 'line'
        
        % Generate adjacency matrix for line
        % N.B. Could do this in one line with spdiags, but this is easier to parse.
        Adj = spdiags(ones(numNodes,1), 1, numNodes,numNodes);
        Adj = Adj + Adj';
        
        % Node locations along a line
        nodePos = [(0:numNodes-1)' zeros(numNodes,1)];
        
    case 'ring'
        
        % For backwards compatibility, make default single-neighbour on
        % ring.
        if numel(graphParams) == 0
            numNeigh = 1;
        else
            numNeigh = graphParams{1};
        end
        
        
        % For backwards compatibility, no requirement to give [xSize,ySize]
        % for plotting.
        if numel(graphParams) == 1
            nodePosFlag = false;
            nodePos = [];
        else
            nodePosFlag = true;
            radius = graphParams{2};
        end
        
        % Generate adjacency matrix for line and then add extra edge
        % N.B. Could do this in one line with spdiags, but this is easier to parse.
        AdjNear = spdiags(ones(numNodes,numNeigh), ...
            1:numNeigh, numNodes,numNodes);
        AdjFar = spdiags(ones(numNodes,numNeigh), ...
            ((numNodes-numNeigh):(numNodes-1)), numNodes,numNodes);
        Adj = AdjNear+AdjFar;
        Adj = Adj + Adj';
        
        if nodePosFlag
            
            % Node locations evenly spaced around a circle (with radius
            % variation for visualisation)
            evenlySpaced = linspace(0,1,numNodes+1)';
            evenlySpaced(end) = [];
            varAmount = 0.05;
            posVar = (1 - varAmount) + 0.05*(-1).^(linspace(1,numNodes,numNodes)');
            nodePos = radius*posVar.*[cos(2*pi*evenlySpaced),sin(2*pi*evenlySpaced)];
            
            figure
            plot(graph(Adj),'XData',nodePos(:,1),'YData',nodePos(:,2));
            
        end
        
    case 'tree'
        
        % For backwards compatibility, no requirement to give [xSize,ySize]
        % for plotting.
        if numel(graphParams) == 0
            nodePosFlag = false;
            nodePos = [];
        else
            nodePosFlag = true;
            xSize = graphParams{1}(1);
            ySize = graphParams{1}(2);
        end
        
        % Generate Pruefer sequence
        pruefer = randsample(numNodes,numNodes-2,true);
        
        % Number of edges to add to each node
        numEdgesToAdd = ones(1,numNodes) + sum(pruefer==(1:numNodes),1);
        
        % Generate empty adjacency matrix
        Adj = zeros(numNodes);
        
        % Pruefer decoding algorithm
        for iSeq = 1:(numNodes-2)
            
            % First node with only one edge to add and its target
            lowestLeaf = find(numEdgesToAdd==1,1);
            prueferVal = pruefer(iSeq);
            
            % Create link between this node and the appropriate number in Pruefer
            % sequence
            Adj(lowestLeaf,prueferVal) = 1;
            
            % Decrease the number of edges needed on each of those nodes
            numEdgesToAdd([lowestLeaf prueferVal]) =...
                numEdgesToAdd([lowestLeaf prueferVal])-1;
        end
        
        % Add final edge
        finalEdge = find(numEdgesToAdd==1);
        Adj(finalEdge(1),finalEdge(2)) = 1;
        
        % Make adjacency matrix symmetric
        Adj = sparse(Adj + Adj');
        
        % Find node positions from Matlab's plot function and rescale them
        if nodePosFlag
            figure
            h = plot(graph(Adj),'Layout','force');
            nodeXPos = h.XData';
            nodeYPos = h.YData';
            
            % Rescale appropriately
            nodeXPos = nodeXPos - min(nodeXPos);
            nodeYPos = nodeYPos - min(nodeYPos);
            nodeXPos = nodeXPos*(xSize-1)/max(nodeXPos)+1;
            nodeYPos = nodeYPos*(ySize-1)/max(nodeYPos)+1;
            
            nodePos = [nodeXPos nodeYPos];
        end

    case 'tree with extras'
        
        % For backwards compatibility, no requirement to give [xSize,ySize]
        % for plotting.
        if numel(graphParams) == 0
            nodePosFlag = false;
            nodePos = [];
        else
            nodePosFlag = true;
            numExtraEdges = graphParams{1};
            xSize = graphParams{2}(1);
            ySize = graphParams{2}(2);
        end
        
        % Generate Pruefer sequence
        pruefer = randsample(numNodes,numNodes-2,true);
        
        % Number of edges to add to each node
        numEdgesToAdd = ones(1,numNodes) + sum(pruefer==(1:numNodes),1);
        
        % Generate empty adjacency matrix
        Adj = zeros(numNodes);
        
        % Pruefer decoding algorithm
        for iSeq = 1:(numNodes-2)
            
            % First node with only one edge to add and its target
            lowestLeaf = find(numEdgesToAdd==1,1);
            prueferVal = pruefer(iSeq);
            
            % Create link between this node and the appropriate number in Pruefer
            % sequence
            Adj(lowestLeaf,prueferVal) = 1;
            
            % Decrease the number of edges needed on each of those nodes
            numEdgesToAdd([lowestLeaf prueferVal]) =...
                numEdgesToAdd([lowestLeaf prueferVal])-1;
        end
        
        % Add final edge
        finalEdge = find(numEdgesToAdd==1);
        Adj(finalEdge(1),finalEdge(2)) = 1;

        % Make adjacency matrix symmetric
        Adj = Adj + Adj';
        
        % Add random edges 
        addedEdges = 0;
        while addedEdges < numExtraEdges
            chosenNodes = ceil(rand(2,1)*numNodes);
            if Adj(chosenNodes(1),chosenNodes(2)) == 0
                Adj(chosenNodes(1),chosenNodes(2)) = 1;
                Adj(chosenNodes(2),chosenNodes(1)) = 1;
                addedEdges = addedEdges+1;
            end
        end
        
        % Store as sparse
        Adj = sparse(Adj);
        
        % Find node positions from Matlab's plot function and rescale them
        if nodePosFlag
            figure
            h = plot(graph(Adj),'Layout','force');
            nodeXPos = h.XData';
            nodeYPos = h.YData';
            
            % Rescale appropriately
            nodeXPos = nodeXPos - min(nodeXPos);
            nodeYPos = nodeYPos - min(nodeYPos);
            nodeXPos = nodeXPos*(xSize-1)/max(nodeXPos)+1;
            nodeYPos = nodeYPos*(ySize-1)/max(nodeYPos)+1;
            
            nodePos = [nodeXPos nodeYPos];
        end
        
    case 'lattice'
        
        % Extract y-dimension 
        xSize = graphParams{1}(1);
        ySize = graphParams{1}(2);
        
        if numNodes ~= xSize * ySize
            error('For lattice type, numNodes and xDim*yDim must be identical');
        end
        
        % Generate grid of connections
        Adj = zeros(numNodes);
        for i=1:numNodes
            % Connections down from everything except last row
            if rem(i,ySize) ~= 0
                Adj(i,i+1) = 1;
            end
            % Connections left from everything except first column
            if i>ySize
                Adj(i,i-ySize) = 1;
            end
        end
        
        Adj = sparse(Adj+Adj');
        
        % Node positions
        nodeIndices = (1:numNodes)';
        nodePos = [floor((nodeIndices-1)/ySize),rem(nodeIndices-1,ySize)]+1;
        
    case 'latent space'
        
        % Extract parameters
        
        % Size of 2D space along one edge. Larger L can decrease the
        % probability of nodes being connected to each other by increasing
        % distance beween them.
        xSize = graphParams{3}(1);
        ySize = graphParams{3}(2);
        
        % Parameter alpha in latent space model (see Hoff, Raftery,
        % Salter-Townshend papers). Larger alpha increases probability of
        % nodes being connected to each other.
        alpha = graphParams{1};
        
        scaling = graphParams{2};
        
        % Node positions
        nodePos = [xSize*rand(numNodes,1) ySize*rand(numNodes,1)]*scaling;
        
        % Distances between nodes
        xDiff = nodePos(:,1) - nodePos(:,1)';
        yDiff = nodePos(:,2) - nodePos(:,2)';
        euclidDist = (xDiff.^2 + yDiff.^2).^(1/2);
        
        % Probability of connection
        probConnect = (1 + exp(-alpha)*exp(euclidDist)).^(-1);
        
        % Delete everything on diagonal and below
        probConnect(tril(true(numNodes))) = 0;
        
        % Generate upper triangular matrix of random numbers with zeros on diagonal
        adjRand = triu(rand(numNodes),1);
        
        % Use this to create adjacency matrix. (Note the need to be careful because
        % of the zero entries in adjRand.
        Adj = adjRand>(1-probConnect);
        Adj = Adj + Adj';
        Adj = sparse(Adj);
        
        nodePos = nodePos/scaling;
        
        % Cosmetic change to fit the original grid
        nodePos(:,1) = nodePos(:,1) - min(nodePos(:,1));
        nodePos(:,2) = nodePos(:,2) - min(nodePos(:,2));
        nodePos(:,1) = nodePos(:,1)*(xSize-1)/max(nodePos(:,1))+1;
        nodePos(:,2) = nodePos(:,2)*(ySize-1)/max(nodePos(:,2))+1;
        
        figure
        plot(graph(Adj),'XData',nodePos(:,1),'YData',nodePos(:,2));
        
        case 'latent gridlike space'
        
        % Extract parameters
        
        % Size of 2D space along one edge. Larger L can decrease the
        % probability of nodes being connected to each other by increasing
        % distance beween them.
        xSize = graphParams{3}(1);
        ySize = graphParams{3}(2);
        
        if numNodes ~= xSize * ySize
            error('For lattice type, numNodes and xDim*yDim must be identical');
        end
        
        % Parameter alpha in latent space model (see Salter-Townshend
        % papers). Larger alpha increases probability of nodes being
        % connected to each other.
        alpha = graphParams{1};
        
        scaling = graphParams{2};
        
        % Node positions (adding an error of less than 1/2 to standard
        % positions
        nodeIndices = (1:numNodes)';
        nodePos = [floor((nodeIndices-1)/ySize),rem(nodeIndices-1,ySize)]+1;
        nodePos = nodePos - 1*(1/2 + rand(numNodes,2));
        nodePos = nodePos*scaling;
                
        % Distances between nodes
        xDiff = nodePos(:,1) - nodePos(:,1)';
        yDiff = nodePos(:,2) - nodePos(:,2)';
        euclidDist = (xDiff.^2 + yDiff.^2).^(1/2);
        
        % Probability of connection
        probConnect = (1 + exp(-alpha)*exp(euclidDist)).^(-1);
        
        % Delete everything on diagonal and below
        probConnect(tril(true(numNodes))) = 0;
        
        % Generate upper triangular matrix of random numbers with zeros on diagonal
        adjRand = triu(rand(numNodes),1);
        
        % Use this to create adjacency matrix. (Note the need to be careful because
        % of the zero entries in adjRand.
        Adj = adjRand>(1-probConnect);
        Adj = Adj + Adj';
        Adj = sparse(Adj);
        
        nodePos = nodePos/scaling;
        
        % Cosmetic change to fit the original grid
        nodePos(:,1) = nodePos(:,1) - min(nodePos(:,1));
        nodePos(:,2) = nodePos(:,2) - min(nodePos(:,2));
        nodePos(:,1) = nodePos(:,1)*(xSize-1)/max(nodePos(:,1))+1;
        nodePos(:,2) = nodePos(:,2)*(ySize-1)/max(nodePos(:,2))+1;
        
        figure
        plot(graph(Adj),'XData',nodePos(:,1),'YData',nodePos(:,2));
        
    otherwise
        
        error('Graph type not recognised')
    
end

%% Use Adj to generate edgeArray

[edgeRows, edgeCols] = find(Adj);
edgeArray = [edgeRows edgeCols];