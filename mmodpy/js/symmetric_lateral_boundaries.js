// Define model
var m = Model.GetFromID(1);

// Determine the boundaries of the soil domain
var edgeFlag = AllocateFlag();
Node.Select(edgeFlag, "Select 4 nodes to define the soil domain lateral extents", m);
var boundaryNodes = Node.GetFlagged(m, edgeFlag);
var xMin = 1e9
var xMax = -1e9
var yMin = 1e9
var yMax = -1e9
for (var i = 0; i<boundaryNodes.length; i++){
    var x_i = boundaryNodes[i].x
    var y_i = boundaryNodes[i].y
    if (x_i < xMin){
        xMin = x_i
        }
    if (y_i < yMin){
        yMin = y_i
        }
    if (x_i > xMax){
        xMax = x_i
        }
    if (y_i > yMax){
        yMax = y_i
        }
    }

// Create arrray of z coord, set number pairs
zArray = new Array

// Set tolerance for z coordinate search
var tol = 0.01

// Loop through all nodes and create sets based on z coordinate
Node.ForEach(m, zSets)
function zSets(currentNode){
	var setTrigger = 0
	// Check if node is part of soil domain boundary and add to zArray if proper set already exists
	for (var i = 0; i<zArray.length; i++){
        if (currentNode.x == xMax || currentNode.x==xMin || currentNode.y==yMax || currentNode.y==yMin){
            if (currentNode.z>=zArray[i][0]-tol && currentNode.z<=zArray[i][0]+tol){
                var currentSet=Set.GetFromID(m, zArray[i][1], Set.NODE)
                currentSet.Add(currentNode.nid);
                setTrigger = 1;
                }
            }
        }

	// Create new set if node is part of soil domain and proper nodal z set does not exist in zArray
	if (setTrigger == 0){
        if (currentNode.x == xMax || currentNode.x==xMin || currentNode.y==yMax || currentNode.y==yMin){
            zArray.push([currentNode.z,Set.NextFreeLabel(m, Set.NODE)]);
            var newSet = new Set(m,Set.NextFreeLabel(m, Set.NODE), Set.NODE);
            newSet.Add(currentNode.nid);
            }
        }
    }

// Create NRBs for all sets in zArray
for (var i = 0; i<zArray.length; i++){
	var nrb = new NodalRigidBody(m, zArray[i][1], NodalRigidBody.NextFreeLabel(m))
    }

Message("Script complete. Check for errors.")