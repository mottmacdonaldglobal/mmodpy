// Define model
var m = Model.GetFromID(1);

// Create arrray of z coord, set number pairs
zArray = new Array

// Set tolerance for z coordinate search
var tol = 0.01

// Loop through all nodes and create sets based on z coordinate
Node.ForEach(m, zSets)
function zSets(currentNode){
	var setTrigger = 0
	// Check if nodal z already exists in the zArray and add if true
	for (var i = 0;i<zArray.length;i++){
		if (currentNode.z>=zArray[i][0]-tol && currentNode.z<=zArray[i][0]+tol){
			var currentSet=Set.GetFromID(m, zArray[i][1], Set.NODE)
			currentSet.Add(currentNode.nid);
			setTrigger = 1;
			}
	}

	// Create new set if nodal z does not exist in zArray
	if (setTrigger == 0){
		zArray.push([currentNode.z,Set.NextFreeLabel(m, Set.NODE)]);
		var newSet = new Set(m,Set.NextFreeLabel(m, Set.NODE), Set.NODE);
		newSet.Add(currentNode.nid);
		}
}

// Create NRBs for all sets in zArray
for (var i = 0;i<zArray.length;i++){
	var nrb = new NodalRigidBody(m, zArray[i][1], NodalRigidBody.NextFreeLabel(m))
}

