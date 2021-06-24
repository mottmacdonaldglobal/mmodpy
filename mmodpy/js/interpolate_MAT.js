//---------------------------------------------------------------
// Developed by Kevin Stanton
//
// This script assigns material cards to elements in a 3D soil
// model based on the material cards defined in the nearest part
// with the same title from an index soil column model. The index
// model may contain any number of soil columns and must be
// opened in Oasys Primer after first opening the soil block 
// model caontianing the full 3D geometry and appropriate part
// titles.
//---------------------------------------------------------------

var soilBlockModel = Model.Select("Select Soil Block Model");
var matIndexModel = Model.Select("Select Material Index Model");
var resultantModel = soilBlockModel;
var parts = Part.GetAll(soilBlockModel);
var units = new Array;
for (var i=0;i<parts.length;i++){
	units.push([i+1,parts[i].heading])
}

var matIndexArray = new Array;

Solid.ForEach(matIndexModel, buildMatIndexArray);

function buildMatIndexArray(s){
	var elNum = s.eid
	var partNum = s.pid
	var partObj = Part.GetFromID(matIndexModel,partNum);
	var matNum = partObj.mid
	var matObj = Material.GetFromID(soilBlockModel,matNum);
	var typeIndex = 99999999
	var secNum = 1
	var title = partObj.heading
	var matAssignedPart
	var currIndex = 0
	var rangeStart = 0
	var rangeEnd = 99999999 
	
	for (var i=0;i<units.length;i++){
		if (title == units[i][1]) {
			rangeStart=1000*(i+1)
			rangeEnd=rangeStart+999
			currIndex=rangeStart
			while (Part.GetFromID(soilBlockModel,currIndex) != null) {currIndex++}			
			if (currIndex>rangeEnd) {Message ("***********ERROR*********")}
			var matAssignedPart = new Part(resultantModel, currIndex, secNum, matNum, title + " (Mat " + matNum + " from Index Model)");
			typeIndex = rangeStart
		}
	}	
	
	var centroidX
	var centroidY
	var centroidZ

	[centroidX,centroidY,centroidZ] = getSolidCentroid(s,matIndexModel)
		matIndexArray.push([elNum,title,matNum,matAssignedPart,centroidX,centroidY,centroidZ,typeIndex])
}

Solid.ForEach(resultantModel,assignMatPart);

function assignMatPart(currentSolid){
	
	var centroidX
	var centroidY
	var centroidZ

	[centroidX,centroidY,centroidZ] = getSolidCentroid(currentSolid,resultantModel)

	var newpartNumtoUpdate = findNearestIndexPart(currentSolid,centroidX,centroidY,centroidZ)

	currentSolid.pid=newpartNumtoUpdate
}

function getSolidCentroid(currentSolid,currentModel){
	
	var numNodes = currentSolid.nodes
	var sumX=0
	var sumY=0
	var sumZ=0
	
	for (var i=1;i<=numNodes;i++) {
		var nodeIndexToExtract = "n" + i.toString()
		var nodeToExtract = currentSolid[nodeIndexToExtract.toString()]
		var currNode = Node.GetFromID(currentModel, nodeToExtract);
		sumX=sumX+currNode.x
		sumY=sumY+currNode.y
		sumZ=sumZ+currNode.z
	}
	
	var centX=sumX/numNodes
	var centY=sumY/numNodes
	var centZ=sumZ/numNodes

	return [centX,centY,centZ]
}

function findNearestIndexPart(currSol,cX,cY,cZ){

	var totDist = 9999999999999
	var newPartNum = 9999999999
	var elNum = currSol.eid
	var partNum = currSol.pid
	var blockPartObj = Part.GetFromID(soilBlockModel,partNum)
	var blockPartTitle = blockPartObj.heading

	for (var i=0;i<matIndexArray.length;i++){
		var valX=matIndexArray[i][4]
		var valY=matIndexArray[i][5]
		var valZ=matIndexArray[i][6]
		
		var distX = cX-valX
		var distY = cY-valY
		var distZ = cZ-valZ
		
		var step_one =Math.pow(distX,2)+Math.pow(distY,2)+Math.pow(distZ,2)
		var calcDist = Math.pow(step_one,0.5)
		
		var partID = matIndexArray[i][3].pid
		var indexPartTitle = matIndexArray[i][1]
		
		if (calcDist<=totDist &&  indexPartTitle==blockPartTitle){
			totDist=calcDist
			newPartNum = partID
		}
	}
return 	newPartNum
}

