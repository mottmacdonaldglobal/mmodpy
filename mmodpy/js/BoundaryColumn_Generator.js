
// User input
//
var tolerance = 1;  //tolerance for nodes on boundary
var col_width_fac = 5;  //plan size of column = width of main mesh * col_width_fac
var col_offset_fac = 3; //column offset from main mesh by col_offset_fac*width of main mesh * col_width_fac

//lysmer input
var lysmer_ro = 1.8;
var lysmer_vs = 1000;
var lysmer_vp = 1600;
var SDOrient = 11;   //start id for orientation vectors for dampers (*DEFINE_SD_ORIENTATION)

//end of user input
//-------------------------------------------------

//
var m = Model.Select("Select model 1");
if (!m)
{
	Error("No model!");
	Exit();
}

var flag = AllocateFlag();
m.ClearFlag(flag);

Solid.Select(flag,"Select solids");
m.PropagateFlag(flag);


// Get all nodes in model into a variable
var allNodes = Node.GetFlagged(m, flag);
var allSolids = Solid.GetAll(m);

// User needs to input the position of the boundaries in the variables below
var xmin = 99999999;
var xmax = -99999999;
var ymin = 99999999
var ymax = -99999999;
var zmin = 99999999;
var zmax = -99999999;

for(i=0;i<allNodes.length;i++)
{
    xmin = Math.min(xmin,allNodes[i].x);
    xmax = Math.max(xmax,allNodes[i].x);
    ymin = Math.min(ymin,allNodes[i].y);
    ymax = Math.max(ymax,allNodes[i].y);
    zmin = Math.min(zmin,allNodes[i].z);
    zmax = Math.max(zmax,allNodes[i].z);    
}

var originx = (xmin + xmax)/2;
var originy = (ymin + ymax)/2;

Message("Centre " + originx + ", " + originy);

if(SDOrient>0)
{
    //create vectors
    var SDOrientX = SDOrient;
    var SDOrientY = SDOrient+1;
    var SDOrientZ = SDOrient+2;
    f = new File("./prtmp_setseg.key", File.WRITE);     // Open temporary file
    f.Writeln("*KEYWORD");                       // Write set keyword data to file
    f.Writeln("*DEFINE_SD_ORIENTATION");
    f.Writeln(SDOrientX + ",0,1.0,0.0,0.0");           
    f.Writeln(SDOrientY + ",0,0.0,1.0,0.0");           
    f.Writeln(SDOrientZ + ",0,0.0,0.0,1.0");           
    f.Writeln("*END");
    f.Close();                                  	    // Close temporary file
    m.Import("prtmp_setseg.key");                       // Import data from file into model
    File.Delete("./prtmp_setseg.key");                  // Delete temporary file

    var pid_lys = Part.NextFreeLabel(m);
    var sid_lys = pid_lys;  //Section.NextFreeLabel(m);
    var mid_lys = pid_lys;  //Material.NextFreeLabel(m);

    var mattype = "*MAT_DAMPER_VISCOUS"
    var mat = new Material(m, mid_lys, mattype);
    mat.SetPropertyByName("DC",  1.0            );

    var sec = new Section(m, sid_lys, Section.DISCRETE, "Lysmer dampers")

    var p = new Part(m, pid_lys, sid_lys, mid_lys, "Lysmer dampers");
    
    var setID = Set.NextFreeLabel(m, Set.NODE);
    var Lys_set = new Set(m, setID, Set.NODE);
    var spc = new Spc(m, setID, 0, 1,1,1,0,0,0,Spc.SET);
    
    var ln = LoadNode.GetAll(m);

}
// This radius will be offset used to position the large soil blocks
var blockwidth = (xmax - xmin) * col_width_fac;
var blockStandOffRadius = col_offset_fac * blockwidth;
var setNum = 0;
var blockNum =0;
var baseNodeSetID = Set.NextFreeLabel(m, Set.NODE);
var BaseNode_Set = new Set(m, baseNodeSetID, Set.NODE);

// Adjust search boundaries by the specified tolerance
xmin = xmin + tolerance;
xmax = xmax - tolerance;
ymin = ymin + tolerance;
ymax = ymax - tolerance;
zmin = zmin + tolerance;
zmax = zmax - tolerance;



// Declare arrays to hold boundary node ID numbers
var xminArr = new Array();
var xmaxArr = new Array();
var yminArr = new Array();
var ymaxArr = new Array();

// variable to hold the max node id - used to dimension array later
var maxNodeID = 0

var BoundNodes = new Array();

var NodSol = new Array();

// Construct an empty array with length large enough to hold all nid pointers
var sparseNodeArr = new Array();

//Populate sparse array with index pointers to actual node array which contains the node objects
for(i=0;i<allNodes.length;i++){
    sparseNodeArr[allNodes[i].nid] = i;
    NodSol[i] = new Array();
}

// Cycle through all nodes, store node IDs which are on the boundaries
allNodes.forEach(function(N) {

    // Keep track of the max nid which has been encountered
    if(N.nid > maxNodeID){maxNodeID = N.nid};

    // Store boundary NIDs in the appropriate array
    if(N.x < xmin){
        xminArr.push(N.nid)
        BoundNodes.push(sparseNodeArr[N.nid]);
    }
    else if(N.x > xmax){
        xmaxArr.push(N.nid)
        BoundNodes.push(sparseNodeArr[N.nid]);
    }
    else if(N.y < ymin){
        yminArr.push(N.nid)
        BoundNodes.push(sparseNodeArr[N.nid]);
    }
    else if(N.y > ymax){
        ymaxArr.push(N.nid)
        BoundNodes.push(sparseNodeArr[N.nid]);
    };
    
}, this)


var BoundNodesChk = new Array();

var boundnodint = new Array();

for (i=0; i<BoundNodes.length; i++)
{
    BoundNodesChk[BoundNodes[i]] = 1;
    
}


Message("finding solids on boundaries");

var SolInt = new Array();

for(i=0;i<allSolids.length;i++)
{
    SolInt[allSolids[i].eid] = i;
    
    check = 0;
    
    intn = sparseNodeArr[allSolids[i].n1];
    if( BoundNodesChk[intn] == 1){NodSol[intn].push(allSolids[i].eid);}
    intn = sparseNodeArr[allSolids[i].n2];
    if( BoundNodesChk[intn] == 1){NodSol[intn].push(allSolids[i].eid);}
    intn = sparseNodeArr[allSolids[i].n3];
    if( BoundNodesChk[intn] == 1){NodSol[intn].push(allSolids[i].eid);}
    intn = sparseNodeArr[allSolids[i].n4];
    if( BoundNodesChk[intn] == 1){NodSol[intn].push(allSolids[i].eid);}
    intn = sparseNodeArr[allSolids[i].n5];
    if( BoundNodesChk[intn] == 1){NodSol[intn].push(allSolids[i].eid);}
    intn = sparseNodeArr[allSolids[i].n6];
    if( BoundNodesChk[intn] == 1){NodSol[intn].push(allSolids[i].eid);}
    intn = sparseNodeArr[allSolids[i].n7];
    if( BoundNodesChk[intn] == 1){NodSol[intn].push(allSolids[i].eid);}
    intn = sparseNodeArr[allSolids[i].n8];
    if( BoundNodesChk[intn] == 1){NodSol[intn].push(allSolids[i].eid);}    
    
}

var breakp;
// Sort X boundary arrays by Y coord smallest to largest
xminArr.sort(function(a,b){return allNodes[sparseNodeArr[a]].y - allNodes[sparseNodeArr[b]].y});
xmaxArr.sort(function(a,b){return allNodes[sparseNodeArr[a]].y - allNodes[sparseNodeArr[b]].y});

// Sort Y boundary arrays by X coord smallest to largest
yminArr.sort(function(a,b){return allNodes[sparseNodeArr[a]].x - allNodes[sparseNodeArr[b]].x});
ymaxArr.sort(function(a,b){return allNodes[sparseNodeArr[a]].x - allNodes[sparseNodeArr[b]].x});


buildBlocks(xminArr);
buildBlocks(xmaxArr); 
buildBlocks(yminArr);
buildBlocks(ymaxArr);
View.Redraw();

var baseNodCount = 0;
//Empty existing array
var finishedNodes = new Array();
var finishedNodes = Node.GetAll(m);
// Create a node set containing all nodes on the base of the model
for(i=0;i<finishedNodes.length;i++)
{
    if(finishedNodes[i].z < zmin){
        BaseNode_Set.Add(finishedNodes[i].nid);
        baseNodCount++;
    }  
}

Message("Node set created containing "+baseNodCount+" nodes...")
Message("Base node set id = " + baseNodeSetID);
Message("Script Finished!!");

//----------------------------------------------------------------------------------------------------------------
//  Function List
//----------------------------------------------------------------------------------------------------------------

function buildBlocks(passedBoundArr){

var tmpMainNodeCol = new Array();
var firstcycle = true;
var inCol = false;
var lastColumn = false;
var sizeOfCol = -1;
// Cycle through a boundary array
var counter = 1;
var complete = false;
var lastCount = 0;
passedBoundArr.push(passedBoundArr[0]);

while(counter < passedBoundArr.length){
    var i = counter;
    // if we are starting a new column of nodes - reset the array to empty by removing all current elements
    if(!inCol && !firstcycle){
        tmpMainNodeCol.splice(0, tmpMainNodeCol.length);
        firstcycle = false;
    }
        // get x and y coord of last node
        var lastNx = allNodes[sparseNodeArr[passedBoundArr[i-1]]].x;
        var lastNy = allNodes[sparseNodeArr[passedBoundArr[i-1]]].y;

        // get x and y coord of current node
        var curNx = allNodes[sparseNodeArr[passedBoundArr[i]]].x;
        var curNy = allNodes[sparseNodeArr[passedBoundArr[i]]].y;
        
        // difference in x and y between last coord and this coord
        var dX = Math.abs(curNx - lastNx);
        var dY = Math.abs(curNy - lastNy);                  

        // check if the nodes are in the same x,y position
        if(dX < tolerance && dY < tolerance){
            // nodes are in same column so add last node to the temp list and set flag to in column
            tmpMainNodeCol.push(allNodes[sparseNodeArr[passedBoundArr[i-1]]].nid);
            inCol = true;
        }
        else if(inCol){
            // Otherwise if we were previouly in column data then add the last node to list
            // but switch flag to false to indicate last column list is complete
            
            tmpMainNodeCol.push(allNodes[sparseNodeArr[passedBoundArr[i-1]]].nid);
            inCol = false;
            // report how many nodes are in the current column
            createBlock(tmpMainNodeCol);
            // reset boolean switch to true indicating we are about to start a new column of nodes
            firstcycle = false;
            var nodesToCut = tmpMainNodeCol.length;
                passedBoundArr.splice(0,nodesToCut)
            
            counter = 0;            // reset counter to 1, ready for the next new cycle
        }
        counter += 1;
    }
    
}

function createBlock(passedNodalCol1){
    blockNum += 1;
    var passedNodalCol  =new Array();
    passedNodalCol = passedNodalCol1;
    //sort passed array by z smallest to largest

    
    var break2;
    passedNodalCol.sort(kdCustom);
    function kdCustom(a,b){
        // Place reference
        var dz = allNodes[sparseNodeArr[a]].z - allNodes[sparseNodeArr[b]].z;
        if (dz < 0) {dz = -1}else if(dz > 0){dz = 1}else {dz = 0};
        return dz;
    }
    

    /*
    // ---------------------------------------------------------------------------------------------
    // create a 2D array of two columns, column 1 = nid and column 2 = z level
    var tmp2DArr = new Array();
    for(i=0;i<passedNodalCol.length;i++){
        //
        var rowOfData = new Array();
        var nid = passedNodalCol[i];
        var zLev = allNodes[sparseNodeArr[nid]].z
        rowOfData.push(nid);
        rowOfData.push(zLev);
        tmp2DArr.push(rowOfData);
    }
    // Sort the 2D array by z column lowest to highest
    tmp2DArr.sort(sortMultDimDes)

    // Andy's sort Function
    function sortMultDimDes(a,b){
        var Sort_Val=1;
        return ((a[Sort_Val] < b[Sort_Val]) ? 1: ((a[Sort_Val] > b[Sort_Val]) ?-1 : 0));
    }


    // overwrite original list of node numbers with node number in the correct order
    for(j=0;j<passedNodalCol.length;j++){passedNodalCol[j] =tmp2DArr[j][0]}
    // ---------------------------------------------------------------------------------------------
    */

    //
    var nodeColx = allNodes[sparseNodeArr[passedNodalCol[0]]].x;
    var nodeColy = allNodes[sparseNodeArr[passedNodalCol[0]]].y;

    // Calculate the corner position of new soil column
    var CornerCoord = new Array();
    CornerCoord = getCentreOfBlock(nodeColx, nodeColy, originx, originy);

    // Create a new array which will be used to store the topology
    var BlockElemTopol = new Array();
    var BNode = new Array();

    // For each node in current passed node column, create 4 nodes at the block position
    for(i=0; i < passedNodalCol.length;i++){
        //
        var currentZ = allNodes[sparseNodeArr[passedNodalCol[i]]].z;
        
        var origModelConnect = allNodes[sparseNodeArr[passedNodalCol[i]]].nid;        // Store node id on original model to be connected to
        BNode.push(origModelConnect);
        
        //
        var hbw = 0.5*blockwidth;
        var n1 = new Node(m, Node.NextFreeLabel(m), CornerCoord[0]-hbw, CornerCoord[1]-hbw, currentZ);
        var n2 = new Node(m, Node.NextFreeLabel(m), CornerCoord[0]+hbw, CornerCoord[1]-hbw, currentZ);
        var n3 = new Node(m, Node.NextFreeLabel(m), CornerCoord[0]+hbw, CornerCoord[1]+hbw, currentZ);
        var n4 = new Node(m, Node.NextFreeLabel(m), CornerCoord[0]-hbw, CornerCoord[1]+hbw, currentZ);

        var surfaceArr = new Array();
        surfaceArr.push(n1.nid);
        surfaceArr.push(n2.nid);
        surfaceArr.push(n3.nid);
        surfaceArr.push(n4.nid);
        BlockElemTopol.push(surfaceArr);
        
        if(i > -1){
            //Create new node set and populate with nodes
            var setID = Set.NextFreeLabel(m, Set.NODE);
            var cnrbNode_Set = new Set(m, setID, Set.NODE);
            cnrbNode_Set.Add(origModelConnect);
            cnrbNode_Set.Add(n1.nid);
            cnrbNode_Set.Add(n2.nid);
            cnrbNode_Set.Add(n3.nid);
            cnrbNode_Set.Add(n4.nid);
            
            var CNRB1 = new NodalRigidBody(m,setID);
            CNRB1.spc = true;
            CNRB1.cmo = 1;
            CNRB1.con1 = 0;
            CNRB1.con2 = 7;
        }
        if(i==0){
            if(SDOrient>0)
            {
                var area_lys = hbw*hbw;
                var area_vs = area_lys*lysmer_ro*lysmer_vs;
                var area_vp = area_lys*lysmer_ro*lysmer_vp;
                var lysn1 = new Node(m, Node.NextFreeLabel(m), CornerCoord[0]-hbw, CornerCoord[1]-hbw, currentZ-1);
                var lysn2 = new Node(m, Node.NextFreeLabel(m), CornerCoord[0]+hbw, CornerCoord[1]-hbw, currentZ-1);
                var lysn3 = new Node(m, Node.NextFreeLabel(m), CornerCoord[0]+hbw, CornerCoord[1]+hbw, currentZ-1);
                var lysn4 = new Node(m, Node.NextFreeLabel(m), CornerCoord[0]-hbw, CornerCoord[1]+hbw, currentZ-1);            
                
                var sp1 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n1.nid, lysn1.nid, SDOrientX, area_vs);
                var sp2 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n1.nid, lysn1.nid, SDOrientY, area_vs);
                var sp3 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n1.nid, lysn1.nid, SDOrientZ, area_vp);
                
                var sp1 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n2.nid, lysn2.nid, SDOrientX, area_vs);
                var sp2 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n2.nid, lysn2.nid, SDOrientY, area_vs);
                var sp3 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n2.nid, lysn2.nid, SDOrientZ, area_vp);

                var sp1 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n3.nid, lysn3.nid, SDOrientX, area_vs);
                var sp2 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n3.nid, lysn3.nid, SDOrientY, area_vs);
                var sp3 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n3.nid, lysn3.nid, SDOrientZ, area_vp);

                var sp1 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n4.nid, lysn4.nid, SDOrientX, area_vs);
                var sp2 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n4.nid, lysn4.nid, SDOrientY, area_vs);
                var sp3 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, n4.nid, lysn4.nid, SDOrientZ, area_vp);            
            
                Lys_set.Add(lysn1.nid);
                Lys_set.Add(lysn2.nid);
                Lys_set.Add(lysn3.nid);
                Lys_set.Add(lysn4.nid);
                
                var lclysX = 0;
                var lclysY = 0;
                var lclysZ = 0;
                
                for(iln=0;iln<ln.length;iln++)
                {
                    if(ln[iln].nid==origModelConnect)
                    {
                        if(ln[iln].dof == 1) {lclysX=ln[iln].lcid;}
                        if(ln[iln].dof == 2) {lclysY=ln[iln].lcid;}
                        if(ln[iln].dof == 3) {lclysZ=ln[iln].lcid;}
                        if(lclysX > 0 && lclysY > 0 && lclysZ >0){break;}
                    }
                }
                
                if (lclysX > 0)
                {
                    var l = new LoadNode(m, LoadNode.POINT, n1.nid, 1, lclysX, area_vs);
                    var l = new LoadNode(m, LoadNode.POINT, n2.nid, 1, lclysX, area_vs);
                    var l = new LoadNode(m, LoadNode.POINT, n3.nid, 1, lclysX, area_vs);
                    var l = new LoadNode(m, LoadNode.POINT, n4.nid, 1, lclysX, area_vs);
                }
                if (lclysY > 0)
                {
                    var l = new LoadNode(m, LoadNode.POINT, n1.nid, 2, lclysY, area_vs);
                    var l = new LoadNode(m, LoadNode.POINT, n2.nid, 2, lclysY, area_vs);
                    var l = new LoadNode(m, LoadNode.POINT, n3.nid, 2, lclysY, area_vs);
                    var l = new LoadNode(m, LoadNode.POINT, n4.nid, 2, lclysY, area_vs);
                }
                if (lclysZ > 0)
                {
                    var l = new LoadNode(m, LoadNode.POINT, n1.nid, 3, lclysZ, area_vp);
                    var l = new LoadNode(m, LoadNode.POINT, n2.nid, 3, lclysZ, area_vp);
                    var l = new LoadNode(m, LoadNode.POINT, n3.nid, 3, lclysZ, area_vp);
                    var l = new LoadNode(m, LoadNode.POINT, n4.nid, 3, lclysZ, area_vp);
                }
            
            }
        }
    }
    
    // create elements from newly created nodes 
    for(j=0;j<(BlockElemTopol.length-1); j++){
        var n1_ID = BlockElemTopol[j][0];
        var n2_ID = BlockElemTopol[j][1];
        var n3_ID = BlockElemTopol[j][2];
        var n4_ID = BlockElemTopol[j][3];

        var n5_ID = BlockElemTopol[j+1][0];
        var n6_ID = BlockElemTopol[j+1][1];
        var n7_ID = BlockElemTopol[j+1][2];
        var n8_ID = BlockElemTopol[j+1][3];
        var curEID = Solid.NextFreeLabel(m);
        
        
        PID = Get_Ele_PID(BNode[j], BNode[j+1]); 
        
        var currentElem = new Solid(m,curEID, PID, n1_ID,n2_ID,n3_ID,n4_ID,n5_ID, n6_ID, n7_ID, n8_ID);
        
    }
    Message("Block_" + blockNum + " complete!!!");

    //----------------------
    //  Internal Functions
    //----------------------
    function getCentreOfBlock(nodeColx, nodeColy, ox, oy){
        //
        
        var dx = nodeColx - ox;
        var dy = nodeColy - oy;

        var angleToNodeCol = Math.atan2(dy, dx);    // Get polar coord angle
        var temptestingAngle = angleToNodeCol/(Math.PI/180);    // Display angle in degrees for testing only

        var blockCorner_Y = oy + Math.sin(angleToNodeCol) * blockStandOffRadius;
        var blockCorner_X = ox + Math.cos(angleToNodeCol) * blockStandOffRadius;
        var newReturnCoord = new Array(blockCorner_X, blockCorner_Y);
        return newReturnCoord;
    }
}



function Get_Ele_PID(nn10, nn20)
{
    var chk1 = 0;
    PID=0;
    
    nn1 = sparseNodeArr[nn10];
    nn2 = sparseNodeArr[nn20];
    
    for(nn=0; nn < NodSol[nn1].length; nn++)
    {
        for(mm=0;mm < NodSol[nn2].length; mm++)
        {
            if(NodSol[nn1][nn] == NodSol[nn2][mm])
            {
                chk1=NodSol[nn1][nn];
                break;
            }
        }
        if(chk1==1) break;
    }
    
    if (chk1 > 0) PID = allSolids[SolInt[chk1]].pid;
    return PID;
}


