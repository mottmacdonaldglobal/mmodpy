//---------------------------------------------------------------
// Developed by Kevin Stanton
//
// This script generates lysmer dampers and load cards for the 
// application of seismic loads via a non-reflecting boundary at 
// the base of an SSI model. Mesh must be rectilinear in plan.
// If rectilinear mesh not possible, try layer_create_full.js.
//---------------------------------------------------------------
//
// User input
//
var tolerance = 1;  //tolerance for nodes on boundary
var col_width_fac = 5;  //plan size of column = width of main mesh * col_width_fac
var col_offset_fac = 3; //column offset from main mesh by col_offset_fac*width of main mesh * col_width_fac

//lysmer input
var lysmer_ro = 4.04;
var lysmer_vs = 1530;
var lysmer_vp = 3000;
var SDOrient = 11;   //start id for orientation vectors for dampers (*DEFINE_SD_ORIENTATION)
var lclysX = 11; //LCIDs for ground motion velocity time histories
var lclysY = 12;
var lclysZ = 13;
		
//end of user input
//-------------------------------------------------

var m = Model.Select("Select model 1");
if (!m){
	Error("No model!");
	Exit();}

// Initialize other variables
var allNodes = Node.GetAll(m);
var baseNodes = new Array;
var xmin = 99999999;
var xmax = -99999999;
var ymin = 99999999
var ymax = -99999999;
var zmin = 99999999;
var zmax = -99999999;

for (var i=0;i<allNodes.length;i++){
    xmin = Math.min(xmin,allNodes[i].x);
    xmax = Math.max(xmax,allNodes[i].x);
    ymin = Math.min(ymin,allNodes[i].y);
    ymax = Math.max(ymax,allNodes[i].y);
    zmin = Math.min(zmin,allNodes[i].z);
    zmax = Math.max(zmax,allNodes[i].z);    
}

for (var i=0;i<allNodes.length;i++){
	if (allNodes[i].z == zmin){
		baseNodes.push(allNodes[i])}  
}

if (SDOrient>0){
    //create vectors
    var SDOrientX = SDOrient;
    var SDOrientY = SDOrient+1;
    var SDOrientZ = SDOrient+2;
    var f = new File("./prtmp_setseg.key", File.WRITE);     // Open temporary file
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
    var spc = new Spc(m, setID, 0, 1,1,1,0,0,0,Spc.SET);}

for (var i=0;i<baseNodes.length;i++){
				
		var area_lys = trib_area(baseNodes[i]);
		var area_vs = area_lys*lysmer_ro*lysmer_vs;
		var area_vp = area_lys*lysmer_ro*lysmer_vp;
		var lysn = new Node(m, Node.NextFreeLabel(m), baseNodes[i].x, baseNodes[i].y, zmin-1);          
		
		var sp1 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, baseNodes[i].nid, lysn.nid, SDOrientX, area_vs);
		var sp2 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, baseNodes[i].nid, lysn.nid, SDOrientY, area_vs);
		var sp3 = new Discrete(m, Discrete.NextFreeLabel(m), pid_lys, baseNodes[i].nid, lysn.nid, SDOrientZ, area_vp);         
	
		Lys_set.Add(lysn.nid);
		
		if (lclysX > 0){
			var l = new LoadNode(m, LoadNode.POINT, baseNodes[i].nid, 1, lclysX, area_vs);}
		if (lclysY > 0){
			var l = new LoadNode(m, LoadNode.POINT, baseNodes[i].nid, 2, lclysY, area_vs);}
		if (lclysZ > 0){
			var l = new LoadNode(m, LoadNode.POINT, baseNodes[i].nid, 3, lclysZ, area_vp);}
	}

function trib_area(n){
	
	var xPos = xmax - n.x
	var xNeg = n.x - xmin
	var yPos = ymax - n.y
	var yNeg = n.y - ymin
	
	for (var t=0;t<baseNodes.length;t++){
		if (baseNodes[t].x > n.x){
			xPos = Math.min(xPos,baseNodes[t].x-n.x)};
		if (baseNodes[t].x < n.x){
			xNeg = Math.min(xNeg,n.x-baseNodes[t].x)};
		if (baseNodes[t].y > n.y){
			yPos = Math.min(yPos,baseNodes[t].y-n.y)};
		if (baseNodes[t].y < n.y){
			yNeg = Math.min(yNeg,n.y-baseNodes[t].y)};
	} 
	
	var xTrib = xPos/2 + xNeg/2
	var yTrib = yPos/2 + yNeg/2
	var tribArea = xTrib*yTrib
	
	return tribArea
}

Message("Finished running create_lysmer_dampers.js");
