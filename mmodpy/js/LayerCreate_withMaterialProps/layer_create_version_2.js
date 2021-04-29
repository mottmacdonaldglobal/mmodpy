// Script to read a csv file to create soil layers     
//
// name:Soil Layers   
// description:Extrude shells to make layered solids, properties from csv file
//   

// DESCRIPTION
// ===========
//
// This script reads a csv file containing data for soil layers. Each row in the
// csv file will create one part with material whose type and properties are 
// given in the csv file. The depth or coordinate of the layer are also given.  
// Layers must be in order, with the uppermost layer first.
//
// The script also extrudes existing shell parts to create the solid elements.
// The shells should be in a horizontal plane at the soil surface. 
// The normal vectors should be facing upwards.
//
// Alternatively, the script can create a column model (1x1 element in plan view) 
// with NRB at each layer if required
//
// Optionally, a node set and *BOUNDARY_PRESC_MOTION or SPC can be created for each
// layer of nodes; in this case the user must first create a set of nodes on
// the shell elements to show which ones should go in the set. The bottom of the
// model may also have a node set and SPC.
//
// At the top of the file, variables can be defined starting with a '$' 
//    e.g.:
//          $freq_max,30
//          $nel_per_wave,5
//          $spc_bottom,7
//          $spc_sides,3
//          $nrb,yes
//          $column,yes
//
// These can be:
//   seg_set_z          Not implemented yet - z-coord at which segment set is created (for LOAD_SEISMIC_SSI)
//   freq_max           max frequency to be considered for wavelength and soil element size (*)
//                             default=30Hz
//   nel_per_wave       number of elements per wavelength (determines soil element size) (*)
//                             default=5;
//      (*) - these two are used only if VS is present in the csv file
//   spc_bottom         SPC code for bottom surface of model (typically 7 (XYZ) if using *LOAD_SEISMIC_SSI)
//   spc_sides          SPC code for sides of the model (requires node set around edge of shells)
//   nrb                Create NRB at every layer (used for 1-element column)
//   bpm_vad    (0,1, or 2) If Boundary_Prescribed_Motion is defined, is it velocity accel or displacement
//   lysmer_ro          Create Lysmer forces and dampers on bottom of model - this parameter gives density
//   lysmer_vs          Shear wave velocity for ditto
//   lysmer_vp          P-wave velocity for ditto
//   lys_lcvx           X-velocity motion vs time curve for Lysmer forces
//   lys_lcvy           Y-velocity motion vs time curve for Lysmer forces
//   lys_lcvz           Z-velocity motion vs time curve for Lysmer forces
//
//
// The first row that doesn't start with a '$' should be the column headers. Each row below this contains
// information to set up one model.  The column headers are:
//
// Header                                  Description
// ------                                  -----------
// thickness                               Thickness of layer                                     
// zcoord                                  Z-coord of top of layer (EITHER depth OR zcoord should be given).
//                                              For first layer, this should be = z of shell elements
// zbot                                    Z-coord of bottom of layer (can be left blank for all except last layer)
// elsize                                  Vertical height of elements (default=layer thickness)
// mtype                                   Material type (supported: 1, 79, 230)
// ro                                      Density
// e                                       Youngs Mod (Mat type 1, 230)            
// pr                                      Poissons Ratio (Mat type 1, 230) 
// k                                       Mat79 bulk modulus
// g                                       Mat79 shear modulus
// lcur                                    Mat79 loadcurve ID
// sclf                                    Mat79 scalefactor
// vs                                      shear wave speed for elsize calculation
// lcbpmx                                  loadcurve for BOUNDARY_PREC_MOTION in x-direction at top of layer
// lcbpmy                                  loadcurve for BOUNDARY_PREC_MOTION in y-direction
// lcbpmz                                  loadcurve for BOUNDARY_PREC_MOTION in z-direction
//                         
//

Message("This is layer_create.js Version 2");

var colname = new Array();
colname[0] = "THICKNESS";
colname[1] = "ZCOORD";
colname[2] = "ELSIZE";
colname[3] = "MTYPE";
colname[4] = "RO";
colname[5] = "E";
colname[6] = "PR";
colname[7] = "K";
colname[8] = "G";
colname[9] = "LCUR";
colname[10] = "SCLF";
colname[11] = "VS";
colname[12] = "LCBPMX";
colname[13] = "LCBPMY";
colname[14] = "LCBPMZ";
colname[15] = "ZBOT";


///////////////////////////////////////////////////////////////////////////////////////
//
//  Initialize

///////////////////////////////////////////////////////////////////////////////////
//  Open the csv file, returning from the function if unable to open it
  var strInputFile = Window.GetFilename("Select file","Select CSV file","csv");
  if (!strInputFile) Exit();
  if ( File.IsFile(strInputFile) )  var fInputFile = new File(strInputFile, File.READ);
  else ErrorTermination("Cannot find or open csv file" + strInputFile);                           

/////////////////////////////////////////////////////////////////////////////////////
//
//  Read data from file

  var strWords = new Array();
  var subWords = new Array();
  var strLine;

  var first = 1;
  var max_param = -1;
  var num_fields = 0;         // Number of fields that a valid line must contain
  var num_layers = 0;         // Number of soil layers
  var column_depth = -1;      // Which column the model name is in
  var column_dir = -1;        // Which column the directory name is in
  var num_param = 0;          // Number of fields in csv file
  var column_param = new Array(); // Which column for each parameter 
                                  // column_param[i] = k means store the i'th column as parameter k
  var param_column = new Array(); // Which parameter is in each column
                                  // param_column[k] = i means the k'th parameter is in the i'th column
  for (var k=0; k<1000; k++) column_param[k] = -1;
  for (var k=0; k<1000; k++) param_column[k] = -1;

  var layer_data = new Array();
  var i_bpm_sides = 0;        // =1 to create BPM on sides
  var seg_set_z = 0.0;        // z at which to create seg set (0 = don't create)
  var nel_per_wave=5;   // number of soil elements per wavelength at max freq
  var freq_max = 25.0;        // max frequency for soil element and wavelength calc
  var spc_bottom = 0;
  var spc_sides = 0;
  var create_nrb = 0;
  var create_column = 0;
  var bpm_vad  = 0;
  var lysmer_ro = 0.0;
  var lysmer_vs = 0.0;
  var lysmer_vp = 0.0;
  var lys_lcvx = 0;
  var lys_lcvy = 0;
  var lys_lcvz = 0;

// Loop over each line
  while ( (strLine = fInputFile.ReadLongLine()) != undefined )
  {
//     Message("Line from csv file: "+strLine);
       strWords = strLine.split(",");     // strWords is an array of strings cut from the line using comma
                                          // as a delimiter

////////////////////////////////////////////////////////////////////////////
// If first character is a '$', read comments
//
      if (strWords[0] != undefined && strWords[0].search(/^\$/) == 0)
      { 

          if (strWords[1] != undefined)
          {

// Check which variable is being set, converting it to lower case
              strVariable = strWords[0].substring(1, strWords[0].length);
              var lowercase = strVariable.toLowerCase();
//              Message("Word0 = " + strWords[0]);
//              Message("Word1 = " + strWords[1]);
	      
              switch(lowercase)
              {
                  case "bpm_layer_sides":      
		  { 
	              	var setting = strWords[1].toLowerCase();
			if (setting=="yes") 
			{
				i_bpm_sides = 1;
				Message("Setting BPM on sides to Yes ");
			}
			else Message("Setting BPM on sides to No  "+" setting="+setting);
			break; 
		  }
                  case "seg_set_z":           
		  { 
			  var seg_set_z = parseFloat(strWords[1]);
			  Message("Found seg_set_z = " + seg_set_z);                 
			  break; 
		  }
                  case "freq_max":           
		  { 
			  var freq_max = parseFloat(strWords[1]);
			  Message("Found freq_max  = " + freq_max );                 
			  break; 
		  }
                  case "nel_per_wave":        
		  { 
			  var nel_per_wave = parseFloat(strWords[1]);
			  Message("Found nel_per_wave = " + nel_per_wave);                 
			  break; 
		  }
                  case "spc_bottom":           
		  { 
			  var spc_bottom = parseFloat(strWords[1]);
			  if (isNaN(spc_bottom)) ErrorTermination("Cannot read line:\n"+strLine);
			  ispc_bottom = Math.round(spc_bottom)
			  Message("Found spc_bottom  = " + ispc_bottom);
			  break; 
		  }
                  case "spc_sides":           
		  { 
			  var spc_sides = parseFloat(strWords[1]);
			  if (isNaN(spc_sides)) ErrorTermination("Cannot read line:\n"+strLine);
			  ispc_sides = Math.round(spc_sides)
			  Message("Found spc_sides  = " + ispc_sides);
			  break; 
		  }
                  case "nrb":           
		  { 
	              	  var setting = strWords[1].toLowerCase();
			  if (setting=="yes") 
			  {
			        var create_nrb = 1;
				Message("Setting NRB to Yes ");
			  }
			  else Message("Setting NRB to No  "+" setting="+setting);
			  break;
		  }
                  case "column":           
		  { 
	              	  var setting = strWords[1].toLowerCase();
			  if (setting=="yes") 
			  {
			        var create_column = 1;
				Message("Setting COLUMN to Yes ");
			  }
			  else Message("Setting COLUMN to No  "+" setting="+setting);
			  break;
		  }
                  case "bpm_vad":           
		  { 
			  var bpm_vad   = parseInt(strWords[1]);
			  Message("Found bpm_vad   = " + bpm_vad  );                 
			  break; 
		  }
		  case "lysmer_ro":
		  {
			  var lysmer_ro = parseFloat(strWords[1]);
			  Message("Found lysmer_ro  = " + lysmer_ro);
			  break;
		  }
		  case "lysmer_vs":
		  {
			  var lysmer_vs = parseFloat(strWords[1]);
			  Message("Found lysmer_vs  = " + lysmer_vs);
			  break;
		  }
		  case "lysmer_vp":
		  {
			  var lysmer_vp = parseFloat(strWords[1]);
			  Message("Found lysmer_vp  = " + lysmer_vp);
			  break;
		  }
                  case "lys_lcvx":          
		  { 
			  var lys_lcvx   = parseInt(strWords[1]);
			  Message("Found lys_lcvx   = " + lys_lcvx  );                 
			  break; 
		  }
                  case "lys_lcvy":          
		  { 
			  var lys_lcvy   = parseInt(strWords[1]);
			  Message("Found lys_lcvy   = " + lys_lcvy  );                 
			  break; 
		  }
                  case "lys_lcvz":          
		  { 
			  var lys_lcvz   = parseInt(strWords[1]);
			  Message("Found lys_lcvz   = " + lys_lcvz  );                 
			  break; 
		  }
		  default: if (lowercase && lowercase !=" ") ErrorTermination("Cannot recognise variable: " + strLine);
              }
          }
      }

////////////////////////////////////////////////////////////////////////////
//    Read column headers
//
      else if (first)    
      {
	      first = 0;
	      for (var i=0; i<strWords.length; i++)
	      {
                      var uppercase = strWords[i].toUpperCase();

		      var found = 0;
		      for (var k=0; k<colname.length; k++)
		      {
			      if (uppercase == colname[k])
			      {
			      		column_param[i] = k;            
			      		param_column[k] = i;            
			      		found = 1; 
			      		num_param++;
					max_param = Math.max(k,max_param);
			      		var prt_i = i+1;
			      		Message("Found header "+strWords[i]+" in column " + prt_i);
			      }
		      }

		      if (!found) 
		      {
       				var answer = Window.Question("Question","Cannot recognise column header: "+strWords[i]
						                       +"\nIgnore this column?",Window.YES|Window.NO);
			      if (answer==Window.NO) ErrorTermination("Unrecognised column");
		      }
	      }

	      Message("Found "+num_param+" columns");
      }

////////////////////////////////////////////////////////////////////////////
//    Read and process a line of data
//
      else if (num_param>1 && strWords.length>=num_param)      //  Normal line of data
      {
	      if (param_column[0]==-1 && param_column[1]==-1) ErrorTermination("Must define either THICKNESS or ZCOORD");             
	      if (param_column[0]!=-1 && param_column[1]!=-1) ErrorTermination("Cannot define both THICKNESS and ZCOORD");             
	      if (param_column[1]!=-1 && param_column[15]==-1) ErrorTermination("If using ZCOORD, must also have ZBOT for bottom layer");

	      var param_list = new Array();            //  Store parameter values
	      for (i=0; i<num_param; i++) param_list[column_param[i]] = parseFloat(strWords[i]);
	      layer_data.push(param_list);

	      num_layers++;
      }
      else
      {
	      if (first) ErrorTermination("Error - Header line not found");
	      ErrorTermination("Error - invalid line: " + strLine);
      }
  }

//  Echo the data to dialog box
Message("num_param = "+num_param);
for (k=0; k<max_param+1; k++) Message("Param "+k+" Column="+param_column[k]);
for (k=0; k<num_param; k++) Message("Column "+k+" Param="+column_param[k]);

for (var i=0; i<=max_param; i++)
{
	if (param_column[i]>-1)
	{
		Message("Layer "+colname[i]);
		for (var k=0; k<num_layers; k++) Message(k+1 + "     " + layer_data[k][i]);
	}
}
fInputFile.Close();

OKToContinue();
	
//   Select or create model;           
if (create_column) var m = new Model();
else var m = Model.Select("Select model");

var flag = AllocateFlag();
m.ClearFlag(flag);

//   Select or create shell parts to extrude
if (create_column)
{
	var n = new Node(m,1,0,0,0)
	var h = new History(m, History.NODE, n.nid, "z=0");
	Message("Created first hist node");
	var n = new Node(m,2,2,0,0)
	var n = new Node(m,3,2,2,0)
	var n = new Node(m,4,0,2,0)
	var p = new Part(m,1,1,1,"Dummy shell");
	p.SetFlag(flag);
	var s = new Shell(m,1,1,1,2,3,4);
	var ns1 = new Set(m,1,Set.NODE,"Side nodes z=0");
	for (var i=1; i<5; i++) ns1.Add(i);
}
else
{
	Part.Select(flag,"Select shell parts to extrude");
}

m.PropagateFlag(flag);

//   Find z-coord of soil surface, max node ID, etc
var n = Node.First(m);
var found = 0;
var nid_next = 1;
var nid_min = 999999999;
var nid_max = -1;
var num_flagged = 0;
while (n)
{
	nid_next = Math.max(nid_next, n.nid+1);
	if (n.Flagged(flag)) 
	{
		found = 1;
		z_surface = n.z;
		nid_min = Math.min(n.nid,nid_min);
		nid_max = Math.max(n.nid,nid_max);
		num_flagged++;
	}
	n = n.Next();
}

//   Decide z-coord of each layer of new nodes
//   Create sub-layers where >1 element per height of layer
var num_sublayers = 0;
var layer_sublayer = new Array();   // gives main layer for each sublayer
var zsublayer = new Array();              // z-coord of layer of nodes
z_prev = z_surface;

for (var i=0; i<num_layers; i++)
{
	if (param_column[0] > -1) zlay = z_prev - layer_data[i][0];
	else if (i<num_layers-1) zlay = layer_data[i+1][1];
	else if (i==num_layers-1) zlay = layer_data[i][15];

	if (param_column[2]>-1 || param_column[11]>-1)   // may be subdivided with user-defined element size
		                 // or by shear wave speed VS (e.g. 5 elements per wavelength at 30Hz)
	{
		dz = Math.abs(zlay - z_prev);
		elsize = dz;
		if (param_column[2]>-1)  
		{
			elsize2 = layer_data[i][2];  // user-defined element size
			if (elsize2>0) elsize = elsize2;
		}
		else 
		{
			Vs = layer_data[i][11];                      // user-input shear wave speed
			if (Math.min(Vs,freq_max,nel_per_wave) > 0)
				elsize = Vs/(freq_max*nel_per_wave);         // nel_per_wave elements per wavelength
		}
		
		Message("elsize="+elsize+"; dz="+dz);
		nsub = Math.round(dz/elsize);                
		nsub = Math.max(1,nsub);
		dzsub = dz/nsub;
		znext = z_prev - dzsub;
		Message("Layer "+i+" has "+nsub + " sublayers");
		for (var k=0; k<nsub; k++)
		{
			layer_sublayer[num_sublayers] = i;
			zsublayer[num_sublayers] = znext;
			Message("...sublayer "+num_sublayers+", z="+znext);
			num_sublayers++;
			znext = znext - dzsub;
		}
	}
	else
	{
		Message("No sublayers - layer "+i+", z="+zlay);
		layer_sublayer[num_sublayers] = i;
		zsublayer[num_sublayers] = zlay;
		num_sublayers++;
	}
	z_prev = zlay;
}


if (!num_flagged) ErrorTerminate("Nothing selected");
if (!found) ErrorTerminate("Cannot find z-coord of soil surface");

var nid_inc = new Array();             // node label inc per layer
ninc_prev = 0;

for (i=0; i<num_sublayers; i++) 
{
	nid_inc[i] = ninc_prev + nid_max + 1 - nid_min;
	ninc_prev = nid_inc[i];
}

//  add another layer of nodes if Lysmer
num_sublayers_n = num_sublayers;        
if (lysmer_ro > 0.0) 
{
	num_sublayers_n = num_sublayers + 1;
	nid_inc[num_sublayers] = nid_inc[num_sublayers-1] + nid_max + 1 - nid_min;
	zsublayer[num_sublayers] = zsublayer[num_sublayers-1] - 1.0;
}

//  create new nodes
for (i=0; i<num_sublayers_n; i++)
{
//	Message("Creating nodes for sublayer "+i+": node ID incr = "+nid_inc[i]);
//	Message("   z = "+zsublayer[i]);
	var n = Node.First(m);
	var first = 1;
	while (n)
	{
		if (n.Flagged(flag))
		{
			var n2 = new Node(m,n.nid+nid_inc[i],n.x,n.y,zsublayer[i]);
			if (first && create_column) var h = new History(m, History.NODE, n.nid+nid_inc[i], "z="+zsublayer[i]);
			first = 0;
		}
		n = n.Next();
	}
}


// check start ID for solids
var solid_next = 1;
var s = Solid.First(m);
while (s)
{
	solid_next = Math.max(solid_next, s.eid);
	s = s.Next();
}
solid_next++;

//  Check start ID for Mat and Part;
var pid_next = 1;
var p = Part.First(m);
var np_flagged = 0;
var pid_min = 999999999;
var pid_max = -1;
while (p)
{
	pid_next = Math.max(pid_next, p.pid);
	if (p.Flagged(flag)) 
	{
		np_flagged++;
		pid_min = Math.min(pid_min,p.pid);
		pid_max = Math.max(pid_max,p.pid);
	}
	p = p.Next();
}

if (!create_column)
{
	var p = Material.First(m);
	while (p)
	{
		pid_next = Math.max(pid_next, p.mid);
		p = p.Next();
	}
	var p =	 Section.First(m);
	while (p)
	{
		pid_next = Math.max(pid_next, p.secid);
		p = p.Next();
	}
	var p =	 Hourglass.First(m);
	while (p)
	{
		pid_next = Math.max(pid_next, p.hgid);
		p = p.Next();
	}
}

pid_next++;

var pid_inc = new Array();            
pid_inc[0]= pid_next - pid_min;

for (i=1; i<num_layers; i++) 
{
	pid_inc[i] = pid_inc[i-1] + pid_max + 1 - pid_min;
}

//  create new part, section, mat, hourglass
//    assume same Material within a layer (not sublayer)
//    assume different parts for each extruded shell part within a layer (not sublayer)
//
var sec = new Section(m, pid_next, Section.SOLID, "Soil")
sec.elform = 1;
var sid = pid_next;

var hg = new Hourglass(m,pid_next);
hg.ihq = 4;
hg.qm = 0.03;
var hgid = pid_next;

for (i=0; i<num_layers; i++)
{
	var mt = Math.round(layer_data[i][3]);                     // create material for layer
	if (mt == 79) var mattype = "*MAT_HYSTERETIC_SOIL";
	else if (mt == 1) var mattype = "*MAT_ELASTIC";
	else if (mt == 230) var mattype = "*MAT_PML_ELASTIC";
	else ErrorTerminate ("Sorry, cannot handle material type "+mt);

	var mat = new Material(m, pid_min+pid_inc[i], mattype);
	mat.SetPropertyByName("RO", layer_data[i][4]);
        var hgidm = 0;

	if (mt==79)
	{
		mat.SetPropertyByName("K0", layer_data[i][7]);
		mat.SetPropertyByName("P0", -1.e10          );
		mat.SetPropertyByName("A0",  1.0            );
		mat.SetPropertyByName("RP",  1.0            );
		mat.SetPropertyByName("LCID", Math.round(layer_data[i][9]));
		mat.SetPropertyByName("SFLC", layer_data[i][10]);
		hgidm = hgid;
	}
	else if (mt==1 || mt==230)
	{
		mat.SetPropertyByName("E", layer_data[i][5]);
		mat.SetPropertyByName("PR", layer_data[i][6]);
		hgidm = 0;
		if (mt==1) hgidm = hgid;
	}


	var psh = Part.First(m);                  // create parts for layer
	while (psh)
	{
		if (psh.Flagged(flag))
		{
			pid_new = psh.pid + pid_inc[i];
			i1 = i+1;
			var partname = "Soil layer "+i1;
			if (np_flagged > 1) partname = partname + " " + psh.heading;
			var p = new Part(m, pid_new, sid, mat.mid, partname);
			p.hgid = hgidm;
		}
		psh = psh.Next();
	}
}

pid_next = pid_max + pid_inc[num_layers-1] + 1;

if (lysmer_ro > 0.0)
{
	var mattype = "*MAT_DAMPER_VISCOUS"
	var mat = new Material(m, pid_next, mattype);
	mat.SetPropertyByName("DC",  1.0            );
	
	var sec = new Section(m, pid_next, Section.DISCRETE, "Lysmer dampers")

	var p = new Part(m, pid_next, pid_next, pid_next, "Lysmer dampers");
	pid_lys = pid_next;
	pid_next++;
}

// Extrude shells to create solids

for (i=0; i<num_sublayers; i++)
{
	inc1 = 0;
	if (i>0) inc1 = nid_inc[i-1];
	inc2 = nid_inc[i] - inc1;
	ilay = layer_sublayer[i];

	var s = Shell.First(m);
	while (s)
	{
		if (s.Flagged(flag))
		{
			var pid = s.pid + pid_inc[ilay];
			n1 = s.n1 + inc1;
			n2 = s.n2 + inc1;
			n3 = s.n3 + inc1;
			n4 = s.n4 + inc1;
			n5 = n1 + inc2;
			n6 = n2 + inc2;
			n7 = n3 + inc2;
			n8 = n4 + inc2;
			var so = new Solid(m,solid_next,pid,n5,n6,n7,n8,n1,n2,n3,n4);
			solid_next++;
		}
		s = s.Next();
	}
}


//   Find node set ID
var set_next = 1;      
var set = Set.First(m, Set.NODE);
while (set)
{
	set_next = Math.max(set_next,set.sid+1);
	set = set.Next();
}

//  Create segment set for *LOAD_SEISMIC_SSI
if (seg_set_z != 0.0)
{
	i_best = -1;
	dz_best = 1e10;
	for (i=0; i<num_sublayers; i++)   // find closest sublayer
	{
		dz = Math.abs(seg_set_z - zsublayer[i]);
		if (dz < dz_best)
		{
			dz_best = dz;
			i_best = i;
		}
	}

	if (i_best==-1) ErrorTerminate("Could not find layer for segment set at z="+seg_set_z);

	var set_s_next = 1;                       // Find ID for segment set
	var set = Set.First(m, Set.SEGMENT);
	while (set)
	{
		set_s_next = Math.max(set_s_next,set.sid+1);
		set = set.Next();
	}

	f = new File("./prtmp_setseg.key", File.WRITE);     // Open temporary file
	f.Writeln("*KEYWORD");                       // Write set keyword data to file
	f.Writeln("*SET_SEGMENT_TITLE");
	f.Writeln("For *LOAD_SEISMIC_SSI");
	f.Writeln(NumberToString(set_s_next,10));      //  Set ID

	var inc1 = nid_inc[i_best];
	var s = Shell.First(m);
	while (s)
	{
		if (s.Flagged(flag))
		{
			n1 = s.n1 + inc1;
			n2 = s.n2 + inc1;
			n3 = s.n3 + inc1;
			n4 = s.n4 + inc1;
			f.Writeln(NumberToString(n1,10) + NumberToString(n2,10) 
				+ NumberToString(n3,10) + NumberToString(n4,10));      //  Segment
		}
		s = s.Next();
	}

	f.Writeln("*END");
	f.Close();                                  	    // Close temporary file
	m.Import("prtmp_setseg.key");                       // Import data from file into model
	File.Delete("./prtmp_setseg.key");                  // Delete temporary file

//   Set of nodes for *LOAD_SEISMIC_SSI
	var ns = new Set(m, set_next, Set.NODE,"Nodes for *LOAD_SEISMIC_SSI");
	var n = Node.First(m);
	while (n)
	{
		if (n.Flagged(flag)) ns.Add(n.nid+inc1);
		n = n.Next();
	}
	set_next++;
}

//   Sets of nodes per layer
sublayer_nset = new Array();
if (spc_sides != 0 || param_column[12]>-1 || param_column[13]>-1 || param_column[14]>-1 || create_nrb)
{
	if (spc_sides != 0) var codes = BreakdownSPC(spc_sides);

//   Select the "master" set on the shells
        if (!create_column)
	{
		Set.Select(Set.NODE,flag,"Select node set on edges of shells for SPC or BPM",m,false);
		var ns0 = Set.First(m, Set.NODE);
		var ns1;
		var count=0;
		while (ns0)
		{
			if (ns0.Flagged(flag)) 
			{
				count++;
				ns1 = ns0;
			}
			ns0 = ns0.Next();
		}
		if (!ns1) ErrorTermination("Did not select a node set?");
		if (count>1) ErrorTermination("Must select only one node set");
		ns1.title = "Side nodes at z=0";
	}

//   Create SPC, BPM for the master node set
	CreateSPC_BPM(ns1.sid, 0);

//   Create new node set at each sublayer of nodes
	for (i=0; i<num_sublayers; i++)
	{
		inc = nid_inc[i];
		var ns = new Set(m, set_next, Set.NODE, "Side nodes, z="+zsublayer[i]);
		ns1.StartSpool();
		while (id = ns1.Spool() ) ns.Add(id+inc);  // copy master set using known node ID increment
	
		ilay = layer_sublayer[i]+1;     // node inc for sublayer i are at bottom of element sublayer
		Message("i="+i+", ilay="+ilay+" num_layers="+num_layers);
		if (layer_data[ilay]) CreateSPC_BPM(set_next,ilay);          //  create SPC and/or BPM
		else if (ilay==num_layers) CreateSPC_BPM(set_next,-1);     //  special treatment for sublayers within
		                                                             //      bottom layer

		set_next++;
	}
	OKToContinue();

	ns1.ClearFlag(flag);     //  So we don't delete the set from the column model

}	

//   Set of nodes at bottom and/or Lysmer treatment
if (spc_bottom != 0 || lysmer_ro > 0.0)
{
	var ibot = num_sublayers_n-1;
	var ns = new Set(m, set_next, Set.NODE,"Nodes on bottom plane");
	var n = Node.First(m);
	while (n)
	{
		if (n.Flagged(flag)) 
		{
			ns.Add(n.nid+nid_inc[ibot]);
			if (lysmer_ro > 0.0) 
			{
				var area = CalcArea(n);
				var area_vs = area*lysmer_ro*lysmer_vs;
				var area_vp = area*lysmer_ro*lysmer_vp;
				var nsp1 = n.nid+nid_inc[ibot-1];
				var nsp2 = n.nid+nid_inc[ibot];
				if (area_vs > 0.0)
				{
					var sp = new Discrete(m, solid_next, pid_lys, nsp1, nsp2, 1, area_vs);
					solid_next++;
					var sp = new Discrete(m, solid_next, pid_lys, nsp1, nsp2, 2, area_vs);
					solid_next++;
					if (lys_lcvx > 0) var l = new LoadNode(m, LoadNode.POINT, nsp1, 1, lys_lcvx, area_vs);
					if (lys_lcvy > 0) var l = new LoadNode(m, LoadNode.POINT, nsp1, 2, lys_lcvy, area_vs);
				}
				if (area_vp > 0.0)
				{
					var sp = new Discrete(m, solid_next, pid_lys, nsp1, nsp2, 3, area_vp);
					solid_next++;
					if (lys_lcvz > 0) var l = new LoadNode(m, LoadNode.POINT, nsp1, 3, lys_lcvz, area_vp);
				}
//				OKToContinue();
			}
		}
		n = n.Next();
	}

	if (lysmer_ro > 0.0) spc_bottom = 7;
	var codes = BreakdownSPC(spc_bottom);
	var spc = new Spc(m, set_next, 0, codes.x,codes.y,codes.z,0,0,0,Spc.SET);
	set_next++;
}	

//  Write the damper orientation vectors for Lysmer
if (lysmer_ro > 0.0)
{
	f = new File("./prtmp_setseg.key", File.WRITE);     // Open temporary file
	f.Writeln("*KEYWORD");                       // Write set keyword data to file
	f.Writeln("*DEFINE_SD_ORIENTATION");
	f.Writeln("         1         0       1.0       0.0       0.0");           
	f.Writeln("         2         0       0.0       1.0       0.0");           
	f.Writeln("         3         0       0.0       0.0       1.0");           
	f.Writeln("*END");
	f.Close();                                  	    // Close temporary file
	m.Import("prtmp_setseg.key");                       // Import data from file into model
	File.Delete("./prtmp_setseg.key");                  // Delete temporary file
}

//  Column model - delete dummy shell
if (create_column) m.DeleteFlagged(flag);

///////////////////////////////////////////////////////////////////////////////////////
//
//  Finish

m.UpdateGraphics();

Termination("Finished!");

/////////////////////////////////////////////////////////////////////////////
//
function Termination(line)
{
	if (fInputFile) fInputFile.Close();
	if (m) m.UpdateGraphics();
	Message(line);
	if (flag) ReturnFlag(flag);
	Exit();
}
/////////////////////////////////////////////////////////////////////////////
//
function ErrorTermination(line)
{
	Window.Error("Error",line,Window.OK);
	Termination(line);
}

///////////////////////////////////////////////////////
//
function OKToContinue()
{
	var answer = Window.Question("Question","OK to continue?", Window.YES | Window.NO | Window.NONMODAL);
	if (answer==Window.NO) Exit();
}

////////////////////////////////////////////////////////
//
function BreakdownSPC(input)
{
	var codes = new Object();
	codes.x = 0;
	codes.y = 0;
	codes.z = 0;
	if (input==1 || input==4 || input==6 || input==7) codes.x = 1;
	if (input==2 || input==4 || input==5 || input==7) codes.y = 1;
	if (input==3 || input==5 || input==6 || input==7) codes.z = 1;
	return codes;
}
////////////////////////////////////////////////////////
//
function CreateSPC_BPM(set_next,ilay)
{
	var bpm_type = PrescribedMotion.SET;
	id = set_next;
	Message("In CreateSPC_BPM, ilay="+ilay+", set_next="+set_next);

	if (create_nrb)                           //   NRB
	{
		var nrb = new Nrb(m, set_next, pid_next);
		nrb.spc = true;
  		nrb.cmo = 1;
  		nrb.con1 = spc_sides;
  		nrb.con2 = 7;             
  		bpm_type = PrescribedMotion.NRBC;
//		bpm_type = PrescribedMotion.RIGID;
		id = pid_next;
		pid_next++;
	}
	if (spc_sides != 0 && !create_nrb)         //  create SPC if required
	{
		var spc = new Spc(m, id, 0, codes.x,codes.y,codes.z,0,0,0,Spc.SET);
	}

	if (ilay<0) return;

	if (param_column[12] > -1)                //   create BPM in X if required - note, same curve for whole layer
	{
		var lcid = layer_data[ilay][12];
		if (lcid>0) var bpm = new PrescribedMotion(m, id, 1, bpm_vad, lcid, bpm_type);              
	}
	if (param_column[13] > -1)                //   BPM in Y
	{
		var lcid = layer_data[ilay][13];
		if (lcid>0) var bpm = new PrescribedMotion(m, id, 2, bpm_vad, lcid, bpm_type);              
	}
	if (param_column[14] > -1)                //   BPM in Z
	{
		var lcid = layer_data[ilay][14];
		if (lcid>0) var bpm = new PrescribedMotion(m, id, 3, bpm_vad, lcid, bpm_type);              
	}
}

//
//
function CalcArea(n)
{
	var sum_a = 0.0;
	var shells = n.GetAttachedShells(false);
	for (var ref1=0; ref1<shells.length; ref1++)
	{
		var shl = shells[ref1];
		if (shl.Flagged(flag))
		{
			var n1sid = shl.n1;
			var n2sid = shl.n2;
			var n3sid = shl.n3;
			if (shl.nodes>3) var n4sid = shl.n4;
			var ns1 = Node.GetFromID(m, n1sid);
			var ns2 = Node.GetFromID(m, n2sid);
			var ns3 = Node.GetFromID(m, n3sid);
			if (shl.nodes>3) var ns4 = Node.GetFromID(m, n4sid);
			dx1 = ns3.x - ns1.x;
			dy1 = ns3.y - ns1.y;
			dz1 = ns3.z - ns1.z;
			ns5 = ns3;
			if (shl.nodes > 3) ns5 = ns4;
			dx2 = ns5.x - ns2.x;
			dy2 = ns5.y - ns2.y;
			dz2 = ns5.z - ns2.z;
			var crprd = dx1*dy2 - dx2*dy1;
			d_area = crprd*0.5/shl.nodes;
//			Message("For node "+n.nid+": Shell " + shl.eid + " has area per node = " + d_area);
			sum_a = sum_a + d_area;
		}
	}
//	Message("Total area for node "+n.nid+" = "+sum_a);
	return sum_a;
}



