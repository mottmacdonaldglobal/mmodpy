// Script to read a csv file to create soil layers     
//
// name:Primer Automesh v9   
// description: create 3D mesh from shells defining foundation.
//
//   

// DESCRIPTION
// ===========
//
// Version 3 of this script has been updated so that triangular elements in the shell mesh 
// do not cause a problem as solids. Previously the nodes were being ordered incorrectly
// when they were extruded.
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
//   Variables for creating strip model
//    column                  Yes, No, Strip, plan. 
//                            strip option cerates a strip model for foundation stiffness assessment.
//                            using strip also overwrites some of the other parameters, e.g. nrb, spc_bottom, etc  
//                            see comments when running script.
//                            foundation is modelled as rigid by applying prescribed displacements.
//    strip_full              1= full strip model. 0 = half model (for vertical loading only). default = 1.
//    found_width             *width of foundation. required variable.
//    strip_width             overall width of strip model. default = 10 x found_width.
//    strip_depth             depth (thickness) of strip model in out of plane direction. default = 1.0m.
//    found_max_ele_width     maximum width of elements in foundation zone and refined zone around foundation. default = found_width/10.
//    mesh_max_ele_width      maximum width of elements outside foundation / refined zone. default = 2*found_max_ele_width.
//    presc_disp_z            value of vertical displacement to apply. if omitted, no presc disp created.
//    presc_disp_x            value of horizontal displacement to apply. if omitted, no presc disp created.
//    disp_time_0             time at which prescribed displacement will begin. default = 0.
//    disp_time_ramp          time over which to ramp the prescribed displacement. default = 2s.
//    term_time               analysis termination time. default = disp_time_0 + disp_time_ramp.
//    curve                   denotes stress strain curve. example shown below.
//    endcurves               required keyword to define end of curves.
//    tstep                   *CONTROL_TIMESTEP values DTINIT,TSSFAC,ISDO,TSLIMT,DT2MS,LCTM,ERODE,MS1ST
//                            note the LCTM should be set to the value of the load curve rather than the LCID
//                            if it is not set to zero, a load curve is created with Y values equal to the value specified here
//    
//   
//        * denotes required variables. if other variables are missing, some assumptions are made. see messages 
//          when running script.
//    curve example - creates curve 10003 and 10004.
//
//    $curve,,,,,,,
//    10003,title1,,,,,
//    0.000001,9.95E-07,,,,,,
//    0.00001,9.57E-06,,,,,,
//    0.00003,2.67E-05,,,,,,
//    0.0001,7.26E-05,,,,,,
//    0.0003,0.000147514,,,,,,
//    0.001,0.000247969,,,,,,
//    0.003,0.000421245,,,,,,
//    0.01,0.000689314,,,,,,
//    0.03,0.000829604,,,,,,
//    0.05,0.000843372,,,,,,
//    $curve,,,,,,,
//    10004,title2,,,,,
//    0.000001,9.95E-07,,,,,,
//    0.00001,9.62E-06,,,,,,
//    0.00003,2.70E-05,,,,,,
//    0.0001,7.51E-05,,,,,,
//    0.0003,0.000156979,,,,,,
//    0.001,0.000285068,,,,,,
//    0.003,0.000481402,,,,,,
//    0.01,0.000706519,,,,,,
//    0.03,0.000833743,,,,,,
//    0.05,0.000866711,,,,,,
//    $endcurves,,,
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
// pid                                     part id for layer
// stype                                   soil type (CLAY, SILTY CLAY, PEAT, SILT, SAND)

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
colname[16] = "PID";
colname[17] = "STYPE";

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//
//  Initialize
//
/////////////////////////////////////////////////////////////////////////////////////
//  Open the csv file, returning from the function if unable to open it
if (arguments[1])
{
      var strInputFile = arguments[1];
}else{
      //var strInputFile = "C:/OneDrive - Arup/Docs/P500/Automesh/column.csv";  
      var strInputFile = Window.GetFilename("Select file","Select CSV file","csv");
}
  
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
  var nel_per_wave=5;         // number of soil elements per wavelength at max freq
  var freq_max = 30.0;        // max frequency for soil element and wavelength calc
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
  var mat_hys_cur = 0;
  var bdamp = 0;
  var write_mod = 0;
  var foundation_depth = 0;
  
  var strip_full = 1;
  var strip_depth = 0;                  // depth of model - out of plane
  var strip_width = 0;                 // width of strip model - in plane
  var found_width = 0;                  // footing half width
  var found_max_ele_width = 0;
  var mesh_max_ele_width = 0;
  var mesh_max_ele_window = 0;
  var mesh_max_ele_width2 = 0;

  var presc_disp_x_sw = 0;
  var presc_disp_z_sw = 0;
  var presc_disp_z = 0;
  var presc_disp_x = 0;
  var disp_time_0 = -1;
  var disp_time_ramp = 0;

  var term_time = 0;  
  
  var tstep = 0;
  var tstepa = new Array(0,0,0,0,0,0,0,0);
  
  var con_sol = 0;
  var con_sol1 = new Array(0,0,0,0,0,0,0,0);
  
  var con_energy = 0;
  var con_energy1 = new Array(0,0,0,0,0,0,0,0); 
//-
  var con_mpp_io_nodump = 0; 
  
  var con_solid = 0;
  var con_solid1 = new Array(0,0,0,0,0,0,0,0); 
  
  var con_shell = 0;
  var con_shell1 = new Array(0,0,0,0,0,0,0,0);   
  
  var con_bulk = 0;
  var con_bulk1 = new Array(0,0,0,0,0,0,0,0); 

  var con_cont = 0;
  var con_cont1 = new Array(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
  
  var data_bin_d3pl = 0;
  var data_bin_d3pl1 = new Array(0,0,0,0,0,0,0,0);
  
  var data_bin_d3th = 0;
  var data_bin_d3th1 = new Array(0,0,0,0,0,0,0,0);  
  
  var data_extent_bin = 0;
  var data_extent_bin1 = Array(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0); 
  
  var data_elout = 0;
  var data_elout1 = new Array(0,0,0,0,0,0,0,0); 
  
  var data_glstat = 0;
  var data_glstat1 = new Array(0,0,0,0,0,0,0,0);   
  
  var data_matsum = 0;
  var data_matsum1 = new Array(0,0,0,0,0,0,0,0);  

  var data_nodfor = 0;
  var data_nodfor1 = new Array(0,0,0,0,0,0,0,0);  

  var data_nodout = 0;
  var data_nodout1 = new Array(0,0,0,0,0,0,0,0);  

  var data_sleout = 0;
  var data_sleout1 = new Array(0,0,0,0,0,0,0,0); 

  var data_curout = 0;
  var data_curout1 = new Array(0,0,0,0,0,0,0,0);  

  
// Create or select model
if (Model.GetAll().length==1) var model = Model.First();
else if (!arguments[1] && Model.GetAll().length > 1) var model = Model.Select("Select model");
else if (arguments[1] && Model.GetAll().length > 1) var model = Model.First();
else var model = new Model();
var curLayer = Include.GetFromID(model, model.layer);

///////////////////////////////////////////////////////////////////////////
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
			  }else if (setting=="strip") 
			  {
			        var create_column = 2;
				Message("Setting COLUMN to strip ");
			  }else if (setting=="plan") 
			  {
			        var create_column = 3;
				Message("Setting COLUMN to plan ");
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
                  case "bdamp":          
		  { 
			  var bdamp   = parseFloat(strWords[1]);
			  Message("Found bdamp   = " + bdamp  );                 
			  break; 
		  }              
                  case "curve":          
		  { 
			  var mat_hys_cur   = 1;
                    var Cur_Array = new Array();
                    var Cur_Cnt = -1;
			  Message("Found curves for mat_hys_soil");                 
			  var Lcount = 0; 
                    var strLine = fInputFile.ReadLongLine();
                    strWords = strLine.split(",");                     
                    var lcid =  0;
                    
                    while (strWords[0].search("endcurves") < 0)
                    {
                          //Message("Line " + strLine);
                          
                          if (strWords[0] != undefined)
                          {
                                if (strWords[0].search(/^\$/) != 0)
                                {
                                      if(Lcount==0)
                                      {
                                            var lcid =  parseInt(strWords[0]);
                                            var cur_title = strWords[1];
                                            var x_fac = 1;
                                            var y_fac = 1;
                                            if(parseFloat(strWords[2])>0) x_fac = parseFloat(strWords[2]);
                                            if(parseFloat(strWords[3])>0) y_fac = parseFloat(strWords[3]);

                                            Cur_Cnt++;
                                            Cur_Array[Cur_Cnt] = new Array();
                                            Cur_Array[Cur_Cnt][0] = lcid;
                                            Cur_Array[Cur_Cnt][1] = cur_title;
                                            Cur_Array[Cur_Cnt][2] = x_fac;
                                            Cur_Array[Cur_Cnt][3] = y_fac;
                                            
                                      }else if(Lcount>0 && Lcount<11){
                                            //l.AddPoint(parseFloat(strWords[0]), parseFloat(strWords[1]));
                                            Cur_Array[Cur_Cnt].push([parseFloat(strWords[0]), parseFloat(strWords[1])]);
                                            //Message(parseFloat(strWords[0]) + ", " + parseFloat(strWords[1]));
                                      }else{
                                            ErrorTermination("More than 10 points found for mat hys curve id " + lcid);
                                      }
                                      Lcount++;
                                      
                                }else{
                                      
                                      strVariable = strWords[0].substring(1, strWords[0].length);
                                      var lowercase = strVariable.toLowerCase();
                                      if(lowercase.search("curve")>=0) Lcount = 0;
                                      //Message("New Curve. Lcount " + Lcount);
                                      //Message(lowercase.search("curves"));
                                }
                                
                          }
                          var strLine = fInputFile.ReadLongLine();
                          strWords = strLine.split(",");                    
                    }
                    
                    
                    break; 
		  }
                  case "strip_depth":          
		  { 
			  var strip_depth   = parseFloat(strWords[1]);
			  Message("Found strip_depth   = " + strip_depth  );             
			  break; 
		  } 
                  case "strip_width":          
		  { 
			  var strip_width   = parseFloat(strWords[1]);
			  Message("Found strip_width   = " + strip_width  );                 
			  break; 
		  } 
                  case "found_width":          
		  { 
			  var found_width   = parseFloat(strWords[1]);
			  Message("Found found_width   = " + found_width  );                 
			  break; 
		  }  
                  case "found_max_ele_width":          
		  { 
			  var found_max_ele_width   = parseFloat(strWords[1]);
			  Message("Found found_max_ele_width   = " + found_max_ele_width  );                 
			  break; 
		  } 
                  case "mesh_max_ele_width":          
		  { 
			  var mesh_max_ele_width   = parseFloat(strWords[1]);
                    if(parseFloat(strWords[2]) > 0) mesh_max_ele_window = parseFloat(strWords[2]);
                    if(parseFloat(strWords[3]) > 0) mesh_max_ele_width2 = parseFloat(strWords[3]);                    
			  Message("Found mesh_max_ele_width   = " + mesh_max_ele_width  );                 
			  break; 
		  }
                  case "presc_disp_x":          
		  { 
			  var presc_disp_x   = parseFloat(strWords[1]);

                    if(presc_disp_x == undefined)
                    {
                          presc_disp_x_sw = 0;
                    }else{
                          presc_disp_x_sw = 1;
                    }
                    
			  Message("Found presc_disp_x   = " + presc_disp_x  );                 
			  break; 
		  }
                  case "presc_disp_z":          
		  { 
			  var presc_disp_z   = parseFloat(strWords[1]);
                    
                    if(presc_disp_z == undefined)
                    {
                          presc_disp_z_sw = 0;
                    }else{
                          presc_disp_z_sw = 1;
                    }                    
			  
                    Message("Found presc_disp_z   = " + presc_disp_z  );                 
			  break; 
		  }              
                  case "disp_time_0":          
		  { 
			  var disp_time_0   = parseFloat(strWords[1]);
			  Message("Found disp_time_0   = " + disp_time_0  );                 
			  break; 
		  } 
                  case "disp_time_ramp":          
		  { 
			  var disp_time_ramp   = parseFloat(strWords[1]);
			  Message("Found disp_time_ramp   = " + disp_time_ramp  );                 
			  break; 
		  } 
                  case "term_time":          
		  { 
			  var term_time   = parseFloat(strWords[1]);
			  Message("Found term_time   = " + term_time  );                 
			  break; 
		  } 
                  case "strip_full":          
		  { 
			  var strip_full   = parseInt(strWords[1]);
			  Message("Found strip_full   = " + strip_full  );                 
			  break; 
		  } 
                  case "write_mod":          
		  { 
			  var write_mod   = 1;
                    var model_fname = model.path + "New_Model.k";
                    if (strWords[1]) model_fname = model.path + strWords[1];
			  Message("Found write_mod   = " + write_mod  );                 
			  break; 
		  } 
                  case "foundation_depth":          
		  { 
			  var foundation_depth   = parseFloat(strWords[1]);
			  Message("Found foundation_depth   = " + foundation_depth);                 
			  break; 
		  }  
                  case "control_timestep":          
		  { 
			  tstep = 1;
                    
                    for (ii=1; ii<strWords.length; ii++)
                    {
                         if(parseFloat(strWords[ii]) != 0  && strWords[ii] != "")  tstepa[ii-1] = parseFloat(strWords[ii]);
                    }
			  Message("Found tstep   = " + tstepa);                 
			  break; 
		  }
                  case "control_solution":          
		  { 
			  con_sol = 1;
                    
                    for (ii=1; ii<strWords.length; ii++)
                    {
                         if(parseFloat(strWords[ii]) != 0  && strWords[ii] != "")  con_sol1[ii-1] = parseFloat(strWords[ii]);
                    }
			  Message("Found con_sol   = " + con_sol);                 
			  break; 
		  }
                  case "control_energy":          
		  { 
			  con_energy = 1;
                    
                    for (ii=1; ii<strWords.length; ii++)
                    {
                         if(parseInt(strWords[ii]) != 0 && strWords[ii] != "")  con_energy1[ii-1] = parseInt(strWords[ii]);
                    }
			  Message("Found con_energy   = " + con_energy);                 
			  break; 
		  }
                  case "control_mpp_io_nodump":          
		  { 
			  con_mpp_io_nodump = 1;
                    
			  Message("Found con_mpp_io_nodump   = " + con_mpp_io_nodump);                 
			  break; 
		  }
                  case "control_solid":          
		  { 
			  con_solid = 1;
                    
                    for (ii=1; ii<strWords.length; ii++)
                    {
                         if(parseInt(strWords[ii]) != 0  && strWords[ii] != "")  con_solid1[ii-1] = parseInt(strWords[ii]);
                    }
                    
                    if(parseFloat(strWords[6]) != 0  && strWords[6] != "")  con_solid1[5] = parseFloat(strWords[6]);
                    
			  Message("Found con_solid   = " + con_solid);                 
			  break; 
		  }
                  case "control_shell":          
		  { 
			  con_shell = 1;
                    
                    for (ii=1; ii<strWords.length; ii++)
                    {
                         if(parseInt(strWords[ii]) != 0  && strWords[ii] != "")  con_shell1[ii-1] = parseInt(strWords[ii]);
                    }
                    
                    if(parseFloat(strWords[1]) != 0   && strWords[1] != "")  con_shell1[0] = parseFloat(strWords[1]);
                    
			  Message("Found con_shell   = " + con_shell);                 
			  break; 
		  }
                  case "control_bulk_viscosity":          
		  { 
			  con_bulk = 1;
                    
                    if(parseFloat(strWords[1]) != 0 && strWords[1] != "")  con_bulk1[0] = parseFloat(strWords[1]);
                    if(parseFloat(strWords[2]) != 0 && strWords[2] != "")  con_bulk1[1] = parseFloat(strWords[2]);
                    if(parseInt(strWords[3]) != 0   && strWords[3] != "")  con_bulk1[2] = parseInt(strWords[3]);
                    if(parseInt(strWords[4]) != 0   && strWords[4] != "")  con_bulk1[3] = parseInt(strWords[4]);
                    
                    
			  Message("Found con_bulk   = " + con_bulk);                 
			  break; 
		  }
                  case "control_contact":          
		  { 
			  con_cont = 1;
                    
                    for (ii=1; ii<strWords.length; ii++)
                    {
                         if(parseInt(strWords[ii]) != 0  && strWords[ii] != "")  con_cont1[ii-1] = parseInt(strWords[ii]);
                    }
                    
                    if(parseFloat(strWords[1]) != 0   && strWords[1] != "")  con_bulk1[0] = parseFloat(strWords[1]);
                    if(parseFloat(strWords[2]) != 0   && strWords[2] != "")  con_bulk1[1] = parseFloat(strWords[2]);
                    if(parseFloat(strWords[13]) != 0   && strWords[13] != "")  con_bulk1[12] = parseFloat(strWords[13]);
                    
			  Message("Found con_cont   = " + con_cont);                 
			  break; 
		  } 
                  case "database_binary_d3plot":          
		  { 
			  data_bin_d3pl = 1;
                    
                    for (ii=1; ii<strWords.length; ii++)
                    {
                         if(parseFloat(strWords[ii]) != 0  && strWords[ii] != "")  data_bin_d3pl1[ii-1] = parseFloat(strWords[ii]);
                    }
                    
			  Message("Found data_bin_d3pl   = " + data_bin_d3pl);                 
			  break; 
		  }
                  case "database_binary_d3thdt":          
		  { 
			  data_bin_d3th = 1;
                    
                    for (ii=1; ii<strWords.length; ii++)
                    {
                         if(parseFloat(strWords[ii]) != 0  && strWords[ii] != "")  data_bin_d3th1[ii-1] = parseFloat(strWords[ii]);
                    }
                    
			  Message("Found data_bin_d3th   = " + data_bin_d3th);                 
			  break; 
		  }
                  case "database_extent_binary":          
		  { 
			  data_extent_bin = 1;
                    
                    for (ii=1; ii<strWords.length; ii++)
                    {
                         if(parseInt(strWords[ii]) != 0  && strWords[ii] != "")  data_extent_bin1[ii-1] = parseInt(strWords[ii]);
                    }
                    
			  Message("Found data_extent_bin   = " + data_extent_bin);                 
			  break; 
		  } 
                  case "database_elout":          
		  { 
			  data_elout = 1;
                    data_elout1[0] = parseFloat(strWords[1]);
                    data_elout1[1] = parseInt(strWords[2]);
                    
			  Message("Found data_elout   = " + data_elout);                 
			  break; 
		  }  
                  case "database_glstat":          
		  { 
			  data_glstat = 1;
                    data_glstat1[0] = parseFloat(strWords[1]);
                    if (strWords[2] != "") data_glstat1[1] = parseInt(strWords[2]);
                    
			  Message("Found data_glstat   = " + data_glstat);                 
			  break; 
		  }  
                  case "database_matsum":          
		  { 
			  data_matsum = 1;
                    data_matsum1[0] = parseFloat(strWords[1]);
                    if (strWords[2] != "") data_matsum1[1] = parseInt(strWords[2]);

                    
			  Message("Found data_matsum   = " + data_matsum);                 
			  break; 
		  }  
                  case "database_nodfor":          
		  { 
			  data_nodfor = 1;
                    data_nodfor1[0] = parseFloat(strWords[1]);
                    if (strWords[2] != "") data_nodfor1[1] = parseInt(strWords[2]);
                    
			  Message("Found data_nodfor   = " + data_nodfor);                 
			  break; 
		  }  
                  case "database_nodout":          
		  { 
			  data_nodout = 1;
                    data_nodout1[0] = parseFloat(strWords[1]);
                    if (strWords[2] != "") data_nodout1[1] = parseInt(strWords[2]);
                    
			  Message("Found data_nodout   = " + data_nodout);                 
			  break; 
		  }  
                  case "database_sleout":          
		  { 
			  data_sleout = 1;
                    data_sleout1[0] = parseFloat(strWords[1]);
                    if (strWords[2] != "") data_sleout1[1] = parseInt(strWords[2]);
                    
			  Message("Found data_sleout   = " + data_sleout);                 
			  break; 
		  }
                  case "database_curvout":          
		  { 
			  data_curout = 1;
                    data_curout1[0] = parseFloat(strWords[1]);
                    if (strWords[2] != "") data_curout1[1] = parseInt(strWords[2]);
                    
			  Message("Found data_sleout   = " + data_sleout);                 
			  break; 
		  }              
		  default: 
              if (lowercase && lowercase !=" " )
              {
                    if(lowercase && lowercase.length >0)
                    {
                    ErrorTermination("Cannot recognise variable: " + strLine);
                    }
              }
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

		      if (!found && strWords[i].length >0) 
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
            if (param_column[17]!=-1) param_list[17] = strWords[param_column[17]];
            
            layer_data.push(param_list);

	      num_layers++;
      }
      else
      {
	      if (first) ErrorTermination("Error - Header line not found");
	      ErrorTermination("Error - invalid line: " + strLine);
      }
  }
  
if (param_column[16]==-1)
{
      var pid_spec = 0;
}else{
      var pid_spec = 1;
}
if (param_column[17]==-1)
{
      var p_desc = 0;
}else{
      var p_desc = 1;
}
////////////////////////////////////////////////////////////////////////////  
// checks on input for strip model
if (create_column ==2)
{
     var all_params_defined = 1;
      
     Message("--------------------------------------------------------------------------------------------------");
     Message("You are creating a strip model. The following parameters have been overwritten...");
     Message("Setting NRB to No");
     create_nrb = 0; 
     Message("Setting spc_bottom to 3");
     spc_bottom = 3;  
     Message("Setting spc_sides to 0");
     spc_sides = 0; 
     //Message("Setting lysmer_ro to 0");
     //lysmer_ro = 0;  
     Message("----------------------------------------");
     
     if(found_width ==0)
     {
           ErrorTermination("Foundation width not defined. Please provide found_width");
     }
     
     
     if(strip_depth ==0)
     {
           Message("strip_depth not defined, setting this value to 1.0m");
           strip_depth = 1.0;
           all_params_defined = 0;
     }

     if(strip_width ==0)
     {
           Message("strip_width not defined, setting this value to 20 x found_width.");
           strip_width = found_width*10.0;
           all_params_defined = 0;
     }
     

     if(found_max_ele_width ==0)
     {
           Message("found_max_ele_width not defined, setting this value to found_width/10");
           found_max_ele_width = found_width/10;
           all_params_defined = 0;
     }

     if(mesh_max_ele_width ==0)
     {
           Message("mesh_max_ele_width not defined, setting this value to 2.0 x found_max_ele_width");
           mesh_max_ele_width = 2*found_max_ele_width;
           all_params_defined = 0;
     }

//     if(presc_disp_z ==0)
//     {
//           Message("presc_disp_z not defined, setting this value to -0.05");
//           presc_disp_z = -0.05;
//           all_params_defined = 0;
//     }  

     if(disp_time_0 < 0)
     {
           Message("disp_time_0 not defined, setting this value to 0.0s");
           disp_time_0 = 0.0;
           all_params_defined = 0;
     }  

     if(disp_time_ramp <= 0)
     {
           Message("disp_time_ramp not defined, setting this value to 2.5s");
           disp_time_ramp = 2.5;
           all_params_defined = 0;
     }       
    
     if(term_time <=0)
     {
           Message("term_time not defined, setting this value to disp_time_0 + disp_time_ramp");
           term_time = disp_time_0 + disp_time_ramp;
           all_params_defined = 0;
     }  

     
     found_width = found_width*0.5;
     strip_width = strip_width*0.5;
     
      if(all_params_defined != 1) Message("Not all parameters have been defined manually. Please review assumptions made above!!!")
   
      //OKToContinue("Please review strip model input messages. OK to continue?");
} else if (create_column ==3)
{
     var all_params_defined = 1;
      
     Message("--------------------------------------------------------------------------------------------------");
     Message("You are creating a 3D model. The following parameters have been overwritten...");
     Message("Setting NRB to No");
     //create_nrb = 0; 
     //Message("Setting spc_bottom to 3");
     spc_bottom = 3;  
     Message("Setting spc_sides to 0");
     spc_sides = 0; 
     //Message("Setting lysmer_ro to 0");
     //lysmer_ro = 0;  
     Message("----------------------------------------");
     
     if(found_width ==0)
     {
           //ErrorTermination("Foundation width not defined. Please provide found_width");
           found_width=1;     //set to 1 to avoid potential problems though this should not be used
     }

     if(strip_width ==0)
     {
           Message("strip_width is the 3D mesh size in X-direction. You will be asked to enter a value.");
           strip_width = Window.GetNumber("Mesh Size - X-Direction?","How big would you like the mesh to be in the X-Direction?");
           all_params_defined = 0;
     }
     
     if(strip_depth ==0)
     {
           Message("strip_depth is the 3D mesh size in Y-direction. You will be asked to enter a value.");
           strip_depth =  Window.GetNumber("Mesh Size - Y-Direction?","How big would you like the mesh to be in the Y-Direction?");
           all_params_defined = 0;
     }     
     

     if(found_max_ele_width ==0)
     {
           //Message("found_max_ele_width not defined, setting this value to found_width/10");
           //found_max_ele_width = found_width/10;
           //all_params_defined = 0;
     }

     if(mesh_max_ele_width ==0)
     {
           Message("mesh_max_ele_width not defined. You will be asked to enter a value.");
           mesh_max_ele_width = Window.GetNumber("Element Size?","Enter max element size?");
           all_params_defined = 0;
     }
     
     if(mesh_max_ele_width2 ==0) mesh_max_ele_width2 = mesh_max_ele_width;

//     if(presc_disp_z ==0)
//     {
//           Message("presc_disp_z not defined, setting this value to -0.05");
//           presc_disp_z = -0.05;
//           all_params_defined = 0;
//     }  

     if(disp_time_0 < 0)
     {
           Message("disp_time_0 not defined, setting this value to 0.0s");
           disp_time_0 = 0.0;
           all_params_defined = 0;
     }  

     if(disp_time_ramp <= 0)
     {
           Message("disp_time_ramp not defined, setting this value to 2s");
           disp_time_ramp = 2;
           all_params_defined = 0;
     }       
    
     if(term_time <=0)
     {
           Message("term_time not defined, setting this value to disp_time_0 + disp_time_ramp");
           term_time = disp_time_0 + disp_time_ramp;
           all_params_defined = 0;
     }  

     
     //found_width = found_width*0.5;
     //strip_width = strip_width*0.5;
     
      if(all_params_defined != 1) Message("Not all parameters have been defined manually. Please review assumptions made above!!!")
   
      //OKToContinue("Please review strip model input messages. OK to continue?");
}
/////////////////////////////////////////////////////////////////////////////////////
// 
//   foundation input
  
  //var T_Step_tssfac = 0.8;
  //var T_Step_dt2msf = 0.0;
  
  //var cont_bulk_TYPE = -1;
  
  //var cont_contact_IGNORE = 1;
  
  //var cont_energy_HGEN = 2;
  //var cont_energy_RWEN = 0;
  //var cont_energy_SLNTEN = 2;
  //var cont_energy_RYLEN = 2;
  
  //var cont_mpp_io_nodump = 1;
  
  //var cont_shell_ESORT = 1;
  
  //var cont_solid_ESORT = 1;
  
  //var database_dt = 5.0E-3;
  //var database_d3dt = 1.0E-1;
  //var database_BINARY = 2;
  
  //var database_extent_NEIPH = 0;
  //var database_extent_NEIPS = 0;
  //var database_extent_MAXINIT = 0;
  //var database_extent_STRFLAG = 1;
  //var database_format_IBINARY = 1;
  
/////////////////////////////////////////////////////////////////////////////////////
  var found_nele = rdup(found_width/found_max_ele_width,0);                     // number of elements in foundation
  var refine_nele = found_nele*2;
  var found_ref_nele = found_nele + refine_nele;
  
  var found_ele_width = found_width/found_nele;
  var refine_width = (found_nele + refine_nele)*found_ele_width;
  
  var nele_mesh = rdup((strip_width - refine_width)/mesh_max_ele_width,0);
  var mesh_ele_width = (strip_width - refine_width)/nele_mesh;
  var tot_nele_width = found_nele + refine_nele + nele_mesh;
  
  //Message(found_nele + ", " + refine_nele + ", " + nele_mesh);
  
  var curve_x0 = disp_time_0;
  var curve_h0 = 0;
  var curve_x1 = disp_time_0 + disp_time_ramp;
  var curve_h1 = 1;
 
///
/////////////////////////////////////////////////////////////////////////////////////
//  Echo the data to dialog box
//Message("num_param = "+num_param);
//for (k=0; k<max_param+1; k++) Message("Param "+k+" Column="+param_column[k]);
//for (k=0; k<num_param; k++) Message("Column "+k+" Param="+column_param[k]);

for (var i=0; i<=max_param; i++)
{
	if (param_column[i]>-1)
	{
		//Message("Layer "+colname[i]);
		//for (var k=0; k<num_layers; k++) Message(k+1 + "     " + layer_data[k][i]);
	}
}
fInputFile.Close();

var flag = AllocateFlag();
model.ClearFlag(flag);

//   Select or create shell parts to extrude
if (create_column==1)
{
	var n = new Node(model,1,0,0,0)
	var h = new History(model, History.NODE, n.nid, "z=0");
	Message("Created first hist node");
	var n = new Node(model,2,2,0,0)
	var n = new Node(model,3,2,2,0)
	var n = new Node(model,4,0,2,0)
	var p = new Part(model,1,1,1,"Dummy shell");
	p.SetFlag(flag);
	var s = new Shell(model,1,1,1,2,3,4);

	var ns1 = new Set(model,1,Set.NODE,"Side nodes z=0");
	for (var i=1; i<=4; i++) ns1.Add(i);
           
}
else if (create_column==2)
{
      
      //var h = new History(model, History.NODE, 1, "z=0");
      var p = new Part(model,1,1,1,"Dummy shell");
	p.SetFlag(flag);
      
      if(strip_full)
      {
            i_0 = -tot_nele_width;
            bound_fac = 2.0;
      }else{
            i_0 = 0;
            bound_fac = 1.0;
      }
      
      var first_n = 1;
      var inode = 1;
      var ishl = 1;
      var x_co = 0;
      var ns0 = new Set(model,1,Set.NODE,"Foundation");
           
      for(i=i_0; i<=tot_nele_width;i++)
      {

            
            if(i >= -found_ref_nele && i < found_ref_nele)
            {
                  dx_co = found_ele_width;
            }
            else
            {
                  dx_co = mesh_ele_width;
            }
            
            if(Math.abs(i)<=found_nele)
            {
                  ns0.Add(inode);
                  ns0.Add(inode+1);
            }
            
            var n = new Node(model,inode,x_co,0,0);
            var n = new Node(model,inode+1,x_co,strip_depth,0);
            x_co = x_co + dx_co;
            if(first_n)
            {
                  first_n = 0;
            }
            else
            {
                  var n1 = inode-2;
                  var n2 = inode;
                  var n3 = inode+1;
                  var n4 = inode-1;
                  var s = new Shell(model,ishl,1,n1,n2,n3,n4);
                  ishl++;
            }
            inode=inode+2;
            //Message("x_co " + x_co + " dx_co " + dx_co);
      }
    
      var ns1 = new Set(model,2,Set.NODE,"Side nodes z=0");
      for (var i=1; i<inode; i++) ns1.Add(i);
      
      var label = Curve.NextFreeLabel(model);
      var label2 = label + 1;

      var curve_func = "STEP(TIME," + curve_x0 + "," + curve_h0 + "," + curve_x1 + "," + curve_h1 + ")";
      var curve_func2 = "0.0";
      //Message("curve function " + curve_func);
      var l = new Curve(Curve.CURVE_FUNCTION, model, label, 0, curve_func);
      var l2 = new Curve(Curve.CURVE_FUNCTION, model, label2, 0, curve_func2);

      if(presc_disp_z_sw)
      {
            var label = PrescribedMotion.NextFreeLabel(model);
            var b = new PrescribedMotion(model, ns0.sid, 3, 2, label, PrescribedMotion.SET,label, "Presc vert disp");
            var b2 = new Spc(model, ns0.sid,0,1,0,0,0,0,0,Spc.SET);
            b.sf = presc_disp_z;
      }
      if(presc_disp_x_sw)
      {
            var label = PrescribedMotion.NextFreeLabel(model);
            var b = new PrescribedMotion(model, ns0.sid, 1, 2, label, PrescribedMotion.SET,label, "Presc horiz disp");
            //var b2 = new Spc(model, ns0.sid,0,0,0,1,0,0,0,Spc.SET);
            b.sf = presc_disp_x;
      }      
      
      
}else if (create_column==3)
{
      
      //var h = new History(model, History.NODE, 1, "z=0");

      var tol = 0.01;   //corner points clsoer than this tolerane will be ignored.

      //start by making sure all shell normals point downwards
      //assume we only have horizontal shells
      var s = Shell.First(model);
      var nvector = s.NormalVector();
      if(nvector[2]>0) s.ReverseNormal();
      var f1 = AllocateFlag();
      model.ClearFlag(f1);
      Shell.FlagAll(model, f1);
      Shell.MakeConsistentNormalsFlagged(model, f1, s.eid);
      ReturnFlag(f1);      

      View.Redraw();

      var shell_pid = Part.NextFreeLabel(model);
      var shell_sid = Section.NextFreeLabel(model);
      var shell_mid = Material.NextFreeLabel(model);

      ///////////////////////////////////////////////////////////////////////////////

      Message("PRIMER-Automesh");

      Graphics.DrawingFunction(DrawEdges);

      //assume square Mesh
      var MeshSiz = strip_width;          //Window.GetNumber("Mesh Size?","How Big Would you like the mesh to be?");
      var MeshSizX = strip_width;          //Window.GetNumber("Mesh Size?","How Big Would you like the mesh to be?");
      var MeshSizY = strip_depth;          //Window.GetNumber("Mesh Size?","How Big Would you like the mesh to be?");
      var elsize = mesh_max_ele_width;    //Window.GetNumber("Element Size?","How Big Would you like the elements to be?");

      // Ask user to pick surface edges

      var outer_nodes = new Array();
      var inner_nodess = new Array();

      for (var outer_pick; (outer_pick = Node.Pick(
                  "Pick corner node outer edge")) != null; )
      {	
            outer_nodes = outer_pick.GetFreeEdgeNodes();
            break; 
      }

      if (outer_nodes.length == 0)
      {
            Graphics.DrawingFunction(null);
            View.Redraw();
            Message("Aborted");
      }

      // Find the characteristic length (cl) as the maximum of the distance 
      // to the adjacent nodes
      // also find unit vector of each edge and length of each edge

      var edge_vecs = new Array();

      var outer_node_cls = new Array();
      {
            var nodes = outer_nodes.slice();       // shallow copy
            nodes.unshift(nodes[nodes.length-1]);  // insert last as first
            nodes.push(nodes[1]);                  // add (original) first as last

            for (var i = 1; i < nodes.length-1; i++)
            {
                  outer_node_cls.push(Math.max(distance(nodes[i-1], nodes[i]),  
                              distance(nodes[i], nodes[i+1])));

                  edge_vecs.push([i,vect(nodes[i], nodes[i+1])]);
                              
            }
      }

      var edge_vecs_s = edge_vecs.slice();       // shallow copy
      edge_vecs_s.unshift(edge_vecs_s[edge_vecs_s.length-1]);  // insert last as first

      //find corners
      var corners = new Array();
      var corner_co = new Array();
      var corner_co_x = new Array();
      var corner_co_y = new Array();
      var corner_co_z = new Array();

      for(i=1;i<edge_vecs_s.length;i++)
      {
            dot_p = DotProd(edge_vecs_s[i-1][1],edge_vecs_s[i][1]);
            if(dot_p <0.95)
            {
                  edge_0 = edge_vecs_s[i][0];
                  nodes[edge_0].Sketch(true);
                  corners.push(edge_0);
                  corner_co.push([nodes[edge_0].x,nodes[edge_0].y]);
                  corner_co_x.push(nodes[edge_0].x);
                  corner_co_y.push(nodes[edge_0].y);
                  corner_co_z.push(nodes[edge_0].z);
            }
      }

      corners.push(corners[0]);

      Message("Found " + (corners.length-1) + " Corners");

      // sort ascending
      corner_co_x.sort(compareNumbersAsc);
      corner_co_y.sort(compareNumbersAsc);
      corner_co_z.sort(compareNumbersAsc);

      // calc boundaries of model
      var min_x = corner_co_x[0];
      var max_x = corner_co_x[corner_co_x.length-1];
      var min_y = corner_co_y[0];
      var max_y = corner_co_y[corner_co_y.length-1];
      var min_z = corner_co_z[0];
      var max_z = corner_co_z[corner_co_z.length-1];

      var Cen_x = 0.5*(min_x+max_x);
      var Cen_y = 0.5*(min_y+max_y);
      var Cen_z = 0.5*(min_z+max_z);
      var shl_z = layer_data[0][1];

      var bound_x_min = Cen_x - 0.5*MeshSizX;
      var bound_x_max = Cen_x + 0.5*MeshSizX;
      var bound_y_min = Cen_y - 0.5*MeshSizY;
      var bound_y_max = Cen_y + 0.5*MeshSizY;

      // check for duplicates within tol
      var corner_co_x1 = new Array();
      var corner_co_y1 = new Array();
      corner_co_x1.push(bound_x_min);
      corner_co_y1.push(bound_y_min);

      if(mesh_max_ele_window>0)
      {
      corner_co_x1.push(min_x-mesh_max_ele_window);
      corner_co_y1.push(min_y-mesh_max_ele_window); 
      }       
      
      Message(corner_co_x[0]);
      
      corner_co_x1.push(corner_co_x[0]);
      corner_co_y1.push(corner_co_y[0]);

      for(i=1; i<corner_co_x.length; i++)
      {
            if(Math.abs(corner_co_x[i-1]-corner_co_x[i])> tol)
            {
                  corner_co_x1.push(corner_co_x[i]);
            }
      }
      for(i=1; i<corner_co_y.length; i++)
      {
            if(Math.abs(corner_co_y[i-1]-corner_co_y[i])> tol)
            {
                  corner_co_y1.push(corner_co_y[i]);
            }
      }

      if(mesh_max_ele_window>0)
      {
      corner_co_x1.push(max_x+mesh_max_ele_window);
      corner_co_y1.push(max_y+mesh_max_ele_window); 
      }       
      
      corner_co_x1.push(bound_x_max);
      corner_co_y1.push(bound_y_max);

      //add additional coordinates between corner points - limited by max element size
      var corner_co_x2 = new Array();
      var corner_co_y2 = new Array();
      corner_co_x2.push(corner_co_x1[0]);
      corner_co_y2.push(corner_co_y1[0]);

      for(i=1; i<corner_co_x1.length; i++)
      {
            var tmpx = new Array();
            if (corner_co_x1[i-1] < min_x-(mesh_max_ele_window+1)) 
            {
                  elsize0 = elsize
                  elsize1 = mesh_max_ele_width2;
                  xoff = 0;
                  fact = -1;
                  lensw = 1;
                  //tmpx.push(corner_co_x1[i]);
                  
            }else if ( corner_co_x1[i] > max_x+(mesh_max_ele_window+1))
            {
                  elsize0 = elsize
                  elsize1 = mesh_max_ele_width2;
                  xoff = 1;
                  fact = 1;
                  lensw = 1;
                  
            }else{
                  elsize0 = elsize ;
                  elsize1 = elsize ;
                  xoff = 1;
                  fact = 1;
                  lensw = 0;                  
            }                  
            
            elsize2 = 0.5*(elsize0+elsize1);
            dist = corner_co_x1[i] - corner_co_x1[i-1];
            Nele = Math.max(rdup(dist/elsize2,0),1);
            if(Nele>1)
            {
                  dist1 = lensw*elsize0 + (1-lensw)*dist/Nele;  //Math.min(dist/Nele,elsize0);
                  if (elsize0==elsize1) Ls = 0;
                  else Ls = (dist-Nele*dist1)/(0.5*(Nele*Nele-Nele));
                  
                  xco = corner_co_x1[i-xoff];
                  Nele=Nele-1;
                  for(ii=0;ii<Nele;ii++)
                  {
                        esiz = dist1 + ii*Ls;
                        xco = xco + fact*esiz;
                        tmpx.push(xco);
                  }
            }
            tmpx.push(corner_co_x1[i]);
            tmpx.sort(compareNumbersAsc);
            for(ii=0;ii<tmpx.length;ii++) corner_co_x2.push(tmpx[ii]);
      }

      for(i=1; i<corner_co_y1.length; i++)
      {
            var tmpx = new Array();
            if (corner_co_y1[i-1] < min_y-(mesh_max_ele_window+1)) 
            {
                  elsize0 = elsize
                  elsize1 = mesh_max_ele_width2;
                  xoff = 0;
                  fact = -1;
                  lensw = 1;
                  //tmpx.push(corner_co_x1[i]);
                  
            }else if ( corner_co_y1[i] > max_y+(mesh_max_ele_window+1))
            {
                  elsize0 = elsize
                  elsize1 = mesh_max_ele_width2;
                  xoff = 1;
                  fact = 1;
                  lensw = 1;
                  
            }else{
                  elsize0 = elsize ;
                  elsize1 = elsize ;
                  xoff = 1;
                  fact = 1;
                  lensw = 0;                  
            }                  
            
            elsize2 = 0.5*(elsize0+elsize1);
            dist = corner_co_y1[i] - corner_co_y1[i-1];
            Nele = Math.max(rdup(dist/elsize2,0),1);
            if(Nele>1)
            {
                  dist1 = lensw*elsize0 + (1-lensw)*dist/Nele;  //Math.min(dist/Nele,elsize0);
                  if (elsize0==elsize1) Ls = 0;
                  else Ls = (dist-Nele*dist1)/(0.5*(Nele*Nele-Nele));
                  
                  xco = corner_co_y1[i-xoff];
                  Nele=Nele-1;
                  for(ii=0;ii<Nele;ii++)
                  {
                        esiz = dist1 + ii*Ls;
                        xco = xco + fact*esiz;
                        tmpx.push(xco);
                  }
            }
            tmpx.push(corner_co_y1[i]);
            tmpx.sort(compareNumbersAsc);
            for(ii=0;ii<tmpx.length;ii++) corner_co_y2.push(tmpx[ii]);            
      }
            
      //create nodes
      var nlist = new Array();
      for(i=0; i<corner_co_x2.length; i++)
      {
            xco = corner_co_x2[i];
            nlist[i] = new Array();
            
            for(j=0; j<corner_co_y2.length; j++)
            {
                  yco = corner_co_y2[j];
                  var label = Node.NextFreeLabel(model);
                  var n = new Node(model, label, xco, yco, shl_z);
                  nlist[i].push(label);
            }
      }

      //create shells
      for(i=1; i<nlist.length; i++)
      {
            for(j=1; j<nlist[i].length; j++)
            {
                  var label = Shell.NextFreeLabel(model);
                  N1 = nlist[i-1][j-1];
                  N2 = nlist[i][j-1];
                  N3 = nlist[i][j];
                  N4 = nlist[i-1][j];
                  
                  var s = new Shell(model, label, shell_pid, N1, N2, N3, N4);

            }
      }

     
      
      var p = new Part(model,shell_pid,shell_sid,shell_mid,"Dummy shell");
	p.SetFlag(flag); 

      Message("Done creating surface")       
             
}
else
{
	Part.Select(flag,"Select shell parts to extrude");
}


if(mat_hys_cur)
{
      for(icur=0; icur<=Cur_Cnt;icur++)
      {
            label = Cur_Array[icur][0];
            cur_tit = Cur_Array[icur][1];
            x_fac = Cur_Array[icur][2];
            y_fac = Cur_Array[icur][3];                  
            var l = new Curve(Curve.CURVE, model, label);
            l.heading = cur_tit;
            l.sfa = x_fac;
            l.sfo = y_fac;                  

            for (ip=4; ip<Cur_Array[icur].length;ip++)
            {
                  l.AddPoint(Cur_Array[icur][ip][0],Cur_Array[icur][ip][1]);
            }
      }
}


model.PropagateFlag(flag);

//   Find z-coord of soil surface, max node ID, etc
//var n = Node.First(model);
var found = 0;
var nid_next = 1;
var nid_min = 999999999;
var nid_max = -1;
var num_flagged = 0;
var n_all = Node.GetAll(model);
var flagged_nodes = new Array();
for (i=0; i<n_all.length; i++)
//while (n)
{
	nid_next = Math.max(nid_next, n_all[i].nid+1);
	if (n_all[i].Flagged(flag)) 
	{
		flagged_nodes.push(i);
            found = 1;
		z_surface = n_all[i].z;
		nid_min = Math.min(n_all[i].nid,nid_min);
		nid_max = Math.max(n_all[i].nid,nid_max);
		num_flagged++;
	}
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
		
		//Message("elsize="+elsize+"; dz="+dz);
		nsub = Math.round(dz/elsize);                
		nsub = Math.max(1,nsub);
		dzsub = dz/nsub;
		znext = z_prev - dzsub;
		//Message("Layer "+i+" has "+nsub + " sublayers");
		for (var k=0; k<nsub; k++)
		{
			layer_sublayer[num_sublayers] = i;
			zsublayer[num_sublayers] = znext;
			//Message("...sublayer "+num_sublayers+", z="+znext);
			num_sublayers++;
			znext = znext - dzsub;
		}
	}
	else
	{
		//Message("No sublayers - layer "+i+", z="+zlay);
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

//find min and max coords
var lay_nod_lst_bound = new Array();
var flagged_nodes_bound = new Array();
var flagged_nodes_corner = new Array();
var min_x = 9999999999;
var max_x = -9999999999;
var min_y = 9999999999;
var max_y = -9999999999;    


for (i=0; i<flagged_nodes.length; i++)
//var n = Node.First(model);
//while (n)
{
      min_x = Math.min(min_x,n_all[flagged_nodes[i]].x);
      max_x = Math.max(max_x,n_all[flagged_nodes[i]].x);
      min_y = Math.min(min_y,n_all[flagged_nodes[i]].y);
      max_y = Math.max(max_y,n_all[flagged_nodes[i]].y);
      flagged_nodes_bound[i] = -1;
      flagged_nodes_corner[i] = -1;
}

var label_s = Set.NextFreeLabel(model, Set.NODE);
var s = new Set(model, label_s, Set.NODE);
s.title = "Bound_Nodes_layer_0";
lay_nod_lst_bound.push(s.sid);

var label_s = Set.NextFreeLabel(model, Set.NODE);
var sl1 = new Set(model, label_s, Set.NODE);
sl1.title = "Corner_1";

var label_s = Set.NextFreeLabel(model, Set.NODE);
var sl2 = new Set(model, label_s, Set.NODE);
sl2.title = "Corner_2";

var label_s = Set.NextFreeLabel(model, Set.NODE);
var sl3 = new Set(model, label_s, Set.NODE);
sl3.title = "Corner_3";

var label_s = Set.NextFreeLabel(model, Set.NODE);
var sl4 = new Set(model, label_s, Set.NODE);
sl4.title = "Corner_4";

for (i=0; i<flagged_nodes.length; i++)
//var n = Node.First(model);
//while (n)
{
      dx0 = Math.abs(min_x-n_all[flagged_nodes[i]].x);
      dx1 = Math.abs(max_x-n_all[flagged_nodes[i]].x);
      dy0 = Math.abs(min_y-n_all[flagged_nodes[i]].y);
      dy1 = Math.abs(max_y-n_all[flagged_nodes[i]].y);
      min_d = Math.min(dx0,dx1,dy0,dy1);
      if(min_d <= 0.001)
      {
            s.Add(n_all[flagged_nodes[i]].nid);
            flagged_nodes_bound[i] = 1;
      }
      if(dx0 <= 0.001 && dy0 <= 0.001)
      {
            sl1.Add(n_all[flagged_nodes[i]].nid);
            flagged_nodes_corner[i] = 1;
      } 
      if(dx0 <= 0.001 && dy1 <= 0.001)
      {
            sl2.Add(n_all[flagged_nodes[i]].nid);
            flagged_nodes_corner[i] = 2;
      }
      if(dx1 <= 0.001 && dy1 <= 0.001)
      {
            sl3.Add(n_all[flagged_nodes[i]].nid);
            flagged_nodes_corner[i] = 3;
      } 
      if(dx1 <= 0.001 && dy0 <= 0.001)
      {
            sl4.Add(n_all[flagged_nodes[i]].nid);
            flagged_nodes_corner[i] = 4;
      }         
}

//  create new nodes

for (i=0; i<num_sublayers_n; i++)
{
      var label_s = Set.NextFreeLabel(model, Set.NODE);
      var s = new Set(model, label_s, Set.NODE);
      s.title = "Bound_Nodes_layer_" + (i+1);
      lay_nod_lst_bound.push(s.sid);
      
	var first = 1;
      for (j=0; j<flagged_nodes.length; j++)
      {
            var n2 = new Node(model,n_all[flagged_nodes[j]].nid+nid_inc[i],n_all[flagged_nodes[j]].x,n_all[flagged_nodes[j]].y,zsublayer[i]);
            if (first && create_column==1) var h = new History(model, History.NODE, n_all[flagged_nodes[j]].nid+nid_inc[i], "z="+zsublayer[i]);
            first = 0;
            
            if(flagged_nodes_bound[j]>0)
            {
                  s.Add(n2.nid);
            }
            if(flagged_nodes_corner[j]==1)
            {
                  sl1.Add(n2.nid);
            }
            if(flagged_nodes_corner[j]==2)
            {
                  sl2.Add(n2.nid);
            } 
            if(flagged_nodes_corner[j]==3)
            {
                  sl3.Add(n2.nid);
            } 
            if(flagged_nodes_corner[j]==4)
            {
                  sl4.Add(n2.nid);
            }             
	}
}


// check start ID for solids
var solid_next = 1;
var s = Solid.First(model);
while (s)
{
	solid_next = Math.max(solid_next, s.eid);
	s = s.Next();
}
solid_next++;

//  Check start ID for Mat and Part;
var pid_next = 1;
var p = Part.First(model);
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
	var p = Material.First(model);
	while (p)
	{
		pid_next = Math.max(pid_next, p.mid);
		p = p.Next();
	}
	var p =	 Section.First(model);
	while (p)
	{
		pid_next = Math.max(pid_next, p.secid);
		p = p.Next();
	}
	var p =	 Hourglass.First(model);
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

if(pid_spec == 1) var sid = layer_data[0][16];
else var sid = pid_next;

var sec = new Section(model, sid, Section.SOLID, "Soil")
sec.elform = 1;

if(pid_spec == 1) var hgid = layer_data[0][16];
else var hgid = pid_next;
var hg = new Hourglass(model,hgid);
hg.ihq = 4;
if (create_column == 2) hg.qm = 1e-5;
else hg.qm = 0.001;

var s_pt_id = Set.NextFreeLabel(model, Set.PART);
var s_pt = new Set(model, s_pt_id, Set.PART);
s_pt.title = "soil_parts"

//Message(pid_spec + "   " + hgid);

for (i=0; i<num_layers; i++)
{
	var mt = Math.round(layer_data[i][3]);                     // create material for layer
	if (mt == 79) var mattype = "*MAT_HYSTERETIC_SOIL";
	else if (mt == 1) var mattype = "*MAT_ELASTIC";
	else if (mt == 230) var mattype = "*MAT_PML_ELASTIC";
	else ErrorTerminate ("Sorry, cannot handle material type "+mt);
      
      var create_part = 1;

      if(pid_spec == 1)
      {
            pid1 = layer_data[i][16];
            var p = Part.GetFromID(model, pid1);

            if(p)
            {
                  create_part = 0;
                  if (mt==79)
                  {
                        var mat = Material.GetFromID(model, pid1);
                        var K0 = mat.GetPropertyByName("K0");
                        var LC = mat.GetPropertyByName("LCID");
                        var SF = mat.GetPropertyByName("SFLC");
                        
                        check_par(i, pid1, "K0", K0, layer_data[i][7]);
                        check_par(i, pid1, "LCID", LC, Math.round(layer_data[i][9]));
                        check_par(i, pid1, "SFLC", SF, layer_data[i][10]);
                  }
                  
            }else{
                  var mat = new Material(model, pid1, mattype);
                  var partname = "Soil layer "+(i+1);
                  if (p_desc) partname = partname + " " + layer_data[i][17];
                  var p = new Part(model, pid1, sid, mat.mid, partname);
                  p.hgid = hgid;
                  mat.title = partname;
                  if (p_desc) 
                  {
                        p.colour = GetColour(layer_data[i][17]);
                        if(mat_hys_cur)
                        {
                              if(c = Curve.GetFromID(model, layer_data[i][9])) c.heading = partname;;   
                        }
                  }
                  s_pt.Add(pid1);
            }
            
      }else{      
            var mat = new Material(model, pid_min+pid_inc[i], mattype);
            var psh = Part.First(model);                  // create parts for layer
            while (psh)
            {
                  if (psh.Flagged(flag))
                  {
                        pid_new = psh.pid + pid_inc[i];
                        i1 = i+1;
                        var partname = "Soil layer "+i1;
                        if (np_flagged > 1) partname = partname + " " + psh.heading;
                        var p = new Part(model, pid_new, sid, mat.mid, partname);
                        p.hgid = hgid;
                        s_pt.Add(pid_new);
                  }
                  psh = psh.Next();
            }            
      }
	      
      if(create_part)
      {
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
      }


}

//pid_next = pid_max + pid_inc[num_layers-1] + 1;
pid_next = Part.NextFreeLabel(model);

if(Math.abs(bdamp) > 0) 
{
      var nar = new Array([s_pt_id,-10],[bdamp,-10])
      var fname = "./bdamp.k"
      var f = new File(fname, File.WRITE);
      f.Writeln("*KEYWORD");
      write_KW2(f, true, "*DAMPING_PART_STIFFNESS_SET", nar);
      f.Writeln("*END");
      f.Close();
      
      if (model.layer > 0)
      {
            var curLayer = Include.GetFromID(model, model.layer);
            model.ImportInclude(fname, curLayer);     // Import data from file into model
            File.Delete(fname);                       // Delete temporary file
            curLayer.MakeCurrentLayer();
      }
      else
      {
            model.Import(fname);                 // Import data from file into model
            File.Delete(fname);                  // Delete temporary file
      }
}


if (lysmer_ro > 0.0)
{
	var mattype = "*MAT_DAMPER_VISCOUS"
	var mat = new Material(model, pid_next, mattype);
	mat.SetPropertyByName("DC",  1.0            );
	
	var sec = new Section(model, pid_next, Section.DISCRETE, "Lysmer dampers")

	var p = new Part(model, pid_next, pid_next, pid_next, "Lysmer dampers");
      p.colour = GetColour('LYSMER');
	pid_lys = pid_next;
	//pid_next++;
      pid_next = Part.NextFreeLabel(model);
}

// Extrude shells to create solids
var s_all = Shell.GetAll(model);
var s_flag = new Array();
for (i=0; i<s_all.length; i++)
{
      if (s_all[i].Flagged(flag))
      {
            s_flag.push(i);
      }
}

for (i=0; i<num_sublayers; i++)
{
	inc1 = 0;
	if (i>0) inc1 = nid_inc[i-1];
	inc2 = nid_inc[i] - inc1;
	ilay = layer_sublayer[i];
      
      for (j=0; j<s_flag.length; j++)
      {

            if(s_all[s_flag[j]].nodes == 3)
            {
                  if(pid_spec == 1)
                  {
                        var pid = layer_data[ilay][16];
                  }else{
                        var pid = s_all[s_flag[j]].pid + pid_inc[ilay];
                  }
                  n4 = s_all[s_flag[j]].n1 + inc1;
                  n3 = s_all[s_flag[j]].n3 + inc1;
                  n6 = s_all[s_flag[j]].n2 + inc1;

                  n1 = n4 + inc2;
                  n2 = n3 + inc2;
                  n5 = n6 + inc2;

                  var so = new Solid(model,solid_next,pid,n1,n2,n3,n4,n5,n6);
                  solid_next++;
            }
            else if(s_all[s_flag[j]].nodes == 4)
            {
                  if(pid_spec == 1)
                  {
                        var pid = layer_data[ilay][16];
                  }else{
                        var pid = s_all[s_flag[j]].pid + pid_inc[ilay];
                  }
                  
                  n5 = s_all[s_flag[j]].n1 + inc1;
                  n6 = s_all[s_flag[j]].n2 + inc1;
                  n7 = s_all[s_flag[j]].n3 + inc1;
                  n8 = s_all[s_flag[j]].n4 + inc1;

                  n1 = n5 + inc2;
                  n2 = n6 + inc2;
                  n3 = n7 + inc2;
                  n4 = n8 + inc2;

                  var so = new Solid(model,solid_next,pid,n1,n2,n3,n4,n5,n6,n7,n8);
                  solid_next++;
            }
      }
	
}

///////////////////////////////////////////
// Remove solids above/within foundations.

var foundation_z = Cen_z;           //nodes[corners[0]].z;
Message("Foundation Z Coord: "+ foundation_z);


// Get solid elements which are above the foundations
var solids = Solid.GetAll(model);
var upper_solids = new Object;
for(i=0; i<solids.length; i++)
{
	var p = Part.GetFromID(model, solids[i].pid)
	if (isStringMatch(p.heading.toUpperCase(), "SOIL"))
	{
		var centroid = SolCoords(solids[i],model);
		  
		// Replace z-coordinate with top surface elevation
		centroid[2] = TopCoord(solids[i],model)

		if(centroid[2] > foundation_z)
		{
			upper_solids[solids[i].eid] = solids[i];
			upper_solids[solids[i].eid].centroid = centroid;
		}
	}
}

var del_flag = AllocateFlag();
for(i in upper_solids)
{
	var cent = upper_solids[i].centroid
	var count = 0; 
	for(j=0; j<=corners.length-2; j++)
	{
		n1 = nodes[corners[j]];
		n2 = nodes[corners[j+1]];

		if(cent[1] > n1.y && cent[1] < n2.y || cent[1] > n2.y && cent[1] < n1.y ) //centroid is between the line ends
		{
			if(cent[0] > n1.x)
			{
				count = count + 1;
			}
		}
	}
	if(count%2 != 0) // even - inside shape
	{
		upper_solids[i].SetFlag(del_flag);
	}
}
model.DeleteFlagged(del_flag);
ReturnFlag(del_flag);

///////////////////////////////////////////

//   Find node set ID
var set_next = 1;      
var set = Set.First(model, Set.NODE);
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
	var set = Set.First(model, Set.SEGMENT);
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
	var s = Shell.First(model);
	while (s)
	{
		if (s.Flagged(flag))
		{
			if(s.nodes == 3)
			{
				n1 = s.n1 + inc1;
				n2 = s.n2 + inc1;
				n3 = s.n3 + inc1;
				f.Writeln(NumberToString(n1,10) + NumberToString(n2,10) 
					+ NumberToString(n3,10) );      //  Segment
			}
			else if(s.nodes == 4)
			{
				n1 = s.n1 + inc1;
				n2 = s.n2 + inc1;
				n3 = s.n3 + inc1;
				n4 = s.n4 + inc1;
				f.Writeln(NumberToString(n1,10) + NumberToString(n2,10) 
					+ NumberToString(n3,10) + NumberToString(n4,10));      //  Segment
			}
		}
		s = s.Next();
	}

	f.Writeln("*END");
      f.Close();                                  	    // Close temporary file
      
      if (model.layer > 0)
      {
            var curLayer = Include.GetFromID(model, model.layer);
            model.ImportInclude("prtmp_setseg.key", curLayer);    // Import data from file into model      
            File.Delete("./prtmp_setseg.key");                    // Delete temporary file
            curLayer.MakeCurrentLayer();
      }
      else
      {
            model.Import("prtmp_setseg.key");               // Import data from file into model
            File.Delete("./prtmp_setseg.key");              // Delete temporary file
      }

//   Set of nodes for *LOAD_SEISMIC_SSI
	var ns = new Set(model, set_next, Set.NODE,"Nodes for *LOAD_SEISMIC_SSI");
	var n = Node.First(model);
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
	if (spc_sides != 0) var codes = BreakdownSPC(Math.abs(spc_sides));

//   Select the "master" set on the shells
        if (!create_column)
	{
		Set.Select(Set.NODE,flag,"Select node set on edges of shells for SPC or BPM",model,false);
		var ns0 = Set.First(model, Set.NODE);
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
            
      //   Create SPC, BPM for the master node set
            CreateSPC_BPM(ns1.sid, 0); 

            //   Create new node set at each sublayer of nodes
            for (i=0; i<num_sublayers; i++)
            {
                  inc = nid_inc[i];
                  var ns = new Set(model, set_next, Set.NODE, "Side nodes, z="+zsublayer[i]);
                  ns1.StartSpool();
                  while (id = ns1.Spool() ) ns.Add(id+inc);  // copy master set using known node ID increment

                  ilay = layer_sublayer[i]+1;     // node inc for sublayer i are at bottom of element sublayer
                  //Message("i="+i+", ilay="+ilay+" num_layers="+num_layers);
                  if (layer_data[ilay]) CreateSPC_BPM(set_next,ilay);          //  create SPC and/or BPM
                  else if (ilay==num_layers) CreateSPC_BPM(set_next,-1);     //  special treatment for sublayers within
                                                                               //      bottom layer

                  set_next++;
            }
            //OKToContinue();

            ns1.ClearFlag(flag);     //  So we don't delete the set from the column model

	}else{
            for (i=0; i<num_sublayers; i++)
            {
                  sid = lay_nod_lst_bound[i];
                  var label_nrb = NodalRigidBody.NextFreeLabel(model);
                  var nrb = new NodalRigidBody(model, sid, label_nrb);
                  nrb.spc = true;
                  nrb.cmo = 1;
                  nrb.con1 = spc_sides;
                  nrb.con2 = 7; 
            }
            
      }
}	

//   Set of nodes at bottom and/or Lysmer treatment
if (spc_bottom != 0 || lysmer_ro > 0.0)
{
	var ibot = num_sublayers_n-1;
	var ns = new Set(model, set_next, Set.NODE,"Nodes on bottom plane");
	var n = Node.First(model);
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
					var sp = new Discrete(model, solid_next, pid_lys, nsp1, nsp2, 1, area_vs);
					solid_next++;
					var sp = new Discrete(model, solid_next, pid_lys, nsp1, nsp2, 2, area_vs);
					solid_next++;
					if (lys_lcvx > 0) var l = new LoadNode(model, LoadNode.POINT, nsp1, 1, lys_lcvx, area_vs);
					if (lys_lcvy > 0) var l = new LoadNode(model, LoadNode.POINT, nsp1, 2, lys_lcvy, area_vs);
				}
				if (area_vp > 0.0)
				{
					var sp = new Discrete(model, solid_next, pid_lys, nsp1, nsp2, 3, area_vp);
					solid_next++;
					if (lys_lcvz > 0) var l = new LoadNode(model, LoadNode.POINT, nsp1, 3, lys_lcvz, area_vp);
				}
//				OKToContinue();
			}
		}
		n = n.Next();
	}

	if (lysmer_ro > 0.0) spc_bottom = 7;
	var codes = BreakdownSPC(spc_bottom);
	var spc = new Spc(model, set_next, 0, codes.x,codes.y,codes.z,0,0,0,Spc.SET);
	set_next++;
}

if (create_column==1)
{     
      //model.control.timestep.exists = true;
      //model.control.timestep.tssfac = T_Step_tssfac; 
      //model.control.timestep.dt2msf = T_Step_dt2msf;
      
      //model.control.energy.exists = true;
      //model.control.energy.hgen = cont_energy_HGEN;
      //model.control.energy.rwen = cont_energy_RWEN;
      //model.control.energy.slnten = cont_energy_SLNTEN;
      //model.control.energy.rylen = cont_energy_RYLEN;

      //if(cont_mpp_io_nodump) model.control.mpp_io_nodump.exists = true;
      
      //model.control.solid.exists = true;
      //model.control.solid.esort = cont_solid_ESORT;
      
      //model.database.elout.exists = true;
      //model.database.glstat.exists = true;
      //model.database.matsum.exists = true;
      //model.database.nodfor.exists = true;
      //model.database.nodout.exists = true;
      //model.database.sleout.exists = true;
      //model.database.binary.d3plot.exists = true;
      //model.database.binary.d3thdt.exists = true;
      //model.database.extent_binary.exists = true;
            
      //model.database.elout.dt = database_dt;
      //model.database.glstat.dt = database_dt;
      //model.database.matsum.dt = database_dt;
      //model.database.nodfor.dt = database_dt;
      //model.database.nodout.dt = database_dt;
      //model.database.sleout.dt = database_dt;
      //model.database.binary.d3plot.dt = database_d3dt;
      //model.database.binary.d3thdt.dt = database_dt;
      
      //model.database.elout.binary = database_BINARY;
      //model.database.glstat.binary = database_BINARY;
      //model.database.matsum.binary = database_BINARY;
      //model.database.nodfor.binary = database_BINARY;
      //model.database.nodout.binary = database_BINARY;
      //model.database.sleout.binary = database_BINARY;
      
      //model.database.extent_binary.neiph = database_extent_NEIPH;
      //model.database.extent_binary.neips = database_extent_NEIPS;
      //model.database.extent_binary.maxint = database_extent_MAXINIT;
      //model.database.extent_binary.strflg = database_extent_STRFLAG;

      //model.database.format.exists = true;
      //model.database.format.ibinary = database_format_IBINARY;      
}

if (create_column==2)
{
      
      var ns_X = new Set(model, set_next, Set.NODE,"Nodes X Plane");
      var ns_Y = new Set(model, set_next+1, Set.NODE,"Nodes Y Plane");
      
   	var n = Node.First(model);
	while (n)
	{  
            if(n.x == 0 || n.x == strip_width*bound_fac)
            {
                  ns_X.Add(n.nid);
            }else{
                  ns_Y.Add(n.nid);
            }
            n = n.Next();
      }
      
      var spc = new Spc(model, ns_X.sid, 0, 1,1,0,0,0,0,Spc.SET);
      var spc = new Spc(model, ns_Y.sid, 0, 0,1,0,0,0,0,Spc.SET);
      
      model.control.termination.exists = true;
      model.control.termination.endtim = term_time;
      
      //model.control.timestep.exists = true;
      //model.control.timestep.tssfac = T_Step_tssfac; 
      //model.control.timestep.dt2msf = T_Step_dt2msf; 

      //model.control.bulk_viscosity.exists = true;
      //model.control.bulk_viscosity.ibq = cont_bulk_TYPE;   

      //model.control.contact.exists = true;
      //model.control.contact.ignore = cont_contact_IGNORE;  
      
      //model.control.energy.exists = true;
      //model.control.energy.hgen = cont_energy_HGEN;
      //model.control.energy.rwen = cont_energy_RWEN;
      //model.control.energy.slnten = cont_energy_SLNTEN;
      //model.control.energy.rylen = cont_energy_RYLEN;

      //if(cont_mpp_io_nodump) model.control.mpp_io_nodump.exists = true;
      
      //model.control.shell.exists = true;
      //model.control.shell.esort = cont_shell_ESORT;
      
      //model.control.solid.exists = true;
      //model.control.solid.esort = cont_solid_ESORT;
      
      //model.database.elout.exists = true;
      //model.database.glstat.exists = true;
      //model.database.matsum.exists = true;
      //model.database.nodfor.exists = true;
      //model.database.nodout.exists = true;
      //model.database.sleout.exists = true;
      //model.database.binary.d3plot.exists = true;
      //model.database.binary.d3thdt.exists = true;
      //model.database.extent_binary.exists = true;
            
      //model.database.elout.dt = database_dt;
      //model.database.glstat.dt = database_dt;
      //model.database.matsum.dt = database_dt;
      //model.database.nodfor.dt = database_dt;
      //model.database.nodout.dt = database_dt;
      //model.database.sleout.dt = database_dt;
      //model.database.binary.d3plot.dt = database_d3dt;
      //model.database.binary.d3thdt.dt = database_dt;
      
      //model.database.elout.binary = database_BINARY;
      //model.database.glstat.binary = database_BINARY;
      //model.database.matsum.binary = database_BINARY;
      //model.database.nodfor.binary = database_BINARY;
      //model.database.nodout.binary = database_BINARY;
      //model.database.sleout.binary = database_BINARY;
      
      //model.database.extent_binary.neiph = database_extent_NEIPH;
      //model.database.extent_binary.neips = database_extent_NEIPS;
      //model.database.extent_binary.maxint = database_extent_MAXINIT;
      //model.database.extent_binary.strflg = database_extent_STRFLAG;

      //model.database.format.exists = true;
      //model.database.format.ibinary = database_format_IBINARY;
     
      //var c = new History(model, History.NODE, 1, "Foundation Node");
      var nfg = new NodalForceGroup(model, ns0.sid);
      if (create_column != 2)
      {
            var fl = AllocateFlag();
            ns0.SetFlag(fl);
            model.PropagateFlag(fl)
            var n = Node.GetFlagged(model,fl)[0];
            model.ClearFlag(fl);
            ReturnFlag(fl);
            var c = new History(model, History.NODE, n.nid, "Foundation Node");
      }
      
}

if(tstep)
{
      Message(tstepa);
      model.control.timestep.exists = true;
      model.control.timestep.tssfac =     tstepa[1];
      model.control.timestep.isdo =       tstepa[2];
      model.control.timestep.tslimt =     tstepa[3]; 
      model.control.timestep.dt2ms =      tstepa[4]; 
      model.control.timestep.erode =      tstepa[6]; 
      model.control.timestep.ms1st =      tstepa[7];

      if (tstepa[5] !=0)
      {
            var label = Curve.NextFreeLabel(model);
            var l = new Curve(Curve.CURVE, model, label);
            l.AddPoint(0, tstepa[5]);
            l.AddPoint(1000, tstepa[5]);
            model.control.timestep.lctm = label;
      }

}

if(con_energy)
{           
      model.control.energy.exists = true;
      model.control.energy.hgen = con_energy1[0];
      model.control.energy.rwen = con_energy1[1];
      model.control.energy.slnten = con_energy1[2];
      model.control.energy.rylen = con_energy1[3];
}      
 
if(con_sol)
{           
      model.control.solution.exists = true;
      model.control.solution.soln = con_sol1[0];
      model.control.solution.nlq = con_sol1[1];
      model.control.solution.isnan = con_sol1[2];
      model.control.solution.lcint = con_sol1[3];
} 
//----
if(con_mpp_io_nodump)
{           
      model.control.mpp_io_nodump.exists = true;
} 

if(con_solid)
{           
      model.control.solid.exists = true;
      model.control.solid.esort = con_solid1[0];
      model.control.solid.fmatrix = con_solid1[1];
      model.control.solid.niptets = con_solid1[2];
      model.control.solid.swlocl = con_solid1[3];
      model.control.solid.psfail = con_solid1[4];
      model.control.solid.t10jtol = con_solid1[5];
      model.control.solid.icohed = con_solid1[6];
} 

if(con_shell)
{           
      model.control.shell.exists = true;
      model.control.shell.wrpang = con_shell1[0];
      model.control.shell.esort = con_shell1[1];
      model.control.shell.irnxx = con_shell1[2];
      model.control.shell.istupd = con_shell1[3];
      model.control.shell.theory = con_shell1[4];
      model.control.shell.bwc = con_shell1[5];
      model.control.shell.miter = con_shell1[6];
      model.control.shell.proj = con_shell1[7];
}
  
if(con_bulk)
{           
      model.control.bulk_viscosity.exists = true;
      model.control.bulk_viscosity.q1 = con_bulk1[0];
      model.control.bulk_viscosity.q2 = con_bulk1[1];
      model.control.bulk_viscosity.ibq = con_bulk1[2];
      model.control.bulk_viscosity.btype = con_bulk1[3];

} 

if(con_cont)
{           
      model.control.contact.exists = true;
      model.control.contact.slsfac = con_cont1[0];
      model.control.contact.rwpnal = con_cont1[1];
      model.control.contact.islchk = con_cont1[2];
      model.control.contact.shlthk = con_cont1[3];
      model.control.contact.penopt = con_cont1[4];
      model.control.contact.thkchg = con_cont1[5];
      model.control.contact.orien = con_cont1[6];
      model.control.contact.enmass = con_cont1[7];
      model.control.contact.usrstr = con_cont1[8];
      model.control.contact.usrfrc = con_cont1[9];
      model.control.contact.nsbcs = con_cont1[10];
      model.control.contact.interm = con_cont1[11];
      model.control.contact.xpene = con_cont1[12];
      model.control.contact.ssthk = con_cont1[13];
      model.control.contact.ecdt = con_cont1[14];
      model.control.contact.tiedprj = con_cont1[15];
}  

if(data_bin_d3pl)
{           
      model.database.binary.d3plot.exists
      model.database.binary.d3plot.dt = data_bin_d3pl1[0];
}

if(data_bin_d3th)
{           
      model.database.binary.d3thdt.exists = true;
      model.database.binary.d3thdt.dt = data_bin_d3th1[0];
}

if(data_extent_bin)
{           
      model.database.extent_binary.exists = true;
      model.database.extent_binary.neiph = data_extent_bin1[0];
      model.database.extent_binary.neips = data_extent_bin1[1];
      model.database.extent_binary.maxint = data_extent_bin1[2];
      model.database.extent_binary.strflg = data_extent_bin1[3];
      model.database.extent_binary.sigflg = data_extent_bin1[4];
      model.database.extent_binary.epsflg = data_extent_bin1[5];
      model.database.extent_binary.rltflg = data_extent_bin1[6];
      model.database.extent_binary.engflg = data_extent_bin1[7];
      model.database.extent_binary.cmpflg = data_extent_bin1[8];
      model.database.extent_binary.ieverp = data_extent_bin1[9];
      model.database.extent_binary.beamip = data_extent_bin1[10];
      model.database.extent_binary.dcomp = data_extent_bin1[11];
      model.database.extent_binary.shge = data_extent_bin1[12];
      model.database.extent_binary.stssz = data_extent_bin1[13];
      model.database.extent_binary.n3thdt = data_extent_bin1[14];
      model.database.extent_binary.ialemat = data_extent_bin1[15];
}
 

if(data_elout)
{           
      model.database.elout.exists = true;
      model.database.elout.dt = data_elout1[0];
      model.database.elout.binary = data_elout1[1];
}

if(data_glstat)
{           
      model.database.glstat.exists = true;
      model.database.glstat.dt = data_glstat1[0];
      model.database.glstat.binary = data_glstat1[1];
} 

if(data_matsum)
{           
      model.database.matsum.exists = true;
      model.database.matsum.dt = data_matsum1[0];
      model.database.matsum.binary = data_matsum1[1];

} 

if(data_nodfor)
{           
      model.database.nodfor.exists = true;
      model.database.nodfor.dt = data_nodfor1[0];
      model.database.nodfor.binary = data_nodfor1[1];
} 

if(data_nodout)
{           
      model.database.nodout.exists = true;
      model.database.nodout.dt = data_nodout1[0];
      model.database.nodout.binary = data_nodout1[1];
} 

if(data_sleout)
{           
      model.database.sleout.exists = true;
      model.database.sleout.dt = data_sleout1[0];
      model.database.sleout.binary = data_sleout1[1];
} 

 if(data_curout)
{           
      model.database.curvout.exists = true;
      model.database.curvout.dt = data_curout1[0];
      model.database.curvout.binary = data_curout1[1];
} 

//  Write the damper orientation vectors for Lysmer
if (lysmer_ro > 0.0)
{
	f = new File("./lysmer_orient.key", File.WRITE);     // Open temporary file
	f.Writeln("*KEYWORD");                       // Write set keyword data to file
	f.Writeln("*DEFINE_SD_ORIENTATION");
	f.Writeln("         1         0       1.0       0.0       0.0");           
	f.Writeln("         2         0       0.0       1.0       0.0");           
	f.Writeln("         3         0       0.0       0.0       1.0");           
	f.Writeln("*END");
      f.Close();                                  	    // Close temporary file

      if (model.layer > 0)
      {
            var curLayer = Include.GetFromID(model, model.layer);
            model.ImportInclude("lysmer_orient.key", curLayer);   // Import data from file into model
            File.Delete("./lysmer_orient.key");                   // Delete temporary file
            curLayer.MakeCurrentLayer();
      }
      else
      {
            model.Import("lysmer_orient.key");                       // Import data from file into model
            File.Delete("./lysmer_orient.key");                  // Delete temporary file
      }
}

//  Column model - delete dummy shell
if (create_column) model.DeleteFlagged(flag);

if(foundation_depth>0)
{
      var fx1 = AllocateFlag();
      model.ClearFlag(fx1); 
      var fx2 = AllocateFlag();
      model.ClearFlag(fx2);
      var GL = 0;
      var FL = GL-foundation_depth;
      
      var ss1 = new Set(model,1,Set.SOLID,"lateral push elements");
      var ss2 = new Set(model,2,Set.SOLID,"vertical push elements");
      
      var sol1 = Solid.First(model);

      while(sol1)
      {
            var SCo = SolCoords(sol1,model);
            if (SCo[0]<= (strip_width + found_width) && SCo[2] >= FL)
            {
                  ss1.Add(sol1.eid);
            }else if(SCo[2] >= FL)
            {
                  ss2.Add(sol1.eid);
            }
            
            sol1=sol1.Next();
      }

      var ss1Flag = AllocateFlag();
      var ss2Flag  = AllocateFlag();
      model.ClearFlag(ss1Flag);
      model.ClearFlag(ss2Flag);
      ss1.SetFlag(ss1Flag);
      ss2.SetFlag(ss2Flag);
      model.PropagateFlag(ss1Flag);
      model.PropagateFlag(ss2Flag);

      if(presc_disp_z)
      {
            var sol_ele = Solid.GetFlagged(model, ss1Flag);
            var pids  = new Array();
            for ( var s = 0; s < sol_ele.length; s++)
            {
                  if (pids.indexOf(sol_ele[s].pid) == -1) pids.push(sol_ele[s].pid);
            }
            
            var narr = new Array();
            ns0.StartSpool();
            while(id = ns0.Spool()) 
            {
                  var node = Node.GetFromID(model, id);
                  for (var iets = 0; iets < pids.length; iets++)
                  {
                        part = Part.GetFromID(model, pids[iets]);
                        node2 = part.ClosestNode(node.x, node.y, (node.z - foundation_depth)  );
                        node2 = Node.GetFromID(model, node2);
                        if ( Math.abs(node2.z - (node.z - foundation_depth)) < 1e-2 ) narr.push(node2.nid);
                  }
            }  
            ns0.Empty();
            for (var iii = 0; iii < narr.length; iii++) 
            {
                  if (!ns0.Contains(narr[iii]) ) ns0.Add(narr[iii]);
            }
      }

      if(presc_disp_x)
      {
            var nodes = Node.GetFlagged(model, ss1Flag);
            ns0.Empty();
            for (var ni = 0; ni < nodes.length; ni++ )
            {
                  if (nodes[ni].Flagged(ss2Flag)) 
                  {
                        if (!ns0.Contains(nodes[ni].nid)) ns0.Add(nodes[ni].nid);
                  }
            }
      }
      

      if (presc_disp_x || presc_disp_z)
      {
            var delflag = AllocateFlag();
            model.ClearFlag(delflag);
            ss1.SetFlag(delflag);
            model.PropagateFlag(delflag);
            model.DeleteFlagged(delflag);
            model.ClearFlag(delflag);
            ReturnFlag(delflag);
      } 
      
      if (presc_disp_z)
      {
            var delflag = AllocateFlag();
            model.ClearFlag(delflag);
            ss2.SetFlag(delflag);
            model.PropagateFlag(delflag);
            model.DeleteFlagged(delflag);
            model.ClearFlag(delflag);
            ReturnFlag(delflag);
      }  

      if (create_column == 2)
      {
            var fl = AllocateFlag();
            ns0.SetFlag(fl);
            model.PropagateFlag(fl)
            var n = Node.GetFlagged(model,fl)[0];
            model.ClearFlag(fl);
            ReturnFlag(fl);
            var c = new History(model, History.NODE, n.nid, "Foundation Node");
      }
}
if(write_mod) model.Write(model_fname);

///////////////////////////////////////////////////////////////////////////////////////
//
//  Finish

model.UpdateGraphics();

Termination("Finished!");

/////////////////////////////////////////////////////////////////////////////
//
function Termination(line)
{
	if (fInputFile) fInputFile.Close();
	if (model) model.UpdateGraphics();
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
// function OKToContinue(ques)
// {
// 	if(ques)
//       {
//             txt = ques;
//       }else{
//             txt = "OK to continue?";
//       }
      
//       var answer = Window.Question("Question",txt, Window.YES | Window.NO | Window.NONMODAL);
// 	if (answer==Window.NO) Exit();
// }

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
		var nrb = new Nrb(model, set_next, pid_next);
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
		var spc = new Spc(model, id, 0, codes.x,codes.y,codes.z,0,0,0,Spc.SET);
	}

	if (ilay<0) return;

	if (param_column[12] > -1)                //   create BPM in X if required - note, same curve for whole layer
	{
		var lcid = layer_data[ilay][12];
		if (lcid>0) var bpm = new PrescribedMotion(model, id, 1, bpm_vad, lcid, bpm_type);              
	}
	if (param_column[13] > -1)                //   BPM in Y
	{
		var lcid = layer_data[ilay][13];
		if (lcid>0) var bpm = new PrescribedMotion(model, id, 2, bpm_vad, lcid, bpm_type);              
	}
	if (param_column[14] > -1)                //   BPM in Z
	{
		var lcid = layer_data[ilay][14];
		if (lcid>0) var bpm = new PrescribedMotion(model, id, 3, bpm_vad, lcid, bpm_type);              
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
			var ns1 = Node.GetFromID(model, n1sid);
			var ns2 = Node.GetFromID(model, n2sid);
			var ns3 = Node.GetFromID(model, n3sid);
			if (shl.nodes>3) var ns4 = Node.GetFromID(model, n4sid);
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

function rd(a,b)	//round function
{
c = Math.pow(10,b);
d = (Math.round(a * c))/c;
return d;
//var str = new String(d);
//return eval(str.substr(0,str.indexOf(".")+b+1));   
}

function rdup(a,b)
{
      c = rd(a,b);
      if(c < a) c = c + Math.pow(10,-b);
      return c;
}
////////////////////////////////////////////////////////////////////////////////
// Global functions


function round(x, d)
{
	var n = Math.pow(10, d);
	return Math.round(x * n) / n;
}

function distance(n1, n2)
{
	return Math.sqrt( Math.pow(n1.x - n2.x, 2) + 
			Math.pow(n1.y - n2.y, 2) +
			Math.pow(n1.z - n2.z, 2) );
}

function DrawEdges()
{
	Graphics.Start();

	Graphics.LineWidth(3);
	Graphics.LineColour(Colour.RED);
	Graphics.LineStyle(Graphics.SOLID_LINE);
    
	if (outer_nodes)
	{
		for (var i = 0; i < outer_nodes.length; i++)
		{
			var j = i + 1;
			if (j == outer_nodes.length)
				j = 0;
			
			Graphics.Line(outer_nodes[i].x, outer_nodes[i].y, outer_nodes[i].z, 
					outer_nodes[j].x, outer_nodes[j].y, outer_nodes[j].z);
		}
	}

	Graphics.LineStyle(Graphics.DOT_LINE);

	if (inner_nodess)
	{
		inner_nodess.forEach(function(inner_nodes)
		{
			for (var i = 0; i < inner_nodes.length; i++)
			{
				var j = i + 1;
				if (j == inner_nodes.length)
					j = 0;
				
				Graphics.Line(inner_nodes[i].x, inner_nodes[i].y, inner_nodes[i].z, 
						inner_nodes[j].x, inner_nodes[j].y, inner_nodes[j].z);
			}
		});
	}

	Graphics.Finish();
}
function sortMultiDimensional(a,b)
{
Sort_Val = 0;
return ((a[Sort_Val] < b[Sort_Val]) ? -1 : ((a[Sort_Val] > b[Sort_Val]) ? 1 : 0));
}

function compareNumbersAsc(a, b) 
{
  return a - b;
}
function compareNumbersDes(a, b) 
{
  return b - a;
}

function CrossProd(a,b,c)
{
   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];
}

function DotProd(a,b)
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

function GetLocCoord(aa,bb,cc,dd,ee)
//(coordarray,localXVect,localYVect,localZVect,local coords)
{
	ee[0] = DotProd(aa, bb);
	ee[1] = DotProd(aa, cc);
	ee[2] = DotProd(aa, dd);
}

function GetGlobCoord(aa,bb,cc,dd,ee)
//(loc coord array,localXVect,localYVect,localZVect,glob coords)
{
	ee[0] = aa[0]*bb[0] + aa[1]*cc[0] + aa[2]*dd[0];
	ee[1] = aa[0]*bb[1] + aa[1]*cc[1] + aa[2]*dd[1];
	ee[2] = aa[0]*bb[2] + aa[1]*cc[2] + aa[2]*dd[2];
}

function UnitVec(a)
{
   length = Math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
   if (length > 0.0)
   {
      a[0] = a[0]/length;
      a[1] = a[1]/length;
      a[2] = a[2]/length;
   }
}

function ErrorTermination(message)
{
	Window.Error("Error",message,Window.OK);
	ReturnFlag(my_flag);
	Exit();
}
function OKToContinue(line)
{
	var answer = Window.Information("Info",line+"\nOK to continue?", Window.YES|Window.NO|Window.NONMODAL);
	if (answer==Window.NO) Exit();
}
function radians(aa)
{
return aa*(3.141592654/180);
}

function degrees(aa)
{
return aa*(180/3.141592654);
}

function Sign(aa)
{
if (aa>=0){
      bb=1;
}else{
      bb=-1;
}
return bb;
}

function vect(n1, n2)
{
	var vec = new Array();
      vec[0] = n2.x - n1.x;
      vec[1] = n2.y - n1.y;
      vec[2] = n2.z - n1.z;
      vec[3] = Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
      UnitVec(vec);                
      return vec;
                  
}

function GetColour(soil)
{
      var colour1 = soil.toUpperCase();
      if(colour1 == "PEAT")                                 var colour2 = Colour.RGB(130,80,0);
      else if(colour1 == "SILT")                            var colour2 = Colour.RGB(225,230,150);
      else if(colour1 == "CLAY")                            var colour2 = Colour.RGB(190,145,0);
      else if(colour1 == "SILTY CLAY")                      var colour2 = Colour.RGB(190,155,50);
      else if(colour1 == "SAND")                            var colour2 = Colour.RGB(255,220,100);
      else if(colour1 == "SILTY SAND")                      var colour2 = Colour.RGB(255,230,150);
      else if(colour1 == "LYSMER")                          var colour2 = Colour.RGB(130,130,130);
      else                                                  var colour2 = Colour.RGB(145,135,125);
      return colour2
}

function write_KW2(f, first, KW, PRINTS)
{
// Writes the *CONSTRAINED_LINEAR_LOCAL cards to file <f>.
// Assumes <f> has been opened for writing.
//
// If <first> is true then we need to write the first card

// Keyword and first card
  if(first)
  {
    f.Writeln(KW);
  }
    mystr = "";
    for (ip = 0; ip<PRINTS.length; ip++)
    {
        
        if(PRINTS[ip][1] < 0)
        {
            str1 = FixCol(PRINTS[ip][0], -PRINTS[ip][1])
        }else{
            num = PRINTS[ip][0];
            if(PRINTS[ip][0] <0)
            {
                str1 = FixCol(num.toFixed(PRINTS[ip][1]-3),PRINTS[ip][1]); 
            }else{
                str1 = FixCol(num.toFixed(PRINTS[ip][1]-2),PRINTS[ip][1]); 
            }
        }
        mystr = mystr + str1;   //FixCol(PRINTS[ip][0], PRINTS[ip][1])
    }
  
// Second card
  f.Writeln(mystr);
}


function FixCol(Num1, Space)
{
      StrNum = Num1.toString();
      NSpace = Space - StrNum.length;
      //StrOut = StrNum;
      if (NSpace>0)
      {
            for (inum = 0; inum <NSpace; inum++){StrNum = " " + StrNum};
      }
      return StrNum
}

function SolCoords(isol,mx)
{
      var nnod = isol.nodes;

      var SC = new Array(0,0,0);
      
      for(ix=1; ix<=nnod; ix++)
      {
            var tag = "n"+(ix);
            var nt = Node.GetFromID(mx, isol[tag]);
            
            SC[0] = SC[0] + nt.x;
            SC[1] = SC[1] + nt.y;
            SC[2] = SC[2] + nt.z;
      }
      for(ix=0; ix<3; ix++) SC[ix] = SC[ix] / nnod;
      
      return SC
}

function TopCoord(isol,mx)
{
      var nnod = isol.nodes;

      var ztop = 0.0;
      
      for(ix=1; ix<=nnod; ix++)
      {
            var tag = "n"+(ix);
            var nt = Node.GetFromID(mx, isol[tag]);
            
            if (ix == 1) ztop = nt.z
            else if (nt.z > ztop) ztop = nt.z
      }
      
      return ztop
}

function check_par(la,pa,vn,va,vb)
{
      diff = Math.abs((va-vb)/va);
      //if(diff > 1e-6)    OKToContinue("Layer " + la + " has PID " + pa + " which has already been created. " + 
      //"\n Variable " + vn + " does not match. " + 
      //"\n Previous value = " + va + "." + 
      //"\n Value for this layer = " + vb + ". Original values will be kept. Please check!!")
      if(diff > 1e-6)    Message("WARNING!!! \n Layer " + la + " has PID " + pa + " which has already been created. " + 
      "\n Variable " + vn + " does not match. " + 
      "\n Previous value = " + va + "." + 
      "\n Value for this layer = " + vb + ". Original values will be kept. Please check!!")      
      
}

function isStringMatch(str, str_to_match)
{
    return (str.indexOf(str_to_match) > -1);
}
