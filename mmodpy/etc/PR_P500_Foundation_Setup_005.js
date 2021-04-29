// Javascript to create foundation side springs
// Date: 28 Aug 2017
// Version: 1.0
// name: H+V Springs
// description: H+V Springs

// find outer corners of the foundation
var outer_nodes = new Array();
var inner_nodess = new Array();

var pick = Window.Message("Corner Node", "Pick Foundation Corner Node");
for (var outer_pick; (outer_pick = Node.Pick(
          "Pick corner node outer edge")) != null; )
{	
    outer_nodes = outer_pick.GetFreeEdgeNodes();
    break; 
}

if (outer_nodes == null)
{
    Message("Aborted. Did you select solids instead of shells?");
    Exit();
}
if (outer_nodes.length == 0)
{
    Message("Aborted");
    Exit();
}

var m = Model.GetFromID(outer_nodes[0].model);

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
for(i=0; i<corners.length; i++)
{
	Message("Corner Nodes: "+ nodes[corners[i]].nid);
}

var foundation_z = nodes[corners[0]].z;


/////////////////////////////////
///////// Extrude Shells ////////
var answerS = Window.Question("Create foundation solids?", "Do you want to create foundation solids?");

if (answerS == Window.YES)
{
      // Get Part ID for grade beams
      var p_gr_id = Window.GetInteger("PID Grade beams", "Enter PID for grade beams", Part.NextFreeLabel(m));
      if( p_gr_id == null) Exit();

      var beam_length = Window.GetNumber("Spring Length", "Enter spring length",0.01);
      if( beam_length == null) Exit();

      var grade_height = Window.GetNumber("Grade Beam Height", "Enter Grade Beam Height",0.5);
      if (grade_height == null) Exit();

      var s_num = Window.GetNumber("Number of solids", "Number of solids to be created through the height of the grade beam",2);
      if (s_num == null) Exit();
}
var height = 0 - nodes[corners[0]].z;

if (answerS == Window.YES)
{
      // calulate layers
      var s_num = Math.round(Math.max(s_num,1));
      var int = grade_height / s_num;
      Message("Foundation Layers: "+s_num);

      var sflag = AllocateFlag();
      var merge_flag = AllocateFlag();
      var select  = Window.Message("Shells Extrude", "Select Shells to Extrude");
      Shell.Select(sflag, 'Select shells to extrude', m, false);

      var s = Shell.GetAll(m);
      for(i=0; i<s.length; i++)
      {
            if( s[i].Flagged(sflag) )
            {
                  var flip = 0;
                  var norm = s[i].NormalVector();
                  if( norm[2] < 0)
                  {
                        s[i].ReverseNormal();
                        flip = 1;	
                  }

                  for(j=1; j<=s_num; j++)
                  {
                        var z1 = (int*(j-1)) - height;
                        var z2 = (int*j) - height;
                        if(s[i].nodes > 3)
                        {
                              var n = new Object;
                              for(k=1; k<=4; k++)
                              {
                                    var nlabel = Node.Last(m).nid+1;
                                    var node = Node.GetFromID(m, s[i]["n"+k]);

                                    // check if node is outer
                                    var check = 0;
                                    for(l=0; l<outer_nodes.length; l++)
                                    {
                                          if(node.nid == outer_nodes[l].nid) check = 1;
                                    }

                                    if(check == 1)
                                    {
                                          var coords = adjust_postion(node, beam_length);
                                    }
                                    else
                                    {
                                          var coords = new Object;
                                          coords.x = node.x;
                                          coords.y = node.y;
                                    }
                                    n[k] = new Node(m, nlabel  ,coords.x, coords.y, z1); 
                                    n[k+4] = new Node(m, nlabel+1, coords.x, coords.y, z2); 
                              }
                              var slabel = Solid.NextFreeLabel(m);      //Solid.Last(m).eid+1;
                              var solid = new Solid(m, slabel, p_gr_id, n[1].nid, n[2].nid, n[3].nid, n[4].nid, n[5].nid, n[6].nid, n[7].nid, n[8].nid);
                        }
                        else
                        {
                              var n = new Object;
                              for(k=1; k<=3; k++)
                              {
                                    var nlabel = Node.Last(m).nid+1;
                                    var node = Node.GetFromID(m, s[i]["n"+k]);

                                    // check if node is outer
                                    var check = 0;
                                    for(l=0; l<=outer_nodes.length; l++)
                                    {
                                          if(node.nid == outer_nodes[l]) check = 1;
                                    }

                                    if(check == 1)
                                    {
                                          var coords = adjust_postion(node, beam_length);
                                    }
                                    else
                                    {
                                          var coords = new Object;
                                          coords.x = node.x;
                                          coords.y = node.y;
                                    }
                                    n[k] = new Node(m, nlabel  , coords.x, coords.y, z1); 
                                    n[k+3] = new Node(m, nlabel+1, coords.x, coords.y, z2); 
                              }
                              var slabel = Solid.NextFreeLabel(m);      //Solid.Last(m).eid+1;
                              var solid = new Solid(m, slabel, p_gr_id, n[1].nid, n[2].nid, n[5].nid, n[4].nid, n[3].nid, n[6].nid);
                        }
                  }
                  if(flip == 1)
                  {
                        s[i].ReverseNormal();
                  }
            }	
      }
      var sec_label = Section.FirstFreeLabel(m);
      var sec = new Section(m, sec_label, Section.SOLID, 'Grade Beam');
      sec.elform = 2;
      var mat_label = Material.FirstFreeLabel(m);
      var mat = new Material(m, mat_label, "NULL");

      var p = Part.GetFromID(m, p_gr_id);
      if (p == null || !p.exists) {var p = new Part(m, p_gr_id, sec_label, mat_label, 'Grade Beam');}
      p.SetFlag(merge_flag);
      m.PropagateFlag(merge_flag);
      var merged = m.MergeNodes(merge_flag, 0.005);
      Message("m"+merged);

      var p_gr = Part.GetFromID(m, p.pid);

      ReturnFlag(merge_flag);
      m.UpdateGraphics();
}

///////////////////////////////////////////////////////
/// Get Horizontal Spring and Grade Beam Parameters ///

var answerV = Window.Question("Create horizontal Springs?", "Do you want to create horizontal springs?");

if (answerV == Window.YES)
{
      // Get Part IDs for horizontal springs
      var pxid = Window.GetInteger("PID Springs PX", "Enter PID for Soil springs horizontal - PX", Part.NextFreeLabel(m));
      if( pxid == null) Exit();

      var nxid = Window.GetInteger("PID Springs NX", "Enter PID for Soil springs horizontal - NX", pxid+1);
      if( nxid == null) Exit();

      var pyid = Window.GetInteger("PID Springs PX", "Enter PID for Soil springs horizontal - PY", nxid+1);
      if( pyid == null) Exit();

      var nyid = Window.GetInteger("PID Springs NY", "Enter PID for Soil springs horizontal - NY", pyid+1);
      if( nyid == null) Exit();

      var rows = Window.GetInteger("Rows", "Enter number of spring rows",1);
      if(rows == null) Exit();
      
      var base_spacing = Window.GetNumber("Spacing", "Enter spring spacing",0.2);
      if( base_spacing == null) Exit();
      
      var answer = Window.Question("Spring Location", "Spring Location Based on Grade Beam Height");
      if(answer == Window.YES) var s_height = grade_height;
      else var s_height = height;

      if (answerS == Window.NO)
      {
            var beam_length = Window.GetNumber("Spring Length", "Enter spring length",0.01);
            if( beam_length == null) Exit();
      }
}

// set up contact spring sets
var set = Set.First(m, Set.NODE);
while(set)
{
      if (set.title == "Soil Contact Set")
      {
            var soil_set = set;
            break;
      }
      set = set.Next();
} 

var set = Set.First(m, Set.NODE);
while(set)
{
      if (set.title == "Foundation Contact Set")
      {
            var foundation_set = set;
            break;
      }
      set = set.Next();
}

var set = Set.First(m, Set.PART);
while(set)
{
      if (set.title == "Foundation Parts Contact Set")
      {
            var foundation_part_set = set;
            break;
      }
      set = set.Next();
}

// Create tied contacts
if(!soil_set || !foundation_set || !foundation_part_set)
{
      // Create new node sets for tied contact
      var set_label = Set.LastFreeLabel(m, Set.NODE);
      var soil_set = new Set(m, set_label, Set.NODE);
      soil_set.title = "Soil Contact Set";

      var set_label = Set.LastFreeLabel(m, Set.NODE);
      var foundation_set = new Set(m, set_label, Set.NODE);
      foundation_set.title = "Foundation Contact Set";

      // select soil parts for the contact
      var pick = Window.Message("Spring Contact", "Select Soil Parts for Spring Contact");
      var pflag = AllocateFlag();
      Part.Select(pflag, 'Select parts', m, false);
      var p = Part.GetFlagged(m, pflag);
      var set_label = Set.LastFreeLabel(m, Set.PART);
      var soil_part_set = new Set(m, set_label, Set.PART);
      soil_part_set.title = "Soil Parts Contact Set";
      for(i=0; i<p.length; i++)
      {
            soil_part_set.Add(p[i].pid);
      }

      var clabel = Contact.LastFreeLabel(m);
      var c1 = new Contact(m, "TIED_NODES_TO_SURFACE", clabel, "Soil to Spring Contact");
      c1.sstyp = 4;
      c1.ssid = soil_set.sid;
      c1.mstyp = 2;
      c1.msid = soil_part_set.sid;
      c1.soft = 1;

      // select foundation parts for the contact
      var pick = Window.Message("Spring Contact", "Select Foundation Parts for Spring Contact");
      var pflag = AllocateFlag();
      Part.Select(pflag, 'Select parts', m, false);
      var p = Part.GetFlagged(m, pflag);
      var set_label = Set.LastFreeLabel(m, Set.PART);
      var foundation_part_set = new Set(m, set_label, Set.PART);
      foundation_part_set.title = "Foundation Parts Contact Set";
      for(i=0; i<p.length; i++)
      {
            foundation_part_set.Add(p[i].pid);
      }

      var clabel = Contact.LastFreeLabel(m);
      var c2 = new Contact(m, "TIED_NODES_TO_SURFACE", clabel, "Foundation to Spring Contact");
      c2.sstyp = 4;
      c2.ssid = foundation_set.sid;
      c2.mstyp = 2;
      c2.msid = foundation_part_set.sid;
      c2.soft = 1;
      c2.offset_flag = 1; // penalty-based to prevent possible clash with NRBs on foundations
}

// Create Horizontal Springs
for(i=0; i<=corners.length-2; i++)
{
      // get corner nodes
      n1 = nodes[corners[i]];
      n2 = nodes[corners[i+1]];
      //Message("N1: "+n1.nid +"  N2: "+n2.nid);

      //calculate length
      var xlen = n2.x - n1.x; 
      var ylen = n2.y - n1.y;
      var len = Math.sqrt(Math.pow(xlen, 2) + Math.pow(ylen, 2));
      //Message("X len: "+xlen+"  Y Len: "+ylen+"  Len: "+len);

      //calculate spacing
      var number = Math.round( len / base_spacing );
      var x_spacing = xlen / number;
      var y_spacing = ylen / number;
      //Message("Num: "+number+"  Space: "+x_spacing+"  "+y_spacing);

      // calculate beam direction and typical shell element size
      var e_size = 0;
      var shell_array = n1.GetAttachedShells(false);
      var c_coords = [0,0];
      // find centroid of shell connected to the corner node
      for(j=0; j<shell_array.length; j++)
      {
            var coords = shell_array[j].IsoparametricToCoords(0, 0);
            c_coords[0] = c_coords[0] + coords[0];
            c_coords[1] = c_coords[1] + coords[1];

            e_size = Math.max(shell_array[j].Length(), e_size);
      }

      //setup third node for direction vector calc
      var n3 = new Object;
      n3.x = c_coords[0]/shell_array.length;
      n3.y = c_coords[1]/shell_array.length;
      n3.z = n1.z;
      var rs_vec = getvector(n1,n2,n3); // direction vector calc (n1 - first corner node, n2 - second corner node, n3 - summed shell centroids)

      if (answerV == Window.YES)
      {
            for(j=1; j<=number; j++) // cycle spring rows and create springs
            {
                  var x_offset = (j*x_spacing)-(x_spacing/2);
                  var y_offset = (j*y_spacing)-(y_spacing/2);
                  var x_coord = n1.x + x_offset;
                  var y_coord = n1.y + y_offset;

                  for(k=1; k<=(2*rows); k=k+2)
                  {
                        var int = s_height/(2*rows);
                        var z_coord = foundation_z + ( k*int );
                        
                        var nlabel = Node.LastFreeLabel(m);
                        var bn1 = new Node(m, nlabel, x_coord, y_coord, z_coord);
                        var nlabel = Node.LastFreeLabel(m);
                        var bn2 = new Node(m, nlabel, x_coord+(beam_length*rs_vec[0]), y_coord+(beam_length*rs_vec[1]), (z_coord+(beam_length*rs_vec[2])) );

                        //check orientation/vector and set part id
                        var beam_vec = new Object;
                        beam_vec[1] = bn2.x - bn1.x; 
                        beam_vec[2] = bn2.y - bn1.y; 
                        beam_vec[3] = bn2.z - bn1.z;
                        var base_vec = new Object;
                        base_vec[1] = 1; 
                        base_vec[2] = 0; 
                        base_vec[3] = 0; 
                        var angle  = dot_product_angle(beam_vec, base_vec);
                        
                        if(angle.deg > -10 && angle.deg < 10) var bpid = pxid; // positive x
                        else if(angle.deg > 170 && angle.deg < 190) var bpid = nxid; // negative x
                        else if(angle.deg > 80 && angle.deg < 100) // positive or negative y
                        {
                              var base_vec = new Object;
                        base_vec[1] = 0; 
                        base_vec[2] = 1; 
                        base_vec[3] = 0; 
                        var angle2  = dot_product_angle(beam_vec, base_vec);
                              
                              if(angle2.deg > -10 && angle2.deg < 10) var bpid = pyid; // positive y
                              else if(angle2.deg > 170 && angle2.deg < 190) var bpid = nyid; // negative y
                        }

                        var blabel = Beam.LastFreeLabel(m);
                        var b = new Beam(m, blabel, bpid, bn1.nid, bn2.nid);
                  
                        // add node to contact sets    
                        soil_set.Add(bn1.nid);
                        foundation_set.Add(bn2.nid);
                  }
            }
      }
}
//
//-----------End of Ian Bruce's script-------------------------------------------------------
//


///////////////////////////////
/// Create Vertical Springs ///

var answer = Window.Question("Create vertical Springs?", "Do you want to create vertical springs?");
if (answer == Window.YES)
{
      Message("Creating vertical springs");
      
      Model.BlankAll();
      foundation_part_set.StartSpool();
      while (pid_unblank = foundation_part_set.Spool())
      {
            part_unblank = Part.GetFromID(m,pid_unblank);
            part_unblank.Unblank();
      }
      View.Redraw();

      var vspr_length = Window.GetNumber("Spring length", "Enter spring length",0.01);
      if( vspr_length == null) Exit();

      //start by making sure all shell normals point downwards
      //assume we only have horizontal shells
      var s = Shell.First(m);
      var nvector = s.NormalVector();
      if(nvector[2]>0) s.ReverseNormal();
      var f1 = AllocateFlag();
      m.ClearFlag(f1);
      Shell.FlagAll(m, f1);
      Shell.MakeConsistentNormalsFlagged(m, f1, s.eid);
      ReturnFlag(f1); 

      // Generate vertical springs
      while(true)
      {
            View.Redraw();
            //var nn1 = nodes[corners[0]];
            //var nn2 = nodes[corners[1]];
            //var nn3 = Node.Pick("Pick node 3");  
            
            if (!pid_inside) var PID = Part.NextFreeLabel(m);
            else var PID = pid_inside;
            
            var pid_inside = Window.GetInteger("PID Inside Springs", "Enter PID for springs on inside of grade beam / footing. Note PID for springs on outside will be incremented by 1",PID);
            if( pid_inside == null) break;

            var pid_edge = pid_inside+1;
            
            if (!vspr_No) {var numspring = 6;}
            else {var numspring = vspr_No;}

            var vspr_No = Window.GetInteger("Number of springs across section", "Enter Number of springs across section", numspring);
            if( vspr_No == null) break;

            if (nele_edge == null) {var numedgespring = 1;}
            else {var numedgespring = nele_edge;}
      
            var nele_edge = Window.GetInteger("Number of Edge Springs", "Enter number of springs to be created at edge of footing", numedgespring);      
            if( nele_edge == null) break;

            var nn1 = Node.Pick("Pick node 1");
            if( nn1 == null) break;
            var nn2 = Node.Pick("Pick node 2");
            if( nn2 == null) break;
            var nn3 = Node.Pick("Pick node 3"); 
            if( nn3 == null) break;      
            
            var vec1 = new Array();
            vec1 = vect(nn1, nn2);
            var vec2 = new Array();
            vec2 = vect(nn1, nn3);
            
            var glob_origin = new Array(nn1.x, nn1.y, nn1.z);
            var glob_origin3 = new Array(nn3.x, nn3.y, nn3.z);
            var loc_origin = new Array();
            var loc_origin3 = new Array();
            
            var vec3 = new Array();
            CrossProd(vec1,vec2,vec3);
            UnitVec(vec3);
            
            CrossProd(vec3,vec1,vec2);     
            
            GetLocCoord(glob_origin,vec1,vec2,vec3,loc_origin);
            GetLocCoord(glob_origin3,vec1,vec2,vec3,loc_origin3);   
            
            vec2[3] = Math.abs(loc_origin3[1]-loc_origin[1]);
            
            n_pt_2 = vspr_No; //Math.round(vec2[3]/vspr_spacing);
            n_pt_2 = Math.max(n_pt_2,1);
            space_pt_2 = vec2[3]/vspr_No;
            
            n_pt_1 = Math.round(vec1[3]/space_pt_2);
            n_pt_1 = Math.max(n_pt_1,1);      
            space_pt_1 = vec1[3]/n_pt_1;   

            var Loc_co = new Array(0,0,loc_origin[2]);
            
            for(ispx=0; ispx<n_pt_1; ispx++)
            {
                  Loc_co[0] = loc_origin[0] + (space_pt_1/2 + ispx*space_pt_1);
                  for(ispy=0; ispy<n_pt_2; ispy++)
                  {
                        Loc_co[1] = loc_origin[1] + (space_pt_2/2 + ispy*space_pt_2);
                        var Glo_co = new Array(0,0,0);
                        GetGlobCoord(Loc_co,vec1,vec2,vec3,Glo_co);
                        
                        var nNo =  Node.NextFreeLabel(m);
                        Nx1 = new Node(m, nNo, Glo_co[0],Glo_co[1], Glo_co[2]-vspr_length);
                        var nNo2 =  Node.NextFreeLabel(m);
                        Nx2 = new Node(m, nNo2, Glo_co[0],Glo_co[1], Glo_co[2]); 
                        
                        var bNo =  Beam.NextFreeLabel(m);
                        if(ispy < nele_edge || ispy > n_pt_2-nele_edge-1){pid = pid_edge}else{ pid = pid_inside};
                        bx = new Beam(m,bNo,pid,nNo,nNo2);
                        
                        soil_set.Add(nNo);
                        foundation_set.Add(nNo2);
                  }
            }      
      }
}
/// remove shell elements
if (answerS == Window.YES)
{
      var answerD = Window.Question("Delete shells?", "Do you want to delete old foundation shells?");
      if (answerD == Window.YES)
      {
            m.PropagateFlag(sflag);
            m.DeleteFlagged(sflag);
            ReturnFlag(sflag);
      }
}
Message("Done");
Exit();

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
function getsolidcentroid(s)
{
	var n_count = s.nodes  // number of node in the element

	var tx = 0;
	var ty = 0;
	var tz = 0;

	for(j=1; j<=n_count; j++)
	{
		var tag = "n"+(j);
		tx = tx + Node.GetFromID(m, s[tag]).x
		ty = ty + Node.GetFromID(m, s[tag]).y
		tz = tz + Node.GetFromID(m, s[tag]).z
	}

	var cen = new Object;
	cen.x = tx / n_count;
	cen.y = ty / n_count;
	cen.z = tz / n_count;
		
	return cen; 
}

function adjust_postion(node, base_length)
{
	var shell_array = node.GetAttachedShells(false);
	var coords = new Object;
	var xvec = 0;
	var yvec = 0;

	for(var i=0; i<shell_array.length; i++) // sum attached shell centroid x/y vector
	{
		var iso = shell_array[i].IsoparametricToCoords(0, 0);
		var x_tmp = Sign(iso[0] - node.x,1e-6);
		var y_tmp = Sign(iso[1] - node.y,1e-6);

		var xvec = xvec + x_tmp;
		var yvec = yvec + y_tmp;
	}

	// change to a unit vector
	var vecmag = Math.sqrt(Math.pow(xvec,2)+Math.pow(yvec,2));
	var xvec = xvec / vecmag;
	var yvec = yvec / vecmag;    

	if(shell_array.length == 2)
	{
		var vlen = base_length;
	}
	else
	{
            var vlen = Math.sqrt(Math.pow(base_length,2)+Math.pow(base_length,2));
	}

	coords.x = node.x + (vlen)* xvec;
	coords.y = node.y + (vlen)* yvec;
	
	return coords
}
function getvector(n1,n2,n3)
{
    n1x = n1.x;
    n1y = n1.y;
    n1z = n1.z;
    n2x = n2.x;
    n2y = n2.y;
    n2z = n2.z;
    n3x = n3.x;
    n3y = n3.y;
    n3z = n3.z;

    var n12 = new Object; // x axis
    n12[1] = (n2x)-(n1x);
    n12[2] = (n2y)-(n1y);
    n12[3] = (n2z)-(n1z);

    var xcos = new Object; // x cosine
    xcos[1] = n12[1] / ( Math.abs(n12[1]) + Math.abs(n12[2]) + Math.abs(n12[3]) ); // X x cosine
    xcos[2] = n12[2] / ( Math.abs(n12[1]) + Math.abs(n12[2]) + Math.abs(n12[3]) ); // X y cosine
    xcos[3] = n12[3] / ( Math.abs(n12[1]) + Math.abs(n12[2]) + Math.abs(n12[3]) ); // X z cosine

    var n13 = new Object; // xz plane
    n13[1] = (n3x)-(n1x);
    n13[2] = (n3y)-(n1y);
    n13[3] = (n3z)-(n1z);

    var na13 = new Object; // temp y cosine
    na13[1] = n13[1] / ( Math.abs(n13[1]) + Math.abs(n13[2]) + Math.abs(n13[3]) ); // Temp Y x cosine
    na13[2] = n13[2] / ( Math.abs(n13[1]) + Math.abs(n13[2]) + Math.abs(n13[3]) ); // Temp Y y cosine
    na13[3] = n13[3] / ( Math.abs(n13[1]) + Math.abs(n13[2]) + Math.abs(n13[3]) ); // Temp Y z cosine

    var ztcos = new Object; // temp z cosine
    ztcos[1] =  (xcos[2] * na13[3]) - (xcos[3] * na13[2]); 
    ztcos[2] =  (xcos[3] * na13[1]) - (xcos[1] * na13[3]);
    ztcos[3] =  (xcos[1] * na13[2]) - (xcos[2] * na13[1]); 

    var zcos = new Object; // z cosine
    zcos[1] = ztcos[1] / ( Math.abs(ztcos[1]) + Math.abs(ztcos[2]) + Math.abs(ztcos[3]) ); // Z x cosine
    zcos[2] = ztcos[2] / ( Math.abs(ztcos[1]) + Math.abs(ztcos[2]) + Math.abs(ztcos[3]) ); // Z y cosine
    zcos[3] = ztcos[3] / ( Math.abs(ztcos[1]) + Math.abs(ztcos[2]) + Math.abs(ztcos[3]) ); // Z z cosine

    var ytcos = new Object; // temp y cosine
    ytcos[1] =  (zcos[2] * xcos[3]) - (zcos[3] * xcos[2]); 
    ytcos[2] =  (zcos[3] * xcos[1]) - (zcos[1] * xcos[3]);
    ytcos[3] =  (zcos[1] * xcos[2]) - (zcos[2] * xcos[1]);

    var ycos = new Object; // y cosine
    ycos[1] = ytcos[1] / ( Math.abs(ytcos[1]) + Math.abs(ytcos[2]) + Math.abs(ytcos[3]) ); // Z x cosine
    ycos[2] = ytcos[2] / ( Math.abs(ytcos[1]) + Math.abs(ytcos[2]) + Math.abs(ytcos[3]) ); // Z y cosine
    ycos[3] = ytcos[3] / ( Math.abs(ytcos[1]) + Math.abs(ytcos[2]) + Math.abs(ytcos[3]) ); // Z z cosine

    // Full Consine Matrix
    //Message(xcos[1]+"  "+xcos[2]+"   "+xcos[3]);
    //Message(ycos[1]+"  "+ycos[2]+"   "+ycos[3]);
    //Message(zcos[1]+"  "+zcos[2]+"   "+zcos[3]);
    return [ ycos[1], ycos[2], ycos[3] ];
}

function distance(n1, n2)
{
	return Math.sqrt( Math.pow(n1.x - n2.x, 2) + 
			Math.pow(n1.y - n2.y, 2) +
			Math.pow(n1.z - n2.z, 2) );
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

function DotProd(a,b)
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

function CrossProd(a,b,c)
{
   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];
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

function dot_product_angle(vdot1, vdot2)
{
	var mag1 = Math.sqrt(Math.pow(vdot1[1],2) + Math.pow(vdot1[2],2) + Math.pow(vdot1[3],2));
	var mag2 = Math.sqrt(Math.pow(vdot2[1],2) + Math.pow(vdot2[2],2) + Math.pow(vdot2[3],2));

	var vector_dot = vdot1[1]*vdot2[1] + vdot1[2]*vdot2[2] + vdot1[3]*vdot2[3]; // dot product

	var angle = new Object;
	angle.rad = Math.acos( vector_dot / (mag1 * mag2) );
	angle.deg = angle.rad / Math.PI*180;

	return angle;
}

function Sign(aa,toli)
{
if (aa >= toli){
      bb=1;
}else if (aa <= -toli){
      bb=-1;
}else{
      bb=0;
}
return bb;
}
///////////////////

