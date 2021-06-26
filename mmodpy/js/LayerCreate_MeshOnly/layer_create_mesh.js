//-------------------------------------------------------------------------------------
//  drag shell elements down to a specifed level
//-------------------------------------------------------------------------------------
//  user input
//

var num_sublayers = 50; //number of elements on drag
var drag_arr = new Array(0,0,-1);
var drag_Len = 100;
var bias_sol = 0;
 

//-------------------------------------------------------------------------------------
var m = Model.Select("Select model");

var drag_vec = new Array(drag_arr[0]*drag_Len, drag_arr[1]*drag_Len, drag_arr[2]*drag_Len);
var drag_unit = new Array(drag_arr[0]*drag_Len, drag_arr[1]*drag_Len, drag_arr[2]*drag_Len);
UnitVec(drag_unit);

var flag = AllocateFlag();
m.ClearFlag(flag);

Part.Select(flag,"Select shell parts to extrude");

m.PropagateFlag(flag);

var s_all = Shell.GetAll(m);
var s_flag = new Array();
for (i=0; i<s_all.length; i++)
{
      if (s_all[i].Flagged(flag))
      {
            s_flag.push(i);
      }
}

// check normals of shells
var sref = s_flag[0];
var nvector = s_all[sref].NormalVector();
var nz = Node.GetFromID(m, s_all[sref].n1);

var dot_p = DotProd(nvector,drag_unit);
if(dot_p > 0) s_all[sref].ReverseNormal();

//if(nz.z > drag_lev)
//{
//      if(nvector[2]<0) s_all[sref].ReverseNormal();
//}else{
//      if(nvector[2]>0) s_all[sref].ReverseNormal();
//}

Shell.MakeConsistentNormalsFlagged(m, flag, s_all[sref].eid);
     
View.Redraw();

var solid_pid = Part.NextFreeLabel(m);
var Solid_sid = Section.NextFreeLabel(m);
var Solid_mid = Material.NextFreeLabel(m);

var n_all = Node.GetAll(m);
var flagged_nodes = new Array();
var num_flagged = 0;

for (i=0; i<n_all.length; i++)
{
	if (n_all[i].Flagged(flag)) 
	{
		flagged_nodes.push(i);
		num_flagged++;
	}
}

if (!num_flagged) ErrorTerminate("Nothing selected");

// create new nodes
var new_nodes = new Array();
var new_node_ref = new Array();
var first = 1;

bia = bias_sol;

if(!bia) bia = 1.0;
if(bia<0) bia = Math.pow(-bia,1/(num_sublayers-1));
var bfac = 0; 
for(j=0; j<num_sublayers; j++)     bfac = bfac + Math.pow(bia,j) ; 

var rat = 0;

for (i=0; i<=num_sublayers; i++)
{
      
      new_nodes[i] = new Array();
      //new_node_ref[i] = new Array();
      if(i>0) rat = rat + Math.pow(bia,(i-1))/bfac;
      rat1 = 1-rat;       
	
      for (j=0; j<flagged_nodes.length; j++)
      {
            if(first)
            {
                  new_nodes[i][j] = n_all[flagged_nodes[j]].nid;
                  new_node_ref[n_all[flagged_nodes[j]].nid] = j;
            }else{
                  var nid = Node.NextFreeLabel(m);
                  new_nodes[i][j] = nid;
                  
                  x_0 = n_all[flagged_nodes[j]].x;
                  y_0 = n_all[flagged_nodes[j]].y;
                  z_0 = n_all[flagged_nodes[j]].z;
                  
                  x_targ = x_0 + drag_vec[0];
                  y_targ = y_0 + drag_vec[1];
                  z_targ = z_0 + drag_vec[2];
                  

                  //rat = i/num_sublayers;
                  //rat1 = 1-rat;
                  
                  x = rat*x_targ + rat1*x_0;
                  y = rat*y_targ + rat1*y_0;
                  z = rat*z_targ + rat1*z_0;
                  
                  var n2 = new Node(m,nid,x,y,z);
            }
	}
      first = 0;
}



// Extrude shells to create solids
for (i=0; i<num_sublayers; i++)
{     
      for (j=0; j<s_flag.length; j++)
      {
            SID = Solid.NextFreeLabel(m);

            if(s_all[s_flag[j]].nodes == 3)
            {
                  
                  n1 = new_nodes[i][new_node_ref[s_all[s_flag[j]].n1]];
                  n2 = new_nodes[i][new_node_ref[s_all[s_flag[j]].n2]];
                  n5 = new_nodes[i][new_node_ref[s_all[s_flag[j]].n3]];

                  n4 = new_nodes[i+1][new_node_ref[s_all[s_flag[j]].n1]];
                  n3 = new_nodes[i+1][new_node_ref[s_all[s_flag[j]].n2]];
                  n6 = new_nodes[i+1][new_node_ref[s_all[s_flag[j]].n3]];

                  var so = new Solid(m,SID,solid_pid,n1,n2,n3,n4,n5,n6);
            }
            else if(s_all[s_flag[j]].nodes == 4)
            {         
                  
                  n1 = new_nodes[i][new_node_ref[s_all[s_flag[j]].n1]];
                  n4 = new_nodes[i][new_node_ref[s_all[s_flag[j]].n2]];
                  n3 = new_nodes[i][new_node_ref[s_all[s_flag[j]].n3]];
                  n2 = new_nodes[i][new_node_ref[s_all[s_flag[j]].n4]];

                  n5 = new_nodes[i+1][new_node_ref[s_all[s_flag[j]].n1]];
                  n8 = new_nodes[i+1][new_node_ref[s_all[s_flag[j]].n2]];
                  n7 = new_nodes[i+1][new_node_ref[s_all[s_flag[j]].n3]];
                  n6 = new_nodes[i+1][new_node_ref[s_all[s_flag[j]].n4]];

                  var so = new Solid(m,SID,solid_pid,n1,n2,n3,n4,n5,n6,n7,n8);
            }
      }
}

//Message("Number of elements" + so.nid + "   " + n1 + "   " + n2 + "   " + n3);

Message("Done");
Exit();

//------------------------------------------------------------------------------------------------
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
