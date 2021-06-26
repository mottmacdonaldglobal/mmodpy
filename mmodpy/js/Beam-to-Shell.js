// Javascript to create shells from beams
// Date: 12 May 2016
// Version: 1.0
// name: Beam to Shell
// description: create shells from beam

var m = Model.GetFromID(1)
var bflag = AllocateFlag();
var pflag = AllocateFlag();

var qplane = Window.Question("Plane Orientation", "YES = XY Plane : NO = XZ Plane");
if(qplane == Window.YES) var plane = "xy";
if(qplane == Window.NO)  var plane = "xz";

Beam.Select(bflag, 'Select beams', m);

var b = Beam.GetAll(m);
len = b.length;

var part_list = new Object;

for(i=0; i<len; i++)
{
	if (b[i].Flagged(bflag) )
	{
		// check to see if part if already created
		var check = 0;
		for(j in part_list)
		{
			if(j == b[i].pid)
			{
				check = part_list[j];
			}
		}

		if(check != 0)
		{
			var pid = check;
		}
		else
		{
		    	var bpart = Part.GetFromID(m, b[i].pid);
    			var bsec = Section.GetFromID(m, bpart.secid);

			// Integrated Beam
			if( bsec.elform == 0 || bsec.elform == 1)
			{
				if( bsec.cst == 0) // rectanglur section
				{
					if(plane == "xy")
					{
    						var t1 = bsec.ts1;
    						var t2 = bsec.tt1;
					}
					else if(plane == "xz")
					{
    						var t2 = bsec.ts1;
    						var t1 = bsec.tt1;
					}
				}
				else // non-rectangular section
				{
					var dim = dim_gui(b[i]);
					if(plane == "xy")
					{
						var t1 = dim[1];
						var t2 = dim[2];
					}
					else if(plane == "xz")
					{
						var t2 = dim[1];
						var t1 = dim[2];
					}
				}
			}
			// Resultant Beam
			else if( bsec.elform == 2)
			{
				if( bsec.stype == "SECTION_11")// rectangular section
				{
    					if(plane == "xy")
					{
    						var t1 = bsec.d1;
    						var t2 = bsec.d2;
					}
					else if(plane == "xz")
					{
    						var t2 = bsec.d1;
    						var t1 = bsec.d2;
					}
				}
				else // non-rectangular section
				{
					var dim = dim_gui(b[i]);
					if(plane == "xy")
					{
						var t1 = dim[1];
						var t2 = dim[2];
					}
					else if(plane == "xz")
					{
						var t2 = dim[1];
						var t1 = dim[2];
					}
				}
			}
			// Unknown section type
			else
			{
				var dim = dim_gui(b[i]);
				if(plane == "xy")
				{
					var t1 = dim[1];
					var t2 = dim[2];
				}
				else if(plane == "xz")
				{
					var t2 = dim[1];
					var t1 = dim[2];
				}
			}
       
			// create part for bolt nulls 
			var sid = Section.Last(m).secid + 1;
			var sec = new Section(m, sid, Section.SHELL, "Beam Shell");
			sec.t1 = t2;
			sec.t2 = t2;
			sec.t3 = t2;
			sec.t4 = t2;

			var pid = Part.Last(m).pid + 1;
			var p = new Part(m, pid, sid, bpart.mid, "Beam Shell");

			part_list[ b[i].pid ] = pid;
		}

		// find beam vector
		n1 = Node.GetFromID(m, b[i].n1);
		n2 = Node.GetFromID(m, b[i].n2);
		n3 = Node.GetFromID(m, b[i].n3);

		n1x = n1.x;
		n1y = n1.y;
		n1z = n1.z;

		n2x = n2.x;
		n2y = n2.y;
		n2z = n2.z;

		n3x = n3.x;
		n3y = n3.y;
		n3z = n3.z;

	    	// Find beam x axis
		var bx_vec = new Object; 
		bx_vec[1] = (n2x)-(n1x);
		bx_vec[2] = (n2y)-(n1y);
		bx_vec[3] = (n2z)-(n1z);

	    	// Find beam z axis
		var n13 = new Object;
		n13[1] = (n3x)-(n1x);
		n13[2] = (n3y)-(n1y);
		n13[3] = (n3z)-(n1z);
		var bz_vec = x_product(bx_vec, n13);

	    	// Find beam y axis
		var by_vec = x_product(bx_vec, bz_vec);

		// create node pattern at N1
		var nid = Node.Last(m).nid + 1;
		var r1 = new Object;
		r1[1] = new Node(m, nid ,n1x, n1y, n1z+t1/2);
		r1[2] = new Node(m, nid+1 ,n1x, n1y, n1z-t1/2);

		// create node pattern at N2
		var nid = Node.Last(m).nid + 1;
		var r2 = new Object;
		r2[1] = new Node(m, nid ,n2x, n2y, n2z+t1/2);
		r2[2] = new Node(m, nid+1 ,n2x, n2y, n2z-t1/2);

		vbase = new Object;
		vbase[1] = 0;
		vbase[2] = 0;
		vbase[3] = 1;

		for(j=1; j<=2; j++) // rotate pattern nodes
		{
			if(plane == "xy")
			{
				vec = by_vec;
			}
			else if(plane == "xz")
			{
				vec = bz_vec;
			}
		    	rotate(vbase, vec, n1, r1[j]); 
		    	rotate(vbase, vec, n2, r2[j]); 
		}

	    	// create shell elements 1
		if (Shell.Last(m) == null) var eid = 1;
		else var eid = Shell.Last(m).eid + 1;

		var s = new Shell(m, eid, pid, r1[1].nid, r2[1].nid, r2[2].nid, r1[2].nid);
	}
}
/////////////////// dim gui
function dim_gui(b)
{
	var dim = new Object;
	
	wa = new Window("Section Dimensions Entry", 0.2, 0.3, 0.5, 0.8);

	cwl = new Widget(wa, Widget.LABEL, 0, 55, 10, 20, "Width");
	cw = new Widget(wa, Widget.TEXTBOX, 5, 50, 25, 35, "1.0");
	cw.background = Widget.DARKBLUE; 
	cw.foreground = Widget.WHITE; 

	ctl = new Widget(wa, Widget.LABEL, 0, 55, 45, 55, "Thickness");
	ct = new Widget(wa, Widget.TEXTBOX, 5, 50, 60, 70, "1.0");
	ct.background = Widget.DARKBLUE; 
	ct.foreground = Widget.WHITE; 

	var ske = new Widget(wa, Widget.BUTTON, 0, 55, 80, 90, "Sketch Beam");
	ske.background = Widget.DARKBLUE; 
	ske.foreground = Widget.WHITE; 
	ske.onClick = sk_clicked;

	var cre = new Widget(wa, Widget.BUTTON, 0, 25, 95, 105, "Ok");
	cre.background = Widget.DARKGREEN; 
	cre.foreground = Widget.WHITE; 
	cre.onClick = cr_clicked;
		
	var can = new Widget(wa, Widget.BUTTON, 30, 55, 95, 105, "Cancel");
	can.background = Widget.DARKRED; 
	can.foreground = Widget.WHITE; 
	can.onClick = ca_clicked;
	
	wa.Show();

	/////////////////////////////////////////
	function ca_clicked()
	{
		Exit();
	}
	/////////////////////////////////////////
	function sk_clicked()
	{
		b.Sketch();
	}
	/////////////////////////////////////////
	function cr_clicked()
	{
		wa.Hide();
		dim[1] = Number(cw.text);	
		dim[2] = Number(ct.text);
	}
	/////////////////////////////////////////
	return dim;
}
/////////////////// rotate 
function rotate(v1, v2, n, rn)
{

vn1 = vec_norm(v1); // normalize vector 1
vn2 = vec_norm(v2); // normalize vector 2

if( (vn1[1] == vn2[1] && vn1[2] == vn2[2] && vn1[3] == vn2[3]) || (vn1[1] == -vn2[1] && vn1[2] == -vn2[2] && vn1[3] == -vn2[3]) ) // vectors are already aligned - don't rotate
{
	Message("Vectors are aligned no rotation needed");
}
else // vectors not aligned - rotate
{
	// offset node to 0,0,0
	rn.x = rn.x - n.x;
	rn.y = rn.y - n.y;
	rn.z = rn.z - n.z;
	
	vc = x_product(vn1,vn2); // cross product of normalized vectors

	vcn = vec_norm(vc); // normalize cross product vector

	v_dot = dot_product(v1,v2); // dot product of basic vectors
	v_angle = Math.acos(v_dot);

	rot_mat = new Object;
	rot_mat[1] = new Object;
	rot_mat[2] = new Object;
	rot_mat[3] = new Object;

	cos = Math.cos(v_angle)
	sin = Math.sin(v_angle)
	P1 = Math.pow(vcn[1],2)
	P2 = Math.pow(vcn[2],2)
	P3 = Math.pow(vcn[3],2)

	rot_mat[1][1] = cos + P1 * (1 - cos) ;
	rot_mat[2][1] = vcn[2]*vcn[1]*(1-cos) + vcn[3]*sin ; 
	rot_mat[3][1] = vcn[3]*vcn[1]*(1-cos) - vcn[2]*sin ; 

	rot_mat[1][2] = vcn[1]*vcn[2]*(1-cos) - vcn[3]*sin ; 
	rot_mat[2][2] = cos + P2 * (1 - cos) ;
	rot_mat[3][2] = vcn[3]*vcn[2]*(1-cos) + vcn[1]*sin ; 

	rot_mat[1][3] = vcn[1]*vcn[3]*(1-cos) + vcn[2]*sin ; 
	rot_mat[2][3] = vcn[2]*vcn[3]*(1-cos) - vcn[1]*sin ; 
	rot_mat[3][3] = cos + P3 * (1 - cos) ;

	//Message(rot_mat[1][1] +"  "+ rot_mat[1][2] +"  "+ rot_mat[1][3] );
	//Message(rot_mat[2][1] +"  "+ rot_mat[2][2] +"  "+ rot_mat[2][3] );
	//Message(rot_mat[3][1] +"  "+ rot_mat[3][2] +"  "+ rot_mat[3][3] );

	x = rn.x;
	y = rn.y;
	z = rn.z;

	rn.x = (x * rot_mat[1][1]) + (y * rot_mat[1][2]) + (z * rot_mat[1][3]);
	rn.y = (x * rot_mat[2][1]) + (y * rot_mat[2][2]) + (z * rot_mat[2][3]);
	rn.z = (x * rot_mat[3][1]) + (y * rot_mat[3][2]) + (z * rot_mat[3][3]);

	// reposition

	rn.x = rn.x + n.x;
	rn.y = rn.y + n.y;
	rn.z = rn.z + n.z;

}
}
///////////////////
/////////////////// cross product 3 axis
function x_product(vcross1, vcross2)
{
	var vc = new Object;
	vc[1] = vcross1[2]*vcross2[3] - vcross1[3]*vcross2[2];
	vc[2] = vcross1[3]*vcross2[1] - vcross1[1]*vcross2[3];
	vc[3] = vcross1[1]*vcross2[2] - vcross1[2]*vcross2[1];
	vc.mag = Math.sqrt(Math.pow(vc[1],2) + Math.pow(vc[2],2) + Math.pow(vc[3],2));
	return vc;
}
///////////////////
/////////////////// dot product 3 axis 
function dot_product(vdot1, vdot2)
{
	var mag1 = Math.sqrt(Math.pow(vdot1[1],2) + Math.pow(vdot1[2],2) + Math.pow(vdot1[3],2));
	var vdot1_1 = vdot1[1]/mag1;
	var vdot1_2 = vdot1[2]/mag1;
	var vdot1_3 = vdot1[3]/mag1;
	
	var mag2 = Math.sqrt(Math.pow(vdot2[1],2) + Math.pow(vdot2[2],2) + Math.pow(vdot2[3],2));
	var vdot2_1 = vdot2[1]/mag2;
	var vdot2_2 = vdot2[2]/mag2;
	var vdot2_3 = vdot2[3]/mag2;

	var vector_dot = vdot1_1*vdot2_1 + vdot1_2*vdot2_2 + vdot1_3*vdot2_3; // dot product
	return vector_dot;
}
///////////////////
/////////////////// vector normailize (3 axis) 
function vec_norm(v)
{
	var mag = Math.sqrt(Math.pow(v[1],2) + Math.pow(v[2],2) + Math.pow(v[3],2));
	var vnorm = new Object;
	vnorm[1] = v[1]/mag;
	vnorm[2] = v[2]/mag;
	vnorm[3] = v[3]/mag;
	return vnorm ;
}
///////////////////
