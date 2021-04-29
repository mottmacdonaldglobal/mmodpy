// JavaScript: triangular_prism_dupnode_fixer.js
// Author: Lauren Biscombe
// Date: 05 April 2015
// Job: 229746 P500
//
//
//
//
// DESCRIPTION
//
// If you have triangular prism-shaped solid elements with duplicate nodes, this will eliminate the duplicate nodes from the 
// solid definition, turning the 8-noded element into a 6-noded element.
// This problem occurs when using the script that extrudes a planar mesh into layers of solid elements; it can't handle triangles.
//
//


Message("Starting \"triangular_prism_dupnode_fixer.js\"...");


// Select a model. Model selected automatically if only one present.
var m = Model.Select("Select a model:"); 
if(m == null) Exit();


// Pop up reminder about memory.
var answer = Window.Warning("Memory", "Make sure you set Memory to at least 500 when running this script! Did you?", Window.YES|Window.NO);
if (answer == Window.NO) Exit();

// Pop up estimate of script runtime.
//var answer2 = Window.Warning("Time estimate", "Divide total number of solid elements by 44,444 for a minutes estimate of script runtime. Proceed?:", Window.YES|Window.NO);
//if (answer2 == Window.NO) Exit();

// Pop up about deleting old elements.
Window.Warning("Post-script fixes", "After running script, please delete duplicate elements with element IDs offset by 1000000 from previous final EID. Do not delete associated nodes, loads, etc.", Window.YES|Window.NO);

// Select set of solid elements with duplicate nodes. This set can be created after checking the model. Right click on the solid element
// error and create SET_SOLID.
var dupFlag = AllocateFlag();
var dups = Set.Select(Set.SOLID, dupFlag, "Select Solid Set containing solids with duplicate nodes:", m, false);
if(dups == null) Exit();

var ssolid = Set.GetAll(m, Set.SOLID);
for(i=0; i<ssolid.length; i++)
{
	if(ssolid[i].Flagged(dupFlag))
	{
		setid = ssolid[i].sid;
		i = ssolid.length;
	}
}
var setsolid = Set.GetFromID(m, setid, Set.SOLID);
var stotal = setsolid.total

// Edit nodes defining each solid in dupSet. To define a 6-noded solid element with 8 nodes, you must have the nodes in the following order:
// n1, n2, n3, n4, n5, n5, n6, n6
var sNextLabel = Solid.NextFreeLabel(m);

var id;
var co1= 0;
var co2= 0;
setsolid.StartSpool();
while (id = setsolid.Spool() )
{
	var s = Solid.GetFromID(m, id);

	var sn1 = s.n1;
	var sn2 = s.n2;
	var sn3 = s.n3;
	var sn5 = s.n5;
	var sn6 = s.n6;
	var sn7 = s.n7;

	var sEID = id;
	s.eid = sNextLabel + 1000000;
	sNextLabel++;

	var sNew = new Solid(m, sEID, s.pid, sn2, sn1, sn5, sn6, sn3, sn7);

	co1 = co1 + 1;
	co2 = co2 + 1;
	if( (co1 /stotal * 100 ) > 5 )
	{
		co1 = 0;
		Message( Math.round(co2 /stotal * 100) +"%Complete");
	}
}
