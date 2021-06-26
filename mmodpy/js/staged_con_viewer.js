// Javascript to view construction stages
// Date: 08 Nov 2018
// Author: Ian Bruce (Arup)
// Version: 4.0
// name: Staged Construction Viewer 
// description: Tool to visualize which parts are in each construction stage
  
var m = Model.GetFromID(1);

var dcs_data = {};
var dscp_data = {};

Message("Processing Construction Stage Data");

// Define_Construction_Stages
var cs = ConstructionStages.GetAll(m);
dcs_co = cs.length;
for(i=1; i<=dcs_co; i++)
{
	dcs_data[i] = new Object;
	dcs_data[i].st = cs[i-1].ats; // stage start time
	dcs_data[i].et = cs[i-1].ate; // stage end time
}

// Define_Staged_Construction_Part
var scp = StagedConstructionPart.GetAll(m);
dscp_co = scp.length;
for(i=1; i<=dscp_co; i++)
{
	dscp_data[i] = new Object;
	
	if(scp[i-1].option == 54) // Part
	{
		dscp_data[i].st = 0; //  not a set id
		dscp_data[i].pt = scp[i-1].id; // part id
		dscp_data[i].ss = scp[i-1].stga; // stage added
		dscp_data[i].es = scp[i-1].stgr; // stage removed
	}
	else if(scp[i-1].option == 48) // Part Set
	{
		dscp_data[i].st = scp[i-1].id; // part set id
		dscp_data[i].pt = 0; //  not a part 
		dscp_data[i].ss = scp[i-1].stga; // stage added
		dscp_data[i].es = scp[i-1].stgr; // stage removed
	}
}

Message("Finished Processing Construction Stage Data");

///// Set the view up for Stage 0
dcs_data[0] = new Object;
dcs_data[0].st = 0.0;
dcs_data[0].et = 0.0;
cstage = 0;

sflag = AllocateFlag();
Part.FlagAll(m, sflag);  // apply flag to all parts

for (j=1; j<=dscp_co; j++) // cycle through *DSCP array
{
	if(dscp_data[j].ss > cstage ) // clear flag from parts with a start stage > 0
	{
		if(dscp_data[j].st == 0 ) // not _SET
		{
			Part.GetFromID(m, dscp_data[j].pt).ClearFlag(sflag);
		}
		else // if _SET cycle through the *Set_Part
		{
			var pid;
			pset = Set.GetFromID(m, dscp_data[j].st, Set.PART);
			pset.StartSpool();
			while (pid = pset.Spool() )
			{
				Part.GetFromID(m, pid).ClearFlag(sflag);
			}
		}
	}
}

Part.BlankAll(m);
Part.UnblankFlagged(m, sflag);
Part.UnflagAll(m, sflag);
View.Redraw();

/////////////////// GUI

var we = new Window("Construction Stages", 0.1, 0.15, 0.8, 0.9 );

var l = new Widget(we, Widget.LABEL, 5, 55, 1, 2, "");

var t1 = new Widget(we, Widget.BUTTON, 10, 28, 5, 15, "<<");
t1.background = Widget.DARKRED; 
t1.foreground = Widget.WHITE;
t1.t = "rwd";
t1.onClick = bt_clicked;

var t2 = new Widget(we, Widget.BUTTON, 32, 50, 5, 15, ">>");
t2.background = Widget.DARKGREEN; 
t2.foreground = Widget.WHITE;
t2.t = "fwd";
t2.onClick = bt_clicked;

var lbl1 = new Widget(we, Widget.LABEL, 10, 45, 18, 23, "Current Stage: "+cstage);
var lbl2 = new Widget(we, Widget.LABEL, 5, 55, 25, 30, "Start Time: "+dcs_data[cstage].st+"  End Time: "+dcs_data[cstage].et);

var exit = new Widget(we, Widget.BUTTON, 10, 50, 34, 45, "Exit");
exit.background = Widget.DARKBLUE; 
exit.foreground = Widget.WHITE; 
exit.onClick = ex_clicked;

we.Show();

////////////////////////////////////////////////////////////////////////////////
function ex_clicked()
{
	Exit();
}
//////////////////////////////////////////////////////////////////////////////
function bt_clicked()
{
	if(this.t == "fwd") cstage = Math.min(dcs_co, cstage + 1); 
	if(this.t == "rwd") cstage = Math.max(0, cstage - 1);
        lbl1.text = "Current Stage: "+cstage;	
        lbl2.text = "Start Time: "+dcs_data[cstage].st+"  End Time: "+dcs_data[cstage].et;	

	Part.FlagAll(m, sflag);   // apply flag to all parts

	for (j=1; j<=dscp_co; j++) // cycle through *DSCP array
	{
		if(dscp_data[j].ss > cstage)  // clear flag from parts with a start stage > current stage
		{
			if(dscp_data[j].st == 0 ) // not _SET
			{
				Part.GetFromID(m, dscp_data[j].pt).ClearFlag(sflag);
			}
			else // if _SET cycle through the *Set_Part
			{
				var pid;
				pset = Set.GetFromID(m, dscp_data[j].st, Set.PART);
				pset.StartSpool();
				while (pid = pset.Spool() )
				{
					Part.GetFromID(m, pid).ClearFlag(sflag);
				}
			}
		}
		if(dscp_data[j].es <= cstage && dscp_data[j].es != 0)  // clear flag from parts with a end stage <= current stage
		{
			if(dscp_data[j].st == 0 ) // not _SET
			{
				Part.GetFromID(m, dscp_data[j].pt).ClearFlag(sflag);
			}
			else // if _SET cycle through the *Set_Part
			{
				var pid;
				pset = Set.GetFromID(m, dscp_data[j].st, Set.PART);
				pset.StartSpool();
				while (pid = pset.Spool() )
				{
					Part.GetFromID(m, pid).ClearFlag(sflag);
				}
			}
		}
	}

	Part.BlankAll(m);
	Part.UnblankFlagged(m, sflag);
	Part.UnflagAll(m, sflag);
	View.Redraw();
}
//////////////////////////////////////////////////////////////////////////////
