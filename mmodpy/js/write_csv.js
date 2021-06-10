// <AUTHOR/> M BOWERS <\AUTHOR>
// This script is a snippet of code which can be dropped into any primer/d3Plot Java to create and write out to a csv file of your choice

var filename = "C:\\output.csv"; //change to whatever you want

var file = new File(filename, File.WRITE); //open file for writing

var data = new Array(); //replace with your array

for (var i=0; i<data.length; i++) //for all cells in data
{
	file.Writeln(data[i]); //Print contents of cell i to file

	/*
	//if you're writing more than one thing on each line, separate with commas like so:
	file.Writeln(data[i] + "," + data[i] + "," + data[i]);
	*/

}

file.Close(); //close file
