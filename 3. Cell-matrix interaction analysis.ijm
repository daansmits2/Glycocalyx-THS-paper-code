close("*");

open();

if(isOpen("Log")){
	close("Log");
}
if(isOpen("ROI Manager")){
	close("ROI Manager");
}
if(isOpen("Results")){
	close("ROI Results");
}



Filebasename = File.nameWithoutExtension;

imagedir = getDirectory("image");


newdir=imagedir+Filebasename;


// Get the directory path from the user
roidir = getDirectory("Select a target folder for the .zip file for rois");
// Split the path by the file separator to get individual directories
splitPath = split(roidir, File.separator);
// Extract the last directory name (short name)
shortName = splitPath[splitPath.length - 1];
// Print the short name

// Split the string by "_"
splitString = split(Filebasename, "_");
// Replace the third part with the shortName
splitString[2] = shortName;
// Join the parts back together with "_"
Filebasename = String.join(splitString, "_");
// Print the modified string
print(Filebasename);





if (File.isDirectory(roidir + Filebasename + "_Analysis")) {
        print("Analysis directory already exists in the selected directory.");
    } else {
        // Create "roidir" directory in the selected directory
        File.makeDirectory(roidir + Filebasename + "_Analysis");
        print("Directory for analysis created in the selected directory.");
    }
    

roidir = roidir + File.separator + Filebasename + "_Analysis" + File.separator;
savepath = roidir + Filebasename;

print(savepath);


rename("image");

chromatic_aberration_correction = getBoolean("Do you want to control for chromatic abberration?");


if(chromatic_aberration_correction==1){
	run("Split Channels");
	selectWindow("C1-image");
	rename("Temp-image");
	slicenr = nSlices;
	setSlice(slicenr);
	newslicenr=slicenr+1;
	run("Add Slice");
	setSlice(1);
	run("Slice Keeper", "first=2 last="+newslicenr+" increment=1");
	rename("C1-image");
	close("Temp-image");
	waitForUser;
	run("Concatenate...", "open image1=C1-image image2=C2-image image3=C3-image");
	run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
	rename("Image");
}




run("Subtract...", "value=10000 stack");

setSlice(nSlices/2);

	Stack.setChannel(1);
run("Enhance Contrast", "saturated=0.35");
resetMinAndMax();
run("Magenta");

	Stack.setChannel(2);
run("Enhance Contrast", "saturated=0.35");
resetMinAndMax();
run("Cyan");
	Stack.setChannel(3);
run("Enhance Contrast", "saturated=0.35");
resetMinAndMax();
run("Yellow");
//run("Channels Tool...");
Property.set("CompositeProjection", "Sum");
Stack.setDisplayMode("composite");

saveAs("Tiff", savepath+"_bgcor.tiff");




rename("image");



	
	title = "Manual ROI specification?";
	  width=1024; height=1024;
	  specifyroi = newArray("No", "Yes");
	  Dialog.create("Do you want to manually specify ROI nr?");
	  Dialog.addChoice("Manually specify ROI:", specifyroi);
	  Dialog.show();
	  specifyroi = Dialog.getChoice();

		
		

	if(specifyroi=="Yes"){
		  roinr=1;
		  title = "Untitled";
		  Dialog.create("RoiNR of Roi to be analyzed?");
		  Dialog.addNumber("RoiNR:", roinr);
		  Dialog.show();
		  roinr = Dialog.getNumber();

	}



// Check if the file exists
if(File.exists(roidir+File.separator+Filebasename+"_Allrois_RoiSet.zip")==0){
	  roinr=1;
}else{
     roiManager("Open",roidir+File.separator+Filebasename+"_Allrois_RoiSet.zip");
     numberofrois = roiManager("count");
     roinr = numberofrois +1;
}






do{	
	
	if(roinr>1){
	run("ROI Manager...");
	roiManager("Open",roidir+File.separator+Filebasename+"_Allrois_RoiSet.zip");
	selectWindow("image");
	run("Select None");
	run("Duplicate...", "duplicate");
	rename("temp");
	roiManager("Show All with labels");
	close("ROI Manager");
	close("temp");
	selectWindow("image");
	}
	
	setTool("rotrect");
	selectWindow("image");
	waitForUser("Draw a square surrounding fiber/fragment of interest. Start with begin of structure");
	
	if(roinr>1){
	run("ROI Manager...");
	roiManager("Open",savepath+"_Allrois_RoiSet.zip");
	selectWindow("image");
	roicount = roiManager("count");
	if(roiManager("count")==roinr){
			roicount = roiManager("count");
			selectWindow("ROI Manager");
			roiManager("select", roicount-1);
			roiManager("delete");
			selectWindow("image");
			roiManager("add");
			roiManager("select", roicount-1);
			roiManager("rename", roinr);
			selectWindow("ROI Manager");
			roiManager("deselect");
			roiManager("Save", savepath+"_Allrois_RoiSet.zip");
			selectWindow("image");
			close("ROI Manager");
		}else{
		roiManager("add");
		roicount = roiManager("count");
		roiManager("select", roicount-1);
		roiManager("rename", roinr);
		selectWindow("ROI Manager");
		roiManager("deselect");
		roiManager("Save", savepath+"_Allrois_RoiSet.zip");
		selectWindow("image");
		close("ROI Manager");
		}
	}
	if(roinr==1){
	selectWindow("image");
	roiManager("add");
	selectWindow("ROI Manager");
	roiManager("select", roinr-1);
	roiManager("rename", roinr);
	roiManager("deselect");
	roiManager("Save", savepath+"_Allrois_RoiSet.zip");
	close("ROI Manager");
	}
	
	run("Flatten");
	saveAs("Tiff", savepath+"_"+roinr+"_flattened.tiff");
	close();
	
	
	selectWindow("image");
	run("Duplicate...", "duplicate");
	
	setTool("polyline");
	selectWindow("image-1");
	

	
	rename("roi");

	
	run("In [+]");
	run("In [+]");
	run("In [+]");
	run("In [+]");
	
	Stack.setChannel(1);
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(2);
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(3);
	run("Enhance Contrast", "saturated=0.35");
	
	Stack.setActiveChannels("101");
	
	saveAs("Tiff", savepath+"_"+roinr+".tiff");
	rename("roi");
	
	
		Stack.setDisplayMode("color");
		Stack.setChannel(2);
		roiwidth = getWidth();
		roiheight = getHeight();
		rectwidth = 20;
		
		nrectangles = floor(roiwidth/rectwidth);

		setTool("rectangle");
		
		x=0;
		y=0;
		nchannels=3;
		
		
	 title = "Single-slice analysis?";
	  width=1024; height=1024;
	  singleslice = newArray("No", "Yes");
	  Dialog.create("Do you want to analyze the image as single slice instead of MAX projection?");
	  Dialog.addChoice("Single slice:", singleslice);
	  Dialog.show();
	  singleslice = Dialog.getChoice();
	  
		
		for (i = 1; i <= nrectangles; i++) {
			run("Roi Defaults...", "color=white stroke=1 group=0");
			selectWindow("roi");
			makeRectangle(x, y, rectwidth, roiheight);
			Stack.setDisplayMode("composite");
			waitForUser("Set slice number so that the collagen fiber or segment of interest is in focus in the ROI");
			currentslice = getSliceNumber();
			currentslice = Math.ceil(currentslice/3);
			
			bottomslice = currentslice -(1);
			topslice = currentslice +(1);

			
			run("Duplicate...", "duplicate");

			rename("temp");
			
			if(singleslice=="Yes"){
				run("Duplicate...", "duplicate slices="+currentslice);
			}else{
				run("Make Substack...", "channels=1-3 slices="+bottomslice+"-"+topslice);
			}
			
			rename(i);
			close("temp");
			
			selectWindow(i);
			
			
			x=x+rectwidth;
			
			
			if(i==2){
				run("Combine...", "stack1=1 stack2=2");
				selectWindow("Combined Stacks");
				rename("Combined");
			
			}
			if(i>2){
				selectWindow("Combined");
				run("Combine...", "stack1=Combined stack2="+i);
				selectWindow("Combined Stacks");
				rename("Combined");
			}
		}
		
		
		Stack.setChannel(1);
		run("Magenta");
		run("Enhance Contrast", "saturated=0.35");
		Stack.setChannel(2);
		run("Enhance Contrast", "saturated=0.35");
		run("Cyan");
		Stack.setChannel(3);
		run("Enhance Contrast", "saturated=0.35");
		run("Yellow");
		
		
		saveAs("Tiff", savepath+"_"+roinr+"_isolatedsection.tiff");
		rename("multistack");
		
		
		
		
		if(singleslice=="No"){
		run("Z Project...", "projection=[Max Intensity]");
		
		saveAs("Tiff", savepath+"_"+roinr+"_isolatedsection_MAX.tiff");
		close("roi");
		rename("roi");
		}
		run("Duplicate...", "duplicate");
		
		
		lineposition = rectwidth-1;
		for (i = 1; i <= nrectangles; i++) {
			
			makeLine(lineposition, 0, lineposition, roiheight);
			Roi.setStrokeColor("white");
			Roi.setStrokeWidth(1);
			run("Draw", "stack");
			lineposition = lineposition + rectwidth;
		}
		
		saveAs("Tiff", savepath+"_"+roinr+"_isolatedsection_MAX_sectioning.tiff");
		rename("roi-1");
		close("roi-1");
		selectWindow("roi");
		
		if(singleslice=="Yes"){
			close("roi");
			selectWindow("multistack");
			rename("roi");
		}
		
linindexboolean = getBoolean("Is there a fiber of interest that you want to retrieve the linearity index of?");

if(linindexboolean==1){
	counter=0;
		do{// do..while loop to create multiple multi-focus images;
			counter=counter+1;
			close("ROI Manager");
			// Prompt the user to draw the segmented line using the "Segmented Line" tool
			run("Select None");
			
			setTool("polyline");
			waitForUser("Draw segmented line from cell-matrix adhesion to end of the fiber");
			roiManager("add");
			
			
			// Get the segmented line coordinates from ROI Manager
			roiManager("Show All");
			roiManager("Select", 0);
			
			
			// Get measurements
			run("Measure");
			totallength = getResult("Length");
			run("Clear Results");
			
			roiManager("Show All with labels");
			
			setTool("line");
			waitForUser("Draw line from the begin to the end of the previously-drawn line");
			roiManager("add");
			
			roiManager("Select", 1);
			// Get measurements
			run("Measure");
			shortestlength = getResult("Length");
			run("Clear Results");
			
			// Calculate linearity index
			linearityindex = shortestlength / totallength;
			
			roiManager("save", savepath + "_LinearityIndexROIs.zip");
			
			// Close ROI Manager
			roiManager("Reset");
			close("ROI Manager");
			
			additionalfiber = getBoolean("Is there another fiber of which linearity index should be determined?");
			
			
			if(counter==1){
				sumlinearity = linearityindex;
			}else{
				sumlinearity = sumlinearity + linearityindex;
			}
			
			
		}while(additionalfiber==1)
			print("ROInr: ", roinr);
			print("Linearity index: ", linearityindex);
			print("Average linearity calculated over: ", counter, " Fibers");
			linearityindex = sumlinearity/counter;

}else{
			linearityindex="NA";
			print("Linearity index: ", linearityindex);
}


// Delete .csv files matching to the currently-analyzed ROInr

// Get a list of files in the target folder
fileList = getFileList(savepath);


for (i = 0; i < fileList.length; i++) {
  currentFile = fileList[i];

  // Check if the file matches the specified pattern
  if (endsWith(Filebasename, "_Roinr_" + roinr)) {
    fileToDelete = targetFolder + currentFile;
    File.delete(fileToDelete);
    print("Deleted file: " + fileToDelete);
  }
}


	if(isOpen("Results")){
		close("Results");
	}
	
	title = "Section Type?";
	  width=1024; height=1024;
	  types = newArray("Relaxed", "Tense", "NA", "Background");
	  Dialog.create("Choose fiber type");
	  Dialog.addChoice("Type:", types);
	  Dialog.show();
	  type = Dialog.getChoice();
		  
		
	title = "Branching Phenotype?";
	  width=1024; height=1024;
	  branch = newArray("Main", "Branch", "NA");
	  Dialog.create("Choose branching phenotype of chosen structure");
	  Dialog.addChoice("Type:", branch);
	  Dialog.show();
	  branch = Dialog.getChoice();
		
		
		
		
	lineoverfiber = getBoolean("Do you want to draw a line over a fiber?");
	startroi = 0;
	
	additionallineroi = "Yes";

 do{//do..while loop for line rois in the particular multi-focus image
 	
 	
	// draw line over a fiber before cell-matrix interaction
	
			if(lineoverfiber==1){
			
					
						Stack.setActiveChannels("111");
						Stack.setDisplayMode("composite");
						setTool("polyline");
						waitForUser("Draw line over the fiber to the cell-matrix interaction");
			
			
					selectWindow("roi");
					Stack.setActiveChannels("111");
			
					roiManager("Add");
					roiManager("Select", startroi);
					Stack.setChannel(2);
					roiManager("Rename", "Col-fibertostartofcell");
					  run("Clear Results");
					  profile = getProfile();
					  for (i=0; i<profile.length; i++)
					      setResult("Value", i, profile[i]);
					  updateResults;
					  selectWindow("Results");
					  
					saveAs("Results",savepath+"_"+roinr+"_"+type+"_"+linearityindex+"_col_foreground_"+branch+"_beforecellstart.csv");
			
					Roi.setStrokeWidth(3);
					
		
			}else{
						Stack.setActiveChannels("111");
						Stack.setDisplayMode("composite");
						setTool("polyline");
						waitForUser("Draw small ROI to indicate beginning");
			
			
					selectWindow("roi");
					Stack.setActiveChannels("111");
					
					roiManager("Add");
					roiManager("Select", startroi);
					Stack.setChannel(2);
					roiManager("Rename", "Col-fibertostartofcell");
					  run("Clear Results");
					  profile = getProfile();
					  for (i=0; i<profile.length; i++)
					      setResult("Value", i, profile[i]);
					  updateResults;
					  selectWindow("Results");
					  
					saveAs("Results",savepath+"_"+roinr+"_"+type+"_"+linearityindex+"_col_foreground_"+branch+"_beforecellstart.csv");
			
					Roi.setStrokeWidth(3);
					
			}


	// draw line from cell-matrix interaction towards the cell body
			
				selectWindow("roi");
				Stack.setActiveChannels("111");
				Stack.setDisplayMode("composite");
				setTool("polyline");
				waitForUser("Draw line that matches collagen fiber or structure of interest");
				
			
			
					
					
					Roi.setStrokeWidth(3);
					roiManager("Add");
					
					roiManager("Set Color", "blue");
					roiManager("Select", startroi+1);
				    Stack.setActiveChannels("111");
				    Stack.setDisplayMode("composite");
				    run("Flatten");
				    run("Flatten");
					selectWindow("roi-1");
					saveAs("Tiff", savepath+"_Roinr_"+roinr+"_lineplot_col_flattened.tiff");
					close("roi-1");
					selectWindow("roi");
					
					
					selectWindow("roi");
					Stack.setActiveChannels("111");
			
					
					Stack.setChannel(2);
					roiManager("Rename", "Col");
					  run("Clear Results");
					  profile = getProfile();
					  for (i=0; i<profile.length; i++)
					      setResult("Value", i, profile[i]);
					  updateResults;
					  selectWindow("Results");
					  
					saveAs("Results",savepath+"_"+roinr+"_"+type+"_"+linearityindex+"_col_foreground_"+branch+"_aftercellstart.csv");
					
				
					
					
					roiManager("Select", startroi+1);
				    Stack.setActiveChannels("111");
				    Stack.setDisplayMode("composite");
				    run("Flatten");
					selectWindow("roi-1");
					saveAs("Tiff", savepath+"_Roinr_"+roinr+"_lineplot_b1int_flattened.tiff");
					close();
					selectWindow("roi");
					
					
					Stack.setActiveChannels("100");
					selectWindow("roi");
					waitForUser("If needed, move roi so that it aligns with the B1 integrin channel");
					roiManager("Add");
					roiManager("Select", startroi+2);
					roiManager("Rename", "B1int");
					Stack.setChannel(1);
					  run("Clear Results");
					  profile = getProfile();
					  for (i=0; i<profile.length; i++)
					      setResult("Value", i, profile[i]);
					  updateResults;
					  selectWindow("Results");
					saveAs("Results",savepath+"_"+roinr+"_"+type+"_"+linearityindex+"_B1int_foreground_"+branch+"_aftercellstart.csv");
					
					
					roiManager("Select", startroi +2);
					Stack.setActiveChannels("111");
				    Stack.setDisplayMode("composite");
				    run("Flatten");
					selectWindow("roi-1");
					saveAs("Tiff", savepath+"_Roinr_"+roinr+"_lineplot_gcx_flattened.tiff");
					rename("flattenedimage");
					close("flattenedimage");
					selectWindow("roi");
					
					if(isOpen("roi-2")){
						close("roi-2");
					}
					
					Stack.setActiveChannels("001");
					selectWindow("roi");
					Stack.setChannel(3);
					run("Clear Results");
					  profile = getProfile();
					  for (i=0; i<profile.length; i++)
					      setResult("Value", i, profile[i]);
					  updateResults;
					  selectWindow("Results");
					waitForUser("If needed, move roi so that it aligns with the gcxchannel");
					roiManager("Add");
					roiManager("Select", startroi+3);
					roiManager("Rename", "gcx");
					
					
					selectWindow("roi");
					Stack.setChannel(3);
					  run("Clear Results");
					  profile = getProfile();
					  for (i=0; i<profile.length; i++)
					      setResult("Value", i, profile[i]);
					  updateResults;
					  selectWindow("Results");
					saveAs("Results",savepath+"_"+roinr+"_"+type+"_"+linearityindex+"_gcx_foreground_"+branch+"_aftercellstart.csv");
					
					
					
					
					
					
					selectWindow("roi");
					
						
					Roi.setStrokeWidth(3);
					
					
				//draw roi for background measurement
					Stack.setActiveChannels("101");
					waitForUser("Draw line for to measure the background of B1int and gcx channels");
					
					
					
					roiManager("Add");
					roiManager("Select", startroi+4);
					roiManager("Rename", "B1int_gcx_background");
					
					Stack.setActiveChannels("100");
					Stack.setChannel(1);
					  run("Clear Results");
					  profile = getProfile();
					  for (i=0; i<profile.length; i++)
					      setResult("Value", i, profile[i]);
					  updateResults;
					  selectWindow("Results");
					saveAs("Results",savepath+"_"+roinr+"_"+type+"_"+linearityindex+"_B1int"+"_background_"+branch+"_aftercellstart.csv");
					
					Stack.setActiveChannels("001");
					roiManager("Select", startroi+4);
					Stack.setChannel(3);
					  run("Clear Results");
					  profile = getProfile();
					  for (i=0; i<profile.length; i++)
					      setResult("Value", i, profile[i]);
					  updateResults;
					  selectWindow("Results");
					saveAs("Results",savepath+"_"+roinr+"_"+type+"_"+linearityindex+"_gcx"+"_background_"+branch+"_aftercellstart.csv");
					
					
					
					roiManager("Save", savepath+"_"+roinr+"_"+type+".zip");
					

					close("ROI manager");
					
	
					
					title = "Additionallineroi?";
				  width=1024; height=1024;
				  additionallineroi = newArray("Yes", "No");
				  Dialog.create("Draw additional line ROI on current image?");
				  Dialog.addChoice("Type:", additionallineroi);
				  Dialog.show();
				  additionallineroi = Dialog.getChoice();	
					
				roinr = roinr +1;
			}while(additionallineroi=="Yes");
				
				
					close("roi");
					
	
	selectWindow("image");	
	newsegment = getBoolean("Analyze new segment?");
}while(newsegment==1);


close("Results");
close("image");
close("Results");

close("*");