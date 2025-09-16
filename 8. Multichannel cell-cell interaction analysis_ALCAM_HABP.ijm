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

if (indexOf(Filebasename, "ALCAM")!=1 || indexOf(Filebasename, "HABP") != -1) {
    // Code to execute only if "ALCAM-" is found in Filebasename
    print("ALCAM found in filename: " + Filebasename);
    // your code here


chromatic_aberration_correction = getBoolean("Control for chromatic abberation by shifting the \n far red channel (ALCAM) one slice down?");


if(chromatic_aberration_correction==1){
	slicenr = nSlices;
	setSlice(slicenr);
	newslicenr=slicenr+1;
	run("Add Slice");
	setSlice(1);
	run("Slice Keeper", "first=2 last="+newslicenr+" increment=1");
}



}


newdir=imagedir+Filebasename;


roidir = getDirectory("Select a target folder for the .zip file for rois");


if (File.isDirectory(roidir + Filebasename + "_Analysis")) {
        print("Analysis directory already exists in the selected directory.");
    } else {
        // Create "roidir" directory in the selected directory
        File.makeDirectory(roidir + Filebasename + "_Analysis");
        print("Directory for analysis created in the selected directory.");
    }
    

roidir = roidir + File.separator + Filebasename + "_Analysis" + File.separator;
savepath = roidir + Filebasename;


  title = "BGcor";
  Dialog.create("Background correction?");
  Dialog.addChoice("Airyscan Background subtraction?", newArray("Yes", "No"));
  Dialog.show();
  bgcor = Dialog.getChoice();

if(bgcor=="Yes"){
run("Subtract...", "value=10000 stack");
}
saveAs("Tiff", savepath+"_bgcor.tiff");




rename("image");
run("Grays");"


  roinr=1;
  title = "Untitled";
  Dialog.create("RoiNR of Roi to be analyzed?");
  Dialog.addNumber("RoiNR:", roinr);
  Dialog.show();
  roinr = Dialog.getNumber();


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
	
 	Dialog.create("RoiNR of Roi to be analyzed?");
 	Dialog.setInsets(200,200,200);
	items = newArray("No", "Yes");
	Dialog.addChoice("Open existing roi?", items);
	Dialog.show();
  	openexistingroi = Dialog.getChoice();
  	
  	
  	setTool("rotrect");
	selectWindow("image");
	setSlice(nSlices/2);
	run("Enhance Contrast", "saturated=0.35");
  	if(openexistingroi=="Yes"){
  		path = File.openDialog("Choose an image file");
  		roiManager("Open", path);
  		roiManager("select", 0);
  		waitForUser("Select correct ROI and slice. Click OK when done");
  		close("ROI Manager");
  	}
	if(openexistingroi=="No"){
	waitForUser("Draw a square surrounding fiber/fragment of interest. Start with begin of structure");
	}
	
	
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
	
	
	
	selectWindow("image");
	run("Duplicate...", "duplicate");
	
	setTool("polyline");
	selectWindow("image-1");
	
	
	rename("roi");

	
	
	saveAs("Tiff", savepath+"_Roinr_"+roinr+".tiff");
	rename("roi");
	
		roiwidth = getWidth();
		roiheight = getHeight();
		rectwidth = 20;
		
		nrectangles = floor(roiwidth/rectwidth);

		setTool("rectangle");
		
		x=0;
		y=0;
		nchannels=1;
		
		//Projection thickness (odd number)
		Projectionthickness = 3;
		zslicedistance = (Projectionthickness-1)/2;
		
		
		for (i = 1; i <= nrectangles; i++) {
			run("Roi Defaults...", "color=white stroke=1 group=0");
			selectWindow("roi");
			makeRectangle(x, y, rectwidth, roiheight);
			
			run("Enhance Contrast", "saturated=0.35");
			//waitForUser("Set slice number so that the collagen fiber or segment of interest is in focus in the ROI");
			currentslice = getSliceNumber();
			print("Selected slice: ", currentslice);
			
			bottomslice = currentslice -(zslicedistance);
			topslice = currentslice +(zslicedistance);
			
			print("Slice range based on selection: ",bottomslice, "-",topslice);
			
			File.saveString(currentslice, savepath+"_Roinr_"+roinr+"_"+"Midslice_"+currentslice+".txt");
			
			run("Duplicate...", "duplicate");

			rename("temp");
			run("Make Substack...", "channels=1-"+nchannels+" "+ "slices="+bottomslice+"-"+topslice);
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
		
		
		saveAs("Tiff", savepath+"_Roinr_"+roinr+"_isolatedsection.tiff");
		rename("roi");
		
		
		
		saveAs("Tiff", savepath+"_Roinr_"+roinr+"_isolatedsection_Section.tiff");
		close("roi");
		rename("roi");
		
		
		
		run("Duplicate...", "duplicate");


// do.. while loop for number of ROIs to be generated
do{
	run("Select None");
// do while loop for linethickness
linethickness =1;
print("linethickness: ", linethickness);

print(getTitle());
		if(linethickness==1){
		lineposition = rectwidth-1;
			for (i = 1; i <= nrectangles; i++) {
				
				makeLine(lineposition, 0, lineposition, roiheight);
				Roi.setStrokeColor("gray");
				Roi.setStrokeWidth(rectwidth);
				run("Draw", "stack");
				lineposition = lineposition + rectwidth;
			}

		saveAs("Tiff", savepath+"_Roinr_"+roinr+"_isolatedsection_Section_linesectioning.tiff");
		rename("roi-1");
		close("roi-1");
		}
		selectWindow("roi");
		
		run("Enhance Contrast", "saturated=0.35");

		counter=0;
		newline=1;
	

	

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





 numlinesizes=10;
 
 

title = "Exclusively single-cell segment or including cell-cell contact?";
	  width=1024; height=1024;
	  ref = newArray("Contact", "Control");
	  Dialog.create("Choose reference of analysis");
	  Dialog.setInsets(200,200,200);
	  Dialog.addChoice("Type:", ref);
	  Dialog.show();
	  ref = Dialog.getChoice();	
	  
	  if(ref=="Contact"){
	  	numofmeasurments=3;
	  }
	  if(ref=="Control"){
	  	numofmeasurments=2;
	  }



for(l=1; l<=numofmeasurments; l++){
	
	 title = "Type of segment that will be drawn?";
				  width=1024; height=1024;
				  segment = newArray("Single-cell segment", "Cell-cell-contact", "Background");
				  Dialog.create("Choose type of segment that will be analysed next");
				  Dialog.setInsets(200,200,200);
				  Dialog.addChoice("Type:", segment);
				  Dialog.show();
				  segment = Dialog.getChoice();	
	
setTool("polyline");
setSlice(nSlices-1);
waitForUser("Draw line over the chosen segment. For single-cell, draw over the transition until the start of a fixed cell-cell contact region");


   labelArray=newArray();
			    
			    
			   

if(linethickness==1){
		
 		Dialog.create("Open existing line profile?");
 		Dialog.setInsets(200,200,200);
		items = newArray("No", "Yes");
		Dialog.addChoice("Open existing line profile ROI?", items);
		Dialog.show();
  		openexistingline = Dialog.getChoice();
  		
  		
  		title = "CellNR? where lineplot will start?";
	  width=1024; height=1024;
	  types = newArray("1", "2");
	  Dialog.create("Choose CellNR type");
	  Dialog.setInsets(200,200,200);
	  Dialog.addChoice("Type:", types);
	  Dialog.show();
	  type = Dialog.getChoice();
		  
  		
  		setSlice(Projectionthickness-1);
  		
		if(openexistingline=="Yes"){
  		path = File.openDialog("Choose an image file");
  		roiManager("Open", path);
  		roiManager("select", 0);
  		roiManager("delete");
  		close("ROI Manager");
  		waitForUser("OK like this?");
	  	}
	  	
	  	
	  	
  	
	if(isOpen("Results")){
		close("Results");
	}
	drawline = 1;
		  
	}





	for(linethickness=1; linethickness<=numlinesizes; linethickness++){		
	
	
	

		  
	if(drawline==1){
		print("Line OKAY");
		Roi.setStrokeWidth(linethickness);
	}
		
		
		
		selectWindow("roi");

		
			// Initialize an array to store the mean profile
			meanProfile = newArray(nSlices);
			
			// Loop through each slice and generate line profiles
			for (slice=1; slice<=nSlices; slice++) {
			    selectImage("roi");
			    setSlice(slice);
			
			    profile = getProfile();
			    
			   if(slice==1){
				     getSelectionCoordinates(x, y);
				     for (g=0; g<x.length; g++){
				     	if(g==x.length-1){
				         print("Please continue the plot profile from coordinates: "+x[g]+" "+y[g]);
				     	}
				     }
			   }
			    
			    
			    // Add the current profile values to the meanProfile array
			    for (i = 0; i < profile.length; i++) {
			        if (slice == 1) {
			            meanProfile[i] = profile[i];
			        } else {
			            meanProfile[i] += profile[i];
			        }
			    }
			    
			 
			    
			    
			    // Display the mean profile values in the Results window
			for (i = 0; i < meanProfile.length; i++) {
			    setResult("Value", i, meanProfile[i]);
			    labelArray[i] = segment;
			}
			
			
			
			updateResults();
			selectWindow("Results");
			
			
			   
			}
			
			// Calculate the mean profile by dividing each value by the number of slices
			for (i = 0; i < meanProfile.length; i++) {
			    meanProfile[i] /= nSlices;
			}
			
			// Display the mean profile values in the Results window
			for (i = 0; i < meanProfile.length; i++) {
			    setResult("Value", i, meanProfile[i]);
			}
			
		roiManager("add");
		roiManager("Show All");
		run("Select None");
		
		updateResults();
		selectWindow("Results");
		
				  
				  
		  
		saveAs("Results",savepath+"_"+roinr+"_"+type+"_"+ref+"_"+segment+"_"+linethickness+"_gcx_cell-cell-contact.csv");
		close("Results");
		
		if(linethickness==1){
			
			roiManager("Select", 0); //Delete the first ROI
			roiManager("Save", savepath+"_"+roinr+"_"+type+"_"+ref+"_"+segment+"_Lineprofile.roi");
		}
		roiManager("Select", 0); //Delete the first ROI
		roiManager("delete");
		
		rename("roi");
	
		
		
		if(drawline!=1){
			print("Line not OKAY");
			close("roi");
			File.delete(savepath+roinr+".tiff");
			run("ROI Manager...");
			roiManager("Open",savepath+"_Allrois_RoiSet.zip");
			roiManager("select", roinr-1);
			roiManager("Delete");
			roiManager("Save", savepath+"_Allrois_RoiSet.zip");
			close("ROI Manager");
		}
		
		print("linethickness: ", linethickness);
	}
	
	
	
}
		roinr = roinr +1;
		roiManager("Save", savepath+"_"+roinr+"_"+type+".zip");
		close("ROI Manager");
		newline = getBoolean("Analyze new segment?");
		
		roinr=roinr+1;
			
	}while(newline==1);
	
		roinr=roinr+1;
		
		roiManager("Save", savepath+"_"+roinr+"_"+type+".zip");
		close("ROI Manager");

	close("roi");
	selectWindow("image");	
	newsegment = getBoolean("Perform new ROI analysis?");
}while(newsegment==1);


close("Results");
close("image");
close("Results");