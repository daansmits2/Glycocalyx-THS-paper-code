

close("*");

open();





linethickness=20;

imagename = File.nameWithoutExtension;
path = File.directory;

setTool("rectangle");
waitForUser("Draw rectangle to measure background");

roiManager("add");

roiManager("save", path + imagename + "_BgRoi.roi" );



run("Measure");

roiManager("delete");

backgroundMean = getResult("Mean", nResults()-1);
	
			close("Results");
			
			
			title = "SegmentRegions";
			  Dialog.create("Regionto be analyzed?");
			  Dialog.addChoice("Segment:", newArray("Single-cell-membrane", "Transitionzone", "Cell-cell-contact", "Granule"));
			  Dialog.show();
			  region = Dialog.getChoice();
			  
			  
			  title = "CellNR";
			  Dialog.create("CellNR?");
			  Dialog.addChoice("CellNR:", newArray("1", "2"));
			  Dialog.show();
			  cellnr = Dialog.getChoice();
			
if(region=="Granule"){
	iterations =1;
}else{
	iterations = 1;
}
			
			
			
for(s = 1; s<=iterations; s++){
			

			
			
			
			  title = "Segment";
			  Dialog.create("Subsegment to be analyzed?");
			  Dialog.addChoice("Subsegment:", newArray("Unbleached_Start_Single-cell-segment", "Bleached_Mid_Transitionzone", "Unbleached_Start_Cell-cell contact", 
			  "Bleached_Mid_Single-cell-segment", "Bleached_Start_Cell-cell-contact", "Unbleached_End_Single-cell-segment", "Unbleached_End_Cell-cell-contact", "Bleached_Mid_Granule", "Bleached_Mid_Cell-cell-contact"));
			  Dialog.show();
			  segment = Dialog.getChoice();
			
			
		
			
			setTool("polyline");
			waitForUser("Draw line over the chosen segment");
			
			run("Line Width...", "line="+linethickness);
			roiManager("add");
			roiManager("save", path + imagename +"_"+region+"_"+segment+ "_LineRoi.roi" );
			roiManager("delete");
			
			
			run("Set Measurements...", "mean integrated median display redirect=None decimal=3");
			
			// Clear the results window
			run("Clear Results");
			// Get the title of the image for context
			imageName = getTitle();
			
			
			for(n = 1; n <= nSlices; n++){
			    setSlice(n);  // Set the current slice
			   
			    time =  Stack.getFrameInterval();
			    getPixelSize(unit, pw, ph, pd);


			    profile = getProfile();  // Get the profile for the current slice
			    
			    // Display the profile in the Results window
			    for (i = 0; i < profile.length; i++) {
			        row = i + (n-1)*profile.length;  // Calculate row index
			        setResult("X", row, i);  // X coordinate (pixel position along the line)
			        setResult("Intensity", row, profile[i]);  // Intensity value
			        setResult("Background", row, backgroundMean);  // Add the background mean intensity to the row
			        setResult("Frame", row, n);  // Slice number
			        
			        setResult("Region", row, region);  
			        setResult("Segment", row, segment);  
			        setResult("Image", row, imageName);  // Image name
			        setResult("CellNR", row, cellnr);  // Image name
			        setResult("Pixelsize_micron", row, pw);
			        setResult("Time", row, time);
			    }
			}



saveAs("Results", path + imagename +"_"+region+ "_"+segment+ "_Profiles.csv");
	
			close("Results");
}