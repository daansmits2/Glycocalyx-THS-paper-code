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
File.makeDirectory(newdir + "_SNR_analysis");

savepath = newdir + "_SNR_analysis"+ File.separator+Filebasename;

if(File.exists(savepath+"_"+"AllROIs_SNR.zip")){
roiManager("Open", savepath+"_"+"AllROIs_SNR.zip");
};


bgsub=getBoolean("Background subtract?");


if(bgsub==1){
run("Subtract...", "value=10000 stack");
}else{
}

setSlice((nSlices/2)+1);
run("Enhance Contrast", "saturated=0.35");

saveAs("Tiff", savepath+"_bgcor.tiff");




rename("image");
	
	
	
	run("ROI Manager...");
	selectWindow("image");
	run("Select None");
	run("Duplicate...", "duplicate");
	rename("temp");
	roiManager("Show All with labels");
	
	close("temp");
	
	selectWindow("image");
	
	highres=getBoolean("high res image?");
	
	if(highres==1){
	setTool("line");
	}else{
		setTool("polyline");
	}
	
	if(highres==1){
	zslicesfrombottom =15;
	}else{
	zslicesfrombottom=2;
	}
	
	waitForUser("Set slice to ventral membrane");
  
  
  	bottom =getSliceNumber();
  	
	//waitForUser("Set slice to dorsal membrane");
	
	//top = getSliceNumber();
	
  slicenr = bottom + zslicesfrombottom;

  
  run("Set Slice...",  "slice="+slicenr);

	// Slice number is just an estimate. Please edit when the slice is not completely perpendicular.
		
	do{
	
	selectWindow("image");
	
	waitForUser("Draw line over GCX from extracellular -> intracellular");
	
	if(highres==1){
	Roi.setStrokeWidth(1);
	}else{
	Roi.setStrokeWidth(0);
	}
	
	roiManager("add");
	
	roinr = roiManager("count");
	
	roinr = roinr-1;
	
		roiManager("Select", roinr);
		  run("Clear Results");
		  profile = getProfile();
		  for (i=0; i<profile.length; i++)
		      setResult("Value", i, profile[i]);
		  updateResults;
		  selectWindow("Results");
		  
		saveAs("Results",savepath+"_"+roinr+"_SNR.csv");
		
	newsegment = getBoolean("Analyze new segment?");
}while(newsegment==1);

close("Results");	


roiManager("save", savepath+"_"+"AllROIs_SNR.zip");
	
	close("ROI Manager");
	close("*");
	
	
	
	