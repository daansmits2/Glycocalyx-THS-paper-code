// This macro batch processes all the files in a folder and any
// subfolders in that folder.

   requires("1.33s"); 
   dir = getDirectory("Choose a Directory ");
   setBatchMode(true);
   count = 0;
   countFiles(dir);
   n = 0;
   processFiles(dir);
   //print(count+" files processed");
   
   function countFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/"))
              countFiles(""+dir+list[i]);
          else
              count++;
      }
  }

   function processFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/"))
              processFiles(""+dir+list[i]);
          else {
             showProgress(n++, count);
             path = dir+list[i];
             processFile(path);
          }
      }
  }

  function processFile(path) {
       if (endsWith(path, ".tif")) {
       	
       	
       	
       	close("*");
       	  
       	if(isOpen("Raw Results")){
       		close("Raw Results");
       	};
       	if(isOpen("Endresults")){
       		close("Endresults");
       	};
      	if(isOpen("ROI Manager")){
       		close("ROI Manager");
       	};
       	
           open(path);
			run("Clear Results");
			//print(File.nameWithoutExtension);
			run("Set Measurements...", "mean display redirect=None decimal=3");
			file_basename = File.nameWithoutExtension;
			print(file_basename);
			
			requires("1.32f");
		  title = getTitle;
		  width = getWidth;
		  height = getHeight;
		  depth = nSlices;
		  getPixelSize(unit, pw, ph, pd);
		
			pxsize = 1/pw;			
			
			directory = getDirectory("image");
			
			File.makeDirectory(directory + File.separator + file_basename);
			subdir = directory + file_basename;
			par=File.getParent(directory);
			
			print(subdir);
			print(par);
			dir_noext = File.getNameWithoutExtension(directory);
			
	
			roiManager("reset");
			run("Clear Results");
			
			rename("image");
			selectWindow("image");
				
			 GCAMP="C1";
			 DNAdamage="C2";
			 GCX="C3";
			 BF="C4";
			 thres="Huang";
			 thres2="Triangle";

			zstack = nSlices>4;
			
				if(zstack==1){
				run("Z Project...", "projection=[Max Intensity]");
				};
				
				run("Split Channels");
				if(zstack==1){
				selectWindow(DNAdamage+"-MAX_image");
				}else{
				selectWindow(DNAdamage+"-image");
				};
				run("Enhance Contrast", "saturated=0.35");
				
				run("Duplicate...", " ");
				
				run("Subtract Background...", "rolling=100 stack");
				
				run("Enhance Contrast", "saturated=0.35");
				
				run("Gaussian Blur...", "sigma=10");
				
			
				setAutoThreshold(thres+" dark");
				
				
				
				l=lengthOf(file_basename);
				l=l-6;
				suffix = substring(file_basename, l);
				print(suffix);
	
				setOption("BlackBackground", true);
				run("Convert to Mask");	
				//run("Erode");
				run("Dilate");
				run("Erode");
				run("Erode");
				run("Watershed");

				//segmenting nuclei
				run("Analyze Particles...", "size=300-30000 pixel exclude add");
				

				//select original DNAdamage window and plot ROIs on it to measure the number of foci
				
				if(zstack==1){
				selectWindow(GCAMP+"-MAX_image");
				}else{
				selectWindow(GCAMP+"-image");
				};
				
				run("Enhance Contrast", "saturated=0.35");
				roiManager("Show All");
				
				wait(500);
				
		
				run("Flatten");
				
				save(subdir+File.separator+"segmentationoverlay"+ ".png");
				save(directory+file_basename+"_segmentationoverlay"+ ".png");
				
				
				
				//run("8-bit");
				numROIs = roiManager("count");
				roiManager("deselect");
				
				roiManager("save", subdir+File.separator+file_basename+".zip");
				
				for(z=0; z<numROIs;z++) {
					
					// loop through ROIs
					//print(file_basename + "_"+z);
					
				
				if(zstack==1){
				selectWindow(GCAMP+"-MAX_image");
				}else{
				selectWindow(GCAMP+"-image");
				}
				
				roiManager("Select", z);
				run("Set Measurements...", "mean redirect=None decimal=3");
				run("Measure");
				
				
				
						if(zstack==1){
						selectWindow(GCAMP+"-MAX_image");
						}else{
							selectWindow(GCAMP+"-image");
						};
				};
				selectWindow("Results");
saveAs("Results", par + File.separator + dir_noext + suffix+"_gcamp_results.csv");
			};
close("*");
close("Results");
close("Threshold");
close("B&C");
close("ROI Manager");
      };

