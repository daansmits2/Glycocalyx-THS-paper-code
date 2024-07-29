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
           open(path);
			run("Clear Results");
			//print(File.nameWithoutExtension);
			run("Set Measurements...", "mean display redirect=None decimal=3");
			file_basename = File.nameWithoutExtension;
			print(file_basename);
			
			directory = getDirectory("image");
			//setBatchMode(true);
			setBatchMode(false);
			roiManager("reset");
			run("Clear Results");
			//roiManager("delete");
			rename("image");
			selectWindow("image");
				
				run("Duplicate...", "title=[53bp1] duplicate channels=2");
				rename(file_basename);
				run("Duplicate...", "title=[thresholded]");
				
				//Nuclear segmentation based on 53BP1
				selectWindow("thresholded");
				//run("Subtract Background...", "rolling=1000 slice");
				run("Gaussian Blur...", "sigma=3");
				run("Threshold...");
				setAutoThreshold("Triangle dark");
				waitForUser("Set the correct threshold");
				if (getBoolean("Satisfied with the threshold?", "Yes", "No")==1) {
				run("Make Binary");
				run("Watershed");
				run("Analyze Particles...", "size=1000-100000 pixel circularity=0.05-1.00 exclude add");

						
				//select original 53BP1 window and plot ROIs on it to measure the number of foci
				selectWindow(file_basename);
				run("Enhance Contrast", "saturated=0.35");
				run("Enhance Contrast", "saturated=0.35");
				roiManager("Show All");
				run("Gaussian Blur...", "sigma=2");
				//run("8-bit");
				numROIs = roiManager("count");
				for(z=0; z<numROIs;z++) {// loop through ROIs
					run("Set Measurements...", "mean min  redirect=None decimal=3");
					roiManager("Select", z);
					run("Find Maxima...", "prominence=25 output=Count");
				}
			}
saveAs("Results", directory + File.separator + file_basename + "_foci_results.csv");
close("*");
close("Results");
close("Threshold");
close("B&C");
close("ROI Manager");
      }
  }
