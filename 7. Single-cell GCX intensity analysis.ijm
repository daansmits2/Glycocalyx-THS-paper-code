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
       if (endsWith(path, ".czi")) {
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
			print(subdir);
			
			setBatchMode(true);
			roiManager("reset");
			run("Clear Results");
			
			rename("image");
			selectWindow("image");
				
				
				run("Z Project...", "projection=[Max Intensity]");

				run("Split Channels");
				selectWindow("C3-MAX_image");
				
				run("Gaussian Blur...", "sigma=2");
					setAutoThreshold("Otsu dark");
	
				setOption("BlackBackground", true);
				run("Convert to Mask");	
				run("Erode");

				//segmenting nuclei
				run("Analyze Particles...", "size=10-Infinity pixel exclude add");

				//select original 53BP1 window and plot ROIs on it to measure the number of foci
				selectWindow("C1-MAX_image");
				run("Enhance Contrast", "saturated=0.35");
				roiManager("Show All");
				//run("8-bit");
				numROIs = roiManager("count");
				for(z=0; z<numROIs;z++) {// loop through ROIs
					print(file_basename + "_"+z);
					run("Set Measurements...", "center integrated redirect=None decimal=3");
					roiManager("Select", z);
				
				selectWindow("C1-MAX_image");
				
				run("Measure");
							x = getResult("XM");
							xcor = x*pxsize;
							y = getResult("YM");
							ycor = y*pxsize;
				run("Clear Results");
				selectWindow("C1-MAX_image");
				rectsize= 35;
				makeRectangle(xcor-(rectsize/2), ycor-(rectsize/2), rectsize, rectsize);
				run("Duplicate...", "title="+z);
				
				run("In [+]");				
				run("In [+]");				
				run("In [+]");
				run("In [+]");
					
				run("Duplicate...", "title=threshold");
				selectWindow("threshold");
				run("In [+]");				
				run("In [+]");				
				run("In [+]");
				run("In [+]");
				
				
				run("Median...", "radius=2");
				
				setAutoThreshold("Otsu dark");
				run("Convert to Mask");
				
				//segmenting GCX signal
				run("Create Selection");
				
				roiManager("Add");
				
				selectWindow(z);
				
				roiManager("Select", numROIs);
				
				
				
				save(subdir+File.separator+z+".tif");
				rename(z);
				run("Set Measurements...", "mean redirect=None decimal=3");
				
				
				run("Measure");
				
				
				roiManager("Select", numROIs);
				roiManager("Delete");
				close("threshold");
				
					nR = nResults;
					label = newArray(nR);
					x1 = newArray(nR);
					y1 = newArray(nR);
					x2 = newArray(nR);
					y2 = newArray(nR);
					
					// Grab the old results
						x1 = getResult("Mean", 0);

					// Rename the old table
					IJ.renameResults("Raw Results");
					
					// Make the new table
					if(z>0){
					selectWindow("Endresults");
					IJ.renameResults("Results");
					}
						setResult("Mean", z, x1);
					updateResults();
					IJ.renameResults("Endresults");
					
					selectWindow("Raw Results");
					IJ.renameResults("Results");
					close("Results");
					
					close(z);
				
				selectWindow("C1-MAX_image");
				}
				selectWindow("Endresults");
saveAs("Results", directory + File.separator + file_basename + "_gcx_results.csv");
			};
close("*");
close("Results");
close("Threshold");
close("B&C");
close("ROI Manager");
      };
  
