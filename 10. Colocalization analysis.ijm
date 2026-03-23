
  requires("1.33s"); 
  
  close("Results");
   dir = getDirectory("Choose a Directory");


autothrA = "Otsu";
autothrB = "Li";
ratio_dbcoths = 1;
thra = 300;
thrb = 50;
rep= 0;

// for MV3 exp
// thra = 80 (ConA)
//thra = 80;
// thrb = 30 (THS) , 100 (DBCO)
//thrb = 30;
// 100/ 30 = 3.33
//ratio_dbcoths = 3.33;

manual = true;
aspecificbackground = 0;
skip_files_with = "Ac4+"


setBatchMode(true);
checkimages = getBoolean("check images?");
if(checkimages ==1){
	setBatchMode(false);
}else{
	setBatchMode(true);
}

   parent= File.getParent(dir);
   print(parent);
   count = 0;
   countFiles(dir);
   n = 0;
   
   processFiles(dir);
   
   
   print(count+" files processed");
   
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

  function processFile(path){ // analysis function

       if (endsWith(path, "bgcor.tiff") || endsWith(path, "bgcor.tiff") ) {
		roiManager("reset");


			//setBatchMode(true);
			
			close("*");
			
			ch1="ConA";
			ch2="Gcx";
			ch3="DAPI";
			
			
			
			open(path);
			
			Stack.setChannel(2);
			run("Subtract...", "value="+aspecificbackground);
			
			
			print("\\Clear");
			//run("Clear Results");
			//print(File.nameWithoutExtension);
			file_basename = File.nameWithoutExtension;
			print(file_basename);
			image = File.name;
			
			if(nSlices>2){
			
			run("Slice Remover", "first=3 last=3 increment=1");
			};
			
			
			imagedirectory = File.directory();
			print(imagedirectory);
			
			
			ev = startsWith(file_basename, skip_files_with);
			print(ev);
			if(ev==1){
				continue;
			}
			file_basename = File.nameWithoutExtension;
			
				if(indexOf(file_basename, "DBCO") != -1){
				    print("Correction for background ratio-wise performed");
				    thrb = thrb * ratio_dbcoths;
				}else{
				    print("No ratio-wise background correction");
				}


			
					
			rename("image");
			
			Property.set("CompositeProjection", "Sum");
			Stack.setDisplayMode("composite");
			
			Stack.setActiveChannels("110");
			
			 roi_path = imagedirectory + file_basename + "_BgRoi.roi";
    
			    // Check if the ROI file already exists
			    if (File.exists(roi_path)) {
			    	
			        print("ROI already exists: " + roi_path);
			             roiManager("Open", roi_path);
        roiManager("Select", 0);
			    } else {
			    	
			setTool("rectangle");
			waitForUser("Draw rectangle to measure background");
			        roiManager("add");
			        roiManager("save", roi_path);
			        print("ROI saved: " + roi_path);
			    };
			    
			 
			Stack.setChannel(1);
			getStatistics(area, mean, min, max, std, histogram);
			backgroundMean_1 = mean;
			print(backgroundMean_1);
			run("Select None");
			selectWindow("image");
			run("Subtract...", "value="+backgroundMean_1);
			Stack.setChannel(2);
			roiManager("select", 0);
			getStatistics(area, mean, min, max, std, histogram);
			backgroundMean_2 = mean;
			print(backgroundMean_2);
			run("Select None");
			run("Subtract...", "value="+backgroundMean_2);
			selectWindow("image");
			roiManager("delete");
			
			
			
				  requires("1.32f");
				  title = getTitle;
				  width = getWidth;
				  height = getHeight;
				  depth = nSlices;
				  getPixelSize(unit, pw, ph, pd);
				
				  print("Title: " + title);
				  print("Size: " + width*pw+"x"+height*ph+"x"+depth*pd+" " + unit);
				  if (unit!="pixel" || pd!=1) {
				      print("Pixel Size: "+pw+"x"+ph+"x"+pd + " " + unit);
				      if (pw==ph)
				          print("Resolution: "+1/pw+" pixels per "+unit);
				      else {
				          print("X Resolution: "+1/pw+" pixels per "+unit);
				          print("Y Resolution: "+1/ph+" pixels per "+unit);
				      }
				  }
			
			run("Set Measurements...", "mean standard modal min redirect=None decimal=3");
		
			if(rep==1){				
				file_basename = File.nameWithoutExtension; // Get base name of the file
				file_basename = replace(file_basename, "Ac4+", "Ac4ManNAz");
				file_basename = replace(file_basename, "Ac4ManNAz+", "Ac4ManNAz");
				file_basename = replace(file_basename, "+", "");
			}
				print(file_basename);

				imagedirectory = File.directory;           // Get directory of the file
				saveAs("tiff", imagedirectory + file_basename + ".tiff");
				rename("image");
				cellcount = 1;  
				
				numpattern = "\\d+";                       // Regex for one or more digits
				searchString = file_basename + "_" + numpattern + "_CellRoi.roi"; // Combine basename and pattern
				files = getFileList(imagedirectory);       // Get list of all files in the directory
				
				for (i = 0; i < files.length; i++) {
				    if (files[i].matches(searchString)) {  // Use the dynamic regex to match
				        print(files[i] + " matches the pattern.");
				        fullFilePath = imagedirectory + files[i]; // Construct the full file path
				        roiManager("open",fullFilePath);               // Open the matched file
				    } else {
				        print(files[i] + " does not match.");
				    }
				}
				nRois = roiManager("count");	
				
	if(nRois==0){   // only proceed to define ROIs if there are 0 rois previously defined
					
				
					    do{   //proceed with defining cell ROIs
						
					    setBatchMode(false);
					    setTool("rectangle");
					    selectImage("image");
					   	
						setTool("polygon");
						waitForUser("Put cell of interest in ROI");		
					   
					
					    roiManager("add");
					    print("Drawing ROI for cellnr: ", cellcount);
					    cellpath = imagedirectory + file_basename +"_"+cellcount;
					    cellpath_save  = cellpath + "_CellRoi.roi";
					    
					    rois = roiManager('count');
					    roiManager("select", rois-1);
						roiManager("save", cellpath_save);
					    
					    run("Select None");

					 	
						proceed = getBoolean("analyze another cell?");
					 	cellcount ++;
					 	
					   }while(proceed==true);	//do..while loop for defining individual cells
					    
			}else{
			checkimages =0;	 // if loop for if no Rois can be found from previous quantifications
			}
			n = roiManager('count');
			roiManager("reset");
			
			for (i = 1; i <= n; i++) {  //start the ROI loop (individual cell analysis) 
    		
			print(n, " Rois of whole cells inside this image were discovered");
    		print("CellNR: ",i);
    		print("Path: ", imagedirectory + file_basename + "_" + i+"_CellRoi.roi");
    		roiManager("open", imagedirectory + file_basename + "_"+i+"_CellRoi.roi");
		    	selectWindow("image");
		    	roiManager("select", 0);
		    	run("Duplicate...", "duplicate");
		    	run("Clear Outside");
		    	run("Select None");
		    	save(imagedirectory + file_basename + "_" + i);
    	
    	
			 intracellular_roi_path = imagedirectory +file_basename+"_"+ i + "_IntracellularRoi.roi";
    
			    // Check if the ROI file already exists
			    if (File.exists(intracellular_roi_path)) {
			    	
			        print("ROI already exists: " + intracellular_roi_path);
			       roiManager("Open", intracellular_roi_path);
			     
        			roiManager("Select", 1);
			    }else{
				setTool("polygon");		
				
				
				Stack.setChannel(1);
				waitForUser("Encircle intracellular part");
			   
				
				roiManager("add");
				nRois = roiManager("count");
				roiManager("select", nRois-1);
				roiManager("save", intracellular_roi_path);
				
			    };
			   
				roiManager("select", 1);
				
				
				
				run("Clear", "slice");
				Stack.setChannel(2);
				
				run("Clear", "slice");
				
				roiManager("select", 1);
				roiManager("delete");
			
				run("Select None");
			
		        
				run("Convert to Mask", "method=Huang background=Dark calculate black create");
				selectWindow("MASK_image-1");
				run("Create Selection");
				roiManager("add");
				selectWindow("image-1");
				newroi = roiManager('count');
				roiManager("select", newroi-1);
				
				Stack.setChannel(1);
			    getStatistics(area, mean, min, max, std, histogram);
			    max_1 = max;
			    mean_1 = mean;
			    min_1 = min;
			    print(max_1, mean_1, min_1);
			    Stack.setChannel(2);
			    getStatistics(area, mean, min, max, std, histogram);
			    max_2 = max;
			    mean_2 = mean;
			    min_2 = min;
			    
			    print(max_2, mean_2, min_2);
			    run("Select None");
			    
				roiManager("select", newroi-1);
				
			    roiManager("delete");
			    selectWindow("image-1");
			    run("Select None");
			    close("MASK_image-1");
			    
			    
				

			    
							
				run("Set Measurements...", "mean standard modal min redirect=None decimal=3");
				
				roiManager("reset");
				
				Stack.setChannel(1);
				resetMinAndMax;
				
				//max_1 = 5500;
				//max_2 = 7000;
				//setMinAndMax(min_1, max_1);
				
				run("Enhance Contrast", "saturated=1");
				Stack.setChannel(2);
				//resetMinAndMax;
				//setMinAndMax(min_2, max_2);
				
				run("Enhance Contrast", "saturated=1");
				
				
				saveAs("tiff", imagedirectory + file_basename + "_"+i+ "_image_membraneonly.tiff");
				
				rename("image-1");
				
				if(manual==true){
				run("BIOP JACoP", "channel_a=1 channel_b=2 threshold_for_channel_a=[Use Manual Threshold Below] threshold_for_channel_b=[Use Manual Threshold Below] manual_threshold_a="+thra+" manual_threshold_b="+thrb+" get_pearsons get_manders get_overlap get_fluorogram costes_block_size=5 costes_number_of_shuffling=100");
				}else{
				run("BIOP JACoP", "channel_a=1 channel_b=2 threshold_for_channel_a="+autothrA+" threshold_for_channel_b="+autothrB+" manual_threshold_a="+thra+" manual_threshold_b="+thrb+" get_pearsons crop_rois get_manders get_overlap get_fluorogram costes_block_size=5 costes_number_of_shuffling=100");	
				};
				saveAs("tiff", imagedirectory + file_basename + "_"+i+ "_colocalization.tiff");
				
  
			selectWindow("Results");
			saveAs("Results", imagedirectory+"Coloc_results.csv");
			
			
						if(checkimages ==1){
					waitForUser("Check image");
				};// closes the if check images option;
		
				  // setBatchMode(true);
				 setResult("Name", nResults-1, file_basename);
				setResult("Max_A", nResults-1, max_1);
				
				setResult("Max_B", nResults-1, max_2);
				
				setResult("Mean_A", nResults-1, mean_1);
				
				setResult("Mean_B", nResults-1, mean_2);
				
				setResult("Cellcount", nResults-1, intracellular_roi_path);
				
				  
				  setResult("PixelSize",nResults-1, pw);
				  
				updateResults();
				close;
				run("Select None");
				
						
		
	
			//close("Log");
			close("Threshold");
			close("B&C");
			close("Channels");
			close("image-1");
		
		}; //closes the ROI loop;
		
		
  };  // closes the scanning in folder loop;
  
  }; // finishes the analysis function
  
  close("*");