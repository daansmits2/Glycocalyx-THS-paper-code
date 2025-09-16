
// Aligns two channels laterally, which is later useful for colocalization analysis. This macro is aimed to process single-slice Airyscan Processed images which has at least two channels.


   dir = getDirectory("Choose a Directory");
   parent= File.getParent(dir);
   print(parent);
   count = 0;
   countFiles(dir);
   n = 0;
   
   registration = getBoolean("perform registration?");
	if(registration!=1){
		print("No registration performed");
	}else{
		print("Registration performed");
	}
   
   
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

  function processFile(path){



       if (endsWith(path, "Processing.czi")) {



setBatchMode(true);
close("*");

ch1="GCX";
ch2="ConA";
ch3="DAPI";

open(path);
print("\\Clear");
run("Clear Results");
//print(File.nameWithoutExtension);
file_basename = File.nameWithoutExtension;
print(file_basename);


imagedirectory = File.directory();
print(imagedirectory);


run("Subtract...", "value=10000 stack");
Stack.setChannel(1);
run("Enhance Contrast", "saturated=0.35");
Stack.setChannel(2);
run("Enhance Contrast", "saturated=0.35");
Stack.setChannel(3);
run("Enhance Contrast", "saturated=0.35");


saveAs("tiff", imagedirectory + file_basename + "_bgcor.tiff");

			if(registration==1){
			
			
						roiManager("reset");
						run("Clear Results");
						//roiManager("delete");
						rename("ch1");
							
										
							run("Split Channels");
							selectWindow("C1-ch1");
							rename("Glycocalyx_probe");
							selectWindow("C2-ch1");
							rename("Con_A");
							
							selectWindow("C3-ch1");
							rename(ch3);
							run("16-bit");
							
							//setBatchMode(true);
							run("bUnwarpJ", "source_image=Con_A target_image=Glycocalyx_probe registration=Accurate image_subsample_factor=0 initial_deformation=[Very Coarse] final_deformation=Fine divergence_weight=0 curl_weight=0 landmark_weight=0 image_weight=1 consistency_weight=10 stop_threshold=0.01 verbose");
							selectWindow("Registered Target Image");
							run("Slice Remover", "first=3 last=5 increment=1");
							run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
							
							//setBatchMode(true);
							Stack.setChannel(1);
							resetMinAndMax();
							run("16-bit");
							run("Make Composite", "display=Composite");
							Stack.setChannel(1);
							run("Enhance Contrast", "saturated=0.35");
							//run("Apply LUT");
							run("Yellow");
							Stack.setChannel(2);
							resetMinAndMax();
							run("16-bit");
							run("Enhance Contrast", "saturated=0.35");
							run("Magenta");
							
							
							
							
							saveAs("tiff", imagedirectory + file_basename + "_bgcor_registered.tiff");
				}
			}
				
  }