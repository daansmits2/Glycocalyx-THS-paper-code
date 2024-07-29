
// Aligns two channels laterally, which is later useful for colocalization analysis. This macro is aimed to process single-slice Airyscan Processed images which has at least two channels.

close("*");

ch1="GCX"
ch2="ConA"
ch3="DAPI"

open();
print("\\Clear");
run("Clear Results");
//print(File.nameWithoutExtension);
file_basename = File.nameWithoutExtension;
print(file_basename);


imagedirectory = File.directory();
print(imagedirectory);



run("Subtract...", "value=10100 stack");
Stack.setChannel(1);
run("Enhance Contrast", "saturated=0.35");
Stack.setChannel(2);
run("Enhance Contrast", "saturated=0.35");
Stack.setChannel(3);
run("Enhance Contrast", "saturated=0.35");


saveAs("tiff", imagedirectory + file_basename + "_bgcor.tiff");


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
				
				setBatchMode(true);
				run("bUnwarpJ", "source_image=Con_A target_image=Glycocalyx_probe registration=Accurate image_subsample_factor=0 initial_deformation=[Very Coarse] final_deformation=Fine divergence_weight=0 curl_weight=0 landmark_weight=0 image_weight=1 consistency_weight=10 stop_threshold=0.01");
				selectWindow("Registered Target Image");
				run("Slice Remover", "first=3 last=3 increment=2");
				run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
				
				setBatchMode(true);
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