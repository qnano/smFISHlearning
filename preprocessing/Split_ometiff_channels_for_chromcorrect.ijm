/* Seperate calibration images from data images in multichannel ometiffs
 *   5 November 2019
 *
 *  --USAGE--
 *   -specify indir
 *   -check that channels are setup properly
 *   -call from the command line:
 *     > fiji --headless --console -macro ~/src/FIJI_macros/Split_ometiff_channels_for_chromcorrect.ijm
 *
 */


setBatchMode(true);
run("ImageJ2...", "scijavaio=true loglevel=INFO");

// Get data dir and create output dirs
indir = "/usr/people/bioc1301/data/AdultBrain_smFISH_MASTER/250220_280220_smFISH_MBONs/20200331_repeat/";
//indir = getDirectory("choose image directory");
list = getFileList(indir);
File.makeDirectory(indir + "image");
File.makeDirectory(indir + "cal");

// set up loop to get image files
for (i=0; i<list.length; i++) { 
	if(endsWith(list[i],".r3d")){
		showProgress(i+1, list.length);
		print("processing ... "+i+1+"/"+list.length+"\n         "+list[i]);
		image_name = list[i];
		path=indir+list[i];
		
		// Open file
    		open(path);
		//run("Bio-Formats", "open=path color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		original = getTitle();
		run("Duplicate...", "duplicate channels=3-4");
		cal = getTitle();
		selectWindow(cal);
	
		// specify unit, which was lost in Huygens processing
		run("Properties...", "unit=micron");
		
		// save images
		saveAs("Tiff", indir+"cal/"+list[i]);
		close();

		selectWindow(original);
		run("Duplicate...", "duplicate channels=1-2");
		image = getTitle();
		selectWindow(image);

                // specify unit, which was lost in Huygens processing
                run("Properties...", "unit=micron");

		// save images
		saveAs("Tiff", indir+"image/"+list[i]);
		close();
		close();
	}
}



setBatchMode(false);
run("Quit");
