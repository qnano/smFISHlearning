/* Template to extract individual image series from .mvd2 files and save them as .tifs
 * 28 February 2020
 * --USAGE---
 * fiji --headless --console -macro ~/src/FIJI_macros/convert_mvd2_to_tif.ijm
 */

setBatchMode(true);
run("Bio-Formats Macro Extensions");

// Select .mvd file
//id = File.openDialog("Choose a file");
//id = "/Users/joshtitlow/tmp/transfer/060320_smFISH_MBONs/060320_smFISH_MBONs.mvd2";
id = "/Users/joshtitlow/tmp/smFISH_MS/20161210_DlgYFP_smFISH_2Channel/20161210_DlgYFP_smFISH_2Channel.mvd2";
print("Image path: " + id);

// Select output directory
//savedir = getDirectory("Choose a Storage Directory");
File.makeDirectory("ometiffs");
savedir = "/Users/joshtitlow/tmp/smFISH_MS/20161210_DlgYFP_smFISH_2Channel/ometiffs/";

// Determine the number of series in the file
Ext.setId(id);
Ext.getSeriesCount(seriesCount);
print("Series count: " + seriesCount);

// Process each image series individually
	//for (s=0; s<seriesCount; s++) {
	for (s=0; s<2; s++) {
		showProgress(s+1, seriesCount);
		print("processing ... " + "series" +s);

		// Open file
		//run("Bio-Formats", "open=id autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s);
    	run("Bio-Formats", "open=id color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s);
    	//open("/Users/joshtitlow/tmp/smFISH_MS/20161210_DlgYFP_smFISH_2Channel/20161210_DlgYFP_smFISH_2Channel.mvd2");
    	title = getTitle();
		savename = savedir + substring(title, 39) + ".ome.tiff";
		print (savename);
		//saveAs("Tiff", savename);
		run("OME-TIFF...", "save=" + savename + " export compression=Uncompressed");
		close();
	}

setBatchMode(false);
