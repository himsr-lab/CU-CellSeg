![header](https://user-images.githubusercontent.com/19319377/116188917-678dcd00-a6e5-11eb-91ca-9f0418843dff.png)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4599644.svg)](https://doi.org/10.5281/zenodo.4599644)
# CU-CellSeg
## ImageJ2 macro for the cell segmentation of multi-channel images

### Segmentation method
The CU-CellSeg macro implements a "classic" cell segmentation method with a marker-based cell expansion: In short, pixel probability maps for nucleus channels (mandatory) and cell matrix channels (optional) are thresholded and segmented by ImageJ's [watershed algorithm](https://imagej.nih.gov/ij/docs/guide/146-29.html#sub:Watershed) (nuclei) and by the ["Find maxima..." function](https://imagej.nih.gov/ij/docs/guide/146-29.html#sub:Find-Maxima) (cell matrix) to create individual cellular compartments. Typical nucleus channels would be markers for DAPI, dsDNA, or histone. Cell matrix channels could be cytoplasm markers (beta-tubulin, keratin, vimentin, ...) or membrane markers (V-ATPase, HLA class 1, CD8, ...), or both. If no cell matrix channel is specified (default), cell outlines are generated from the expansion of the nuclei outlines by a fixed radius. Both cell matrix options generate non-overlapping cell outlines.

![nuclei](https://user-images.githubusercontent.com/19319377/116176320-34404380-a6cf-11eb-998f-4f9d501c8398.png) ![matrix](https://user-images.githubusercontent.com/19319377/116176328-373b3400-a6cf-11eb-9588-298a12cf4f00.png)
**Figure 1: Segmentation example.** Detail from the center of the composite image (top). Left side: Grayscale dsDNA channel with nuclei outlines (blue). Right side: Grayscale beta-tubulin channel with cell matrix outlines (red).

### Software documentation
The documentation of our macros is in the corresponding source code: You can view the source code on GitHub by following the links to the macros.

### Software requirements
CU-CellSeg requires a recent version of the [Fiji](https://fiji.sc/) image processing package:
* ImageJ2 app (>= 1.53e)
* Bio-Formats plugin (>= 6.4.0)
* Trainable Weka Segmentation plugin (>= 3.2.24)

Any multi-channel image that can be imported with the Bio-Formats plugin can be processed by CU-CellSeg. However, you might have to adjust the [`suffixes`](https://github.com/christianrickert/CU-CellSeg/blob/1469060902e4277354b93b9649c065e6e55826b2/CU-CellSeg.ijm#L69) variable to select the file extensions for your specific instrument. In addition, if the metadata extraction and therefore the slice labeling fails, you will have to refer to individual channels explicitly by slice number, e.g. `"1", "14"` instead of using the convenient pattern matching available for labeled channels, e.g. `"beta-tubulin", "dsDNA"`.

### Example files
The [example folder](https://github.com/christianrickert/CU-CellSeg/tree/main/example) contains a single [MIBIscope™ image](https://github.com/christianrickert/CU-CellSeg/blob/main/example/20200109_3232_Run-16_FOV1_Final_3232_Top_R3C1_Tonsil.tiff?raw=true) (800x800 px) and two corresponding pixel-classifier models - trained with dsDNA ([nu.model](https://github.com/christianrickert/CU-CellSeg/blob/main/example/nu.model?raw=true)) and beta-tubulin ([ce.model](https://github.com/christianrickert/CU-CellSeg/blob/main/example/ce.model?raw=true)), respectively.
Running CU-CellSeg with the default [`userThresholds`](https://github.com/christianrickert/CU-CellSeg/blob/639e2b154d8455282ee769e337859a0713fe247f/CU-CellSeg.ijm#L75) and with the [`cellMatrixChannels`](https://github.com/christianrickert/CU-CellSeg/blob/639e2b154d8455282ee769e337859a0713fe247f/CU-CellSeg.ijm#L78) variable set to `newArray("beta-tubulin")`, should yield results identical to the data in the [results subfolder](https://github.com/christianrickert/CU-CellSeg/tree/main/example/20200109_3232_Run-16_FOV1_Final_3232_Top_R3C1_Tonsil). However, the measurements table (400 MiB) has been compressed with 7-zip before uploading.

CU-CellSeg produces distinct result files for every multi-channel image in the batch:
* `*.csv` - [measurements table](https://github.com/christianrickert/CU-CellSeg/blob/main/example/20200109_3232_Run-16_FOV1_Final_3232_Top_R3C1_Tonsil/20200109_3232_Run-16_FOV1_Final_3232_Top_R3C1_Tonsil.7z?raw=true)
* `*.tif` - [overview image](https://github.com/christianrickert/CU-CellSeg/blob/main/example/20200109_3232_Run-16_FOV1_Final_3232_Top_R3C1_Tonsil/20200109_3232_Run-16_FOV1_Final_3232_Top_R3C1_Tonsil.tif?raw=true)
* `*.txt` - [log file](https://github.com/christianrickert/CU-CellSeg/blob/main/example/20200109_3232_Run-16_FOV1_Final_3232_Top_R3C1_Tonsil/20200109_3232_Run-16_FOV1_Final_3232_Top_R3C1_Tonsil.txt)
* `*.zip` - [rois](https://github.com/christianrickert/CU-CellSeg/blob/main/example/20200109_3232_Run-16_FOV1_Final_3232_Top_R3C1_Tonsil/20200109_3232_Run-16_FOV1_Final_3232_Top_R3C1_Tonsil.zip?raw=true)

### Copyright notices
The [MIBIscope™ image](https://mibi-share.ionpath.com/tracker/overlay/sets/16/116) was kindly provided by Ionpath for demonstration purposes.
