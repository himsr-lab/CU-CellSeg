/*  Copyright 2021 Regents of the University of Colorado
 *
 *  This file is part of CU-CellSeg.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Author:       Christian Rickert <christian.rickert@cuanschutz.edu>
 *  Group:        Human Immune Monitoring Shared Resource (HIMSR)
 *                University of Colorado, Anschutz Medical Campus
 *
 *  Title:        CU-CellSeg
 *  Summary:      ImageJ2 macro for the cell segmentation of multi-channel images
 *
 *  DOI:          https://doi.org/10.5281/zenodo.4599644
 *  URL:          https://github.com/christianrickert/CU-CellSeg/
 *
 *  Description:
 *
 *  Creates normalized projections for the selected 'regionChannels' and classifies
 *  these into the 'regionClasses' using the "Trainable Weka Segmentation" plugin.
 *  The corresponding region labels can be median-filtered to suppress smaller
 *  regions while preserving the general dimensions (no blurring).
 *  The resulting segmentation maps will be converted into "ROI Manager" entries.
 *  The macro saves a log file, the binary segmentation maps, and the regions of
 *  interest into the regions folder of each input image.
 *  The traces used for the training and the model used for the classification
 *  are stored with the input images.
 *
 *  Dependencies:
 *
 *  CU-CellSeg requires a recent version of CU-MacroLibrary to be installed:
 *  https://github.com/christianrickert/CU-MacroLibrary/
 */

print("\\Clear");

run("Bio-Formats Macro Extensions");

batchMode = true;
medianFilter = 5;  // radius [units]
regionChannels = newArray("autofluorescence", "ck", "dapi");
regionClasses = newArray("Tumor (positive)", "Stroma (negative)", "Glass (neutral)");
regionFolder = "regions";
scalingFactor = 0.25;  // increase classification speed by downscaling
suffixes = newArray(".tif", ".tiff");
files = getFilesInFolder("Select the first TIFF of your dataset", suffixes);
targetName = "tu_st_gl";  // class label and file output
versionString = "v1.00 (2022-02-11)";
processFolder(files);

/*
 *  Loop
 */

// Function to process files with matching suffixes from a folder
function processFolder(files)
{
  files_length = files.length;
  for ( i = 0; i < files_length; ++i )
  {
    processFile(files[i]);
  }
}

function processFile(file)
{
  // clear previous run
  initializeRun(versionString);
  clearAllSelections();

  // print current file name
  print("\n*** Processing file ***");
  print("\t" + file);

  // read image file
  filePath = File.getDirectory(file);
  fileName = File.getName(file);
  fileSlices = readImage(file);
  fileTitle = getTitle();

  // get calibration data from image file
  pixelCalibration = getPixelCalibration();
  if ( pixelCalibration[0] == "pixels" )  // no conversion required
    toPixels = 1.0;
  else  // calculate conversion factor
    toPixels = Math.round(1.0 / 0.5 * (pixelCalibration[1] + pixelCalibration[2]));
  print("\tCalibration: " + toPixels + " pixel per " + pixelCalibration[0]);

  // create region projection
  setBatchMode(batchMode);  // avoid screen flickering during stack preparation
  projectedRegion = projectStack(fileTitle, fileSlices, regionChannels, targetName);
  if ( scalingFactor != 1.0 )  // (optional) rescale image file
  {
    print("\tScaling factor: " + scalingFactor);
    projectedRegion = projectedRegion + "->(*f)";
    getDimensions(width, height, channels, slices, frames);
    run("Scale...", "x=" + v2p(scalingFactor) + " y=" + v2p(scalingFactor) +
                   " z=1.0 width=" + v2p(width) + " height=" + v2p(height) + " depth=" + v2p(slices) +
                   " interpolation=Bicubic average process create" +
                   " title=" + v2p(projectedRegion));
  }
  setBatchMode("exit and display");  // batch mode incompatible with Trainable Weka Segmentation plugin

  // create region classification
  classifiedRegion = classifyImage(projectedRegion, targetName, filePath);

  // batch mode safe from here
  toggleBatchMode(batchMode, false);

  // (optional) rescale image file back
  if ( scalingFactor != 1.0 )
  {
    invertedScalingFactor = 1.0 / scalingFactor;
    print("\tScaling factor: " + invertedScalingFactor);
    classifiedRegion = classifiedRegion + "->(*1/f)";
    getDimensions(width, height, channels, slices, frames);
    run("Scale...", "x=" + v2p(invertedScalingFactor) + " y=" + v2p(invertedScalingFactor) +
                   " z=1.0 width=" + v2p(width) + " height=" + v2p(height) + " depth=" + v2p(slices) +
                   " interpolation=Bicubic average process create" +
                   " title=" + v2p(classifiedRegion));
  }

  // Median-filter probability map
  if ( medianFilter > 0 )
  {
    print("\tMedian filter: " + medianFilter + " " + pixelCalibration[0]);
    medianFilter = Math.round(toPixels * medianFilter);  // radius, from units to pixels
    run("Median...", "radius=" + v2p(medianFilter));  // preserve edges
  }


  // get regions of interest (tissue segmentation)
   run("Grays");  // change LUT to grayscale
   run("Add...", "value=1");  // detect all regions (first region with value of zero)
   updateDisplayRange(NaN, NaN);
   getRoisFromMasks("region", true);

   rois = roiManager("count");
   for ( i = 0; i < rois; ++i )
   {
     renameRegion(i, toString(i + 1) + delimiter + regionClasses[i]);
   }
   roiManager("Deselect");
   roiManager("Show All");

  // save and clear run, free memory
  finalizeRun(filePath, fileName);
  freeMemory();

  // restore user interface and display results
  toggleBatchMode(batchMode, false);
}

/*
 *  Functions
 */

// Function to classify images for more robust segmentation results
function classifyImage(image, target, path)
{
  print("\n*** Classifying " + target + " image ***");

  selectWindow(image);
  setMetadata("Label", image);  // Weka displays label in its window
  normalizePixelValues();  // normalize for stable classification results
  output = runWekaClassifier(image, target, path);
  return output;
}

// Function to finalize a segmentation run by saving all results
function finalizeRun(path, name)
{
  print("\n*** Saving results to files ***");

  label = File.getNameWithoutExtension(name);
  directory = path + File.separator + label;  // create folder for result files
  File.makeDirectory(directory);
  folder = directory + File.separator + regionFolder;  // create subfolder for result files
  File.makeDirectory(folder);
  result = folder + File.separator + targetName;  // generic full result path (without extension)
  zipFile = result + ".zip";
  print("\tWriting: " + zipFile);
  waitForFileDeletion(zipFile);
  roiManager("save", zipFile);
  close("\\Others");  // close all images except for the selected image
  tifFile = result + ".tif";
  print("\tWriting: " + tifFile);
  waitForFileDeletion(tifFile);
  run("Bio-Formats Exporter", "save=[" + tifFile + "] export compression=zlib");
  Ext.close();  // close active Bio-Formats dataset
  txtFile = result + ".txt";
  print("\tWriting: " + txtFile);
  waitForFileDeletion(txtFile);
  printDateTimeStamp();
  saveLogFile(txtFile);
}

// Function to create a projected image from an image stack
function projectStack(image, labels, channels, target)
{
  // The slices used for the projection are identified by matching the slice labels
  // with the user-defined channel list. We're then creating a copy of the user-defined
  // channel or a projection of the list user-defined channels. Before the projection,
  // each slice is normalized by its median value to balance all pixel intensities.
  print("\n*** Projecting " + target + " image ***");
  output = target + "->proj";
  channelsLength = channels.length;
  channelMatches = 0;
  labelsLength = labels.length;
  stackSelection = "";

  stackSelection = getSlicesFromLabels(labels, channels);
  selectionStack = makeSubstack(image, stackSelection, "stack-" + target);
  selectWindow(selectionStack);

  slices = nSlices;
  for (i = 1; i <= slices; ++i)
  {
    setSlice(i);
    normalizePixelValues();  // normalize for balanced projection results
  }

  if ( slices > 1 )
    run("Z Project...", "projection=[Sum Slices]");  // project stack to image
  renameImage("", output);
  close("stack-*");  // close selection stack
  print("\tChannels:");
  Array.print(stackSelection);
  return output;
}

// Function to train and run a Weka classification
function runWekaClassifier(image, target, path)
{
  // The tricky part was to wait until the Trainable Weka Classifier plugin has
  // started or until the computation of the probability maps was completed,
  // since these variable times are highly system-dependent. This problem was
  // solved by frequently checking for the currently selected image name:
  // By default, newly opened or created images get the focus in ImageJ2.
  output = image + "->reg";
  classifier = path + target +".model";
  data = path + target +".arff";
  title =   "Finalize classifier training";
  message = "Draw selections in the Weka window and\n" +
            "assign these to their respective classes:\n" +
            "Train the classifier to update the results.\n \n" +
            "The macro will save the new classifier,\n" +
            "upon confirming this dialog with OK,\n" +
            "but stop execution with Cancel.";

  run("Trainable Weka Segmentation");  // start the Trainable Weka Segmentatio plugin
  waitForWindow("Trainable Weka Segmentation");  // title contains changing version number
  call("trainableSegmentation.Weka_Segmentation.setFeature", "Entropy=true");
  call("trainableSegmentation.Weka_Segmentation.changeClassName", "0", regionClasses[0]);
  call("trainableSegmentation.Weka_Segmentation.changeClassName", "1", regionClasses[1]);
  if ( regionClasses.length > 2 )
    call("trainableSegmentation.Weka_Segmentation.createNewClass", regionClasses[2]);
  if ( !File.exists(classifier) )  // classifier missing in folder
  {
    print("\tNo classifier file in dataset folder...");
    waitForUser(title, message);
    call("trainableSegmentation.Weka_Segmentation.saveClassifier", classifier);
    call("trainableSegmentation.Weka_Segmentation.saveData", data);
  }
  print("\tUsing classifier found in dataset folder...");
  call("trainableSegmentation.Weka_Segmentation.loadClassifier", classifier);
  call("trainableSegmentation.Weka_Segmentation.getResult");  // get result label image
  waitForWindow("Classified image");  // computation time machine-dependent
  close("Trainable Weka Segmentation*");  // title changes with version
  renameImage("", output);

  return output;
}
