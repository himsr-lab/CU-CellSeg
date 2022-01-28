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
 *  The corresponding probability maps are then thresholded by the user and the
 *  binary segmentation maps transfered into "ROI Manager" entries.
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
regionChannels = newArray("keratin");
regionClasses = newArray("Tumor (positive)", "Stroma (negative)", "Glass (neutral)");
regionFolder = "regions";
userThresholds = newArray(false, 0.5, 1e30);  // default values
targetName = "tu_st_gl";  // class label and file output
suffixes = newArray(".tif", ".tiff");
files = getFilesInFolder("Select the first TIFF of your dataset", suffixes);
versionString = "v1.00 (2022-01-28)";
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

  // create donor classification
  projectedRegion = projectStack(fileTitle, fileSlices, regionChannels, targetName);
  classifiedRegion = classifyImage(projectedRegion, targetName, filePath);

  // projection and pixel classification incompatible with batch mode, safe from here
  toggleBatchMode(batchMode, false);

  // create binary map with theshold values upon request
  setUserThresholds(userThresholds);
  if ( userThresholds[1] != -1e30 || userThresholds[2] != 1e30 )
  {
    setOption("BlackBackground", true);  // don't invert LUT
    run("Convert to Mask", "method=Default background=Dark black");

    for (i = 0; i < nSlices; ++i)  // thresholding runs on stack, rescaling not
    {
      setSlice(i + 1);
      rescalePixelValues(NaN, NaN, 0, 1);
    }

  }

  // get regions of interest (tissue segmentation)
  regionClassesLength = regionClasses.length;
  for (i = 0; i < regionClassesLength; ++i)
  {
    setSlice(i + 1);
    getRoisFromMasks(regionClasses[i], true);
  }

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

// Function reads user-defined threshold values
function getUserThresholds(thresholds)
{
  title =   "Finalize thresholds limits";
  message = "Set the limits in the Threshold window and\n" +
            "cover the background with the blue mask:\n" +
            "Leave the white sample regions unmasked.\n \n" +
            "The macro will apply the new thresholds,\n" +
            "upon confirming this dialog with OK,\n" +
            "but stop execution with Cancel.";

  toggleBatchMode(batchMode, true);  // stay in batch mode, but show current image
  run("Threshold...");
  call("ij.plugin.frame.ThresholdAdjuster.setMethod", "Default");  // preset Window defaults
  call("ij.plugin.frame.ThresholdAdjuster.setMode", "Over/Under");
  waitForUser(title, message);
  getThreshold(thresholds[1], thresholds[2]);
  run("Close");
  toggleBatchMode(batchMode, true);  // hide image and (re-)enter batch mode
}

// Function to create a projected image from an image stack
function projectStack(image, slices, channels, target)
{
  // The slices used for the projection are identified by matching the slice labels
  // with the user-defined channel list. We're then creating a copy of the user-defined
  // channel or a projection of the list user-defined channels. Before the projection,
  // each slice is normalized by its median value to balance all pixel intensities.
  print("\n*** Projecting " + target + " image ***");
  output = target + "->proj";
  channelsLength = channels.length;
  channelMatches = 0;
  slicesLength = slices.length;
  stackSelection = newArray(0);

  selectWindow(image);

  // find matching channels with slice labels
  for (i = 1; i <= slicesLength; ++i)  // iterate through slices
  {

    for (j = 0; j < channelsLength; ++j)  // match slice names with channels
    {
      slice = toLowerCase(slices[i - 1]);  // label pattern: "#" or "name (channel/mass)"
      if ( slice == channels[j] ||
          slice.contains(toLowerCase(channels[j])  + " ") )  // matching pattern: "name "
      {
        stackSelection = Array.concat(stackSelection, i);
        channelMatches += 1;
      }
    }

  }

  // extract matching slice labels for projection
  run("Duplicate...", "title=" + v2p("stack-" + target) + " duplicate");

  for (k = slicesLength; k > 0; --k)  // reverse iterate through slices
  {
    setSlice(k);
    if ( isInArray(stackSelection, k) )  // keep
      normalizePixelValues();  // normalize for balanced projection results
    else if ( nSlices > 1 )  // remove
      run("Delete Slice");
  }

  if ( nSlices > 1 )
  {
    run("Z Project...", "projection=[Sum Slices]");  // project stack to image
    close("stack-*");  // close projection stack
  }
  renameImage("", output);
  print("\tTarget: \"" + target + "\", Slices:");
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
  output = image + "->prob";
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
  call("trainableSegmentation.Weka_Segmentation.getProbability");
  waitForWindow("Probability maps");  // computation time machine-dependent
  close("Trainable Weka Segmentation*");  // title changes with version
  renameImage("", output);

//  while ( nSlices() > 1 )  // use only first probability map
//  {
//    setSlice(nSlices());
//    run("Delete Slice");
//  }

  return output;
}

// Function reads user-defined threshold values
function setUserThresholds(thresholds)
{
  setThreshold(thresholds[1], thresholds[2]);  // preset default values
  if ( thresholds[0] == false )  // check for custom values
    getUserThresholds(thresholds);
  setThreshold(thresholds[1], thresholds[2]);
  thresholds[0] = true;
  print("\tThresholds: " + thresholds[1] + " (lower), " + thresholds[2] + " (upper)");
}
