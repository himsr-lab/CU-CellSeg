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
 *
 *  Removes a background (donor) contribution in multi-channel images from
 *  sample (receptor) channels by scaling each sample channel pixel with the
 *  inverse of its background probability. The probability map is generated
 *  from a corresponding "Trainable Weka Segmentation" plugin model.
 *  The background-corrected receptor value R*[x,y] for each given pixel R[x,y]
 *  can be calculated using the following formula:
 *
 *  R*[x,y] = (1 - (d * p(D[x,y])) * (R[x,y] - o)			(1)
 *
 *  The donor factor 'd' is a normalization factor that can be used to increase
 *  or decrease the background correction scaling. The scaling is based on the
 *  background probability 'p', derived from the donor channel 'D[x,y]'.
 *  In addition to the scaling, the receptor values can be corrected by applying
 *  the receptor offset 'o', which is applied before the scaling.
 *  As an alterntive to the background correction scaling, the donor probability
 *  map can be binarized by thresholding:
 *  In this case, all receptor values that fall within the background mask will
 *  be set to zero and all receptor values that fall within the sample mask
 *  will only be corrected by the offset 'o'.
 *
 *  Dependencies:
 *
 *  CU-CellSeg requires a recent version of CU-MacroLibrary to be installed:
 *  https://github.com/christianrickert/CU-MacroLibrary/
 *
 *  Version:
 *
 *  v1.00 (2021-10-26)
 */

print("\\Clear");

run("Bio-Formats Macro Extensions");

batchMode = true;
donorChannels = newArray("gold");
donorFactor = 1.0;
receptorChannels = newArray(0);  // optional, all if not specified
receptorOffset = 0.0;
userThresholds = newArray(false, -1e30, 1e30);  // default values
targetNames = newArray("do");  // class label and file output
suffixes = newArray(".tif", ".tiff");
files = getFilesInFolder("Select the first TIFF of your dataset", suffixes);
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
  print("\\Clear");
  run("Close All");

  // print current file name
  print("\n*** Processing file ***");
  print("\t" + file);

  // read image file
  filePath = File.getDirectory(file);
  fileName = File.getName(file);
  fileSlices = readImage(file);
  fileTitle = getTitle();

  // create donor classification
  projectedDonor = projectStack(fileTitle, fileSlices, donorChannels, targetNames[0]);
  classifiedDonor = classifyImage(projectedDonor, targetNames[0], filePath);

  // projection and pixel classification incompatible with batch mode, safe from here
  toggleBatchMode(batchMode, false);

  // create binary map with theshold values upon request
  setUserThresholds(userThresholds);
  if ( userThresholds[1] != -1e30 || userThresholds[2] != 1e30 )
  {
    setOption("BlackBackground", true);  // don't invert LUT
    run("Convert to Mask", "method=Default background=Dark black");
    rescalePixelValues(NaN, NaN, 0, 1);
  }

  // apply donor factor to donor classification
  run("Multiply...", "value=" + v2p(donorFactor));

  // apply classification or binary map to recipient channels
  toggleBatchMode(batchMode, true);  // make probability map availabe for 'imageCalculator'
  receptorChannelsLength = receptorChannels.length;
  if ( receptorChannelsLength == 0 )  // apply to all channels
    imageCalculator("Multiply 32-bit stack", fileName, classifiedDonor);
  else  // apply to selected channels only
  {
    fileSlicesLength = fileSlices.length;
    for (i = 1; i <= fileSlicesLength; ++i)  // iterate through slices
    {

      for (j = 0; j < receptorChannelsLength; ++j)  // match slice names with channels
      {
        slice = toLowerCase(fileSlices[i - 1]);  // label pattern: "#" or "name (channel/mass)"
        if ( slice == receptorChannels[j] ||
             slice.contains(toLowerCase(receptorChannels[j])  + " ") )  // matching pattern: "name "
        {
          selectWindow(fileName);  // needs to be here: timing issue with TWS image release
          setSlice(i);
          run("Add...", "value=" + v2p(receptorOffset) + " slice");
          imageCalculator("Multiply 32-bit", fileName, classifiedDonor);
        }
      }

    }

  }
  toggleBatchMode(batchMode, true);
  setSlice(nSlices());  // move channel slider to the right

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
  directory = path + File.separator + "no-bg";  // create subfolder for result files
  File.makeDirectory(directory);
  result = directory + File.separator + label;  // generic full result path (without extension)
  selectWindow(name);
  close("\\Others");  // close all images except for the selected image
  tifFile = result + ".tif";
  print("\tWriting: " + tifFile);
  waitForFileDeletion(tifFile);
  run("Bio-Formats Exporter", "save=[" + tifFile + "] export compression=LZW");
  Ext.close();  // close active Bio-Formats dataset
  printDateTimeStamp();
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
  stackSelection = "";

  selectWindow(image);

  for (i = 1; i <= slicesLength; ++i)  // iterate through slices
  {

    for (j = 0; j < channelsLength; ++j)  // match slice names with channels
    {
      slice = toLowerCase(slices[i - 1]);  // label pattern: "#" or "name (channel/mass)"
      if ( slice == channels[j] ||
          slice.contains(toLowerCase(channels[j])  + " ") )  // matching pattern: "name "
      {
        if ( stackSelection.length > 0 )  // append slices
          stackSelection = stackSelection + ",";
        stackSelection = stackSelection + toString(i);
        channelMatches += 1;
      }
    }

  }

  if ( channelMatches <= 1 )  // copy slice from stack
  {
    if ( channelMatches == 1 )  // select matching channel
      setSlice(stackSelection);
    run("Duplicate...", "title=slice-" + target);
  }
  else if ( channelMatches >= 2 ) // stack matching channels
  {
    run("Make Substack...", "channels=" + v2p(stackSelection));
    renameImage("", "stack-" + target);

    for (i = 0; i < channelMatches; ++i)
    {
      setSlice(i + 1);
      normalizePixelValues();  // normalize for balanced projection results
    }

    run("Z Project...", "projection=[Sum Slices]");  // project stack to image
    close("stack-*");  // close projection stack
  }
  renameImage("", output);
  print("\tChannels: \"" + stackSelection + "\" (" + target + ")");
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
  call("trainableSegmentation.Weka_Segmentation.changeClassName", "0", "sample (receptor)");
  call("trainableSegmentation.Weka_Segmentation.changeClassName", "1", "background (donor)");
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

  while ( nSlices() > 1 )  // use only first probability map
  {
    setSlice(nSlices());
    run("Delete Slice");
  }

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
