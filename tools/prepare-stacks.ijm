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
 *  Extracts, projects, and normalizes channels from muli-channel images -
 *  single channels are saved in their corresponding channel folders and
 *  projected channels can be saved in the 'projectionTarget' folder, if the
 *  'projectionChannels' array contains at least one channel name or slice number.
 *  The output images can then be stacked by loading the images into ImageJ2 and
 *  calling the "Images to Stack" function. Load the resulting stack with the
 *  "Trainable Weka Segmentation" plugin to achieve stable segmentation results.
 *
 *  Dependencies:
 *
 *  CU-CellSeg requires a recent version of CU-MacroLibrary to be installed:
 *  https://github.com/christianrickert/CU-MacroLibrary/
 *
 *  Version:
 *
 *  v1.00 (2022-01-28)
 */

print("\\Clear");

run("Bio-Formats Macro Extensions");

projectionTarget = "nuclei";
projectionChannels = newArray("dapi", "dsdna");
suffixes = newArray(".tif", ".tiff");
files = getFilesInFolder("Select the first TIFF of your dataset", suffixes);
processFolder(files);

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
  toggleBatchMode(true, false);

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

  if ( projectionChannels.length > 0 )
  {
    // create projections from multiple channels
    projectedChannels = projectStack(fileTitle, fileSlices, projectionChannels, projectionTarget);
    selectWindow(projectedChannels);
    normalizePixelValues();  // normalize for stable classification results
    finalizeRun(filePath, fileName, projectionTarget);
  }

  // create images from individual channels
  fileSlicesLength = fileSlices.length;
  for ( i = 0; i < fileSlicesLength; ++i )
  {
    extractedChannel = projectStack(fileTitle, fileSlices, newArray(fileSlices[i]), fileSlices[i]);
    selectWindow(extractedChannel);
    finalizeRun(filePath, fileName, fileSlices[i]);
  }

  // close active Bio-Formats dataset
  Ext.close();
  close("*");

  toggleBatchMode(true, false);
}

// Function to finalize a segmentation run by saving all results
function finalizeRun(path, name, channel)
{
  label = File.getNameWithoutExtension(name);
  directory = path + File.separator + channel;  // create subfolder for result files
  File.makeDirectory(directory);
  result = directory + File.separator + label;  // generic full result path (without extension)
  tifFile = result + ".tif";
  print("\tWriting: " + tifFile);
  waitForFileDeletion(tifFile);
  run("Bio-Formats Exporter", "save=[" + tifFile + "] export compression=zlib");
  close();
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
