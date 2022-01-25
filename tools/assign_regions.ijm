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
 *  First loads all regions of interest with tissue segments from the regions folder and
 *  then loads all regions of interest with cell segmentation data from hte cells folder.
 *  Using the Roi Manager, the macro calculates the relative overlap (percentage) between
 *  the tissue segments and the cellular compartments. The results are stored in a table
 *  and finally exported to a comma-separated value file in the regions folder.
 *  Please organize your data in the following file system structure:
 *
 *  ./image.tiff
 *    ./image/
 *      ./cells/
 *        *.roi
 *        *.zip
 *      ./regions/
 *        *.roi
 *        *.zip
 *
 *  Dependencies:
 *
 *  CU-CellSeg requires a recent version of CU-MacroLibrary to be installed:
 *  https://github.com/christianrickert/CU-MacroLibrary/
 */

print("\\Clear");
run("ROI Manager...");  // start before batch mode
batchMode = true;
cellFolder = "cells";
fileSuffixes = newArray(".tif", ".tiff");
files = getFilesInFolder("Select the first TIFF of your dataset", fileSuffixes);
regionFolder = "regions";
roiSuffixes = newArray(".roi", ".zip");
versionString = "v1.00 (2021-01-25)";
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
  // speed up processing
  toggleBatchMode(batchMode, false);

  // prepare next run
  initializeRun(versionString);

  // print current file name
  print("\n*** Processing file ***");
  print("\t" + file);

  // read image file
  filePath = File.getDirectory(file);
  fileName = File.getName(file);
  fileSlices = readImage(file);
  fileLabel = File.getNameWithoutExtension(fileName);

  // load regions of interest (tissue segments)
  print("\n*** Processing regions ***");
  regionPath = filePath + File.separator + fileLabel + File.separator + regionFolder;
  print("\t" + regionPath);
  objects = getFileList(regionPath);  // files and folders
  objectsLength = objects.length;
  for (i = 0; i < objectsLength; ++i)  // tissue segments in folder
  {
    objectPath = regionPath + File.separator + objects[i];
    if( endsWithEither(objectPath, roiSuffixes) )
      roiManager("open", objectPath);
  }

  // store regions of interest (tissue segments)
  regionIndices = newArray();
  regionNames = newArray();

  rois = roiManager("count");
  for ( i = 0; i < rois; ++i )
  {
    regionName = getRegionName(i);
    regionIndices = Array.concat(regionIndices, i);
    regionNames = Array.concat(regionNames, regionName);
  }
  regionNamesLength = regionNames.length;

  // load regions of interest (cell segments)
  cellPath = filePath + File.separator + fileLabel + File.separator + cellFolder;
  print("\t" + cellPath);
  objects = getFileList(cellPath);  // files and folders
  objectsLength = objects.length;
  for (i = 0; i < objectsLength; ++i)  // cell segments in folder
  {
    objectPath = cellPath + File.separator + objects[i];
    if( endsWithEither(objectPath, roiSuffixes) )
      roiManager("open", objectPath);
  }

  // set up results table
  tableName = regionFolder;
  Table.create(tableName);  // creates and resets a table
  regionOverlap = -1;  // value not determined yet
  roiHeader = "Label";  // consistent with Results table

  // write matches of rois (cell segments) and regions (tissue segments) to table
  rois = roiManager("count");
  for ( roi = regionNamesLength; roi < rois; ++roi )
  {
    rowIndex = roi - regionNamesLength;
    roiLabel = fileName + ":" + getRegionName(roi);

    for ( reg = 0; reg < regionNamesLength; ++reg )
    {
      Table.set(roiHeader, rowIndex, roiLabel, tableName);
      if ( isInBounds(reg, roi) )  // fast approximation
        regionOverlap = getRegionOverlap(reg, roi);
      else
        regionOverlap = 0;
      Table.set(regionNames[reg], rowIndex, regionOverlap, tableName);
    }
  }

  // save table and clear run, free memory
  print("\n*** Saving result to file ***");
  csvFile = regionPath + File.separator + regionFolder + ".csv";
  waitForFileDeletion(csvFile);
  Table.save(csvFile);
  clearAllSelections();
  printDateTimeStamp();
  freeMemory();

  // restore user interface and display results
  toggleBatchMode(batchMode, false);
}
