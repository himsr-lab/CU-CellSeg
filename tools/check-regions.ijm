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
 *  Creates the region label image in the "Trainable Weka Segmentation" plugin,
 *  median-filters the result for the currently selected slice and transfers
 *  the corresponding regions to the ROI Manager.
 *  The resulting regions of interest are displayed as red region outlines on
 *  top of the initial image slice and can be analyzed in the "ROI Manager".
 *
 *  Dependencies:
 *
 *  CU-CellSeg requires a recent version of CU-MacroLibrary to be installed:
 *  https://github.com/christianrickert/CU-MacroLibrary/
 *
 *  Version:
 *
 *  v1.00 (2022-02-11)
 */

medianFilter = 10;  // radius [units]
regionClasses = newArray("Tumor (positive)", "Stroma (negative)", "Glass (neutral)");
scalingFactor = 0.25;

run("ROI Manager...");  // start before batch mode
deleteAllRegions();
setBatchMode(true);
focusWindow("Trainable Weka Segmentation");
Stack.getPosition(channel, slice, frame);

// get calibration data from image file
pixelCalibration = getPixelCalibration();
if ( pixelCalibration[0] == "pixels" )  // no conversion required
  toPixels = 1.0;
else  // calculate conversion factor
  toPixels = Math.round(1.0 / 0.5 * (pixelCalibration[1] + pixelCalibration[2]));

// get the result label image
call("trainableSegmentation.Weka_Segmentation.getResult");
waitForWindow("Classified image");  // computation time machine-dependent
setSlice(slice);
run("Duplicate...", "title=Classified");  // selected slice only

// (optional) rescale image file back
if ( scalingFactor != 1.0 )
{
  invertedScalingFactor = 1.0 / scalingFactor;
  getDimensions(width, height, channels, slices, frames);
  run("Scale...", "x=" + v2p(invertedScalingFactor) + " y=" + v2p(invertedScalingFactor) +
                 " z=1.0 width=" + v2p(width) + " height=" + v2p(height) + " depth=" + v2p(slices) +
                 " interpolation=Bicubic average process create" +
                 " title=" + v2p(getTitle() + "->(*1/f)"));
}

// Median-filter probability map
if ( medianFilter > 0 )
{
  medianFilter = Math.round(toPixels * medianFilter);  // radius, from units to pixels
  run("Median...", "radius=" + v2p(medianFilter));  // preserve edges
}

if ( scalingFactor != 1.0 )  // (optional) rescale image file
{
  getDimensions(width, height, channels, slices, frames);
  run("Scale...", "x=" + v2p(scalingFactor) + " y=" + v2p(scalingFactor) +
                 " z=1.0 width=" + v2p(width) + " height=" + v2p(height) + " depth=" + v2p(slices) +
                 " interpolation=Bicubic average process create" +
                 " title=" + v2p(getTitle() + "->(*f)"));
}

// get regions of interest (tissue segmentation)
run("Grays");  // change LUT to grayscale
run("Add...", "value=1");  // detect all regions (first region with value of zero)
updateDisplayRange(NaN, NaN);
getRoisFromMasks("region", true);

// rename regions of interest with class labels
rois = roiManager("count");
for ( i = 0; i < rois; ++i )
{
  renameRegion(i, toString(i + 1) + delimiter + regionClasses[i]);
}
close("Classified*");
setBatchMode(false);
roiManager("show all without labels")
