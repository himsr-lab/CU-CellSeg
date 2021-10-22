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
 *  Creates the probability map in the "Trainable Weka Segmentation" plugin and
 *  and multiplies it with the 'donorFactor'. In a second step, the currently
 *  selected channel of the last image in the "Window" list is corrected by the
 *  'receptorOffset' and then multiplied with the scaled probablity map.
 *  The result is a background-corrected view of the currently selected channel
 *  of the last-opened multi-channel image.
 *
 *  Dependencies:
 *
 *  CU-CellSeg requires a recent version of CU-MacroLibrary to be installed:
 *  https://github.com/christianrickert/CU-MacroLibrary/
 *
 *  Version:
 *
 *  v1.00 (2021-10-22)
 */

donorFactor = 1.0;
receptorOffset = 0.0;
userThresholds = newArray(-1e30, 1e30);  // default values

setBatchMode(true);
close("Background-corrected");
focusWindow("Trainable Weka Segmentation");
call("trainableSegmentation.Weka_Segmentation.getProbability");
waitForWindow("Probability maps");
selectWindow("Probability maps");
run("Duplicate...", "title=Probability duplicate channels=1");
setBatchMode("hide");
if ( userThresholds[0] != -1e30 || userThresholds[1] != 1e30 )
{
  setThreshold(userThresholds[0], userThresholds[1]);
  setOption("BlackBackground", true);  // don't invert LUT
  run("Convert to Mask", "method=Default background=Dark black");
  rescalePixelValues(NaN, NaN, 0, 1);
}
run("Multiply...", "value=" + v2p(donorFactor));
selectWindow(images[imagesLength - 1]);  // select last image
slice = getSliceNumber();
label = getInfo("slice.label");
run("Duplicate...", "title=Background-corrected channels=" + v2p(slice));
run("Add...", "value=" + v2p(receptorOffset) + " slice");
imageCalculator("Multiply 32-bit", "Background-corrected", "Probability");
close("Probability*");
setBatchMode(false);
