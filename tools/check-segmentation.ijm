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
 *  Creates the probability map in the "Trainable Weka Segmentation" plugin,
 *  thresholds the resulting image with the user-specified limits and finally
 *  segments the thresholded image with the "Analyze Particles" function.
 *  The resulting regions of interest are displayed as red cell outlines on
 *  top of the initial image slice and can be analyzed in the "ROI Manager".
 *  The algorithm and settings used are identical to the ones used in CU-CellSeg.
 *
 *  Dependencies:
 *
 *  CU-CellSeg requires a recent version of CU-MacroLibrary to be installed:
 *  https://github.com/christianrickert/CU-MacroLibrary/
 *
 *  Version:
 *
 *  2021-10-08
 */

nucleiFilling = false;  // fixes nuclei with darker centers
userThreshold = newArray(0.75, 1e30);  // upper and lower limit

run("ROI Manager...");  // start before batch mode
setBatchMode(true);
Stack.getPosition(channel, slice, frame);
call("trainableSegmentation.Weka_Segmentation.getProbability");
waitForWindow("Probability maps");
selectWindow("Probability maps");
run("Duplicate...", "title=Probability duplicate channels=1");
setSlice(slice);
setBatchMode("hide");
setThreshold(userThreshold[0], userThreshold[1]);
setOption("BlackBackground", true);
run("Convert to Mask", "method=Otsu background=Dark black");
if ( nucleiFilling )
	run("Analyze Particles...", "size=10-Infinity pixel show=[Masks] include in_situ slice");
run("Watershed", "slice");
run("Analyze Particles...", "size=10-Infinity pixel show=Masks clear add in_situ slice");
close("Probability*");
setBatchMode(false);
resetMinAndMax();
roiManager("show all without labels")
