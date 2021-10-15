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
 *  Removes the background (donor) contribution in multi-channel images from
 *  sample (recipient) channels by scaling each sample channel pixel with the
 *  inverse of its background probability. The probability map is generated from
 *  a corresponding "Trainable Weka Segmentation" plugin model.
 *  The background-corrected pixel value R*[x,y] can be calculated using
 *  the following formula:
 *      
 *  R*[x,y] = f / p(D[x,y] * R[x,y]					(1)
 *
 *  where R[x,y] is the initial pixel value at locaction [x,y] in the recipient
 *  channel and p(D[x,y]) is the background probability as calculated for the
 *  corresponding pixel in the donor channel. The factor f is a fixed scaling
 *  factor that can be used to increase or decrease the background correction.
 *  
 *  Dependencies:
 *
 *  CU-CellSeg requires a recent version of CU-MacroLibrary to be installed:
 *  https://github.com/christianrickert/CU-MacroLibrary/
 *
 *  Version:
 *
 *  v1.00 (2021-10-15)
 */

