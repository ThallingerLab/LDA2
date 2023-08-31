/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */ 

package at.tugraz.genome.lda.analysis;

import java.io.File;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.vos.AutoAnalyteAddVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public interface HeatMapClickListener extends SampleLookup
{
  public boolean heatMapClicked(String experimentName, String resultFilePath,  String moleculeName, String rt, boolean showMSn);
  
  /**
   * showing a bar chart when a specific analyte was clicked
   * @param moleculeName the name of the analyte
   * @param groupName the name of the analyte class
   * @param maxIsotopes the highest isotope number allowed for the quantity
   * @param rtGrouped true when the peaks were grouped by a certain retention time
   * @param settingVO value object specifying the type of displayed values
   * @param prefUnit the preferred unit magnifier, e.g p for pico
   * @param unit a description of the physical unit
   * @return true when the execution was successful
   */
  public boolean analyteClicked(String moleculeName,String groupName, int maxIsotopes, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  
  /**
   * showing a bar chart when a specific analyte was clicked
   * @param moleculeName the name of the analyte
   * @param groupName the name of the analyte class
   * @param maxIsotopes the highest isotope number allowed for the quantity
   * @param rtGrouped true when the peaks were grouped by a certain retention time
   * @param settingVO value object specifying the type of displayed values
   * @param prefUnit the preferred unit magnifier, e.g p for pico 
   * @param unit a description of the physical unit
   * @return true when the execution was successful
   */
  public boolean analyteGroupClicked(String moleculeName,String groupName, int maxIsotopes, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  
  /**
   * showing a bar chart when a specific experiment was clicked
   * @param experimentName the name of the experiment in the heat map
   * @param groupName the name of the analyte class
   * @param maxIsotopes the highest isotope number allowed for the quantity
   * @param rtGrouped true when the peaks were grouped by a certain retention time
   * @param settingVO value object specifying the type of displayed values
   * @param prefUnit the preferred unit magnifier, e.g p for pico
   * @param unit a description of the physical unit
   * @return true when the execution was successful
   */
  public boolean experimentClicked(String experimentName,String groupName, int maxIsotopes, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  
  /**
   * showing a bar chart when a specific sample group was clicked
   * @param experimentGroupName the name of the experiment in the heat map
   * @param groupName the name of the analyte class
   * @param maxIsotopes the highest isotope number allowed for the quantity
   * @param rtGrouped true when the peaks were grouped by a certain retention time 
   * @param settingVO value object specifying the type of displayed values
   * @param prefUnit the preferred unit magnifier, e.g p for pico
   * @param unit a description of the physical unit
   * @return true when the execution was successful
   */
  public boolean experimentGroupClicked(String experimentGroupName,String groupName, int maxIsotopes, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  
  /**
   * showing a bar chart when several analytes were selected
   * @param moleculeNames names of selected analytes
   * @param groupName the name of the analyte class
   * @param maxIsotopes the highest isotope number allowed for the quantity
   * @param rtGrouped true when the peaks were grouped by a certain retention time
   * @param settingVO value object specifying the type of displayed values
   * @param prefUnit the preferred unit magnifier, e.g p for pico
   * @param unit a description of the physical unit
   * @return true when the execution was successful
   */
  public boolean combinedAnalyteSelected(Vector<String> moleculeNames, String groupName, int maxIsotopes, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  
  /**
   * showing a bar chart when several analytes were selected for the grouped view
   * @param moleculeNames names of selected analytes
   * @param groupName the name of the analyte class
   * @param maxIsotopes the highest isotope number allowed for the quantity
   * @param rtGrouped true when the peaks were grouped by a certain retention time
   * @param settingVO value object specifying the type of displayed values
   * @param prefUnit the preferred unit magnifier, e.g p for pico
   * @param unit a description of the physical unit
   * @return true when the execution was successful
   */
  public boolean combinedAnalyteGroupSelected(Vector<String> moleculeNames, String groupName, int maxIsotopes, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  
  public void changeISStatus(String groupName, boolean isGrouped, boolean value);
  public void changeESStatus(String groupName, boolean isGrouped, boolean value);
  public void changeIsotopesUsed(String groupName, boolean isGrouped, int value);
  public void eliminateDoublePeaks(String groupName, String analyteName, String absFilePathStartExp, Vector<String> selectedMods, Vector<String> foundUpdateables);
  public void eliminateAnalyteEverywhere(String groupName, Hashtable<String,String> selectedAnalytes, Vector<String> selectedMods, Vector<String> foundUpdateables);
  
  /**
   * for automatically adding analytes by selecting template probes in the heat map
   * @param groupName the analyte class name
   * @param analyteNames a sorted list of analyte names
   * @param absFilePathStartExps a sorted list of experiments corresponding to the analyte names (they must be in the same sequence as the analyte names)
   * @param selectedMods the modifications of the analytes that have to be searched for
   * @param updateableAndAnalyteBefore information about where in the original file the new analyte has to be added; key: analyte name; value: vector with the information for each experiment where the analyte should be searched for
   * @param maxIsotopes the highest number of isotopes to look for each analyte  (they must be in the same sequence as the analyte names)
   * @param exactProbePosition true when the exact positions of the peak should be taken; false when an automated method shall be used
   */
  public void addAnalytesEverywhereAtPosition(String groupName, Vector<String> analyteNames, Vector<String> absFilePathStartExps, Vector<String> selectedMods,
      Hashtable<String,Vector<AutoAnalyteAddVO>> updateableAndAnalyteBefore, Vector<Integer> maxIsotopes, boolean exactProbePosition);
  
  /**
   * exports the contents of the heat map to mzTab
   * @param exportFile the file where mzTab shall be written
   * @param speciesType which species shall be exported - for details see LipidomicsConstants.EXPORT_ANALYTE_TYPE
   * @param exportDoubleBondPositions true when double bond positions shall be exported
   */
  public void exportMzTab(File exportFile, short speciesType, boolean exportDoubleBondPositions);
  public void exportRdb(File exportFile);
  public void exportMaf(File exportFile);
  
  /**
   * callback to show the settings dialog for exporting
   * @param grouped is it the export settings for the grouped heat map
   */
  public void showExportSettingsDialog(boolean grouped);
  
}
