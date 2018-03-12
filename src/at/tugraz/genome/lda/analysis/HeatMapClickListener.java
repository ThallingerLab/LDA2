/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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
  public boolean heatMapClicked(String experimentName, String resultFilePath,  String moleculeName);
  public boolean analyteClicked(String moleculeName,String groupName, int maxIsotopes, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  public boolean analyteGroupClicked(String moleculeName,String groupName, int maxIsotopes, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  public boolean experimentClicked(String experimentName,String groupName, int maxIsotopes, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  public boolean experimentGroupClicked(String experimentGroupName,String groupName, int maxIsotopes, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  public boolean combinedAnalyteSelected(Vector<String> moleculeNames, String groupName, int maxIsotopes, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  public boolean combinedAnalyteGroupSelected(Vector<String> moleculeNames, String groupName, int maxIsotopes, ResultDisplaySettingsVO settingVO, String prefUnit, String unit);
  
  public void changeISStatus(String groupName, boolean isGrouped, boolean value);
  public void changeESStatus(String groupName, boolean isGrouped, boolean value);
  public void changeIsotopesUsed(String groupName, boolean isGrouped, int value);
  public void eliminateDoublePeaks(String groupName, String analyteName, String absFilePathStartExp, Vector<String> selectedMods, Vector<String> foundUpdateables);
  public void eliminateAnalyteEverywhere(String groupName, Hashtable<String,String> selectedAnalytes, Vector<String> selectedMods, Vector<String> foundUpdateables);
  public void addAnalyteEverywhereAtPosition(String groupName, String analyteName, String absFilePathStartExp, Vector<String> selectedMods, Vector<AutoAnalyteAddVO> updateableAndAnalyteBefore, int maxIsotope, boolean exactProbePosition);
  
  /**
   * exports the contents of the heat map to mzTab
   * @param exportFile the file where mzTab shall be written
   * @param speciesType which species shall be exported - for details see LipidomicsConstants.EXPORT_ANALYTE_TYPE
   */
  public void exportMzTab(File exportFile, short speciesType);
  public void exportRdb(File exportFile);
  public void exportMaf(File exportFile);
  
  /**
   * callback to show the settings dialog for exporting
   * @param grouped is it the export settings for the grouped heat map
   */
  public void showExportSettingsDialog(boolean grouped);
}
