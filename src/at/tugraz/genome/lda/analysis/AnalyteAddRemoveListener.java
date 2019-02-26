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

import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.vos.AddAnalyteVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public interface AnalyteAddRemoveListener
{
  /**
   * adds an analyte to the currently open Excel file
   * @param position the position in class list where the analyte should be added
   * @param analyteDescrVO the required information about the new analyte
   */
  public void addAnalyte(int position, AddAnalyteVO analyteDescrVO);
  public void removeAnalyte(int[] indices);
  public LipidParameterSet getAnalyteInTableAtPosition(int position);
  public void showMs2(int position);
  public void listSelectionChanged(int leadIndex);
  public void changeListSorting(int sortType);
  public void newRule(int position);
  /**
   * re-evaluates the assigned MSn identification
   * according to the currently set rules
   * @param position the position in the JTable
   */
  public void recalculateMSn (int position);
  
  /**
   * pops up a menu for Rt editing
   * @param position the position in the JTable
   */
  public void editRt(int position);
  
  /**
   * starts a new MS1 viewer
   * @param position the position in the displayed table showing a list of species
   */
  public void initANewViewer(int position);
  
  /**
   * 
   * @return true when MS2 view is displayed, false when MS1 view is displayed
   */
  public boolean isMS2Showing();

}
