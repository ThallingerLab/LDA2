/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.swing;

import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;

/**
 * class holding the table data for Display Results table
 * @author Juergen Hartler
 *
 */
public class LipidomicsTableModel extends DefaultTableModel implements TableModel 
{

  private static final long serialVersionUID = 4607090326610106898L;
  
  /** name of analyte is in first column */
  public final static int COLUMN_NAME = 0;
  /** area of analyte is in second column */
  public final static int COLUMN_AREA = 1;

  /** for displaying name string - key: rowIndex; value: name to display*/
  private Hashtable<Integer,String> rowToName_;
  /** for getting the area from MSn object - key: rowIndex; value name to use to get the area*/
  private Hashtable<Integer,String> areaLookup_;
  /** returns the LipidParameterSet by rowIndex - key: rowIndex; value: the LipidParameterSet */
  private Hashtable<Integer,LipidParameterSet> rowToParam_;
  /** sorting changes the order - this returns the original index of the element - key: rowIndex; value: original index*/
  private Hashtable<Integer,Integer> rowToOriginal_;
  /** to avoid displaying a lipid multiple times with different levels of identification (with and without omega double bond position) */
  
  private Hashtable<Integer,Boolean> rowToOmegaEvidenceConflict_;
  
  private Set<String> uniqueIdentifiers_;
  
  private Hashtable<Integer,String> rowToMolecularSpeciesHumanReadable_;
  /** should the data be displayed in MSn style */
  private boolean showMSn_;
  /** should assigned omega double bond positions be displayed */
  private boolean showOmega_;
  /** display the modification of the found hit */
  private boolean showMod_;
  /** do these LipidParameterSets contain any OH information*/
  private boolean hasOh_;
  
  /**
   * constructor providing all available information
   * @param params the params in sorted order
   * @param paramsOriginal the params in their original order (as in the Excel file)
   * @param showMSn should the data be displayed in MSn style
   * @param showMod display the modification of the found hit
   * @param boolean1 
   */
  public LipidomicsTableModel(Vector<LipidParameterSet> params, Vector<LipidParameterSet> paramsOriginal, boolean showMSn, boolean showOmega, boolean showMod){
    showMSn_ = showMSn;
    showOmega_ = showOmega;
    showMod_ = showMod;
    hasOh_ = false;
    rowToName_ = new Hashtable<Integer,String>();
    rowToParam_ = new Hashtable<Integer,LipidParameterSet>();
    rowToOriginal_ = new Hashtable<Integer,Integer>();
    areaLookup_ = new Hashtable<Integer,String>();
    rowToOmegaEvidenceConflict_ = new Hashtable<Integer,Boolean>();
    uniqueIdentifiers_ = new HashSet<String>();
    rowToMolecularSpeciesHumanReadable_ = new Hashtable<Integer,String>();
    int rowCount = 0;
    for (LipidParameterSet param : params){
      if (param.getOhNumber()>LipidomicsConstants.EXCEL_NO_OH_INFO) {
        hasOh_ = true;
      }
      
      boolean hasRequestedOmegaData = hasRequestedOmegaData(showOmega, param);
      boolean hasRequestedMSnData = hasRequestedMSnData(showMSn, param);
      if (hasRequestedOmegaData || hasRequestedMSnData) {
        if (hasRequestedOmegaData) {
          rowCount = showRequestedOmegaData(param, paramsOriginal, rowCount);}
        if (hasRequestedMSnData) {
          rowCount = showRequestedMSnData(param, paramsOriginal, rowCount);}
      } else {
        rowCount = showRequestedParams(param, paramsOriginal, rowCount);}
    }
  }
  
  private boolean hasRequestedMSnData(boolean showMSn, LipidParameterSet param) {
    if (showMSn && param instanceof LipidomicsMSnSet && (((LipidomicsMSnSet)param).getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED)) {
      return true;
    }
    return false;
  }
  
  private boolean hasRequestedOmegaData(boolean showOmega, LipidParameterSet param) {
    if (showOmega && param.hasOmegaInformation()) {
      Vector<DoubleBondPositionVO> assignedLabeledHits = StaticUtils.getAssignedDoubleBondPositions(param.getOmegaInformation());
      if (!assignedLabeledHits.isEmpty()) {
        return true;
      }
    }
    return false;
  }
  
  private int showRequestedOmegaData(LipidParameterSet param, Vector<LipidParameterSet> paramsOriginal, int rowCount) {
    Vector<DoubleBondPositionVO> omegaInfo = param.getOmegaInformation();
    Set<String> molecularSpeciesSet = StaticUtils.getMolecularSpeciesSet(omegaInfo);
    Iterator<String> it = molecularSpeciesSet.iterator();
    while(it.hasNext()) {
      String molecularSpecies = it.next();
      Vector<DoubleBondPositionVO> labeledChainCombiVOsWithUnlabeledAssignment = StaticUtils.getDoubleBondAssignmentsOfMolecularSpecies(omegaInfo, molecularSpecies);
      Vector<DoubleBondPositionVO> assignedHits = StaticUtils.getAssignedDoubleBondPositions(labeledChainCombiVOsWithUnlabeledAssignment);
      if (!assignedHits.isEmpty()) {
        for (DoubleBondPositionVO labeledChainCombiVO : assignedHits) {
          String doubleBondPositionsHumanReadable = labeledChainCombiVO.getDoubleBondPositionsHumanReadable();
          rowToMolecularSpeciesHumanReadable_.put(rowCount, molecularSpecies);
          uniqueIdentifiers_.add(molecularSpecies + "_"+param.getRt());
          areaLookup_.put(rowCount,molecularSpecies);
          if (param.getRt()!=null && param.getRt().length()>0)
            doubleBondPositionsHumanReadable += "_"+param.getRt();
          String paramName = getLipidParamsDisplayString(param,doubleBondPositionsHumanReadable);
          rowToName_.put(rowCount, paramName);
          rowToParam_.put(rowCount, param); 
          for (int i=0; i!=paramsOriginal.size();i++){
            if (getLipidParamsDisplayString(param,param.getNameString()).equalsIgnoreCase(getLipidParamsDisplayString(paramsOriginal.get(i),paramsOriginal.get(i).getNameString()))){
              rowToOriginal_.put(rowCount,i);
              break;
            }
          }
          rowCount++;
        }
      } else {
        if (!showMSn_) {
          // if some molecular species have an ambiguous assignment
          areaLookup_.put(rowCount, molecularSpecies);
          rowToMolecularSpeciesHumanReadable_.put(rowCount, molecularSpecies);
          rowCount = showRequestedParams(param, paramsOriginal, rowCount);
        }
      }
    }
    return rowCount;
  }
  
  private int showRequestedMSnData(LipidParameterSet param, Vector<LipidParameterSet> paramsOriginal, int rowCount) {
    LipidomicsMSnSet msnSet = (LipidomicsMSnSet)param;
    //TODO: ensure using only the new nomenclature here makes sense.
    Vector<String> detected = msnSet.getMSnIdentificationNamesWithSNPositions();
    for (String name : detected)
    {
    	if (uniqueIdentifiers_.contains(name + "_"+param.getRt())) {continue;}
      areaLookup_.put(rowCount, name);
      
      Vector<DoubleBondPositionVO> doubleBondPositionVO = StaticUtils.getDoubleBondAssignmentsOfMolecularSpecies(param.getOmegaInformation(), name);
      if (StaticUtils.getHighAccuracyDoubleBondPositions(doubleBondPositionVO).size() >1) {
        rowToOmegaEvidenceConflict_.put(rowCount, true);
      }
      
      rowToMolecularSpeciesHumanReadable_.put(rowCount, name);
      if (param.getRt()!=null && param.getRt().length()>0)
      	name += "_"+param.getRt();
      String paramName = getLipidParamsDisplayString(param,name);
      rowToName_.put(rowCount, paramName);
      rowToParam_.put(rowCount, param); 
      for (int i=0; i!=paramsOriginal.size();i++){
        if (getLipidParamsDisplayString(param,param.getNameString()).equalsIgnoreCase(getLipidParamsDisplayString(paramsOriginal.get(i),paramsOriginal.get(i).getNameString()))){
          rowToOriginal_.put(rowCount,i);
          break;
        }
      }
      rowCount++;
    }
    return rowCount;
  }
  
  private int showRequestedParams(LipidParameterSet param, Vector<LipidParameterSet> paramsOriginal, int rowCount) {
    Vector<DoubleBondPositionVO> highAccuracyHits = StaticUtils.getHighAccuracyDoubleBondPositions(param.getOmegaInformation());
    if (highAccuracyHits.size() >1) {
      rowToOmegaEvidenceConflict_.put(rowCount, true);
    }
    String paramName = getLipidParamsDisplayString(param,param.getNameString());
    rowToName_.put(rowCount, paramName);
    rowToParam_.put(rowCount, param); 
    for (int i=0; i!=paramsOriginal.size();i++){
      if (paramName.equalsIgnoreCase(getLipidParamsDisplayString(paramsOriginal.get(i),paramsOriginal.get(i).getNameString()))){
        rowToOriginal_.put(rowCount,i);
        break;
      }
    }
    rowCount++;
    return rowCount;
  }

  public Object getValueAt(int rowIndex, int columnIndex)
  {
    String paramName = rowToName_.get(rowIndex);
    if (columnIndex==COLUMN_NAME) return paramName;
    else if (columnIndex==COLUMN_AREA){
      LipidParameterSet param = rowToParam_.get(rowIndex);
      if ((hasRequestedMSnData(showMSn_, param) || hasRequestedOmegaData(showOmega_, param)) && param instanceof LipidomicsMSnSet){
        LipidomicsMSnSet msnSet = (LipidomicsMSnSet)param;
        return String.valueOf((float)msnSet.getRelativeIntensity(areaLookup_.get(rowIndex))*param.getArea());
      }else{
        return String.valueOf(param.getArea());
      }
    }
    return null;
  }
  
//  public String getSumLipidNameAt(int rowIndex){
//    LipidParameterSet param = rowToParam_.get(rowIndex);
//    return getLipidParamsDisplayString(param,param.getNameString());
//  }
  
  public String getDisplayedNameAt(int rowIndex)
  {
  	return rowToName_.get(rowIndex);
  }
  
  /**
   * returns a String representation in the way the data shall be displayed in the table
   * @param params the parameter set to be displayed
   * @param nameString the base string to be displayed
   * @return a String representation in the way the data shall be displayed in the table
   */
  private String getLipidParamsDisplayString(LipidParameterSet params, String nameString){
    String displayString = new String(nameString);
    if (showMod_) displayString+="_"+params.getModificationName();
    return displayString;
  }
  
  public int getRowCount(){
    if (rowToName_==null) return 0;
    return rowToName_.size();
  }

  public int getColumnCount(){
    return 2;
  }
  
  public String getColumnName(int columnIndex){
    if (columnIndex==COLUMN_NAME) return "Name";
    else if (columnIndex==COLUMN_AREA) return "Area";
    return null;
  }
  
  /**
   * is there any MSn information available about this identification
   * @param rowIndex row to display
   * @return true if the element contains MSn evidence
   */
  public boolean hasMS2Evidence(int rowIndex){
    LipidParameterSet param = rowToParam_.get(rowIndex);
    if (StaticUtils.checkMS2Evidence(param)==StaticUtils.MS2_FULL) return true;
    return false;
  }
  
  public boolean hasOmegaEvidenceConflict(int rowIndex) {
    if (showOmega_ && rowToOmegaEvidenceConflict_.containsKey(rowIndex)) {
      return rowToOmegaEvidenceConflict_.get(rowIndex);
    }
    return false;
  }
  
  /**
   * sorting changes the order of the identifications
   * @return returns the original index of the element - key: rowIndex; value: original index
   */
  public Hashtable<Integer,Integer> getPositionToOriginal(){
    return rowToOriginal_;
  }

  /**
   * returns the MSn identification name of certain row in the table, if present
   * @param rowIndex row index in the table
   * @return MSn identification name of certain row in the table, if present
   */
  public String getMSnIdentificationName(int rowIndex){
    if (areaLookup_!=null && areaLookup_.containsKey(rowIndex)) return areaLookup_.get(rowIndex);
    return null;
  }
  
  /** Returns the result names */
  public Hashtable<Integer,String> getRowToName()
  {
    return this.rowToName_;
  }
  
  public Hashtable<Integer,String> getRowToMolecularSpeciesHumanReadable(){
    return this.rowToMolecularSpeciesHumanReadable_;
  }
  
  /**
   * 
   * @param rowIndex row to display
   * @return if there was an MS1 split according to MSn intensities
   */
  public boolean isSplitInstance(int rowIndex){
    LipidParameterSet param = rowToParam_.get(rowIndex);
    if (StaticUtils.checkMS2Evidence(param)==StaticUtils.SPLIT) return true;
    return false;
  }

  /**
   * 
   * @param rowIndex row to display
   * @return if there was an percental split according to MSn intensities
   */
  public boolean isPercentalSplitInstance(int rowIndex){
    LipidParameterSet param = rowToParam_.get(rowIndex);
    if (StaticUtils.checkMS2Evidence(param)==StaticUtils.PERCENTAL_SPLIT) return true;
    return false;
  }

  /**
   * 
   * @return is show MSn selected
   */
  public boolean isShowMSn()
  {
    return showMSn_;
  }
  
  /**
   * 
   * @return is show omega selected
   */
  public boolean isShowOmega(){
    return showOmega_;
  }
  
  /**
   * 
   * @return true when OH information is present in this table
   */
  public boolean hasOhInfo() {
    return hasOh_;
  }
  
  
}
