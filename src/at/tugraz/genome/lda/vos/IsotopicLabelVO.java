/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2021 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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
package at.tugraz.genome.lda.vos;

import java.util.Hashtable;
import java.util.LinkedHashMap;

import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class IsotopicLabelVO
{
  /** the prefix that is used to indicate an isotopic label*/
  protected String labelId_;
  /** the label elements including its amount for each label indicator*/
  protected Hashtable<String,Integer> labelElements_;
  /** the shift in retention time caused by the isotopic label*/
  protected Float rtShift_;
  /** the used prefixes for this label*/
  protected LinkedHashMap<String,Integer> prefixes_;
  /** the omega position this label stands for*/  
  protected int omegaPosition_;
  
  /**
   * constructor for initializing the VO 
   * @param labelId prefix that is used to indicate an isotopic label
   * @param labelElements label elements including its amount for each label indicator
   * @param the prefixes that were used for this isotopic label
   */
  public IsotopicLabelVO(String labelId, Hashtable<String,Integer> labelElements, LinkedHashMap<String,Integer> prefixes) {
    this.labelId_ = labelId;
    this.labelElements_ = labelElements;
    this.prefixes_ = prefixes;
    rtShift_ = null;
    omegaPosition_ = -1;
  }
  
  /**
   * constructor for the values entered in the GUI
   * @param labelId prefix that is used to indicate an isotopic label
   * @param omegaPosition the omega position this label is indicative for
   * @param labelElements label elements including its amount for each label indicator
   * @param rtShift the retention time shift the label is causing
   * @param prefixes the used prefixes (might be empty)
   */
  public IsotopicLabelVO(String labelId, int omegaPosition, Hashtable<String,Integer> labelElements, float rtShift, LinkedHashMap<String,Integer> prefixes) {
    this(labelId,labelElements,prefixes);
    omegaPosition_ = omegaPosition;
    rtShift_ = rtShift;
  }
  

  /**
   * constructor for cloning an existing VO 
   * @param vo the value object to be cloned
   */
  public IsotopicLabelVO(IsotopicLabelVO vo) {
    this.labelId_ = vo.labelId_;
    this.labelElements_ = vo.labelElements_;
    this.prefixes_ = vo.prefixes_;
    this.rtShift_ = vo.getRtShift();
    omegaPosition_ = vo.omegaPosition_;
  }
  

  /**
   * getter for prefix that is used to indicate an isotopic label
   * @return prefix that is used to indicate an isotopic label
   */
  public String getLabelId()
  {
    return labelId_;
  }

  
  /**
   * setter for prefix that is used to indicate an isotopic label
   * @param labelId prefix that is used to indicate an isotopic label
   */
  public void setLabelId(String labelId)
  {
    this.labelId_ = labelId;
  }

  
  /**
   * getter for label elements including its amount for each label indicator
   * @return label elements including its amount for each label indicator
   */
  public Hashtable<String,Integer> getLabelElements()
  {
    return labelElements_;
  }

  
  /**
   * setter for label elements including its amount for each label indicator
   * @param labelElements label elements including its amount for each label indicator
   */
  public void setLabelElements(Hashtable<String,Integer> labelElements)
  {
    this.labelElements_ = labelElements;
  }

  /**
   * getter for the shift in retention time caused by the isotopic label
   * @return the shift in retention time caused by the isotopic label
   */
  public Float getRtShift()
  {
    return rtShift_;
  }

  
  /**
   * setter for the shift in retention time caused by the isotopic label
   * @param rtShift shift in retention time caused by the isotopic label
   */
  public void setRtShift(Float rtShift)
  {
    this.rtShift_ = rtShift;
  }
  
  
  /**
   * getter for the omega position this label stands for
   * @return omega position
   */
  public int getOmegaPosition()
  {
    return omegaPosition_;
  }

  
  /**
   * setter for the omega position this label stands for
   * @param omegaPosition the omega position
   */
  public void setOmegaPosition_(int omegaPosition)
  {
    this.omegaPosition_ = omegaPosition;
  }

  
  /**
   * 
   * @return the used prefixes for this label
   */
  public LinkedHashMap<String,Integer> getPrefixes()
  {
    return prefixes_;
  }

  
  /**
   * compares whether two value objects of this class are the same
   * @param other the other value object to compare to
   * @return true if both contain the same label and the same label elements
   */
  public boolean isEqual(IsotopicLabelVO other) {
    if (!labelId_.contentEquals(other.labelId_))
      return false;
    if (labelElements_.size()!=other.labelElements_.size())
      return false;
    for (String element : this.labelElements_.keySet()) {
      if (!other.labelElements_.containsKey(element))
        return false;
      if (this.labelElements_.get(element).intValue()!=other.labelElements_.get(element).intValue())
        return false;
    }
    return true;
  }
  
  
  public String toString() {
    return labelId_+": n-"+String.valueOf(omegaPosition_)+"; "+StaticUtils.getFormulaInHillNotation_PlusFirst(this.labelElements_, true)+"; RT: "+rtShift_.toString()+"; "+prefixes_.keySet();
  }
  
}
