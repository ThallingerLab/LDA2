/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2023 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.target;

import java.util.Hashtable;

import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * Class containing information about an isotope label
 * 
 * @author Leonida M. Lamp
 *
 */
public class IsotopeLabelVO
{
	/** the label elements including its amount for each label indicator*/
  protected Hashtable<String,Integer> labelElements_;
  /** the omega double bond position this label stands for*/  
  protected int omegaPosition_;
  /** the prefix that is used to indicate an isotope label*/
  protected String labelId_;
  
  /**
   * 
   * @param labelElements
   * @param omegaPosition
   * @param labelId
   */
  public IsotopeLabelVO(
  		Hashtable<String,Integer> labelElements, int omegaPosition, String labelId)
  {
  	this.labelElements_ = labelElements;
  	this.omegaPosition_ = omegaPosition;
  	this.labelId_ = labelId;
  }
  

  /**
   * getter for prefix that is used to indicate an isotope label
   * @return prefix that is used to indicate an isotope label
   */
  public String getLabelId()
  {
    return labelId_;
  }

  
  /**
   * setter for prefix that is used to indicate an isotope label
   * @param labelId prefix that is used to indicate an isotope label
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
  public void setOmegaPosition(int omegaPosition)
  {
    this.omegaPosition_ = omegaPosition;
  }

  
  /**
   * compares whether two value objects of this class are the same
   * @param other the other value object to compare to
   * @return true if both contain the same label and the same label elements
   */
  public boolean isEqual(IsotopeLabelVO other) {
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
    return labelId_+": n-"+String.valueOf(omegaPosition_)+"; "+StaticUtils.getFormulaInHillNotation_PlusFirst(this.labelElements_, true);
  }
  
}
