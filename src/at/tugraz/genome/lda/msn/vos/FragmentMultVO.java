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

package at.tugraz.genome.lda.msn.vos;

import java.util.Objects;

/**
 * Class containing a part of an equation, namely, the fragment name
 * the multiplication factor, if the value shall be added or subtracted, and, if appropriate, the affected position
 * 
 * @author Juergen Hartler
 *
 */
public class FragmentMultVO
{
  /** the name of the fragment */
  private String fragmentName_;
  /** the type of the fragment none, acyl, alkyl or lcb*/
  private short type_;
  /** the multiplication factor */
  private String multFactor_;
  /** shall the value be added or subtracted*/
  private boolean positive_;
  /** if appropriate, the affected position */
  private int position_;
  
  
  /**
   * Constructor requiring all information for the VO; if there is no position appropriate - use 0 or -1
   * @param fragmentName the name of the fragment
   * @param type the chain type (acyl/alkyl/alkenyl/lcb) as defined in LipidomicsConstants
   * @param multFactor multiplication factor
   * @param positive shall the value be added or subtracted (true for adding)
   * @param position if appropriate, the affected position - use 0 or -1
   */
  public FragmentMultVO(String fragmentName, short type, String multFactor,
      Boolean positive, int position)
  {
    super();
    this.fragmentName_ = fragmentName;
    this.type_ = type;
    this.multFactor_ = multFactor;
    this.positive_ = positive;
    this.position_ = position;
  }


  /**
   * 
   * @return the name of the fragment
   */
  public String getFragmentName()
  {
    return fragmentName_;
  }
  
  /**
   * 
   * @return the type of the fragment
   */
  public short getFragmentType() {
    return type_;
  }

  /**
   * 
   * @return multiplication factor value
   */
  public double getMultFactor()
  {
    return Double.parseDouble(multFactor_);
  }

  /**
   * 
   * @return multiplication factor as String
   */
  public String getMultFactorAsString()
  {
    return multFactor_;
  }

  /**
   * 
   * @return shall the value be added or subtracted (true for adding)
   */
  public boolean isPositive()
  {
    return positive_;
  }


  /**
   * 
   * @return the affected position (use if appropriate)
   */
  public int getPosition()
  {
    return position_;
  }
  
  public String toString(){
    String value = "";
    if (positive_) value += "+";
    else value += "-";
    value += fragmentName_;
    if (position_>0)
      value += "["+String.valueOf(position_)+"]";
    value += "*"+multFactor_;
    return value;
  }


  @Override
  public boolean equals(Object obj)
  {
    if (this == obj)
      return true;
    if (obj == null)
      return false;
    if (getClass() != obj.getClass())
      return false;
    FragmentMultVO other = (FragmentMultVO) obj;
    return Objects.equals(fragmentName_, other.fragmentName_)
        && Objects.equals(multFactor_, other.multFactor_)
        && position_ == other.position_ && positive_ == other.positive_
        && type_ == other.type_;
  }
  
}
