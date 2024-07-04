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

import java.util.Hashtable;
import java.util.Objects;
import java.util.Vector;

import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * Object containing all information on one side of the comparator (can be used to calculate the total value)
 * 
 * @author Juergen Hartler
 *
 */
public class ExpressionForComparisonVO
{
  /** the global multiplier of an equation (multiplier outside of the brackets) */
  private String globalMultiplier_;
  /** the fragments and their multiplication factors */
  private Vector<FragmentMultVO> fragments_;
  
  /**
   * Constructor requires the found fragments and its multiplication factors, and the global multiplier
   * @param fragments the fragments and their multiplication factors 
   * @param globalMultiplier the global multiplier of an equation (multiplier outside of the brackets)
   */
  public ExpressionForComparisonVO(Vector<FragmentMultVO> fragments, String globalMultiplier)
  {
    this.fragments_ = fragments;
    this.globalMultiplier_ = globalMultiplier;
  }


  /**
   * 
   * @return the global multiplier value
   */
  public double getGlobalMultiplier()
  {
    return Double.parseDouble(globalMultiplier_);
  }

  /**
   * 
   * @return the global multiplier as String
   */
  public String getGlobalMultiplierString()
  {
    return globalMultiplier_;
  }

  /**
   * 
   * @return the fragments and their multiplication factors
   */
  public Vector<FragmentMultVO> getFragments()
  {
    return fragments_;
  }
  
  /**
   * 
   * @param found the found peak identifications; key: fragment name; value: an object containing the peak identification
   * @param basePeak the area of the base peak (if it is required in this equation)
   * @return the calcualted value for the expression
   */
  public double evaluateExpression(Hashtable<String,CgProbe> found, Float basePeak){
    double sum = 0d;
    for (FragmentMultVO multVO : fragments_){
      double value = 0d;
      if (multVO.getFragmentName().equalsIgnoreCase(IntensityRuleVO.BASEPEAK_NAME)) value = basePeak;
      else if (found!=null && found.containsKey(multVO.getFragmentName())) value = found.get(multVO.getFragmentName()).Area;
      value = value*multVO.getMultFactor();
      if (!multVO.isPositive()) value = value*-1d;
      sum += value;
    }
    sum = sum*Double.parseDouble(globalMultiplier_);
    return sum;
  }
  
  /**
   * 
   * @return the position assignment for a fragment
   */
  public int getPosition(){
    int position = 0;
    for (FragmentMultVO frag : fragments_){
      if (!frag.getFragmentName().equalsIgnoreCase(IntensityRuleVO.BASEPEAK_NAME)){
        position = frag.getPosition();
        break;
      }
    }
    return position;
  }
  
  /**
   * checks if this VO makes use of a fragment with this name
   * @param fragmentName the name of the fragment
   * @return true if the fragment is used
   */
  public boolean containsFragment(String fragmentName){
    boolean fragmentThere = false;
    for (FragmentMultVO fragment : fragments_){
      if (fragment.getFragmentName().equalsIgnoreCase(fragmentName)){
        fragmentThere = true;
        break;
      }
    }
    return fragmentThere;
  }
  
  /**
   * 
   * @return true if there any absolute relationships (comparison to the base peak) are involved
   */
  public boolean isAbsoluteComparison(){
    if (Double.parseDouble(globalMultiplier_)<=0d) return false;
    boolean isAbsolute = false;
    for (FragmentMultVO frag : fragments_){
      if (frag.getFragmentName().equalsIgnoreCase(IntensityRuleVO.BASEPEAK_NAME) && frag.getMultFactor()>0){
        isAbsolute = true;
        break;
      }
    }
    return isAbsolute;
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
    ExpressionForComparisonVO other = (ExpressionForComparisonVO) obj;
    return Objects.equals(fragments_, other.fragments_)
        && Objects.equals(globalMultiplier_, other.globalMultiplier_);
  }
  
}
