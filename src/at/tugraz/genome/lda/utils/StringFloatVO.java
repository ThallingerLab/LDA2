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

package at.tugraz.genome.lda.utils;

/**
 * this value is for sorting Float values stored as String, to retain its original input value 
 * @author Juergen Hartler
 *
 */
public class StringFloatVO
{
  
  /** the String representation of the value*/
  String stringValue_;
  /** the float representation of the value*/
  float floatValue_;
  
  /**
   * @param value the String representation of the value
   */
  public StringFloatVO(String value){
    this.stringValue_ = value;
    this.floatValue_ = Float.parseFloat(value);
  }

  /**
   * @return the String representation of the value
   */
  public String getStringValue()
  {
    return stringValue_;
  }

  /**
   * @return the float representation of the value - for sorting
   */
  public float getFloatValue()
  {
    return floatValue_;
  }
  
  
}
