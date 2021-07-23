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

/**
 * Key/Value pair-class that holds an integer as key, and a double as value
 * @author Juergen Hartler
 *
 */
public class DoubleIntegerVO
{

  /** the key in integer format*/
  private Integer key_;
  /** the value in double format*/
  private Double value_;
  
  /**
   * 
   * @param key key in integer format
   * @param value value in double format
   */
  public DoubleIntegerVO(Integer key, Double value)
  {
    super();
    this.key_ = key;
    this.value_ = value;
  }

  /**
   * @return key in integer format
   */
  public Integer getKey()
  {
    return key_;
  }

  /**
   * @return value in double format
   */
  public Double getValue()
  {
    return value_;
  }

  /**
   * sets the key for this key/value pair VO
   * @param key key in integer format
   */
  public void setKey(Integer key)
  {
    this.key_ = key;
  }

  /**
   * sets the value for this key/value pair VO
   * @param value value in double format
   */
  public void setValue(Double value)
  {
    this.value_ = value;
  }

}
