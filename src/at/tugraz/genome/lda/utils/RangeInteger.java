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

package at.tugraz.genome.lda.utils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class RangeInteger
{
  /** begin value of the range */
  private int start_;
  /** end value of the range */
  private int stop_;
  
  /**
   * constructor of the RangeInteger
   * @param start start value of the range
   * @param stop stop value of the range
   */
  public RangeInteger(int start, int stop)
  {
    this.start_ = start;
    this.stop_ = stop;
  }

  /**
   * checks if a value is inside a range (equals counts to be inside the range)
   * @param value value to be checked
   * @return true if the value is inside the range
   */
  public boolean insideRange(int value){
    if (start_<=value && value<=stop_) return true;
    else return false;
  }

  public int getStart()
  {
    return start_;
  }

  public int getStop()
  {
    return stop_;
  }
  
  /**
   * extends the range to another one (only extension is possible, no reduction),
   * and it verifies that it does not extend a maximum range
   * @param other the range to extend this range
   * @param maxRange the maximum allowed range
   */
  public void extendToOtherRanges(RangeInteger other, RangeInteger maxRange){
    if (other.start_<start_) start_ = other.start_;
    if (start_<maxRange.start_) start_ = maxRange.start_;
    if (other.stop_>stop_) stop_ = other.stop_;
    if (stop_>maxRange.stop_) stop_ = maxRange.stop_;
  }
}
