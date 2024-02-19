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

import java.util.Vector;

/**
 * for specifying ranges and checking if values are within the specified range
 * @author Juergen Hartler
 *
 */
public class Range
{
  /** begin value of the range */
  private float start_;
  /** end value of the range */
  private float stop_;
  
  /**
   * constructor of the Range
   * @param start start value of the range
   * @param stop stop value of the range
   */
  public Range(float start, float stop)
  {
    this.start_ = start;
    this.stop_ = stop;
  }

  /**
   * checks if a value is inside a range
   * @param value value to be checked
   * @return true if the value is inside the range
   */
  public boolean insideRange(double value){
    if (start_<value && value<stop_) return true;
    else return false;
  }

  public float getStart()
  {
    return start_;
  }

  public float getStop()
  {
    return stop_;
  }
  
  public String toString(){
    return start_/60f+"-"+stop_/60f;
  }
  
  public boolean overlap(Range vo){
    // this is a complete coverage
 //   if ((start_<vo.getStart()&&vo.getStop()<stop_) || (vo.getStart()<start_&&stop_<vo.getStop())||
    // these are overlapping regions    
    if  (start_<vo.getStart()&&stop_>vo.getStart() || vo.getStart()<start_&&vo.getStop()>start_){
      return true;
    }else{
      return false;
    }
  }
  
  public void combine(Range vo){
    if (vo.getStart()<start_) start_ = vo.getStart();
    if (vo.getStop()>stop_) stop_ = vo.getStop();
  }
  
  public static Vector<Range> reduce(Range vo, Range reducer){
    Vector<Range> ranges = new Vector<Range>();
    if (reducer.getStart()<vo.getStart() && vo.getStop()<reducer.getStop()) return ranges;
    if (!vo.overlap(reducer)){
      ranges.add(vo);
    }else{
      if (vo.getStart()<reducer.getStart()){
        ranges.add(new Range(vo.getStart(),reducer.getStart()));
      }
      if (vo.getStop()>reducer.getStop()){
        ranges.add(new Range(reducer.getStop(),vo.getStop()));
      }
    }
    return ranges;
  }
}
