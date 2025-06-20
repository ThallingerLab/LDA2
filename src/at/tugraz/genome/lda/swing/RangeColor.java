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

import java.awt.Color;

import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * checks if value is in certain range, and stores a name,
 * a color and a peak identification for this region
 * 
 * @author Juergen Hartler
 *
 */
public class RangeColor extends Range
{
  
  /** the name of the identified peak */
  private String name_;
  /** the color that shall be applied for this peak */
  private Color color_;
  /** the peak identification */
  private CgProbe probe_;
  
  /**
   * constructor providing all of the necessary values
   * @param name name of the identified peak
   * @param color color that shall be applied for this peak
   * @param probe the peak identification
   * @param start start of identification
   * @param stop stop of identification
   */
  public RangeColor (String name, Color color, CgProbe probe, float start, float stop){
    super(start,stop);
    this.name_ = name;
    this.color_ = color;
    this.probe_ = probe;
  }

  /**
   * 
   * @return name of the identified peak
   */
  public String getName()
  {
    return name_;
  }

  /**
   * 
   * @return color that shall be applied for this peak
   */
  public Color getColor()
  {
    return color_;
  }
  
  /**
   * 
   * @return the peak identification
   */
  public CgProbe getProbe(){
    return probe_;
  }
}
