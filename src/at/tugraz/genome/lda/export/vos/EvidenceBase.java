/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.tugraz.genome.lda.export.vos;

/**
 * This class holds basic information that is common to the individual evidence value objects
 * 
 * @author Juergen Hartler
 *
 */
public class EvidenceBase
{
  
  /** the name of the modification*/
  private String modification_;
  /** the charge of the hit*/
  private int charge_;
  
  
  /**
   * constructor for basic evidence value object (a class holding basic information that is common to the individual evidence value objects)
   * @param modification name of the modification
   * @param charge charge of the hit
   */
  public EvidenceBase(String modification, int charge)
  {
    this.modification_ = modification;
    this.charge_ = charge;
  }

  /**
   * 
   * @return name of the modification
   */
  public String getModification()
  {
    return modification_;
  }

  
  /**
   * 
   * @return charge of the hit
   */
  public int getCharge()
  {
    return charge_;
  }
  
}
