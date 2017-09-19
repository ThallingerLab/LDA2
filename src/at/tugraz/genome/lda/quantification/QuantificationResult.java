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

package at.tugraz.genome.lda.quantification;

import java.util.Hashtable;
import java.util.Map;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class QuantificationResult
{
  /** the result identifications */
  private Hashtable<String,Vector<LipidParameterSet>> identifications_;
  /** constants holding the quantitation parameters */
  private LipidomicsConstants constants_;
  /** the levels of identification of each lipid class */
  Map<String,Integer> msLevels_;
  
  /**
   * constructor providing all the necessary information - changes later on are not possible
   * @param identifications the MSn identifications
   * @param constants quantification settings that where used
   * @param msLevels class specific identification MS levels
   */
  public QuantificationResult (Hashtable<String,Vector<LipidParameterSet>> identifications, LipidomicsConstants constants,
      Map<String,Integer> msLevels){
    identifications_ = identifications;
    constants_ = constants;
    msLevels_ = msLevels;
  }

  /**
   * 
   * @return the MS identification hash - key: class
   */
  public Hashtable<String,Vector<LipidParameterSet>> getIdentifications()
  {
    return identifications_;
  }

  /**
   * 
   * @return quantification settings that where used
   */
  public LipidomicsConstants getConstants()
  {
    return constants_;
  }

  /**
   * 
   * @return class specific identification MS levels
   */
  public Map<String,Integer> getMsLevels()
  {
    return msLevels_;
  }
  
  
  
  
}
