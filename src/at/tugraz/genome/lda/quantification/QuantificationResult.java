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
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;

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
  private Map<String,Integer> msLevels_;
  /** the character encoding of the number of hydroxylation sites for the FA*/
  private HydroxyEncoding faHydroxyEncoding_;
  /** the character encoding of the number of hydroxylation sites for the LCB*/
  private HydroxyEncoding lcbHydroxyEncoding_;

  
  /**
   * constructor providing all the necessary information - changes later on are not possible
   * @param identifications the MSn identifications
   * @param constants quantification settings that where used
   * @param msLevels class specific identification MS levels
   * @param faHydroxyEncoding the character encoding of the number of hydroxylation sites for the FA
   * @param lcbHydroxyEncoding the character encoding of the number of hydroxylation sites for the LCB
   */
  public QuantificationResult (Hashtable<String,Vector<LipidParameterSet>> identifications, LipidomicsConstants constants,
      Map<String,Integer> msLevels, HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding){
    identifications_ = identifications;
    constants_ = constants;
    msLevels_ = msLevels;
    faHydroxyEncoding_ = faHydroxyEncoding;
    lcbHydroxyEncoding_ = lcbHydroxyEncoding;
  }
  
  /**
   * constructor to create a deep copy
   * @param that
   */
  public QuantificationResult(QuantificationResult that) {
    this(that.getIdentifications(), that.getConstants(), that.getMsLevels(), that.getFaHydroxyEncoding(), that.getLcbHydroxyEncoding());
  }
  
  /**
   * 
   * @return the MS identification hash => key: class
   */
  public Hashtable<String,Vector<LipidParameterSet>> getIdentifications()
  {
    return identifications_;
  }
  
  /**
   * set the Identifications (needed for OxSheets)
   */
  public void setIdentifications(Hashtable<String,Vector<LipidParameterSet>> identifications)
  {
    this.identifications_ = identifications;
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

  /**
   * 
   * @return the character encoding of the number of hydroxylation sites for the FA
   */
  public HydroxyEncoding getFaHydroxyEncoding()
  {
    return faHydroxyEncoding_;
  }

  /**
   * 
   * @return the character encoding of the number of hydroxylation sites for the LCB
   */
  public HydroxyEncoding getLcbHydroxyEncoding()
  {
    return lcbHydroxyEncoding_;
  }
  
  
  
  
}
