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

import at.tugraz.genome.lda.LipidomicsConstants;

/**
 * value object containing necessare information about a potential MSn fragment itself
 * @author Juergen Hartler
 *
 */
public class FragmentVO
{
  // the name of the fragment
  private String name_;
  // the chemical formula of the fragment
  private String formula_;
  // the charge state of the fragment
  private int charge_;
  // the MSn level at which the fragment is observed
  private int msLevel_;
  // is the presence of the fragment mandatory for this rule - different options possible according to the MANDATORY_... specifications in FragmentRuleVO
  private short mandatory_;
  // the m/z value of the fragment
  private double mass_;
  /** if the fragment is a chain - the type of the chain*/
  private short chainType_;
  
  /**
   * all of the required information for the VO has to be provided in the constructor
   * @param name name of the fragment
   * @param mass  m/z value of the fragment
   * @param formula chemical formula of the fragment
   * @param charge charge state of the fragment
   * @param msLevel MSn level at which the fragment is observed
   * @param mandatory is the presence of the fragment mandatory for this rule
   * @param fromOtherSpecies does this fragment originate from another species (isobar)
   */
  public FragmentVO(String name, double mass, String formula, int charge, int msLevel,
      short mandatory)
  {
    this(name, mass, formula, charge, msLevel, mandatory, LipidomicsConstants.CHAIN_TYPE_FA_ACYL);
  }
    
 /**
   * all of the required information for the VO has to be provided in the constructor
   * @param name name of the fragment
   * @param mass  m/z value of the fragment
   * @param formula chemical formula of the fragment
   * @param charge charge state of the fragment
   * @param msLevel MSn level at which the fragment is observed
   * @param mandatory is the presence of the fragment mandatory for this rule
   * @param fromOtherSpecies does this fragment originate from another species (isobar)
   * @param chainType the type of the fatty acid chain ACYL_CHAIN_/ALKYL_CHAIN/ALKENYL_CHAIN
  */
  public FragmentVO(String name, double mass, String formula, int charge, int msLevel,
        short mandatory, short chainType)
    {    
    this.name_ = name;
    this.mass_ = mass;
    this.formula_ = formula;
    this.charge_ = charge;
    this.msLevel_ = msLevel;
    this.mandatory_ = mandatory;
    this.chainType_ = chainType;
  }

  /**
   * 
   * @return name of the fragment
   */
  public String getName()
  {
    return name_;
  }

  /**
   * 
   * @return chemical formula of the fragment
   */
  public String getFormula()
  {
    return formula_;
  }

  /**
   * 
   * @return charge state of the fragment
   */
  public int getCharge()
  {
    return charge_;
  }

  /**
   * 
   * @return MSn level at which the fragment is observed
   */
  public int getMsLevel()
  {
    return msLevel_;
  }

  /**
   * 
   * @return is the presence of the fragment mandatory for this rule
   */
  public short isMandatory()
  {
    return mandatory_;
  }

  /**
   * 
   * @return m/z value of the fragment
   */
  public double getMass()
  {
    return mass_;
  }
  
  
  
  public short getChainType()
  {
    return chainType_;
  }
  
  
  /**
   * @return String representative of the VO - for debugging purposes
   */
  public String toString(){
    String toPrint =   "Name: "+name_;
    toPrint += "; Mass: "+mass_;
    toPrint += "; Formula: "+formula_;
    toPrint += "; Charge: "+charge_;
    toPrint += "; MS: "+msLevel_;
    toPrint += "; Mand: "+mandatory_;
    return toPrint;
  }
  
}
