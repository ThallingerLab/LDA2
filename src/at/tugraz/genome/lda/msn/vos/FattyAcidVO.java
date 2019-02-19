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

package at.tugraz.genome.lda.msn.vos;

import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * value object holding information about a lipid acid chain
 * i.e. # C atoms, # double bonds, mass value, chemical formula
 * @author Juergen Hartler
 *
 */
public class FattyAcidVO
{
  /** the type of chain: LipidomicsConstants.CHAIN_TYPE_FA or LipidomicsConstants.CHAIN_TYPE_LCB*/
  private short chainType_;
  // prefix to separate fatty acid
  private String prefix_;
  // number of C atoms
  private int cAtoms_;
  // number of double bonds
  private int doubleBonds_;
  // the mass value
  private double mass_;
  // chemical formula
  private String formula_;
  /** the number of OH groups*/
  private int ohNumber_;
  
  /**
   * all the information has to be provided in the constructor
   * @param chainType the type of chain: LipidomicsConstants.CHAIN_TYPE_FA or LipidomicsConstants.CHAIN_TYPE_LCB
   * @param cAtoms number of C atoms
   * @param doubleBonds number of double bonds
   * @param ohNumber
   * @param mass mass value
   * @param formula chemical formula
   */
  public FattyAcidVO(short chainType, String prefix, int cAtoms, int doubleBonds, int ohNumber, double mass, String formula){
    this.chainType_ = chainType;
    this.prefix_ = prefix;
    this.cAtoms_ = cAtoms;
    this.doubleBonds_ = doubleBonds;
    this.ohNumber_ = ohNumber;
    this.mass_ = mass;
    this.formula_ = formula;
  }
  
  
  
  /**
   * 
   * @return the type of chain: LipidomicsConstants.CHAIN_TYPE_FA or LipidomicsConstants.CHAIN_TYPE_LCB
   */
  public short getChainType()
  {
    return chainType_;
  }

 
  /** 
   * allows for chain type correction when a general parser was used
   * @param chainType the new chain type
   */
  public void correctChainType_(short chainType)
  {
    this.chainType_ = chainType;
  }



  /**
   * 
   * @return the prefix
   */
  public String getPrefix()
  {
    return prefix_;
  }

  /**
   * 
   * @return number of C atoms
   */
  public int getcAtoms()
  {
    return cAtoms_;
  }

  /**
   * 
   * @return number of double bonds
   */
  public int getDoubleBonds()
  {
    return doubleBonds_;
  }
  
  
  /**
   * 
   * @return the number of hydroxylation groups
   */
  public int getOhNumber()
  {
    return ohNumber_;
  }


  /**
   * 
   * @return mass value
   */
  public double getMass()
  {
    return mass_;
  }
  
  /**
   * 
   * @return chemical formula
   */
  public String getFormula()
  {
    return formula_;
  }

  
  /**
   * @return name consisting of number of C atoms : number of double bonds
   */
  public String getCarbonDbsId(){
    return StaticUtils.generateLipidNameString(prefix_+String.valueOf(cAtoms_),doubleBonds_);
  }
  
  /**
   * 
   * @return an unique identifier for the chain that includes the chain type, the number of OHs, the #C-atoms, the #double bonds, and the isotope prefix
   */
  public String getChainId() {
    return StaticUtils.encodeLipidNameForCreatingCombis(this,true);
  }
  
  /**
   * 
   * @return the sames as getChainId(), but without the isotope prefix
   */
  public String getChainIdWOPrefix() {
    return StaticUtils.encodeLipidNameForCreatingCombis(this,false);    
  }
  
  /**
   * @return String representing the values of the value object
   */
  public String toString(){
    return ("Name: "+getChainId()+" Mass: "+mass_+" Formula: "+formula_);
  }
  
}
