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

import java.util.Comparator;
import java.util.Objects;

import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * value object holding information about a lipid acid chain
 * i.e. # C atoms, # double bonds, mass value, chemical formula
 * @author Juergen Hartler
 *
 */
public class FattyAcidVO implements Comparable<FattyAcidVO>
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
  /** omegaPosition*/
  private int omegaPosition_;
    // oxidation state
  private String oxState_;
  
  /**
   * all the information has to be provided in the constructor
   * @param chainType the type of chain: LipidomicsConstants.CHAIN_TYPE_FA or LipidomicsConstants.CHAIN_TYPE_LCB
   * @param cAtoms number of C atoms
   * @param doubleBonds number of double bonds
   * @param ohNumber
   * @param mass mass value
   * @param formula chemical formula
   */
  public FattyAcidVO(short chainType, String prefix, int cAtoms, int doubleBonds, int ohNumber, double mass, String formula, String oxState){
    this.chainType_ = chainType;
    this.prefix_ = prefix;
    this.cAtoms_ = cAtoms;
    this.doubleBonds_ = doubleBonds;
    this.ohNumber_ = ohNumber;
    this.mass_ = mass;
    this.formula_ = formula;
    this.omegaPosition_ = -1;
    this.oxState_ = oxState;
  }
  
  public FattyAcidVO(FattyAcidVO other) {
    this(other.getChainType(), other.getPrefix(), other.getcAtoms(), other.getDoubleBonds(), other.getOhNumber(), other.getMass(), other.getFormula(), other.getOxState());
    this.setOmegaPosition(other.getOmegaPosition());
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
  public void correctChainType(short chainType)
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
   * @param prefix the prefix
   */
  public void setPrefix(String prefix)
  {
    this.prefix_ = prefix;
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
   * 
   * @param formula the chemical formula
   */
  public void setFormula(String formula)
  {
    this.formula_ = formula;
  }
  
  /**
   * 
   * @return the omega position of the last double bond
   */
  public int getOmegaPosition()
  {
    return omegaPosition_;
  }


  /**
   * sets the omega position of the last double bond
   * @param omegaPosition the omega position
   */
  public void setOmegaPosition(int omegaPosition)
  {
    this.omegaPosition_ = omegaPosition;
  }

  /**
   * 
   * @return oxidation state
   */
  public String getOxState()
  {
    return oxState_;
  }

  /**
   * @return name consisting of number of C atoms : number of double bonds ;oxState (n- C=C position)
   */
  public String getCarbonDbsId()
  {
    return StaticUtils.generateLipidNameString(prefix_+String.valueOf(cAtoms_),doubleBonds_,omegaPosition_,oxState_);
  }
  
  /**
   * @return an unique identifier for the chain that includes the chain type, the number of OHs, the #C-atoms, the #double bonds, omega position and the isotope prefix
   */
  public String getChainId() {
  	return getChainIdDetailed(true, true);
  }
  
  /**
   * Get a chainId specifying the level of detail the id should include.
   * @param includePrefix
   * @param includeOmegaPosition
   * @return
   */
  public String getChainIdDetailed(boolean includePrefix, boolean includeOmegaPosition) {
    return StaticUtils.encodeLipidNameForCreatingCombis(this,includePrefix,includeOmegaPosition);    
  }
  
  /**
   * @return String representing the values of the value object
   */
  public String toString(){
    return ("Name: "+getChainId()+" Mass: "+mass_+" Formula: "+formula_);
  }
  
  
  @Override
  public boolean equals(Object obj)
  {
    if (this == obj)
      return true;
    if (obj == null)
      return false;
    if (getClass() != obj.getClass())
      return false;
    FattyAcidVO other = (FattyAcidVO) obj;
    return cAtoms_ == other.cAtoms_ && chainType_ == other.chainType_
        && doubleBonds_ == other.doubleBonds_
        && Objects.equals(formula_, other.formula_)
        && Double.doubleToLongBits(mass_) == Double
            .doubleToLongBits(other.mass_)
        && ohNumber_ == other.ohNumber_
        && omegaPosition_ == other.omegaPosition_
        && Objects.equals(prefix_, other.prefix_);
  }
  
  /**
   * @param other
   * @return true if all fields except the omega position are identical
   */
  public boolean equalsNotConsideringOmegaPosition(FattyAcidVO other)
  {
  	return cAtoms_ == other.cAtoms_ && chainType_ == other.chainType_
        && doubleBonds_ == other.doubleBonds_
        && Objects.equals(formula_, other.formula_)
        && Double.doubleToLongBits(mass_) == Double
            .doubleToLongBits(other.mass_)
        && ohNumber_ == other.ohNumber_
        && Objects.equals(prefix_, other.prefix_);
  }

	@Override
	/**
	 * Compares the fields of this class.
	 * In hierarchical order, these fields are compared until there is a difference:
	 * chainType_, cAtoms_, doubleBonds_, ohNumber_, omegaPosition_, prefix_.
	 * @param other
	 * @return
	 */
	public int compareTo(FattyAcidVO other)
	{
		 return Comparator.comparing(FattyAcidVO::getChainType).reversed()
							.thenComparing(FattyAcidVO::getcAtoms)
							.thenComparing(FattyAcidVO::getDoubleBonds)
							.thenComparing(FattyAcidVO::getOhNumber)
							.thenComparing(FattyAcidVO::getOmegaPosition)
							.thenComparing(FattyAcidVO::getPrefix)
							.compare(this, other);
	}
  
}
