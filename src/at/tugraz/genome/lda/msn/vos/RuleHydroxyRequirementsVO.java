/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2019 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

/**
 * value object holding mandatory settings specific to hydroxylation
 * @author Juergen Hartler
 *
 */
public class RuleHydroxyRequirementsVO
{
  /** the OH number*/
  short oh_;
  /** the type of chain: $CHAIN, $ALKYLCHAIN, $ALKENYLCHAIN or $LCB*/
  short chainType_;
  /** whether it is mandatory and to what extend*/
  short mandatory_;
  
  /**
   * 
   * @param oh
   * @param chainType
   * @param mandatory
   */
  public RuleHydroxyRequirementsVO(short oh, short chainType, short mandatory)
  {
    this.oh_ = oh;
    this.chainType_ = chainType;
    this.mandatory_ = mandatory;
  }

  /**
   * 
   * @return OH number
   */
  public short getOh()
  {
    return oh_;
  }

  /** type of chain: $CHAIN, $ALKYLCHAIN, $ALKENYLCHAIN or $LCB
   * 
   * @return 
   */
  public short getChainType()
  {
    return chainType_;
  }
  
  /**
   * sets the chain Type
   * @param chainType type of chain: $CHAIN, $ALKYLCHAIN, $ALKENYLCHAIN or $LCB
   */
  public void setChainType(short chainType)
  {
    this.chainType_ = chainType;
  }

  /**
   * 
   * @return whether it is mandatory and to what extend
   */
  public short getMandatory()
  {
    return mandatory_;
  }

  
  /**
   * sets the mandatory level of this fragment
   * @param mandatory value for mandatory level 
   */
  public void setMandatory(short mandatory)
  {
    this.mandatory_ = mandatory;
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
    RuleHydroxyRequirementsVO other = (RuleHydroxyRequirementsVO) obj;
    return chainType_ == other.chainType_ && mandatory_ == other.mandatory_
        && oh_ == other.oh_;
  }
  
}
