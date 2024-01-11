/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2019 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;

/**
 * object holding all the settings specific to hydroxylation whether rules are mandatory 
 * @author Juergen Hartler
 *
 */
public class RuleHydroxyRequirementSet
{
  
  /** Vector containing the various hydroxylation requirements for one intensity rule*/
  private Vector<RuleHydroxyRequirementsVO> requirements_;
  
  
  /**
   * constructor requiring the various hydroxylation requirements for one intensity rule
   * @param requirements hydroxylation requirements for one intensity rule
   */
  public RuleHydroxyRequirementSet(Vector<RuleHydroxyRequirementsVO> requirements) {
    this.requirements_ = requirements;
  }
  
  
  /**
   * sets the hydroxylation requirements that are undefined to the mandatory parameter
   * @param mandatory the mandatory parameter for undefined hydroxylation requirements
   */
  public void setUndefinedMandatorySettings(short mandatory) {
    for (RuleHydroxyRequirementsVO req : requirements_) {
      if (req.getMandatory()==FragmentRuleVO.MANDATORY_UNDEFINED)
        req.setMandatory(mandatory);
    }
  }
  
  /**
   * checks whether this OH number is present in the requirement settings - if not this rule does not apply for this hydroxylation
   * @param ohNumber the OH number to check
   * @return true when this OH number is present in the requirement settings
   */
  public boolean hasEntry(short ohNumber) {
    for (RuleHydroxyRequirementsVO req : requirements_) {
      if (req.getOh()==ohNumber)
        return true;
    }
    return false;
  }
  
  
  /**
   * get mandatory level
   * @param chainType the type of chain
   * @param ohNumber the oh number
   * @return mandatory level
   */
  public short getMandatory(short chainType, short ohNumber) {
    short mand = FragmentRuleVO.MANDATORY_UNDEFINED;
    Vector<RuleHydroxyRequirementsVO> reqs = getEntry(ohNumber);
    for (RuleHydroxyRequirementsVO req : reqs) {
      if (req.getChainType()==chainType) {
        mand = req.getMandatory();
      }
    }
    return mand;
  }
  
  
  /**
   * returns the hydroxylation requirements for one OH number
   * @param ohNumber the hydroxylation requirements to look for
   * @return list of hydroxylation requirements for one OH number
   */
  public Vector<RuleHydroxyRequirementsVO> getEntry(short ohNumber) {
    Vector<RuleHydroxyRequirementsVO> reqs = new Vector<RuleHydroxyRequirementsVO>();
    for (RuleHydroxyRequirementsVO req : requirements_) {
      if (req.getOh()==ohNumber)
        reqs.add(req);
    }    
    return reqs;
  }
  
  
  /**
   * for rules containing the same chain types: if no chain type is set, set it to ruleChainType
   * for rules containing different chain types: for $LCB chains, set the chain type to $LCB where no chain type is set; for all others: throw an error
   * @param ruleChainType the chain type of the rule
   * @param allFrags the fragments of the intensity rule
   * @param lineNumber the line number of the rule set - for throwing errors
   * @param headSection is this the head section
   * @throws RulesException thrown when there is something wrong with the rule
   */
  public void checkAndCorrectChainTypes(short ruleChainType, Vector<FragmentMultVO> allFrags, int lineNumber, boolean headSection) throws RulesException {
    if (headSection)
      return;
    if (ruleChainType!=IntensityRuleVO.DIFF_CHAIN_TYPES) {
      //there is only one type of chain possible - if no chain type set -> correct - if a chain type is set, it must be this one, otherwise throw an error
      for (RuleHydroxyRequirementsVO vo : requirements_) {
        if (vo.getChainType()==LipidomicsConstants.CHAIN_TYPE_NO_CHAIN)
          vo.setChainType(ruleChainType);
        else if (vo.getChainType()!=ruleChainType)
          throw new RulesException("The set chain type for the \""+FragRuleParser.FRAGMENT_HYDROXY+"\" setting is not possible for this intensity rule! Error at line number "+lineNumber+"!");
      }
    }else {
      boolean hasLcbChain = false;
      Set<Short> availableChainTypes = new HashSet<Short>();
      for (FragmentMultVO frag : allFrags) {
        availableChainTypes.add(frag.getFragmentType());
      }
      if (availableChainTypes.contains(LipidomicsConstants.CHAIN_TYPE_LCB))
        hasLcbChain = true;
      //this is a rule consisting of different chain types
      for (RuleHydroxyRequirementsVO vo : requirements_){
        if (vo.chainType_==LipidomicsConstants.CHAIN_TYPE_NO_CHAIN) {
          //a type CHAIN_TYPE_NO_CHAIN is only allowed when $LCB chains are present
          if (hasLcbChain) vo.setChainType(LipidomicsConstants.CHAIN_TYPE_LCB);
          else throw new RulesException("For \""+FragRuleParser.FRAGMENT_HYDROXY+"\", when there are more fragment types are present, but no "+FragmentRuleVO.LCB_NAME+" is set, the chain type as to be specified (consult manual)! Error at line: "+lineNumber+"!");
        } else {
          if (!availableChainTypes.contains(vo.getChainType()))
            throw new RulesException("The set chain type for the \""+FragRuleParser.FRAGMENT_HYDROXY+"\" setting is not possible for this intensity rule! Error at line number "+lineNumber+"!");
        }
      }
    }
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
    RuleHydroxyRequirementSet other = (RuleHydroxyRequirementSet) obj;
    return Objects.equals(requirements_, other.requirements_);
  }
  
}
