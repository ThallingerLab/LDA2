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

import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MSnDebugVO
{
  
  /** the reason for the discard is unknown */
  public final static int UNKNOWN_REASON = 0;
  /** there is no peak at this m/z value */
  public final static int NO_PEAK_THERE = 1;
  /** the peak is lower than the base peak intensity cutoff */
  public final static int BELOW_BASE_PEAK_CUTOFF = 2;
  
  /** the found chain does not have any matching partner chains to match the total
   *  number of carbon atoms and double bonds */
  public final static int NO_CHAIN_COMINATION_POSSIBLE = 3;
  
  /** the found chain combination is of minor intensity - below the defined threshold */
  public final static int COMBINATION_LOWER_CHAIN_CUTOFF = 4;
  
  /** identifier to be used for removal due to the lack of matching partner chains */
  public final static String ID_COMBI_REMOVE = "COMBINATION";
  
  /** the overall status of the identification */
  private int status_;
  
  /** contains reasons for the discard of head group fragments
   *  key of the hash-table is the fragment name, and
   *  the value contains the reason for the removal*/
  private Hashtable<String,Integer> headGroupStatuses_;
  
  /** contains reasons for the discard of head group fragments due to intensity rules
   * key is the rule identifier, and value is the rule itself */
  private Hashtable<String,IntensityRuleVO> violatedHeadRules_;
  
  /** contains reasons for the discard of potential fatty acid chains
   *  the first key is the name of the fatty acid chain
   *  the second key is either the name of the fragment or the rule identifier of an intensity rule
   *  the value is either the reason for the discard (Integer), or the intensity rule itself */
  private Hashtable<String,Hashtable<String,Object>> violatedChainRules_;
  
  /** contains reasons why a chain combination was discarded */
  private Hashtable<String,Integer> violatedCombinations_;
  
  /** contains the position rules that where not fulfilled by valid chain combinations
   *  the first key is the identifier for the chain combination
   *  the second key is the rule identifier, and
   *  the value is the intensity rule itself */
  private Hashtable<String,Hashtable<String,IntensityRuleVO>> unfulfilledPositionRules_;
  
  /** contains the position assignments that contradict to one another 
   *  the key is the identifier for the chain combinations
   *  the first vector is the sequence of unfulfilled chain combinations for this chain combination
   *  the second vector contains the two contradicting positions rules with their respective assignments*/
  private Hashtable<String,Vector<Vector<IntensityRuleVO>>> contradictingPositionRules_;
  
  /** is the spectrum coverage fulfilled*/
  private boolean spectrumCoverageFulfilled_;
  
  /**
   * constructor: initializes all the necessary hashes
   */
  public MSnDebugVO() {
    headGroupStatuses_ = new Hashtable<String,Integer>();
    violatedHeadRules_ = new Hashtable<String,IntensityRuleVO>();
    violatedChainRules_ = new Hashtable<String,Hashtable<String,Object>>();
    violatedCombinations_ = new Hashtable<String,Integer>();
    spectrumCoverageFulfilled_ = false;
    unfulfilledPositionRules_ = new Hashtable<String,Hashtable<String,IntensityRuleVO>>();
    contradictingPositionRules_ = new Hashtable<String,Vector<Vector<IntensityRuleVO>>>();
  }

  /**
   * 
   * @return the overall status of the identification (same as in LipidomicsMSnSet)
   */
  public int getStatus()
  {
    return status_;
  }

  /**
   * sets the overall status of the identification (same as in LipidomicsMSnSet)
   * @param status the status of the identification
   */
  public void setStatus(int status)
  {
    this.status_ = status;
  }
  
  /**
   * adds a head group fragment that will be discarded and the reason for the discard
   * @param name name of the head group fragment
   * @param status reason for the discard
   */
  public void addDiscardedHeadGroupFragment(String name, int status){
    headGroupStatuses_.put(name, status);
  }
  
  /**
   *
   * @return the discarded head group fragments and the reasons for the discard (key: fragment name; value: reason for the removal)
   */
  public Hashtable<String,Integer> getDiscardedHeadGroupFragments(){
    return headGroupStatuses_;
  }
  
  /**
   * adds violated head intensity rules
   * @param ruleVO the violated intensity rule 
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public void addViolatedHeadRule(IntensityRuleVO ruleVO) throws LipidCombinameEncodingException {
    violatedHeadRules_.put(ruleVO.getReadableRuleInterpretation(Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding()), ruleVO);
  }

  /**
   * 
   * @return the violated head rules (key: rule identifier; value: the rule itself)
   */
  public Hashtable<String,IntensityRuleVO> getViolatedHeadRules()
  {
    return violatedHeadRules_;
  }
  
  /**
   * adds a discarded chain fragment and the reason for the discard
   * @param faName the name of the concerned fatty acid
   * @param fragmentName the name of the fragment
   * @param status reason for the discard
   */
  public void addViolatedChainFragment(String faName, String fragmentName, int status){
    Hashtable<String,Object> violatedRules = new Hashtable<String,Object>();
    if (violatedChainRules_.containsKey(faName)) violatedRules = violatedChainRules_.get(faName);
    violatedRules.put(fragmentName, status);
    violatedChainRules_.put(faName, violatedRules);
  }
  
  /**
   * adds a violated chain intensity rule
   * @param faName the name of the concerned fatty acid
   * @param ruleVO the violated intensity rule 
   */
  public void addViolatedChainRule(String faName, IntensityRuleVO ruleVO) {
    Hashtable<String,Object> violatedRules = new Hashtable<String,Object>();
    if (violatedChainRules_.containsKey(faName)) violatedRules = violatedChainRules_.get(faName);
    violatedRules.put(ruleVO.getRuleIdentifier(), ruleVO);
    violatedChainRules_.put(faName, violatedRules);
  }
  
  /**
   * 
   * @param combiName adds a chain combination that has been discarded
   * @param status reason for the discard
   */
  public void addViolatedCombinations(String combiName, int status){
    violatedCombinations_.put(combiName, status);
  }
  
  /**
   * 
   * @return the discarded chain fragments and chain intensity rules (1st key: fatty acid name; 2nd key: fragment name or rule identifier value: status (Integer) or rule (InntensityRuleVO))
   */
  public Hashtable<String,Hashtable<String,Object>> getViolatedChainRules(){
    return violatedChainRules_;
  }

  /**
   * 
   * @return the discarded chain combinations (key: combination identifier; value: reason for discard)
   */
  public Hashtable<String,Integer> getViolatedCombinations()
  {
    return violatedCombinations_;
  }

  /**
   * 
   * @return if the hit was discarded due to the spectrum coverage
   */
  public boolean isSpectrumCoverageFulfilled()
  {
    return spectrumCoverageFulfilled_;
  }

  /**
   * defines if the spectrum coverage is fulfilled
   * @param spectrumCoverageFulfilled is the spectrum coverage fulfilled?
   */
  public void setSpectrumCoverageFulfilled(boolean spectrumCoverageFulfilled)
  {
    this.spectrumCoverageFulfilled_ = spectrumCoverageFulfilled;
  }
  
  /**
   * adds rules where there was no position assignment possible
   * @param combiName name of the fatty acid chain combination
   * @param ruleVO the unfulfilled intensity rule for positional assignment
   */
  public void addUnfulfilledPositionRule(String combiName, IntensityRuleVO ruleVO){
    Hashtable<String,IntensityRuleVO> unfulfilled = new Hashtable<String,IntensityRuleVO>();
    if (unfulfilledPositionRules_.containsKey(combiName)) unfulfilled = unfulfilledPositionRules_.get(combiName);
    unfulfilled.put(ruleVO.getRuleIdentifier(), ruleVO);
    unfulfilledPositionRules_.put(combiName, unfulfilled);
  }
  
  /**
   * adds position rules that lead to contradictory results
   * @param combiName the name of the affected fatty acid chain combination
   * @param firstRule one of the rules leading to contradictory results
   * @param secondRule the other rule that leads to contradictory results
   */
  public void addContradictingPositionRules(String combiName, IntensityRuleVO firstRule, IntensityRuleVO secondRule){
    Vector<Vector<IntensityRuleVO>> rules = new Vector<Vector<IntensityRuleVO>>();
    if (contradictingPositionRules_.containsKey(combiName)) rules = contradictingPositionRules_.get(combiName);
    Vector<IntensityRuleVO> rulePair = new Vector<IntensityRuleVO>();
    rulePair.add(firstRule);
    rulePair.add(secondRule);
    rules.add(rulePair);
    contradictingPositionRules_.put(combiName,rules);
  }

  /**
   * 
   * @return rules leading to contradictory results (key: combination identifier; 1st vector: enumeration of rules; 2nd vector: contradicting rules)
   */
  public Hashtable<String,Vector<Vector<IntensityRuleVO>>> getContradictingPositionRules()
  {
    return contradictingPositionRules_;
  }

  /**
   * 
   * @return rules where there was no position assignment possible (1st key: combination identifier; 2nd key: rule identifier; value: rule)
   */
  public Hashtable<String,Hashtable<String,IntensityRuleVO>> getUnfulfilledPositionRules()
  {
    return unfulfilledPositionRules_;
  }
  
  
}
