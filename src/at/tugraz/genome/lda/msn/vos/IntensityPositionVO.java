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

import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * Value object containing all necessary information for a rule for intensity comparison
 * plus the two fatty acid that where compared by position finding
 * @author Juergen Hartler
 *
 */
public class IntensityPositionVO extends IntensityRuleVO
{
  private FattyAcidVO biggerFA_;
  private FattyAcidVO smallerFA_;
  
  private boolean negated_;
  private int derivedPosition_;
  
  /**
   * constructor containing the original rule plus the two fatty acids
   * @param rule rule that was applied
   * @param biggerFA fatty acid that was verified at at the greater part of the comparator
   * @param smallerFA fatty acid that was verified at at the lesser part of the comparator
   */
  public IntensityPositionVO(IntensityRuleVO rule, FattyAcidVO biggerFA, FattyAcidVO smallerFA){
    super(rule);
    this.biggerFA_ = biggerFA;
    this.smallerFA_ = smallerFA;
    negated_ = false;
    derivedPosition_ = -1;
  }

  /**
   * @return fatty acid that was verified at at the greater part of the comparator
   */
  public FattyAcidVO getBiggerFA()
  {
    return biggerFA_;
  }

  /**
   * @return fatty acid that was verified at at the lesser part of the comparator
   */
  public FattyAcidVO getSmallerFA()
  {
    return smallerFA_;
  }
  
  /**
   * 
   * @return the human readable String representation for storing in Excel
   */
  public String getReadableRuleInterpretation() {
    Hashtable<String,String> originalToReplacementBigger = new Hashtable<String,String>();
    for (String name : getBiggerNonBasePeakNames()){
      String nameString = name+"["+getBiggerPosition()+"]";
      String value = StaticUtils.getChainFragmentDisplayName(name,biggerFA_.getCarbonDbsId())+"["+getBiggerPosition()+"]";
      originalToReplacementBigger.put(nameString, value);
    }
    Hashtable<String,String> originalToReplacementSmaller = new Hashtable<String,String>();
    for (String name : getSmallerNonBasePeakNames()){
      String nameString = name+"["+getSmallerPosition()+"]";
      String value = StaticUtils.getChainFragmentDisplayName(name,smallerFA_.getCarbonDbsId())+"["+getSmallerPosition()+"]";
      originalToReplacementSmaller.put(nameString, value);
    }
    String[] biggerSmaller = splitToBiggerAndSmallerPart(equation_);
    String biggerPart = biggerSmaller[0];
    String smallerPart = biggerSmaller[1];
    biggerPart = replaceParts(biggerPart, originalToReplacementBigger);
    smallerPart = replaceParts(smallerPart, originalToReplacementSmaller);
    String returnString = "";
    if (equation_.indexOf(">")!=-1){
      returnString = biggerPart+">"+smallerPart;
    } else if (equation_.indexOf("<")!=-1){
      returnString = smallerPart+"<"+biggerPart;
    }
    if (negated_) returnString = biggerFA_.getCarbonDbsId()+"["+derivedPosition_+"] BECAUSE NOT FULFILLED: "+returnString;
    return returnString;
  }
  
  /**
   * checks the objects for equality
   */
  public boolean equals(IntensityRuleVO other){
    if (!super.equals(other)) return false;
    if (!(other instanceof IntensityPositionVO)) return false;
    IntensityPositionVO chain = (IntensityPositionVO)other;
    if (!biggerFA_.getChainId().equalsIgnoreCase(chain.biggerFA_.getChainId())) return false;
    if (!smallerFA_.getChainId().equalsIgnoreCase(chain.smallerFA_.getChainId())) return false;
    return true;
  }
  
  /**
   * gets the identified position for the fatty acid
   * @param fa fatty acid for which the position has to be returned
   * @return the position of the fatty acid
   */
  public int getPositionByFA(String fa){
    if (negated_) return derivedPosition_;
    int position=-1;
    if (biggerFA_.getChainId().equalsIgnoreCase(fa)) position = getBiggerPosition();
    else if (smallerFA_.getChainId().equalsIgnoreCase(fa)) position = getSmallerPosition();
    return position;
  }

  /**
   * inverse method of getReadableRuleInterpretation() - fetches the fatty acids from this interpretation
   * @param storedRule the readable rule interpretation
   * @param ruleVO rule that is required to form an IntensityPositionVO
   * @param chainFragments the found chain fragments with its areas; first key: fatty acid name; second key: name of the fragment; value: CgProbe peak identification object
   * @param missed fragments that were not found in the equation (required for reading of results, since for rules containing "+", not all of the fragments have to be found)
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (not containing type and OH number) cannot be decoded
   */
  public static IntensityPositionVO getFattyAcidsFromReadableRule(String storedRule, IntensityRuleVO ruleVO,
      Hashtable<String,Hashtable<String,CgProbe>> chainFragments, Hashtable<String,String> missed) throws LipidCombinameEncodingException{
    String rule = new String(storedRule);
    int derivedPosition = -1;
    if (rule.indexOf("BECAUSE NOT FULFILLED: ")!=-1){
      String derivedString = rule.substring(0,rule.indexOf("BECAUSE NOT FULFILLED: "));
      derivedString = derivedString.substring(derivedString.indexOf("[")+1,derivedString.indexOf("]"));
      derivedPosition = Integer.parseInt(derivedString);
      rule = rule.substring(rule.indexOf("BECAUSE NOT FULFILLED: ")+"BECAUSE NOT FULFILLED: ".length());
    }
    String ruleBigger = "";
    String ruleSmaller = "";
    if (rule.indexOf(">")!=-1){
      ruleBigger = rule.substring(0,rule.indexOf(">"));
      ruleSmaller = rule.substring(rule.indexOf(">")+1);
    }else if (rule.indexOf("<")!=-1){
      ruleSmaller = rule.substring(0,rule.indexOf("<"));
      ruleBigger = rule.substring(rule.indexOf("<")+1);
    }
    int biggerIndex = ruleBigger.indexOf("["+ruleVO.getBiggerPosition()+"]");
    int smallerIndex = ruleSmaller.indexOf("["+ruleVO.getSmallerPosition()+"]");
    String biggerFA = ruleBigger.substring(0,biggerIndex);
    String smallerFA = ruleSmaller.substring(0,smallerIndex);
    Vector<String> biggerNbpNamesWoPos = ruleVO.getBiggerNonBasePeakNames();
    Vector<String> smallerNbpNamesWoPos = ruleVO.getSmallerNonBasePeakNames();
    biggerFA = IntensityChainVO.extractFANames(biggerFA, biggerNbpNamesWoPos);
    smallerFA = IntensityChainVO.extractFANames(smallerFA, smallerNbpNamesWoPos);
    FattyAcidVO biggerVO = null;
    FattyAcidVO smallerVO = null;
    Vector<String> biggerNbpNamesWiPos = new Vector<String>();
    Vector<String> smallerNbpNamesWiPos = new Vector<String>();
    for (String frag : biggerNbpNamesWoPos) biggerNbpNamesWiPos.add(frag+"["+ruleVO.getBiggerPosition()+"]");
    for (String frag : smallerNbpNamesWoPos) smallerNbpNamesWiPos.add(frag+"["+ruleVO.getSmallerPosition()+"]");    
    if (biggerFA!=null)
      biggerVO = selectChainFromFoundFragments(biggerFA,biggerNbpNamesWiPos,chainFragments,missed);
    else
      biggerVO = selectChainFromFoundFragments(smallerFA,smallerNbpNamesWiPos,chainFragments,missed);    
    if (smallerFA!=null) 
      smallerVO = selectChainFromFoundFragments(smallerFA,smallerNbpNamesWiPos,chainFragments,missed);
    else
      smallerVO = selectChainFromFoundFragments(biggerFA,biggerNbpNamesWiPos,chainFragments,missed);
    IntensityPositionVO posVO = new IntensityPositionVO(ruleVO,biggerVO,smallerVO);
    if (derivedPosition>-1) posVO = createNegatedVO(posVO, derivedPosition);
    return posVO;
  }
  
  public String getRuleValueInterpretation(Hashtable<String,Float> values, Float basePeak) throws RulesException{
    String returnString = super.getRuleValueInterpretation(values, basePeak);
    if (negated_) returnString = "NOT FULFILLED: "+returnString;
    return returnString;
  }
  
  public static IntensityPositionVO createNegatedVO(IntensityPositionVO posVO, int derivedPosition){
    IntensityPositionVO vo = new IntensityPositionVO(posVO,posVO.getBiggerFA(),posVO.getSmallerFA());
    vo.negated_ = true;
    vo.derivedPosition_  = derivedPosition;
    return vo;
  }

  /**
   * @return true if it is a negated rule
   */
  public boolean isNegated()
  {
    return negated_;
  }

  
  
}
