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
import java.util.Objects;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.ShortStringVO;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * Value object containing all necessary information for a rule for intensity comparison
 * plus the two fatty acid that where compared by position finding
 * @author Juergen Hartler
 *
 */
public class IntensityPositionVO extends IntensityRuleVO
{
  /** the chains at the 'greater than' side of this rule*/
  private Hashtable<String,FattyAcidVO> biggerChains_;
  /** the chains at the 'smaller than' side of this rule*/
  private Hashtable<String,FattyAcidVO> smallerChains_;
  /** does this object hold hydroxylation information*/
  private boolean hasOhInfo_;

  private boolean negated_;
  private int derivedPosition_;
  
  /** consists the bigger expression of the equation of missed fragments only*/
  private boolean biggerOnlyMissed_;
  /** consists the smaller expression of the equation of missed fragments only*/
  private boolean smallerOnlyMissed_;
  
  /**
   * 
   * @param rule the corresponding intensity rule VO
   * @param biggerFA the FA on the 'greater than' side of this rule
   * @param smallerFA the FA on the 'smaller than' side of this rule
   * @param hasOhInfo does this object hold hydroxylation information
   * @param biggerOnlyMissed  consists the bigger expression of the equation of missed fragments only
   * @param smallerOnlyMissed consists the smaller expression of the equation of missed fragments only
   */
  public IntensityPositionVO(IntensityRuleVO rule, FattyAcidVO biggerFA, FattyAcidVO smallerFA, boolean hasOhInfo, boolean biggerOnlyMissed, boolean smallerOnlyMissed){
    super(rule);
    this.hasOhInfo_ = hasOhInfo;
    biggerChains_ = new Hashtable<String,FattyAcidVO>();
    smallerChains_ = new Hashtable<String,FattyAcidVO>();
    for (FragmentMultVO frag : biggerExpression_.getFragments()) {
      if (frag.getFragmentType()!=LipidomicsConstants.CHAIN_TYPE_NO_CHAIN)
        biggerChains_.put(frag.getFragmentName(), biggerFA);
    }
    for (FragmentMultVO frag : smallerExpression_.getFragments()) {
      if (frag.getFragmentType()!=LipidomicsConstants.CHAIN_TYPE_NO_CHAIN)
        smallerChains_.put(frag.getFragmentName(), smallerFA);
    }
    negated_ = false;
    derivedPosition_ = -1;
    biggerOnlyMissed_ = false;
    smallerOnlyMissed_ = false;
    biggerOnlyMissed_ = biggerOnlyMissed;
    smallerOnlyMissed_ = smallerOnlyMissed;
  }


  /**
   * @return fatty acid that was verified at at the greater part of the comparator
   */
  public FattyAcidVO getBiggerFA()
  {
    if (biggerChains_!=null && biggerChains_.size()>0)
      return biggerChains_.values().iterator().next();
    else
      return null;
  }

  /**
   * @return fatty acid that was verified at at the lesser part of the comparator
   */
  public FattyAcidVO getSmallerFA()
  {
    if (smallerChains_!=null && smallerChains_.size()>0)
      return smallerChains_.values().iterator().next();
    else
      return null;
  }
  
  /**
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @return the human readable String representation for storing in Excel
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public String getReadableRuleInterpretation(HydroxyEncoding faEncoding, HydroxyEncoding lcbEncoding) throws LipidCombinameEncodingException {
    Hashtable<String,String> originalToReplacementBigger = new Hashtable<String,String>();
    for (String name : getBiggerNonBasePeakNames().keySet()){
      String nameString = name+"["+getBiggerPosition()+"]";
      String value = StaticUtils.getChainFragmentDisplayName(name,StaticUtils.getHumanReadableChainName(this.biggerChains_.get(name), faEncoding,
          lcbEncoding,hasOhInfo_,StaticUtils.isAnAlex123Fragment(name)))+"["+getBiggerPosition()+"]";
      originalToReplacementBigger.put(nameString, value);
    }
    Hashtable<String,String> originalToReplacementSmaller = new Hashtable<String,String>();
    for (String name : getSmallerNonBasePeakNames().keySet()){
      String nameString = name+"["+getSmallerPosition()+"]";
      String value = StaticUtils.getChainFragmentDisplayName(name,StaticUtils.getHumanReadableChainName(this.smallerChains_.get(name), faEncoding,
          lcbEncoding,hasOhInfo_,StaticUtils.isAnAlex123Fragment(name)))+"["+getSmallerPosition()+"]";
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
    if (negated_) returnString = getBiggerFA().getCarbonDbsId()+"["+derivedPosition_+"] BECAUSE NOT FULFILLED: "+returnString;
    return returnString;
  }
  
  /**
   * checks the objects for equality
   */
  public boolean equals(IntensityRuleVO other){
    if (!super.equals(other)) return false;
    if (!(other instanceof IntensityPositionVO)) return false;
    IntensityPositionVO chain = (IntensityPositionVO)other;
    if ((getBiggerFA()==null && chain.getBiggerFA()!=null)||(chain.getBiggerFA()==null && getBiggerFA()!=null))
      return false;
    if (chain.getBiggerFA()!=null && getBiggerFA()!=null) {
      if (!getBiggerFA().getChainId().equalsIgnoreCase(chain.getBiggerFA().getChainId())) return false;
    }
    
    if ((getSmallerFA()==null && chain.getSmallerFA()!=null)||(chain.getSmallerFA()==null && getSmallerFA()!=null))
      return false;
    if (chain.getSmallerFA()!=null && getSmallerFA()!=null) {
      if (!getSmallerFA().getChainId().equalsIgnoreCase(chain.getSmallerFA().getChainId())) return false;
    }
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
    FattyAcidVO faVO;
    if ((faVO = getBiggerFA())!=null && !this.biggerOnlyMissed_ && faVO.getChainId().equalsIgnoreCase(fa)) position = getBiggerPosition();
    else if ((faVO = getSmallerFA())!=null && !this.smallerOnlyMissed_ && faVO.getChainId().equalsIgnoreCase(fa)) position = getSmallerPosition();
    return position;
  }

  /**
   * inverse method of getReadableRuleInterpretation() - fetches the fatty acids from this interpretation
   * @param storedRule the readable rule interpretation
   * @param ruleVO rule that is required to form an IntensityPositionVO
   * @param chainFragments the found chain fragments with its areas; first key: fatty acid name; second key: name of the fragment; value: CgProbe peak identification object
   * @param missed fragments that were not found in the equation (required for reading of results, since for rules containing "+", not all of the fragments have to be found)
   * @param faHydroxyEncoding the OH encodings of the FA moiety
   * @param lcbHydroxyEncoding the OH encodings of the LCB moiety
   * @param hasOhInfo has this intensity rule chain information - this is important for generating the human readable name
   * @param usedAlexMSnTargets is this result from an ALEX123 target list
   * @return the intensity rule VO including the relevant chains
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (not containing type and OH number) cannot be decoded
   */
  public static IntensityPositionVO getFattyAcidsFromReadableRule(String storedRule, IntensityRuleVO ruleVO,
      Hashtable<String,Hashtable<String,CgProbe>> chainFragments, Hashtable<String,Short> missed, HydroxyEncoding faHydroxyEncoding,
      HydroxyEncoding lcbHydroxyEncoding, boolean hasOhInfo, boolean usedAlexMSnTargets) throws LipidCombinameEncodingException{
    String rule = new String(storedRule);
    boolean biggerOnlyMissed = false;
    boolean smallerOnlyMissed = false;
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
    
    Vector<ShortStringVO> biggerNbpNames =  FragmentRuleVO.getLengthSortedFragmentNames(new Hashtable<String,Short>(),ruleVO.getBiggerNonHeadAndBasePeakNames(),new Hashtable<String,Short>());
    Vector<ShortStringVO> biggerPosNames =  FragmentRuleVO.getLengthSortedFragmentNames(new Hashtable<String,Short>(),getNonHeadAndBasePeakNamesPlusPos(ruleVO.biggerExpression_.getFragments()),new Hashtable<String,Short>());
    //the checkIfOnlyMissedValuesArePresent had been introduced because it was not possible to differentiate alkenylated species with same carbon dbs numbers; e.g. P-PE P-16:0/16:0
    if (checkIfOnlyMissedValuesArePresent(biggerPosNames,missed))
      biggerOnlyMissed = true;
    Vector<ShortStringVO> smallerNbpNames = FragmentRuleVO.getLengthSortedFragmentNames(new Hashtable<String,Short>(),ruleVO.getSmallerNonHeadAndBasePeakNames(),new Hashtable<String,Short>());
    Vector<ShortStringVO> smallerPosNames =  FragmentRuleVO.getLengthSortedFragmentNames(new Hashtable<String,Short>(),getNonHeadAndBasePeakNamesPlusPos(ruleVO.smallerExpression_.getFragments()),new Hashtable<String,Short>());
    //the checkIfOnlyMissedValuesArePresent had been introduced because it was not possible to differentiate alkenylated species with same carbon dbs numbers; e.g. P-PE P-16:0/16:0
    if (checkIfOnlyMissedValuesArePresent(smallerPosNames,missed))
      smallerOnlyMissed = true;
    Hashtable<String,FattyAcidVO> biggerChains = IntensityChainVO.extractFANames(ruleBigger, biggerNbpNames, faHydroxyEncoding, lcbHydroxyEncoding, usedAlexMSnTargets);
    Hashtable<String,FattyAcidVO> smallerChains = IntensityChainVO.extractFANames(ruleSmaller, smallerNbpNames, faHydroxyEncoding, lcbHydroxyEncoding, usedAlexMSnTargets);

    FattyAcidVO biggerVO = null;
    FattyAcidVO smallerVO = null;
    if (biggerChains.size()>0)
      biggerVO = biggerChains.values().iterator().next();
    else
      biggerVO = smallerChains.values().iterator().next();    
    if (smallerChains.size()>0)
      smallerVO = smallerChains.values().iterator().next();
    else
      smallerVO = biggerChains.values().iterator().next();
    IntensityPositionVO posVO = new IntensityPositionVO(ruleVO,biggerVO,smallerVO,hasOhInfo,biggerOnlyMissed,smallerOnlyMissed);
    if (derivedPosition>-1) posVO = createNegatedVO(posVO, derivedPosition);
    return posVO;
  }
  
  /**
   * this method checks whether one side of an equation consists of missed fragments only
   * @param names the names of the fragments in form of ShortStringVO
   * @param missed hash containing the names of missed fragments
   * @return true of this side of the equation consists of missed fragments only
   */
  private static boolean checkIfOnlyMissedValuesArePresent(Vector<ShortStringVO> names, Hashtable<String,Short> missed) {
    boolean onlyMissed = false;
    if (names.size()>0) 
      onlyMissed = true;
    for (ShortStringVO vo : names) {
      if (!missed.containsKey(vo.getKey()))
        onlyMissed = false;
    }
    return onlyMissed;
  }
  
  public String getRuleValueInterpretation(Hashtable<String,Float> values, Float basePeak) throws RulesException{
    String returnString = super.getRuleValueInterpretation(values, basePeak);
    if (negated_) returnString = "NOT FULFILLED: "+returnString;
    return returnString;
  }
  
  public static IntensityPositionVO createNegatedVO(IntensityPositionVO posVO, int derivedPosition){
    IntensityPositionVO vo = new IntensityPositionVO(posVO,posVO.getBiggerFA(),posVO.getSmallerFA(),posVO.hasOhInfo_,posVO.biggerOnlyMissed_,posVO.smallerOnlyMissed_);
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

  public CgProbe checkForFragmentAvailability(String frag, Hashtable<String,CgProbe> headFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments) {
    CgProbe probe = super.checkForFragmentAvailability(frag, headFragments, chainFragments);
    if (probe==null) {
      String chainId = null;
      if (biggerChains_.containsKey(frag))
        chainId = biggerChains_.get(frag).getChainId();
      else if (smallerChains_.containsKey(frag))
        chainId = smallerChains_.get(frag).getChainId();
      if (chainId!=null && chainFragments.containsKey(chainId) && chainFragments.get(chainId).containsKey(frag))
        probe = chainFragments.get(chainId).get(frag);
    }
    return probe;
  }
  
  
  @Override
  public boolean equals(Object obj)
  {
    if (this == obj)
      return true;
    if (!super.equals(obj))
      return false;
    if (getClass() != obj.getClass())
      return false;
    IntensityPositionVO other = (IntensityPositionVO) obj;
    return Objects.equals(biggerChains_, other.biggerChains_)
        && biggerOnlyMissed_ == other.biggerOnlyMissed_
        && derivedPosition_ == other.derivedPosition_
        && hasOhInfo_ == other.hasOhInfo_ && negated_ == other.negated_
        && Objects.equals(smallerChains_, other.smallerChains_)
        && smallerOnlyMissed_ == other.smallerOnlyMissed_;
  }
  
  
  /**
   * returns all fragment names which are whether base peak names nor head group fragments
   * @param parts the fragment parts
   * @return fragment names which are whether base peak names nor head group fragments; key: fragment name; value: fragment type
   */
  private static Hashtable<String,Short> getNonHeadAndBasePeakNamesPlusPos(Vector<FragmentMultVO> parts) {
    Hashtable<String,Short> filteredNames = new Hashtable<String,Short>();
    for (FragmentMultVO multVO : parts){
      if (!multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME) && multVO.getFragmentType()!=LipidomicsConstants.CHAIN_TYPE_NO_CHAIN){
        String name = multVO.getFragmentName();
        if (multVO.getPosition()>0)
          name+="["+multVO.getPosition()+"]";
        filteredNames.put(name, multVO.getFragmentType());
      }
    }
    return filteredNames;
  }  
  
}
