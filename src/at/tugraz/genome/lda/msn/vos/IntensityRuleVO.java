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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.utils.RangeInteger;
import at.tugraz.genome.lda.vos.ShortStringVO;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * Value object containing all necessary information for a rule for intensity comparison
 * @author Juergen Hartler
 *
 */
public class IntensityRuleVO
{
  // identifier for the base peak in the rule file
  public final static String BASEPEAK_NAME = "$BASEPEAK";
  
  // different chain types are involved
  public final static short DIFF_CHAIN_TYPES = -1;
  
  // to which section does this rule belong to [HEAD_SECTION, CHAINS_SECTION, or POSITION_SECTION]
  private int type_;
  // the equation as it was entered in the rules file
  protected String equation_;
  // must this rule be fulfilled (e.g. to accept a head group, to accept chain fragment, or to assign a position)
  protected boolean mandatory_;
  /**is this an OR rule*/
  protected boolean orRule_;
  
  protected ExpressionForComparisonVO biggerExpression_;

  protected ExpressionForComparisonVO smallerExpression_;

  /** which type of fatty acid chains are involved; if different chain types are involved use DIFF_CHAIN_TYPES*/
  protected short chainType_;
  
  /** how many hydroxylations must be present for the detection of this fragment; key: number of hydroxylations; value: mandatory - should be null in case of no OH restrictions*/ 
  protected Hashtable<Short,Short> allowedOHs_;
  /** the available types the fragment originates of, i.e. head, acyl, alkyl, alkenyl, lcb*/
  protected Set<Short> availableFragmentTypes_;

  
  /**
   * all of the required information for the VO has to be provided in the constructor
   * @param type to which section does this rule belong to [HEAD_SECTION, CHAINS_SECTION, or POSITION_SECTION]
   * @param originalEquation the equation in its original format
   * @param biggerExpression object containing equation information of the bigger part of the equation
   * @param smallerExpression object containing equation information of the smaller part of the equation
   * @param orRule is this an OR expression
   */
  public IntensityRuleVO(int type, String originalEquation, ExpressionForComparisonVO biggerExpression, ExpressionForComparisonVO smallerExpression, boolean orRule){
    type_ = type;
    equation_ = originalEquation;
    mandatory_ = false;
    biggerExpression_ = biggerExpression;
    smallerExpression_ = smallerExpression;
    chainType_ = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
    availableFragmentTypes_ = new HashSet<Short>();
    if (biggerExpression_!=null) {
      for (FragmentMultVO frag : biggerExpression_.getFragments())
        availableFragmentTypes_.add(frag.getFragmentType());
    }
    if (smallerExpression_!=null) {
      for (FragmentMultVO frag : smallerExpression_.getFragments())
        availableFragmentTypes_.add(frag.getFragmentType());
    }
    this.orRule_ = orRule;
    checkForDiffChainTypes();
  }
  
  /**
   * creates a clone of the object
   * @param rule the IntensityRuleVO to be cloned
   */
  public IntensityRuleVO(IntensityRuleVO rule){
    this(rule.type_,rule.equation_,rule.biggerExpression_,rule.smallerExpression_,rule.orRule_);
    this.mandatory_ = rule.mandatory_;
    this.chainType_ = rule.chainType_;
  }
  
  /**
   * intensity rules do not require any identifier
   * @return an identifier containing the properties of the value object
   */
  public String getRuleIdentifier(){
    return equation_;
  }

  /**
   * @deprecated
   * @return must this rule be fulfilled (e.g. to accept a head group, to accept chain fragment, or to assign a position)
   */
  public boolean isMandatory()
  {
    return mandatory_;
  }

  
  /**
   * @param ohNumber the number of hydroxylation sites
   * @return must the fragment be present - different options possible according to the MANDATORY_... specifications in this class
   */
  public boolean isMandatory(short ohNumber)
  {
    if (ohNumber==LipidomicsConstants.EXCEL_NO_OH_INFO || this.allowedOHs_==null)  
      return mandatory_;
    return (this.allowedOHs_.get(ohNumber)==FragmentRuleVO.MANDATORY_TRUE);
  }

  /**
   * 
   * @return object containing equation information of the more intense part of the equation
   */
  public ExpressionForComparisonVO getBiggerExpression()
  {
    return biggerExpression_;
  }


  /**
   * 
   * @return object containing equation information of the less intense part of the equation
   */
  public ExpressionForComparisonVO getSmallerExpression()
  {
    return smallerExpression_;
  }
  
  /**
   * 
   * @return the affected position (use if appropriate)
   */
  public int getBiggerPosition(){
    return biggerExpression_.getPosition();
  }
  
  /**
   * 
   * @return the affected position (use if appropriate)
   */
  public int getSmallerPosition(){
    return smallerExpression_.getPosition();
  }

  
  /**
   * @return String representative of the VO - for debugging purposes
   */
  public String toString(){
    return this.getRuleIdentifier();
  }
  
  
  
  /**
   * sets if the rule VO is mandatory
   * @param mandatory
   */
  public void setMandatory(boolean mandatory)
  {
    this.mandatory_ = mandatory;
  }

  /**
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @return the human readable String representation for storing in Excel
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public String getReadableRuleInterpretation(HydroxyEncoding faEncoding, HydroxyEncoding lcbEncoding) throws LipidCombinameEncodingException {
    return new String(equation_);
  }
  
  /**
   * returns a String representing the rule, where fragment place holders are replaced by the found intensities of the fragments
   * @param values contains values for the fragments - key is the fragment name, value is the intensity of the fragment
   * @param basePeak value of the base peak - if it is required for this rule
   * @return a String representing the rule, where fragment place holders are replaced by the found intensities of the fragments
   * @throws RulesException if something does not work
   */
  public String getRuleValueInterpretation(Hashtable<String,Float> values, Float basePeak) throws RulesException{
    String returnString = new String(equation_);
    if (basePeak!=null && this.containsBasePeak()) values.put(BASEPEAK_NAME, basePeak);
    Hashtable<String,Short> head = new Hashtable<String,Short>();
    for (String frag : values.keySet())
      head.put(frag, LipidomicsConstants.CHAIN_TYPE_NO_CHAIN);
    Vector<ShortStringVO> lengthSortedNames = FragmentRuleVO.getLengthSortedFragmentNames(head,new Hashtable<String,Short>());
    for (ShortStringVO name : lengthSortedNames){
      Float value = 0f;
      if (values.containsKey(name.getKey())) value = values.get(name.getKey());
      while (returnString.indexOf(name.getKey())!=-1)
        returnString = returnString.substring(0,returnString.indexOf(name.getKey()))+String.valueOf(String.valueOf(value))+returnString.substring(returnString.indexOf(name.getKey())+name.getKey().length());
    }
    return returnString;
  }
  
  /**
   * returns the value of a base peak
   * @param ruleValueInterpretation a String where the fragment names are replaced by the found peak areas
   * @param values the values of the other found fragments; key: fragment name; value: the area of the found fragment
   * @return value of a base peak
   * @throws RulesException thrown if something does not work
   */
  public float extractBasePeakValue(String ruleValueInterpretation, Hashtable<String,Float> values) throws RulesException{
    float basePeak = 0f;
    String containingBPString = getRuleValueInterpretation(values,null);    
    String beforeBasePeak = containingBPString.substring(0,containingBPString.indexOf(BASEPEAK_NAME));
    String afterBasePeak = containingBPString.substring(containingBPString.indexOf(BASEPEAK_NAME)+BASEPEAK_NAME.length());
    String baseString = null;
    if (afterBasePeak.indexOf(BASEPEAK_NAME)!=-1){
      String comp = "";
      if (afterBasePeak.indexOf(">")!=-1){
        comp = ">";
      } else if (afterBasePeak.indexOf("<")!=-1){
        comp = "<";
      }
      afterBasePeak = afterBasePeak.substring(0,afterBasePeak.indexOf(comp));
      String subValueInterpr = ruleValueInterpretation.substring(0,ruleValueInterpretation.indexOf(comp));
      baseString = subValueInterpr.substring(beforeBasePeak.length(),subValueInterpr.length()-afterBasePeak.length());
    } else {
      baseString = ruleValueInterpretation.substring(beforeBasePeak.length(),ruleValueInterpretation.length()-afterBasePeak.length());
    }
    basePeak = Float.parseFloat(baseString);
    return basePeak;
  }

  /**
   * returns the first found fragment that is not a base peak fragment
   * @param headFragments the found head fragments
   * @param chainFragments the found chain fragments
   * @return the first found fragment that is not a base peak fragment
   */
  public CgProbe getAnyNonBasepeakFragment(Hashtable<String,CgProbe> headFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments) {
    CgProbe probe = null;
    Hashtable<String,Short> nonBpNames = getBiggerNonBasePeakNames();
    for (String nonBpName : nonBpNames.keySet()){
      probe = checkForFragmentAvailability(nonBpName,headFragments,chainFragments);
      if (probe!=null) break;
    }
    if (probe!=null) return probe;
    nonBpNames = getSmallerNonBasePeakNames();
    for (String nonBpName : nonBpNames.keySet()){
      probe = checkForFragmentAvailability(nonBpName,headFragments, chainFragments);
      if (probe!=null) break;
    }
    return probe;
  }
  
  /**
   * checks whether this fragment was detected
   * @param frag the name of the fragment
   * @param headFragments the found head fragments
   * @param chainFragments the found chain fragments
   * @return returns the idendified fragment
   */
  public CgProbe checkForFragmentAvailability(String frag, Hashtable<String,CgProbe> headFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments) {
    if (headFragments.containsKey(frag))
      return headFragments.get(frag);
    else
      return null;
  }
  
  
  /**
   * checks the objects for equality
   * @param other other VO to be checked
   * @return true if equal
   */
  public boolean equals(IntensityRuleVO other){
    if (!equation_.equalsIgnoreCase(other.equation_)) return false;
    return true;
  }

  /**
   * 
   * @return which type of fatty acid chains are involved
   */
  public short getChainType()
  {
    return chainType_;
  }


  /**
   * 
   * @return true if this equation contains a placeholder for the base peak
   */
  public boolean containsBasePeak(){
    boolean containsBasepeak = false;
    for (FragmentMultVO multVO : biggerExpression_.getFragments()){
      if (multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)){
        containsBasepeak = true;
        break;
      }
    }
    if (containsBasepeak) return containsBasepeak;
    for (FragmentMultVO multVO : smallerExpression_.getFragments()){
      if (multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)){
        containsBasepeak = true;
        break;
      }
    }
    return containsBasepeak;
  }
  
  /**
   * 
   * @return all fragment names of the bigger side of the expression that are not the base peak
   */
  public Hashtable<String,Short> getBiggerNonBasePeakNames(){
    Hashtable<String,Short> allNonBasepeakNames = new Hashtable<String,Short>();
    for (FragmentMultVO multVO : biggerExpression_.getFragments()){
      if (!multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)){
        String nonBpName = multVO.getFragmentName();
        allNonBasepeakNames.put(nonBpName, multVO.getFragmentType());
      }
    }
    return allNonBasepeakNames;
  }
    
  
  /**
   * 
   * @return all fragment names of the smaller side of the expression that are not the base peak
   */  
  public Hashtable<String,Short> getSmallerNonBasePeakNames(){
    Hashtable<String,Short> allNonBasepeakNames = new Hashtable<String,Short>();
    for (FragmentMultVO multVO : smallerExpression_.getFragments()){
      if (!multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)){
        String nonBpName = multVO.getFragmentName();
        allNonBasepeakNames.put(nonBpName, multVO.getFragmentType());
      }
    }    
    return allNonBasepeakNames;
  }
  
  /**
   * 
   * @return the fragment names on the greater than end of the equation which are whether base peak names nor head group fragments
   */
  protected Hashtable<String,Short> getBiggerNonHeadAndBasePeakNames() {
    return getNonHeadAndBasePeakNames(biggerExpression_.getFragments());
  }
  
  /**
   * 
   * @return the fragment names on the smaller than end of the equation which are whether base peak names nor head group fragments
   */
  protected Hashtable<String,Short> getSmallerNonHeadAndBasePeakNames() {
    return getNonHeadAndBasePeakNames(smallerExpression_.getFragments());
  }

  /**
   * returns all fragment names which are whether base peak names nor head group fragments
   * @param parts the fragment parts
   * @return fragment names which are whether base peak names nor head group fragments; key: fragment name; value: fragment type
   */
  private Hashtable<String,Short> getNonHeadAndBasePeakNames(Vector<FragmentMultVO> parts) {
    Hashtable<String,Short> filteredNames = new Hashtable<String,Short>();
    for (FragmentMultVO multVO : parts){
      if (!multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME) && multVO.getFragmentType()!=LipidomicsConstants.CHAIN_TYPE_NO_CHAIN){
        String name = multVO.getFragmentName();
        filteredNames.put(name, multVO.getFragmentType());
      }
    }
    return filteredNames;
  }
  

  
  /**
   * 
   * @return the first name of a fragment in the equation that is not the base peak
   */
  public String getAnyNonBasePeakName(){
    String nonBpName = "";
    for (FragmentMultVO multVO : biggerExpression_.getFragments()){
      if (!multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)){
        nonBpName = multVO.getFragmentName();
        break;
      }
    }
    if (nonBpName.length()>0) return nonBpName;
    for (FragmentMultVO multVO : smallerExpression_.getFragments()){
      if (!multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)){
        nonBpName = multVO.getFragmentName();
        break;
      }
    }    
    return nonBpName;
  }
  
  /**
   * checks if there where enough fragments found to evaluate the equation
   * @param found the found fragments; key: fragment name; value the peak identification object
   * @return true if enough fragments are found to evaluate this equation
   */
  public boolean enoughFragmentsForRuleEvaluationFound(Hashtable<String,CgProbe> found){
    boolean biggerFragsFound = areFragmentsFound(biggerExpression_,found);
    boolean smallerFragsFound = areFragmentsFound(smallerExpression_,found);
    return (biggerFragsFound && smallerFragsFound);
  }
  
  /**
   * checks if there where enough fragments found to evaluate the equation
   * @param headGroupFragments the detected head group fragments
   * @param chainFragments the detected chain fragments
   * @param chainsToCheck the chains for which this evaluation has to be performed
   * @return true when there were enough fragments found to evaluate the rule
   */
  public boolean enoughFragmentsForRuleEvaluationFound(Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments,
      Vector<FattyAcidVO> chainsToCheck) {
    Hashtable<String,CgProbe> foundFragments = getFragmentsForExpression(biggerExpression_, headGroupFragments, chainFragments,
        chainsToCheck);
    boolean biggerFragsFound = areFragmentsFound(biggerExpression_,foundFragments);
    foundFragments = getFragmentsForExpression(smallerExpression_, headGroupFragments, chainFragments,
        chainsToCheck);
    boolean smallerFragsFound = areFragmentsFound(smallerExpression_,foundFragments);
    return (biggerFragsFound && smallerFragsFound);
  }
  
  /**
   * returns all of the fragments that belong to a certain equation expression
   * @param expression the equation expression
   * @param headGroupFragments the detected head group fragments
   * @param chainFragments the detected chain fragments
   * @param chainsToCheck the chains for which this evaluation has to be performed
   * @return all of the fragments that belong to a certain equation expression; key: fragment name; value: the detected fragment
   */
  private Hashtable<String,CgProbe> getFragmentsForExpression(ExpressionForComparisonVO expression, Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments,
      Vector<FattyAcidVO> chainsToCheck) {
    Hashtable<String,CgProbe> foundFragments = new Hashtable<String,CgProbe>(headGroupFragments);
    for (FragmentMultVO fragVO : expression.getFragments()) {
      for (FattyAcidVO chain : chainsToCheck) {
        if (chain.getChainType()!=fragVO.getFragmentType() || !chainFragments.containsKey(chain.getChainId()))
          continue;
        foundFragments.putAll(chainFragments.get(chain.getChainId()));
      }
    }
    return foundFragments;
  }
  
  /**
   * checks if the any of the found fragments occur in one part of the equation comparator
   * @param expression object containing all information on one side of the comparator
   * @param found the found fragments; key: fragment name; value the peak identification object
   * @return true if any of the found fragments are in this equation part
   */
  private boolean areFragmentsFound(ExpressionForComparisonVO expression, Hashtable<String,CgProbe> found){
    boolean foundFrags = false;
    for (FragmentMultVO multVO : expression.getFragments()){
      String name = multVO.getFragmentName();
      if (name.equalsIgnoreCase(BASEPEAK_NAME) || (found!=null && found.containsKey(name) && multVO.isPositive())){
        foundFrags = true;
        break;
      }
    }
    return foundFrags;
  }
  
  /**
   * checks if the rule is OK
   * @param found the found fragments; key: fragment name; value the peak identification object
   * @param basePeak the base peak area
   * @return true if the rule is fulfilled
   */
  public boolean isRuleFulfilled(Hashtable<String,CgProbe> found, Float basePeak){
    if (this.orRule_) {
      boolean oneFound = false;
      for (FragmentMultVO mult : this.biggerExpression_.getFragments()) {
        if (found.containsKey(mult.getFragmentName())) {
          oneFound = true;
        }
      }
      return oneFound;
    } else {
      double biggerArea = biggerExpression_.evaluateExpression(found, basePeak);
      double smallerArea = smallerExpression_.evaluateExpression(found, basePeak);
      if (biggerArea>smallerArea) return true;
      else return false;
    }
  }
  
  /**
   * is this intensity rule fulfilled or violated
   * @param headGroupFragments the detected head group fragments
   * @param chainFragments the detected chain fragments
   * @param chainsToCheck the chains for which this evaluation has to be performed
   * @param basePeak intensity of the base peak
   * @return true when this rule is fulfilled
   */
  public boolean isRuleFulfilled(Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments,
      Vector<FattyAcidVO> chainsToCheck, Float basePeak) {
    Hashtable<String,CgProbe> foundFragments = getFragmentsForExpression(biggerExpression_, headGroupFragments, chainFragments,
        chainsToCheck);
    double biggerArea = getBiggerArea(foundFragments,basePeak);
    foundFragments = getFragmentsForExpression(smallerExpression_, headGroupFragments, chainFragments,
        chainsToCheck);
    double smallerArea = getSmallerArea(foundFragments,basePeak);
    if (biggerArea>smallerArea) return true;
    else return false; 
  }
  
  /**
   * checks if the rule contains chains of different types
   */
  public void checkForDiffChainTypes(){
    boolean foundOneType = false;
    chainType_ = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
    for (FragmentMultVO multVO : biggerExpression_.getFragments()){
      if (multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)) continue;
      if (!foundOneType){
        chainType_ = multVO.getFragmentType();
        foundOneType = true;
      } else {
        if (chainType_!=multVO.getFragmentType()) {
          chainType_ = DIFF_CHAIN_TYPES;
          break;
        }
      }
    }
    if (chainType_ == DIFF_CHAIN_TYPES) return;
    for (FragmentMultVO multVO : smallerExpression_.getFragments()){
      if (multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)) continue;
      if (!foundOneType){
        chainType_ = multVO.getFragmentType();
        foundOneType = true;
      } else {
        if (chainType_!=multVO.getFragmentType()) {
          chainType_ = DIFF_CHAIN_TYPES;
          break;
        }
      }
    }
  }
  
  /**
   * 
   * @param found set of found fragments; key: fragment name; value the peak identification object
   * @param basePeak the base peak area
   * @return the evaluated expression at the bigger side of the expression
   */
  public double getBiggerArea(Hashtable<String,CgProbe> found, Float basePeak){
    return biggerExpression_.evaluateExpression(found, basePeak);
  }

  /**
   * 
   * @param found set of found fragments; key: fragment name; value the peak identification object
   * @param basePeak the base peak area
   * @return the evaluated expression at the smaller side of the expression
   */
  public double getSmallerArea(Hashtable<String,CgProbe> found, Float basePeak){
    return smallerExpression_.evaluateExpression(found, basePeak);
  }
  
  /**
   * returns the MSLevel this equation affects
   * @param headGroupFragments the found head group fragments with its areas; key: fragment name; value: CgProbe peak identification object
   * @param chainFragments the found chain fragments with its areas; first key: fatty acid name; second key: name of the fragment; value: CgProbe peak identification object
   * @return
   */
  public int getMSLevel(Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments){
    int msLevel = -1;
    CgProbe probe =  getAnyNonBasepeakFragment(headGroupFragments, chainFragments);
    if (probe!=null)
      msLevel = probe.getMsLevel();
    return msLevel;
  }
  
  /**
   * 
   * @param chainType which type of fatty acid chains are involved
   */
  public void setChainType(short chainType)
  {
    this.chainType_ = chainType;
  }
  
  /**
   * checks if this intensity rule VO contains this fragment
   * @param fragmentName the name of the fragment
   * @return true if this fragment is used in the formula
   */
  public boolean containsFragment(String fragmentName){
    return (biggerExpression_.containsFragment(fragmentName) || smallerExpression_.containsFragment(fragmentName));
  }
  
  /**
   * splits an equation to the parts on the bigger and the smaller part of the equation sign
   * @param rule the whole equation
   * @return String[2]: the first part of the array is the bigger, and the  second part the smaller part of the equation
   */
  protected static String[] splitToBiggerAndSmallerPart(String rule){
    String biggerPart = null;
    String smallerPart = null;
    if (rule.indexOf(">")!=-1){
      biggerPart = rule.substring(0,rule.indexOf(">"));
      smallerPart = rule.substring(rule.indexOf(">")+1);
    } else if (rule.indexOf("<")!=-1){
      smallerPart = rule.substring(0,rule.indexOf("<"));
      biggerPart = rule.substring(rule.indexOf("<")+1);
    }
    String[] biggerSmaller = new String[2];
    biggerSmaller[0] = biggerPart;
    biggerSmaller[1] = smallerPart;
    return biggerSmaller;
  }
  
  /**
   * replaces parts of an equation by another String
   * @param rulePart the part of the equation
   * @param originalToReplacement a hashtable containing the string to be replaced and its replacement
   * @return
   */
  protected String replaceParts(String rulePart, Hashtable<String,String> originalToReplacement){
    Vector<RangeInteger> ranges = new Vector<RangeInteger>();
    Hashtable<Integer,String> replacementIndices = new Hashtable<Integer,String>();
    //TODO: the fragmentNames was artificially created since only the name seems to be of interest - might need to be improved in future
    Hashtable<String,Short> fragmentNames = new Hashtable<String,Short>();
    for (String name : originalToReplacement.keySet())
      fragmentNames.put(name, LipidomicsConstants.CHAIN_TYPE_NO_CHAIN);
    Vector<ShortStringVO> lengthSortedNames = FragmentRuleVO.getLengthSortedFragmentNames(fragmentNames,new Hashtable<String,Short>());
    for (ShortStringVO name : lengthSortedNames){
      int startIndex = getStartIndex(name.getKey(),rulePart,ranges);
      ranges.add(new RangeInteger(startIndex,startIndex+name.getKey().length()));
      replacementIndices.put(startIndex, name.getKey());
    }
    String returnString = new String(rulePart);
    List<Integer> idxSorted = new ArrayList<Integer>(replacementIndices.keySet());
    Collections.sort(idxSorted);
    for (int i=(idxSorted.size()-1); i!=-1; i--){
      int idx = idxSorted.get(i);
      String toReplace =  replacementIndices.get(idx);
      String replacement = originalToReplacement.get(toReplace);
      returnString = returnString.substring(0,idx)+replacement+returnString.substring(idx+toReplace.length());
    }
    return returnString;

  }
  
  /**
   * detects the start index for an replacement, and takes car that already found replacements are not covered
   * @param name the String to be replaced
   * @param equation the equation part
   * @param ranges the ranges that are already covered by other replacements
   * @return
   */
  public static int getStartIndex(String name, String equation, Vector<RangeInteger> ranges){
    String subString = new String(equation);
    int removed = 0;
    int idx;
    int start = -1;
    while ((idx = subString.indexOf(name))>-1){
      start = idx+removed;
      int stop = start+name.length();
      boolean found = false;
      for (RangeInteger range : ranges){
        if (range.insideRange(start) && range.getStop()!=start)
          found = true;
        if (range.insideRange(stop) && range.getStart()!=stop)
          found = true;
        if (found) break;
      }
      if (!found){        
        break;
      } else {
        removed = stop;
        subString = subString.substring(idx+name.length());
      }
    }
    return start;
  }

 /**
  * 
  * @return true if there any absolute relationships (comparison to the base peak) are involved
  */
  public boolean isAbsoluteComparison(){
    return (this.biggerExpression_.isAbsoluteComparison() || this.smallerExpression_.isAbsoluteComparison());
  }
  
  
  /**
   * sets how many hydroxylations must be present for the detection of this fragment; key: number of hydroxylations; value: mandatory - should be null in case of no OH restrictions
   * @param allowedOHs
   */
  public void setAllowedOHs(Hashtable<Short,Short> allowedOHs) {
    this.allowedOHs_ = allowedOHs;
  }
  
  
  /**
   * checks whether this fragment is applicable for this OH configuration
   * @param ohNumber the number of hydroxylation sites
   * @return true when the fragment is applicable
   */
  public boolean hydroxylationValid(short ohNumber) {
    if (ohNumber==LipidomicsConstants.EXCEL_NO_OH_INFO || this.allowedOHs_==null)
      return true;
    return (this.allowedOHs_.containsKey(ohNumber));
  }
  
  /**
   * 
   * @return the available types the fragment originates of, i.e. head, acyl, alkyl, alkenyl, lcb
   */
  public Set<Short> getAvailableTypes(){
    return availableFragmentTypes_;    
  }

  /**
   * 
   * @return whether this is an OR rule
   */
  public boolean isOrRule()
  {
    return orRule_;
  }
  
}
