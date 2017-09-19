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
import java.util.Vector;

import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.utils.RangeInteger;
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
  public final static int DIFF_CHAIN_TYPES = -1;
  
  // to which section does this rule belong to [HEAD_SECTION, CHAINS_SECTION, or POSITION_SECTION]
  private int type_;
  // the equation as it was entered in the rules file
  protected String equation_;
  // must this rule be fulfilled (e.g. to accept a head group, to accept chain fragment, or to assign a position)
  protected boolean mandatory_;
  
  protected ExpressionForComparisonVO biggerExpression_;

  protected ExpressionForComparisonVO smallerExpression_;

  /** which type of fatty acid chains are involved; if different chain types are involved use DIFF_CHAIN_TYPES*/
  protected int chainType_;

  
  /**
   * all of the required information for the VO has to be provided in the constructor
   * @param type to which section does this rule belong to [HEAD_SECTION, CHAINS_SECTION, or POSITION_SECTION]
   * @param originalEquation the equation in its original format
   * @param biggerExpression object containing equation information of the bigger part of the equation
   * @param smallerExpression object containing equation information of the smaller part of the equation
   */
  public IntensityRuleVO(int type, String originalEquation, ExpressionForComparisonVO biggerExpression, ExpressionForComparisonVO smallerExpression){
    type_ = type;
    equation_ = originalEquation;
    mandatory_ = false;
    biggerExpression_ = biggerExpression;
    smallerExpression_ = smallerExpression;
    chainType_ = FragmentRuleVO.ACYL_CHAIN;
  }
  
  /**
   * creates a clone of the object
   * @param rule the IntensityRuleVO to be cloned
   */
  public IntensityRuleVO(IntensityRuleVO rule){
    this(rule.type_,rule.equation_,rule.biggerExpression_,rule.smallerExpression_);
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
   * 
   * @return must this rule be fulfilled (e.g. to accept a head group, to accept chain fragment, or to assign a position)
   */
  public boolean isMandatory()
  {
    return mandatory_;
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
   * 
   * @return the human readable String representation
   */
  public String getReadableRuleInterpretation() {
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
    Vector<String> lengthSortedNames = FragmentRuleVO.getLengthSortedFragmentNames(values.keySet(),new HashSet<String>());
    for (String name : lengthSortedNames){
      Float value = 0f;
      if (values.containsKey(name)) value = values.get(name);
      while (returnString.indexOf(name)!=-1)
        returnString = returnString.substring(0,returnString.indexOf(name))+String.valueOf(String.valueOf(value))+returnString.substring(returnString.indexOf(name)+name.length());
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
   * 
   * @return the fatty acid name of the more intense side of the equation (if appropriate)
   */
  public String getBiggerFA()
  {
    return null;
  }

  /**
   * 
   * @return the fatty acid name of the less intense side of the equation (if appropriate)
   */
  public String getSmallerFA()
  {
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
  public int getChainType()
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
  public Vector<String> getBiggerNonBasePeakNames(){
    Hashtable<String,String> allNonBasepeakNames = new Hashtable<String,String>();
    for (FragmentMultVO multVO : biggerExpression_.getFragments()){
      if (!multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)){
        String nonBpName = multVO.getFragmentName();
        allNonBasepeakNames.put(nonBpName, nonBpName);
      }
    }
    return new Vector<String>(allNonBasepeakNames.values());
  }
  
  /**
   * 
   * @return all fragment names of the smaller side of the expression that are not the base peak
   */  
  public Vector<String> getSmallerNonBasePeakNames(){
    Hashtable<String,String> allNonBasepeakNames = new Hashtable<String,String>();
    for (FragmentMultVO multVO : smallerExpression_.getFragments()){
      if (!multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)){
        String nonBpName = multVO.getFragmentName();
        allNonBasepeakNames.put(nonBpName, nonBpName);
      }
    }    
    return new Vector<String>(allNonBasepeakNames.values());
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
   * checks if there where enough fragments found to evaluate the equation, and if the bigger fragments are out of the first group (the latter is for chain fragments)
   * @param found1 first set of found fragments; key: fragment name; value the peak identification object
   * @param found2 second set of found fragments; key: fragment name; value the peak identification object
   * @return boolean[2]: boolean[0] = true if enough fragments are found to evaluate this equation; boolean[0] = true if the bigger fragments are in the first group
   */
  public boolean[] enoughFragmentsForRuleEvaluationFound(Hashtable<String,CgProbe> found1, Hashtable<String,CgProbe> found2){
    boolean biggerFragsFound = areFragmentsFound(biggerExpression_,found1);
    boolean biggerFirst = false;
    boolean smallerFragsFound = false;
    if (biggerFragsFound){
      smallerFragsFound = areFragmentsFound(smallerExpression_,found2);
      biggerFirst = true;
    } else{
      biggerFragsFound = areFragmentsFound(biggerExpression_,found2);
      smallerFragsFound = areFragmentsFound(smallerExpression_,found1);
    }
    boolean enough = (biggerFragsFound && smallerFragsFound);
    boolean[] result = new boolean[2];
    result[0] = enough;
    result[1] = biggerFirst;
    return result;
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
    double biggerArea = biggerExpression_.evaluateExpression(found, basePeak);
    double smallerArea = smallerExpression_.evaluateExpression(found, basePeak);
    if (biggerArea>smallerArea) return true;
    else return false;
  }
  
 /**
  * checks if the rule is OK (for chain fragment comparison)
  * @param found1 first set of found fragments; key: fragment name; value the peak identification object
  * @param found2 second set of found fragments; key: fragment name; value the peak identification object
  * @param basePeak the base peak area
  * @return true if the rule is fulfilled
  */
  public boolean isRuleFulfilled(Hashtable<String,CgProbe> found1, Hashtable<String,CgProbe> found2, Float basePeak){
    boolean biggerFragsFound = areFragmentsFound(biggerExpression_,found1);
    double biggerArea = 0d;
    double smallerArea = 0d;
    if (biggerFragsFound){
      biggerArea = getBiggerArea(found1, basePeak);
      smallerArea = getSmallerArea(found2, basePeak);
    } else {
      biggerArea = getBiggerArea(found2, basePeak);
      smallerArea = getSmallerArea(found1, basePeak);      
    }    
    if (biggerArea>smallerArea) return true;
    else return false;
  }
  
  /**
   * checks if the rule contains chains of different types
   * @param fragRules the fragmentation rules; key: name of the fragment; value: fragmentation rule
   */
  public void checkForDiffChainTypes(Hashtable<String,FragmentRuleVO> fragRules){
    boolean foundOneType = false;
    chainType_ = FragmentRuleVO.ACYL_CHAIN;
    for (FragmentMultVO multVO : biggerExpression_.getFragments()){
      if (multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)) continue;
      if (!fragRules.containsKey(multVO.getFragmentName())) continue;
      FragmentRuleVO ruleVO = fragRules.get(multVO.getFragmentName());
      if (!foundOneType){
        chainType_ = ruleVO.getChainType();
        foundOneType = true;
      } else {
        if (chainType_!=ruleVO.getChainType()){
          chainType_ = DIFF_CHAIN_TYPES;
          break;
        }
      }
    }
    if (chainType_ == DIFF_CHAIN_TYPES) return;
    for (FragmentMultVO multVO : smallerExpression_.getFragments()){
      if (multVO.getFragmentName().equalsIgnoreCase(BASEPEAK_NAME)) continue;
      if (!fragRules.containsKey(multVO.getFragmentName())) continue;
      FragmentRuleVO ruleVO = fragRules.get(multVO.getFragmentName());
      if (!foundOneType){
        chainType_ = ruleVO.getChainType();
        foundOneType = true;
      } else {
        if (chainType_!=ruleVO.getChainType()){
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
   * @param biggerFA the fatty acid of the more intense part (if appropriate)
   * @param smallerFA the fatty acid of the less intense part (if appropriate)
   * @return
   */
  public int getMSLevel(Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments, String biggerFA, String smallerFA){
    int msLevel = -1;
    for (String name : getBiggerNonBasePeakNames()){
      if (headGroupFragments.containsKey(name)){
        msLevel = headGroupFragments.get(name).getMsLevel();
        break;
      } else if (biggerFA!=null && chainFragments.containsKey(biggerFA) && chainFragments.get(biggerFA).containsKey(name)){
        msLevel = chainFragments.get(biggerFA).get(name).getMsLevel();
        break;
      }
    }
    if (msLevel>-1) return msLevel;
    for (String name : getSmallerNonBasePeakNames()){
      if (headGroupFragments.containsKey(name)){
        msLevel = headGroupFragments.get(name).getMsLevel();
        break;
      } else if (smallerFA!=null && chainFragments.containsKey(smallerFA) && chainFragments.get(smallerFA).containsKey(name)){
        msLevel = chainFragments.get(smallerFA).get(name).getMsLevel();
        break;
      }
    }
    
    return msLevel;
  }
  
  /**
   * 
   * @param chainType which type of fatty acid chains are involved
   */
  public void setChainType(int chainType)
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
    Vector<String> lengthSortedNames = FragmentRuleVO.getLengthSortedFragmentNames(originalToReplacement.keySet(),new HashSet<String>());
    for (String name : lengthSortedNames){
      int startIndex = getStartIndex(name,rulePart,ranges);
      ranges.add(new RangeInteger(startIndex,startIndex+name.length()));
      replacementIndices.put(startIndex, name);
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
  
}
