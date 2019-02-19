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
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * Value object containing all necessary information for a rule for intensity comparison
 * plus the fatty acid the comparision concerned
 * @author Juergen Hartler
 *
 */
public class IntensityChainVO extends IntensityRuleVO
{
  
  /** the fatty acid at the greater than side of this rule*/
  private FattyAcidVO biggerFA_;
  /** the fatty acid at the smaller than side of this rule*/
  private FattyAcidVO smallerFA_;
  

  /**
   * constructor containing the original rule plus the fatty acid
   * @param rule rule that was applied
   * @param fa fatty acid that was verified
   */
  public IntensityChainVO(IntensityRuleVO rule, FattyAcidVO fa){
    this(rule,fa,fa);
  }

  
  /**
   * constructor containing the original rule plus the fatty acid
   * @param rule rule that was applied
   * @param biggerFA fatty acid with higher area that was verified
   * @param smallerFA fatty acid with smaller area that was verified
   */
  public IntensityChainVO(IntensityRuleVO rule, FattyAcidVO biggerFA, FattyAcidVO smallerFA){
    super(rule);
    this.biggerFA_ = biggerFA;
    this.smallerFA_ = smallerFA;
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
      String value = StaticUtils.getChainFragmentDisplayName(name, biggerFA_.getCarbonDbsId());
      originalToReplacementBigger.put(name, value);
    }
    Hashtable<String,String> originalToReplacementSmaller = new Hashtable<String,String>();
    for (String name : getSmallerNonBasePeakNames()){
      String value = StaticUtils.getChainFragmentDisplayName(name, smallerFA_.getCarbonDbsId());
      originalToReplacementSmaller.put(name, value);
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
    return returnString;
  }
  
      
  /**
   * checks the objects for equality
   */
  public boolean equals(IntensityRuleVO other){
    if (!equals(other)) return false;
    if (!(other instanceof IntensityChainVO)) return false;
    IntensityChainVO chain = (IntensityChainVO)other;
    if (!biggerFA_.getChainId().equalsIgnoreCase(chain.biggerFA_.getChainId())) return false;
    if (!smallerFA_.getChainId().equalsIgnoreCase(chain.smallerFA_.getChainId())) return false;
    return true;
  }
  
  /**
   * inverse method of getReadableRuleInterpretation() - fetches the fatty acid from this interpretation
   * @param rule the readable rule interpretation
   * @param ruleVO rule that is required to form an IntensityChainVO
   * @param chainFragments the found chain fragments with its areas; first key: fatty acid name; second key: name of the fragment; value: CgProbe peak identification object
   * @param missed fragments that were not found in the equation (required for reading of results, since for rules containing "+", not all of the fragments have to be found)
   * @return IntensityChainVO which includes the fatty acid
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (not containing type and OH number) cannot be decoded
   */
  public static IntensityChainVO getFattyAcidsFromReadableRule(String rule, IntensityRuleVO ruleVO, Hashtable<String,Hashtable<String,CgProbe>> chainFragments,
      Hashtable<String,String> missed) throws LipidCombinameEncodingException{
    String[] biggerSmaller = splitToBiggerAndSmallerPart(rule);
    String biggerPart = biggerSmaller[0];
    String smallerPart = biggerSmaller[1];
    Vector<String> biggerNbpNames = ruleVO.getBiggerNonBasePeakNames();
    Vector<String> smallerNbpNames = ruleVO.getSmallerNonBasePeakNames();
    String biggerFA = extractFANames(biggerPart, biggerNbpNames);
    String smallerFA = extractFANames(smallerPart, smallerNbpNames);
    FattyAcidVO biggerVO = null;
    FattyAcidVO smallerVO = null;
    if (biggerFA!=null)
      biggerVO = selectChainFromFoundFragments(biggerFA,biggerNbpNames,chainFragments,missed);
    if (smallerFA!=null) 
      smallerVO = selectChainFromFoundFragments(smallerFA,smallerNbpNames,chainFragments,missed);
    if (biggerVO==null)
      biggerVO = smallerVO;
    if (smallerVO==null)
      smallerVO = biggerVO;
    return new IntensityChainVO(ruleVO,biggerVO,smallerVO);
  }
  
  
  /**
   * extracts the names of the fatty acids out of an equation part
   * @param rulePart equation part
   * @param names the possible fragment names
   * @return the name of the fatty acid
   */
  public static String extractFANames(String rulePart, Vector<String> names){
    String fa = null;
    char[] chars = rulePart.toCharArray();
    //this is for Alex123 annotation
    if (rulePart.indexOf("FA ")!=-1 && rulePart.length()>(rulePart.indexOf("FA ")+4) && (Character.isDigit(chars[rulePart.indexOf("FA ")+3])||
        (chars[rulePart.indexOf("FA ")+4])=='-' && (chars[rulePart.indexOf("FA ")+3]=='P'||chars[rulePart.indexOf("FA ")+3]=='O'))){
      fa = extracFANamesFromAlex123Annotation(rulePart);
    } else if (rulePart.indexOf("(")!=-1){
      for (String name : names){
        if (rulePart.indexOf(name)==-1) continue;
        String subPart = rulePart.substring(rulePart.indexOf(name)+name.length());
        fa = subPart.substring(subPart.indexOf("(")+1,subPart.indexOf(")"));
        break;
      }
    }
    return fa;
  }
  
  private static String extracFANamesFromAlex123Annotation(String rulePart){
    int start = rulePart.indexOf("FA ")+3;
    int stop = start;
    boolean isChain = true;
    char[] chars = rulePart.toCharArray();
    if (chars.length>(start+1) && chars[start+1]=='-' && (chars[start]=='O' || chars[start]=='P')){
      stop = stop+2;
    }
    while (isChain && stop<rulePart.length()){
      if (Character.isDigit(chars[stop]) || chars[stop]==':')
        stop++;
      else
        isChain=false;
    }
    String fa = rulePart.substring(start,stop);
    return fa;
  }
  
}
