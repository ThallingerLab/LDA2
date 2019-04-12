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

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.ShortStringVO;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * Value object containing all necessary information for a rule for intensity comparison
 * plus the fatty acid the comparision concerned
 * @author Juergen Hartler
 *
 */
public class IntensityChainVO extends IntensityRuleVO
{
  
  /** the chains at the 'greater than' side of this rule*/
  private Hashtable<String,FattyAcidVO> biggerChains_;
  /** the chains at the 'smaller than' side of this rule*/
  private Hashtable<String,FattyAcidVO> smallerChains_;
  /** does this object hold hydroxylation information*/
  private boolean hasOhInfo_;
  

  /**
   * constructor containing the original rule plus the fatty acid
   * @param rule rule that was applied
   * @chain chain the chain object this rule applies to
   * @param hasOhInfo has this intensity rule chain information - this is important for generating the human readable name
   */
  public IntensityChainVO(IntensityRuleVO rule, FattyAcidVO chain, boolean hasOhInfo){
    this(rule,hasOhInfo);
    Vector<FattyAcidVO> chains = new Vector<FattyAcidVO>();
    chains.add(chain);
    assignChainsToExpressions(chains);
  }
  
  /**
   * common constructor of other options
   * @param rule rule that was applied
   * @param hasOhInfo has this intensity rule chain information - this is important for generating the human readable name
   */
  private IntensityChainVO(IntensityRuleVO rule, boolean hasOhInfo) {
    super(rule);
    this.hasOhInfo_ = hasOhInfo;
  }

  
  /**
   * constructor containing the original rule plus the participating chain objects
   * @param rule intensity rule that was applied
   * @param chains chain objects this rule applies to
   * @param hasOhInfo has this intensity rule chain information - this is important for generating the human readable name
   */
  public IntensityChainVO(IntensityRuleVO rule, Vector<FattyAcidVO> chains, boolean hasOhInfo){
    this(rule,hasOhInfo);
    assignChainsToExpressions(chains);
  }
  
  
  /**
   * constructor containing the original rule plus the fatty acids 
   * @param rule rule that was applied
   * @param biggerChains the chains at the 'greater than' side of this rule
   * @param smallerChains chains at the 'smaller than' side of this rule
   * @param hasOhInfo does this object hold hydroxylation information
   */
  private IntensityChainVO(IntensityRuleVO rule, Hashtable<String,FattyAcidVO> biggerChains,
      Hashtable<String,FattyAcidVO> smallerChains, boolean hasOhInfo){
    this(rule,hasOhInfo);
    this.biggerChains_ = biggerChains;
    this.smallerChains_ = smallerChains;
  }
  
  
  /**
   * assigns the chain objects to the fragmentation rules they affect
   * @param chains the available chain objects
   */
  private void assignChainsToExpressions(Vector<FattyAcidVO> chains) {
    biggerChains_ = new Hashtable<String,FattyAcidVO>();
    smallerChains_ = new Hashtable<String,FattyAcidVO>();
    for (FragmentMultVO frag : biggerExpression_.getFragments()) {
      for (FattyAcidVO chain : chains) {
        if (frag.getFragmentType()==chain.getChainType()) {
          biggerChains_.put(frag.getFragmentName(), chain);
          break;
        }
      }
    }
    for (FragmentMultVO frag : smallerExpression_.getFragments()) {
      for (FattyAcidVO chain : chains) {
        if (frag.getFragmentType()==chain.getChainType()) {
          smallerChains_.put(frag.getFragmentName(), chain);
          break;
        }
      }
    }    
  }

  
  /**
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @return the human readable String representation for storing in Excel
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public String getReadableRuleInterpretation(HydroxyEncoding faEncoding, HydroxyEncoding lcbEncoding) throws LipidCombinameEncodingException {
    Hashtable<String,String> originalToReplacementBigger = new Hashtable<String,String>();
    for (String name : getBiggerNonHeadAndBasePeakNames().keySet()){
      String value = StaticUtils.getChainFragmentDisplayName(name,StaticUtils.getHumanReadableChainName(this.biggerChains_.get(name),faEncoding,
          lcbEncoding,hasOhInfo_,StaticUtils.isAnAlex123Fragment(name)));
      originalToReplacementBigger.put(name, value);
    }
    Hashtable<String,String> originalToReplacementSmaller = new Hashtable<String,String>();
    for (String name : getSmallerNonHeadAndBasePeakNames().keySet()){
      String value = StaticUtils.getChainFragmentDisplayName(name,StaticUtils.getHumanReadableChainName(this.smallerChains_.get(name),faEncoding,
          lcbEncoding,hasOhInfo_,StaticUtils.isAnAlex123Fragment(name)));
      originalToReplacementSmaller.put(name, value);
    }
    String biggerPart = "";
    String smallerPart = "";
    if (orRule_) {
      biggerPart = equation_;
    }else {
      String[] biggerSmaller = splitToBiggerAndSmallerPart(equation_);
      biggerPart = biggerSmaller[0];
      smallerPart = biggerSmaller[1];
    }
    biggerPart = replaceParts(biggerPart, originalToReplacementBigger);
    smallerPart = replaceParts(smallerPart, originalToReplacementSmaller);
    String returnString = "";
    if (equation_.indexOf(">")!=-1){
      returnString = biggerPart+">"+smallerPart;
    } else if (equation_.indexOf("<")!=-1){
      returnString = smallerPart+"<"+biggerPart;
    }else if (equation_.indexOf("|")!=-1){
      returnString = biggerPart;
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
    if (chain.biggerChains_.size()!=biggerChains_.size()) return false;
    for (String frag : biggerChains_.keySet()) {
      if ((!chain.biggerChains_.containsKey(frag)) || chain.biggerChains_.get(frag)!=biggerChains_.get(frag))
        return false;
    }    
    if (chain.smallerChains_.size()!=smallerChains_.size()) return false;
    for (String frag : smallerChains_.keySet()) {
      if ((!chain.smallerChains_.containsKey(frag)) || chain.smallerChains_.get(frag)!=smallerChains_.get(frag))
        return false;
    }    
    return true;
  }
  
  /**
   * inverse method of getReadableRuleInterpretation() - fetches the fatty acid from this interpretation
   * @param rule the readable rule interpretation
   * @param ruleVO rule that is required to form an IntensityChainVO
   * @param chainFragments the found chain fragments with its areas; first key: fatty acid name; second key: name of the fragment; value: CgProbe peak identification object
   * @param missed fragments that were not found in the equation (required for reading of results, since for rules containing "+", not all of the fragments have to be found)
   * @param faHydroxyEncoding the OH encodings of the FA moiety
   * @param lcbHydroxyEncoding the OH encodings of the LCB moiety
   * @param hasOhInfo does this object hold hydroxylation information
   * @return IntensityChainVO which includes the fatty acid
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (not containing type and OH number) cannot be decoded
   */
  public static IntensityChainVO getFattyAcidsFromReadableRule(String rule, IntensityRuleVO ruleVO, Hashtable<String,Hashtable<String,CgProbe>> chainFragments,
      Hashtable<String,Short> missed, HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding, boolean hasOhInfo) throws LipidCombinameEncodingException{
    String biggerPart = "";
    String smallerPart = "";
    boolean orRule = rule.indexOf("|")!=-1 && rule.indexOf(">")==-1 && rule.indexOf("<")==-1;
    if (orRule) {
      biggerPart = rule;
    }else {
      String[] biggerSmaller = splitToBiggerAndSmallerPart(rule);
      biggerPart = biggerSmaller[0];
      smallerPart = biggerSmaller[1];
    }
    Vector<ShortStringVO> biggerNbpNames =  FragmentRuleVO.getLengthSortedFragmentNames(new Hashtable<String,Short>(),ruleVO.getBiggerNonHeadAndBasePeakNames(),missed);
    Vector<ShortStringVO> smallerNbpNames = FragmentRuleVO.getLengthSortedFragmentNames(new Hashtable<String,Short>(),ruleVO.getSmallerNonHeadAndBasePeakNames(),missed);
    Hashtable<String,FattyAcidVO> biggerChains = extractFANames(biggerPart, biggerNbpNames, faHydroxyEncoding, lcbHydroxyEncoding);
    Hashtable<String,FattyAcidVO> smallerChains = extractFANames(smallerPart, smallerNbpNames, faHydroxyEncoding, lcbHydroxyEncoding);
    return new IntensityChainVO(ruleVO,biggerChains, smallerChains, hasOhInfo);
  }
  
  
  /**
   * extracts the names of the fatty acids out of an equation part
   * @param ruleOriginal the smaller or bigger part of the equation
   * @param names the possible fragment names
   * @param faHydroxyEncoding the character encoding of the number of hydroxylation sites for the FA
   * @param lcbHydroxyEncoding the character encoding of the number of hydroxylation sites for the LCB
   * @return the names of the chains - key: fragment name; value: decoded chain object
   * @throws LipidCombinameEncodingException thrown when a lipid id (containing type and OH number) cannot be decoded 
   */
  public static Hashtable<String,FattyAcidVO> extractFANames(String ruleOriginal, Vector<ShortStringVO> names, HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding) throws LipidCombinameEncodingException {
    Hashtable<String,FattyAcidVO> fas = new Hashtable<String,FattyAcidVO>();
    String rulePart = new String(ruleOriginal);
    char[] chars = rulePart.toCharArray();
    if ((rulePart.indexOf("FA ")!=-1 && rulePart.length()>(rulePart.indexOf("FA ")+4) && (Character.isDigit(chars[rulePart.indexOf("FA ")+3])||
        (chars[rulePart.indexOf("FA ")+4])=='-' && (chars[rulePart.indexOf("FA ")+3]=='P'||chars[rulePart.indexOf("FA ")+3]=='O'))) ||
        (rulePart.indexOf("LCB ")!=-1 && rulePart.length()>(rulePart.indexOf("LCB ")+5) && Character.isDigit(chars[rulePart.indexOf("LCB ")+4]))){
      fas = extracFANamesFromAlex123Annotation(rulePart,names,faHydroxyEncoding,lcbHydroxyEncoding);
    } else if (rulePart.indexOf("(")!=-1){
      for (ShortStringVO nameVO : names){
        String name = nameVO.getKey();
        if (rulePart.indexOf(name)==-1) continue;
        if ((rulePart.indexOf(name)+name.length()+2)>rulePart.length() || !rulePart.substring(rulePart.indexOf(name)+name.length()).startsWith("("))
          continue;
        String subPart = rulePart.substring(rulePart.indexOf(name)+name.length());
        fas.put(nameVO.getKey(),StaticUtils.decodeHumanReadableChain(subPart.substring(subPart.indexOf("(")+1,subPart.indexOf(")")),faHydroxyEncoding,lcbHydroxyEncoding, false));
        //the next two lines are for backward compatibility - I am not sure whether I can keep it this way
        if (nameVO.getValue()!=LipidomicsConstants.CHAIN_TYPE_MISSED)
          fas.get(nameVO.getKey()).correctChainType(nameVO.getValue());
        rulePart = rulePart.substring(0,rulePart.indexOf(name))+rulePart.substring(rulePart.indexOf(name)+subPart.indexOf(")")+1);
      }
    }
    return fas;
  }
  
  /**
   * extracts the chain names of an Alex123 encoded intensity rule part
   * @param rulePart the the smaller or bigger part of the equation
   * @param names the possible fragment names
   * @param faHydroxyEncoding the character encoding of the number of hydroxylation sites for the FA
   * @param lcbHydroxyEncoding the character encoding of the number of hydroxylation sites for the LCB
   * @return the names of the chains - key: fragment name; value: decoded chain object
   * @throws LipidCombinameEncodingException thrown when a lipid id (containing type and OH number) cannot be decoded
   */
  private static Hashtable<String,FattyAcidVO> extracFANamesFromAlex123Annotation(String rulePart, Vector<ShortStringVO> names, HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding) throws LipidCombinameEncodingException{
    Hashtable<String,FattyAcidVO> fas = new Hashtable<String,FattyAcidVO>();
    String whole = new String(rulePart);
    String part;
    String chainString;
    String oneChain;
    for (ShortStringVO name : names) {
      //System.out.println(name.getKey()+": "+name.getValue());
      if (name.getValue()==LipidomicsConstants.CHAIN_TYPE_NO_CHAIN)
        continue;
      String[] split = new String[2];
      if (name.getKey().indexOf("(")==-1)
        throw new LipidCombinameEncodingException("The fragment: "+name.getKey()+" is not an Alex123 compatible chain part!");
      split[0] = name.getKey().substring(0,name.getKey().indexOf("("));
      split[1] = name.getKey().substring(name.getKey().indexOf("(")+1);
      part = new String(whole);
      int startIndex = 0;
      int endIndex = 0;
      chainString = null;
      oneChain = null;
      while (part.length()>name.getKey().length() && part.indexOf(split[0])!=-1) {
        startIndex += part.indexOf(split[0]);
        part = part.substring(startIndex);
        if (part.indexOf("(")==-1)
          break;
        if (part.substring(part.indexOf("(")+1).startsWith(split[1])) {
          if (part.substring(part.indexOf("(")).indexOf(")")==-1)
            throw new LipidCombinameEncodingException("The rule part: "+rulePart+"cannot be decoded to the corresponding chains");
          endIndex = startIndex+part.indexOf(")")+1;
          if (part.toCharArray().length>(endIndex-startIndex+1)&&part.toCharArray()[endIndex+1]=='[') {
            endIndex = startIndex+part.indexOf("]")+1;
          }
          chainString = whole.substring(startIndex, endIndex); 
          oneChain = part.substring(split[0].length(), part.indexOf("(")).trim();
          whole = whole.substring(0,startIndex)+whole.substring(endIndex);
          break;
        }else {
          part = part.substring(startIndex);
        }
      }
      if (chainString!=null) {
        int prefixStop = 0;
        char[] chainChars = oneChain.toCharArray();
        for (int i=0; i!=chainChars.length; i++) {
          if (Character.isDigit(chainChars[i])) break;
          else prefixStop++;
        }
        String prefix = oneChain.substring(0,prefixStop);
        oneChain = oneChain.substring(prefixStop);
        int oh = 0;
        if (oneChain.indexOf(LipidomicsConstants.ALEX_OH_SEPARATOR)!=-1) {
          try {
            oh = Integer.parseInt(oneChain.substring(oneChain.indexOf(LipidomicsConstants.ALEX_OH_SEPARATOR)+1));
            oneChain = oneChain.substring(0,oneChain.indexOf(LipidomicsConstants.ALEX_OH_SEPARATOR));
          }catch (NumberFormatException nfx) {
            throw new LipidCombinameEncodingException("The rule part: "+rulePart+"cannot be decoded to the corresponding chains: there are problems with the OH encoding!");
          }
        }
        try {
          int[] cAndDbs = StaticUtils.parseCAndDbsFromChainId(oneChain);
          fas.put(name.getKey(), new FattyAcidVO(name.getValue(),prefix,cAndDbs[0],cAndDbs[1],oh,-1, null));
        }
        catch (Exception e) {
          e.printStackTrace();
          throw new LipidCombinameEncodingException("The rule part: "+rulePart+"cannot be decoded to the corresponding chains: there are problems decoding the number of carbon atoms and double bonds!");
        } 
      }
    }
    return fas;
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
  
  
  /**
   *
   * @return all chain objects that participated in this equation (does not return missed chain objects)
   */
  public Vector<FattyAcidVO> getParticipatingChains(){
    Hashtable<String,FattyAcidVO> hash = new Hashtable<String,FattyAcidVO>();
    if (biggerChains_!=null) {
      for (FattyAcidVO chain : this.biggerChains_.values()) {
        if (!hash.containsKey(chain.getChainId()))
          hash.put(chain.getChainId(), chain);
      }
    }
    if (smallerChains_!=null) {
      for (FattyAcidVO chain : this.smallerChains_.values()) {
        if (!hash.containsKey(chain.getChainId()))
          hash.put(chain.getChainId(), chain);
      }
    }
    return new Vector<FattyAcidVO>(hash.values());
  }
}
