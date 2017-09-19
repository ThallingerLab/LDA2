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

package at.tugraz.genome.lda.msn;

import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.msn.vos.IntensityChainVO;
import at.tugraz.genome.lda.msn.vos.IntensityPositionVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * VO inheriting from LipidParameterSet which contains information about MS1 evidence
 * the additional information of this VO reflects MSn evidence
 * @author Juergen Hartler
 *
 */
public class LipidomicsMSnSet extends LipidParameterSet
{
  /** possible statuses of which MSn information can be provided */
  public final static int NO_MSN_PRESENT = 0;
  public final static int DISCARD_HIT = 1;
  public final static int HEAD_GROUP_DETECTED = 2;
  public final static int FRAGMENTS_DETECTED = 3;
  public final static int POSITION_DETECTED = 4;
  
  /** the achieved status by MSn checks */
  private int status_;
  /** the allowed m/z tolerance for the fragment detection */
  private float mzTolerance_;
  /** the detected head group fragments, the key is the fragment name */
  private Hashtable<String,CgProbe> headGroupFragments_;
  /** the intensity rules that were fulfilled by the head group fragments */
  private Hashtable<String,IntensityRuleVO> headIntensityRules_;
  /** the detected chain fragments, the key is the fragment name */
  private Hashtable<String,Hashtable<String,CgProbe>> chainFragments_;
  /** the intensity rules that were fulfilled by the chain fragments */
  private Hashtable<String,Hashtable<String,IntensityChainVO>> chainIntensityRules_;
  /** list of valid chain combinations */
  private Vector<String> validChainCombinations_;
  /** at how many positions the fatty acids may be assigned */
  private int numberOfPositions_;
  /** the intensities of the base peaks */
  private Hashtable<Integer,Float> basePeakValues_;
  
  
  /** position definitions for fatty acid combinations - first key: id of the chain combination; second hash: lookup from the position of the fatty acid in the combination key to the assigned position*/
  private Hashtable<String,Hashtable<Integer,Integer>> positionDefinition_;
  /** 
   * hash containing the fulfilled MSn evidence for a position assignment
   * first key: id of the chain combination
   * second key: position where evidence is provided
   * values: rules that are fulfilled by MSn evidence
  */
  private Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> positionEvidence_;
  
  /**
   * all the required information for this VO has to be provided in the constructor - later adaptions are not possible
   * @param set the MS1 evidence
   * @param status indicator for possible lipid identification
   * @param mzTolerance the m/z range for MSn identification
   * @param headGroupFragments the identified head group fragments
   * @param headIntensityRules the head intensity rules that where fulfilled
   * @param chainFragments the identified chain fragments
   * @param chainIntensityRules the chain intensity rules that where fulfilled
   * @param validChainCombinations list of valid chain combinations
   * @param positionDefinition position definitions for fatty acid combinations - first key: id of the chain combination; second hash: lookup from the position of the fatty acid in the combination key to the assigned position
   * @param positionEvidence hash containing the fulfilled MSn evidence for a position assignment - first key: id of the chain combination; second key: position where evidence is provided; values: rules that are fulfilled by MSn evidence
   * @param numberOfPositions at how many positions the fatty acids may be assigned
   * @param basePeakValues found values for the base peak
   */
  public LipidomicsMSnSet(LipidParameterSet set, int status, float mzTolerance, Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,IntensityRuleVO> headIntensityRules,
      Hashtable<String,Hashtable<String,CgProbe>> chainFragments, Hashtable<String,Hashtable<String,IntensityChainVO>> chainIntensityRules,
      Vector<String> validChainCombinations, Hashtable<String,Hashtable<Integer,Integer>> positionDefinition,
      Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> positionEvidence, int numberOfPositions, 
      Hashtable<Integer,Float> basePeakValues){
    super(set);
    status_ = status;
    mzTolerance_ = mzTolerance;
    headGroupFragments_ = new Hashtable<String,CgProbe>(headGroupFragments);
    headIntensityRules_ =  new Hashtable<String,IntensityRuleVO>(headIntensityRules);
    chainFragments_ = new Hashtable<String,Hashtable<String,CgProbe>>();
    for (String fa : chainFragments.keySet()){
      chainFragments_.put(fa, new Hashtable<String,CgProbe>(chainFragments.get(fa)));
    }
    chainIntensityRules_ = new Hashtable<String,Hashtable<String,IntensityChainVO>>();
    for (String key : chainIntensityRules.keySet()){
      chainIntensityRules_.put(key, new Hashtable<String,IntensityChainVO>(chainIntensityRules.get(key)));
    }
    validChainCombinations_ = new Vector<String>(validChainCombinations);
    positionDefinition_ = new Hashtable<String,Hashtable<Integer,Integer>>();
    for (String key : positionDefinition.keySet()){
      positionDefinition_.put(key, new Hashtable<Integer,Integer>(positionDefinition.get(key)));
    }
    positionEvidence_ = new Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>();
    for (String key1 : positionEvidence.keySet()){
      Hashtable<Integer,Vector<IntensityPositionVO>> subTable = new  Hashtable<Integer,Vector<IntensityPositionVO>>();
      for (Integer key2 : positionEvidence.get(key1).keySet()){
        subTable.put(key2, new Vector<IntensityPositionVO>(positionEvidence.get(key1).get(key2)));
      }
      positionEvidence_.put(key1, subTable);
    }
    numberOfPositions_ = numberOfPositions;
    basePeakValues_ = new  Hashtable<Integer,Float>(basePeakValues);
  }
  
  public LipidomicsMSnSet(LipidomicsMSnSet set){
    this(set,set.status_,set.mzTolerance_,set.headGroupFragments_,set.headIntensityRules_,set.chainFragments_,set.chainIntensityRules_,
        set.validChainCombinations_,set.positionDefinition_,set.positionEvidence_,set.numberOfPositions_,set.basePeakValues_);
  }

  /**
   * 
   * @return detected head group fragments, the key is the fragment name
   */
  public Hashtable<String,CgProbe> getHeadGroupFragments()
  {
    return headGroupFragments_;
  }

  /**
   * 
   * @return intensity rules that were fulfilled by the head group fragments
   */
  public Hashtable<String,IntensityRuleVO> getHeadIntensityRules()
  {
    return headIntensityRules_;
  }
  /**
   * 
   * @return detected chain fragments, the key is the fragment name
   */
  public Hashtable<String,Hashtable<String,CgProbe>> getChainFragments()
  {
    return chainFragments_;
  }

  /**
   * 
   * @return intensity rules that were fulfilled by the chain fragments
   */
  public Hashtable<String,Hashtable<String,IntensityChainVO>> getChainIntensityRules()
  {
    return chainIntensityRules_;
  }
  /**
   * 
   * @return position definitions for fatty acid combinations - first key: id of the chain combination; second hash: lookup from the position of the fatty acid in the combination key to the assigned position
   */
  public Hashtable<String,Hashtable<Integer,Integer>> getPositionDefinition()
  {
    return positionDefinition_;
  }
  /**
   * 
   * @return  hash containing the fulfilled MSn evidence for a position assignment - first key: id of the chain combination - second key: position where evidence is provided - values: rules that are fulfilled by MSn evidence
   */
  public Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> getPositionEvidence()
  {
    return positionEvidence_;
  }
  
  /**
   * 
   * @return the lipid names for the identified lipid species based on MSn evidence - the entry may be a String, or a Vector<String> if several position specific assignments are possible
   */
  public Vector<Object> getMSnIdentificationNames(){
    Vector<Object> names = new Vector<Object>();
    if (status_==NO_MSN_PRESENT || status_==HEAD_GROUP_DETECTED){
      names.add(getNameStringWithoutRt());
    } else if (status_==FRAGMENTS_DETECTED || status_==POSITION_DETECTED){     
      for (String combiName : validChainCombinations_){
        if (status_==POSITION_DETECTED && positionDefinition_.containsKey(combiName)){
          Vector<String> posNames = getPositionSpecificCombiNames(combiName);
          if (posNames.size()==1) names.add(posNames.get(0));
          else names.add(getPositionSpecificCombiNames(combiName));
        }else{
          names.add(combiName);
        }
      }
    }
    return names;
  }
  
  /**
   * for position evidence, returns all the position combinations that are possible with this rules
   * @param combiName the name of the FA combination without position assignments
   * @return
   */
  private Vector<String> getPositionSpecificCombiNames(String combiName){
    Hashtable<Integer,Integer> positions = positionDefinition_.get(combiName);
    String[] fas =  getFAsFromCombiName(combiName);
    Vector<String> unassignedFAs = new Vector<String>();
    Hashtable<Integer,String> definedPositions = new Hashtable<Integer,String>();
    for (int i=0;i!=fas.length;i++){
      if (!positions.containsKey(i) || positions.get(i)==-1) unassignedFAs.add(fas[i]);
      else definedPositions.put(positions.get(i), fas[i]);
    }
    return getPermutedChainPositionNames(definedPositions, unassignedFAs);
  }
  
  /**
   * returns various names if only part of the chains could be correctly assigned
   * @param definedPositions positions where a definite assignment is given
   * @param unassignedFAs fatty acids where no assignment could be made
   * @return the lipid names that are possible with these assignments
   */
  private Vector<String> getPermutedChainPositionNames(Hashtable<Integer,String> definedPositions, Vector<String> unassignedFAs){
    Vector<Hashtable<Integer,String>> possibleCombis = new Vector<Hashtable<Integer,String>>();
    Vector<String> names = new Vector<String>();
    if (unassignedFAs.size()==0){
      possibleCombis.add(definedPositions);
    }else{
      Hashtable<String,String> combis = new Hashtable<String,String>();
      Vector<String> unassigned = new Vector<String>(unassignedFAs);
      for (int i=(definedPositions.size()+unassignedFAs.size()); i<numberOfPositions_;i++){
        unassigned.add("-");
      }
      Vector<String> permutedNames = FragmentCalculator.getPermutedNames(unassigned);
      for (String permutedName : permutedNames){
        if (!combis.containsKey(permutedName)) combis.put(permutedName, permutedName);
      }
      for (String permutedName : combis.keySet()){
        String[] parts =  getFAsFromCombiName(permutedName);
        int permutedPos = 0;
        Hashtable<Integer,String> combiHash = new Hashtable<Integer,String>();
        for (int i=0; i!=(definedPositions.size()+unassigned.size());i++){
          String fa = null;
          if (definedPositions.containsKey(i)) fa = definedPositions.get(i);
          else{
            fa = parts[permutedPos];
            permutedPos++;

          }
          combiHash.put(i, fa);
        }
        possibleCombis.add(combiHash);
      }
    }
    for (Hashtable<Integer,String> combi : possibleCombis){
      String name = "";
      for (int i=0;i!=this.numberOfPositions_;i++){
        if (i!=0) name+= "/";
        if (combi.containsKey(i))
          name += combi.get(i);
        else
          name += "-";
      }
      names.add(name);
    }
    return names;
  }
  
  /**
   * 
   * @return the relative contribution of one fatty acid combination to the total area detected in MS1 - key is the fatty acid combination id
   */
  public Hashtable<String,Double> getChainCombinationRelativeAreas(){
    Hashtable<String,Double> relativeIntensityOfCombination = new Hashtable<String,Double>();
    Hashtable<String,Integer> faChainOccurrences = new Hashtable<String,Integer>();
    for (String combi : validChainCombinations_){
      String[] combiFAs = getFAsFromCombiName(combi);
      for (String combiFA : combiFAs){
        int count = 0;
        if (faChainOccurrences.containsKey(combiFA)) count = faChainOccurrences.get(combiFA);
        count++;
        faChainOccurrences.put(combiFA,count);
      }
    }
    // calculate a total area for each chain found
    Hashtable<String,Double> areas = new Hashtable<String,Double>();
    for (String faName : chainFragments_.keySet()){
      areas.put(faName, getAreaOfOneChain(faName));
    }
    Hashtable<String,Double> combiAreas = new Hashtable<String,Double>();
    double totalArea = 0d;
    for (String combi : validChainCombinations_){
      double area = 0d;
      String[] combiFAs = getFAsFromCombiName(combi);
      for (String combiFA : combiFAs){
        if (areas.containsKey(combiFA))
          area += areas.get(combiFA)/((double)faChainOccurrences.get(combiFA));
      }
      combiAreas.put(combi, area);
      totalArea += area;
    }
    for (String combi : combiAreas.keySet()) relativeIntensityOfCombination.put(combi, combiAreas.get(combi)/totalArea);
    return relativeIntensityOfCombination;
  }
  
  /**
   * for splitting of quantified area - returns relative contribution of identification [0...1] 
   * @param fullName the name of the identified lipid species
   * @return relative contribution of this species [0...1]
   */
  public double getRelativeIntensity(String fullName){
    Hashtable<String,Double> relAreas = getChainCombinationRelativeAreas();
    double area = 0;
    if (!relAreas.containsKey(fullName)){
      for (String combiName : this.validChainCombinations_){
        if (!positionDefinition_.containsKey(combiName)) continue;
        Vector<String> posCombiNames = getPositionSpecificCombiNames(combiName);
        for (String posCombiName : posCombiNames){
          if (posCombiName.equalsIgnoreCase(fullName)) return relAreas.get(combiName);
        }
      }
    } else {
      area = relAreas.get(fullName);
    }
    return area;
  }
  
  /**
   * general method for splitting the id of fatty acid combinations to the individual lipid species
   * @param combiName id of fatty acid combinations
   * @return individual lipid species names
   */
  public static String[] getFAsFromCombiName(String combiName){
    return FragmentCalculator.getFAsFromCombiName(combiName);
  }

  private double getAreaOfOneChain(String faName){
    double area = 0f;
    for (CgProbe probe : chainFragments_.get(faName).values()){
      double oneArea = (double)probe.Area;
      area += oneArea;
    }
    return area;
  }

  /**
   * 
   * @return status of which MSn information can be provided
   */
  public int getStatus()
  {
    return status_;
  }
  
  /**
   * 
   * @return the allowed m/z tolerance for the fragment detection
   */
  public float getMSnMzTolerance(){
    return mzTolerance_;
  }

  /**
   * searches for a fragment for this rule and returns the area of the base peak
   * @param rule the VO containing the rule parameters
   * @return base peak area
   */
  public float getBasePeak(IntensityRuleVO rule)
  {
    CgProbe probe = null;
    Vector<String> nonBpNames = rule.getBiggerNonBasePeakNames();
    for (String nonBpName : nonBpNames){
      probe = getFragment(nonBpName,rule.getBiggerFA(),headGroupFragments_,chainFragments_);
      if (probe!=null) break;
    }
    if (probe!=null) return getBasePeak(probe.getMsLevel());

    nonBpNames = rule.getSmallerNonBasePeakNames();
    for (String nonBpName : nonBpNames){
      probe = getFragment(nonBpName,rule.getSmallerFA(),headGroupFragments_,chainFragments_);
      if (probe!=null) break;
    }   
    return getBasePeak(probe.getMsLevel());
  }

  
  /**
   * 
   * @param msLevel at which MSn do I want to look for the base peak
   * @return the value of the base peak
   */
  public float getBasePeak(int msLevel)
  {
    float result = Float.NaN;
    if (basePeakValues_.containsKey(msLevel)) result = basePeakValues_.get(msLevel);
    return result;
  }
  
  /**
   * searches for fragment areas of an IntensityRuleVO
   * @param ruleVO the VO containing the rule parameters
   * @return hashtable containing the fragment areas; key: name of the fragment; value: area
   */
  public Hashtable<String,Float> getFragmentAreas(IntensityRuleVO ruleVO) {
    return getFragmentAreas(ruleVO,headGroupFragments_,chainFragments_);
  }
  
  /**
   * searches for fragment areas of an IntensityRuleVO
   * @param ruleVO the VO containing the rule parameters
   * @param headGroupFragments the found head group fragments with its areas; key: fragment name; value: CgProbe peak identification object
   * @param chainFragments the found chain fragments with its areas; first key: fatty acid name; second key: name of the fragment; value: CgProbe peak identification object
   * @return hashtable containing the fragment areas; key: name of the fragment; value: area
   */
  public static Hashtable<String,Float> getFragmentAreas(IntensityRuleVO ruleVO, Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments) {
    Hashtable<String,Float> areas = new Hashtable<String,Float>();
    Vector<String> biggerNames =  ruleVO.getBiggerNonBasePeakNames();
    String biggerFA = ruleVO.getBiggerFA();
    for (String biggerName : biggerNames){
      Float biggerArea = getFragmentArea(biggerName,biggerFA,headGroupFragments,chainFragments);
      if (biggerArea == null) continue;
      String identifier = biggerName;
      if (ruleVO.getBiggerPosition()>0) identifier += "["+String.valueOf(ruleVO.getBiggerPosition())+"]";
      areas.put(identifier, biggerArea);
    }
    Vector<String> smallerNames =  ruleVO.getSmallerNonBasePeakNames();
    String smallerFA = ruleVO.getSmallerFA();
    for (String smallerName : smallerNames){
      Float smallerArea = getFragmentArea(smallerName,smallerFA,headGroupFragments,chainFragments);
      if (smallerArea == null) continue;
      String identifier = smallerName;
      if (ruleVO.getSmallerPosition()>0) identifier += "["+String.valueOf(ruleVO.getSmallerPosition())+"]";
      areas.put(identifier, smallerArea);
    }
    return areas;
  }

  /**
   * returns the area for one fragment
   * @param name name of the fragment
   * @param faName fatty acid name - required for chain fragments
   * @param headGroupFragments the found head group fragments with its areas; key: fragment name; value: CgProbe peak identification object
   * @param chainFragments the found chain fragments with its areas; first key: fatty acid name; second key: name of the fragment; value: CgProbe peak identification object
   * @return area of the fragment
   */
  private static Float getFragmentArea(String name, String faName, Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments){
    Float area = null;
    if (name.equalsIgnoreCase(IntensityRuleVO.BASEPEAK_NAME)) return area;
    CgProbe fragment = getFragment(name, faName, headGroupFragments, chainFragments);
    if (fragment!=null) area = fragment.Area;
    else area = -1f;
    return area;    
  }
  
  /**
   * returns the CgProbe peak identification object for one fragment
   * @param name name of the wanted fragment
   * @param faName fatty acid name - required for chain fragments
   * @param headGroupFragments the found head group fragments with its areas; key: fragment name; value: CgProbe peak identification object
   * @param chainFragments the found chain fragments with its areas; first key: fatty acid name; second key: name of the fragment; value: CgProbe peak identification object
   * @return CgProbe peak identification object for one fragment
   */
  private static CgProbe getFragment(String name, String faName, Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments){
    CgProbe fragment = null;
    if (headGroupFragments!=null && headGroupFragments.size()>0 && headGroupFragments.containsKey(name)) fragment = headGroupFragments.get(name); 
    else if (faName!=null && faName.length()>0 && chainFragments!=null){
      if (chainFragments.containsKey(faName) && chainFragments.get(faName).containsKey(name)){
        fragment = chainFragments.get(faName).get(name);
      }else if (chainFragments.containsKey(FragmentCalculator.ALKYL_PREFIX.concat(faName)) && chainFragments.get(FragmentCalculator.ALKYL_PREFIX.concat(faName)).containsKey(name)){
        fragment = chainFragments.get(FragmentCalculator.ALKYL_PREFIX.concat(faName)).get(name);
      }else if (chainFragments.containsKey(FragmentCalculator.ALKENYL_PREFIX.concat(faName)) && chainFragments.get(FragmentCalculator.ALKENYL_PREFIX.concat(faName)).containsKey(name)){
        fragment = chainFragments.get(FragmentCalculator.ALKENYL_PREFIX.concat(faName)).get(name);
      }
    }
    return fragment;
  }
  
  /**
   * returns vector of unique position rule information so that it can be written - each entry is one unique rule
   * @param combiName
   * @return
   */
  public Vector<IntensityRuleVO> getFAsInSequenceAsInRule(String combiName){
    Vector<IntensityRuleVO> result = new Vector<IntensityRuleVO>();
    Hashtable<Integer,Vector<IntensityPositionVO>> posRules = positionEvidence_.get(combiName);
    Hashtable<Integer,Integer> posDef = positionDefinition_.get(combiName);
    if (status_ < POSITION_DETECTED || posDef==null) return result;
    String[] fas = getFAsFromCombiName(combiName);
    Hashtable<String,IntensityRuleVO> usedRules = new Hashtable<String,IntensityRuleVO>();
    for (int i=0;i!=fas.length;i++){
      if (!posDef.containsKey(i) || posDef.get(i)==-1 || posRules.size()==0) continue;
      int position = posDef.get(i);
      Vector<IntensityPositionVO> rules = posRules.get((position+1));
      for (IntensityRuleVO vo : rules){
        if (usedRules.containsKey(vo.getRuleIdentifier()) && usedRules.get(vo.getRuleIdentifier()).equals(vo)) continue;
        result.add(vo);
        usedRules.put(vo.getRuleIdentifier(), vo);
      }
    }
    return result;
  }
}
