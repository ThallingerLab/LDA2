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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.IntensityChainVO;
import at.tugraz.genome.lda.msn.vos.IntensityPositionVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleStringVO;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.voutils.GeneralComparator;
import javafx.util.Pair;

/**
 * TODO: remove ambig pos.
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
  /** the retention times of the used MSn spectra; the first key is the msLevel; second key scan number; value: retention time*/
  private Hashtable<Integer,LinkedHashMap<Integer,Float>> msnRetentionTimes_;
  /** the MS-levels used for each molecular species; key: molecular species; value MS-levels used*/
  private Hashtable<String,Set<Integer>> msLevels_;
  
  /** identifier for msLevels_ for head group only*/
  public final static String MSLEVEL_HEAD_IDENTIFIER = "head";
  /** identifier for msLevels_ for all available molecular species (in a cumulative species export)*/
  public final static String MSLEVEL_ALL_IDENTIFIER = "all";

  
  /** position definitions for fatty acid combinations - first key: id of the chain combination; second hash: lookup from the position of the fatty acid in the combination key to the assigned position*/
  private Hashtable<String,Hashtable<Integer,Integer>> positionDefinition_;
  /** 
   * hash containing the fulfilled MSn evidence for a position assignment
   * first key: id of the chain combination
   * second key: position where evidence is provided
   * values: rules that are fulfilled by MSn evidence
  */
  private Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> positionEvidence_;
  
  /** lookup table from human readable lipid molecular species names to the ones used in the hash tables*/
  private LinkedHashMap<String,String> nameLookupHumReadableToPositionInsensitve_;
  /** lookup table from human readable lipid molecular species names to the ones used in the hash tables for new nomenclature*/
  private LinkedHashMap<String,String> nameLookupPositionSnNomenclature_;
  /** lookup table from human readable lipid molecular species names to the ones used in the hash tables*/
  private LinkedHashMap<String,String> nameLookupPositionInsensitve_;
  /** lipid molecular species identifications where an ambiguous position assignment was made*/
  private Hashtable<String,Vector<String>> ambiguousPositionIdentifications_;
  /** the relative (percentual) value of a species contribution_;*/
  private Hashtable<String,Double> relativeIntensityOfCombination_;
  /** lookup table from encoded ones to the human readable ones*/
  private Hashtable<String,String> chainNameLookupHumanReadable_;
  /** a list of the involved FA objects for this identification*/
  private Hashtable<String,FattyAcidVO> involvedFAs_;
  /** the FA encoding */
  private HydroxyEncoding faEncoding_;
  /** the LCB encoding */
  private HydroxyEncoding lcbEncoding_;
  
  
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
   * @param relativeIntensityOfCombination the relative share on the MS1 intensity of each chain combination - kind of percentual value
   * @param positionDefinition position definitions for fatty acid combinations - first key: id of the chain combination; second hash: lookup from the position of the fatty acid in the combination key to the assigned position
   * @param positionEvidence hash containing the fulfilled MSn evidence for a position assignment - first key: id of the chain combination; second key: position where evidence is provided; values: rules that are fulfilled by MSn evidence
   * @param numberOfPositions at how many positions the fatty acids may be assigned
   * @param basePeakValues found values for the base peak
   * @param msnRetentionTimes the retention times of the used MSn spectra; the first key is the msLevel; second key scan number; value: retention time
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private LipidomicsMSnSet(LipidParameterSet set, int status, float mzTolerance, Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,IntensityRuleVO> headIntensityRules,
      Hashtable<String,Hashtable<String,CgProbe>> chainFragments, Hashtable<String,Hashtable<String,IntensityChainVO>> chainIntensityRules,
      Vector<String> validChainCombinations, Hashtable<String,Double> relativeIntensityOfCombination, Hashtable<String,Hashtable<Integer,Integer>> positionDefinition,
      Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> positionEvidence, int numberOfPositions, 
      Hashtable<Integer,Float> basePeakValues, Hashtable<Integer,LinkedHashMap<Integer,Float>> msnRetentionTimes) throws LipidCombinameEncodingException{
    super(set);
    status_ = status;
    mzTolerance_ = mzTolerance;
    msLevels_ = new Hashtable<String,Set<Integer>>();
    involvedFAs_ = new Hashtable<String,FattyAcidVO>();
    Set<Integer> headLevel = new HashSet<Integer>();
    headGroupFragments_ = new Hashtable<String,CgProbe>(headGroupFragments);
    for (CgProbe probe : headGroupFragments_.values())
      headLevel.add(probe.getMsLevel());
    msLevels_.put(MSLEVEL_HEAD_IDENTIFIER, headLevel);
    Set<Integer> allLevel = new HashSet<Integer>(headLevel);
    Hashtable<String,Set<Integer>> faLevels = new Hashtable<String,Set<Integer>>();
    headIntensityRules_ =  new Hashtable<String,IntensityRuleVO>(headIntensityRules);
    chainFragments_ = new Hashtable<String,Hashtable<String,CgProbe>>();
    for (String fa : chainFragments.keySet()){
      chainFragments_.put(fa, new Hashtable<String,CgProbe>(chainFragments.get(fa)));
      Set<Integer> faSpecificLevel = new HashSet<Integer>();
      for (CgProbe probe : chainFragments.get(fa).values())
        faSpecificLevel.add(probe.getMsLevel());
      allLevel.addAll(faSpecificLevel);
      faLevels.put(fa, faSpecificLevel);
    }
    chainIntensityRules_ = new Hashtable<String,Hashtable<String,IntensityChainVO>>();
    for (String key : chainIntensityRules.keySet()){
      chainIntensityRules_.put(key, new Hashtable<String,IntensityChainVO>(chainIntensityRules.get(key)));
    }
    validChainCombinations_ = new Vector<String>(validChainCombinations);
    for (String combi : validChainCombinations_){
      Set<Integer> combiLevel = new HashSet<Integer>();
      for (FattyAcidVO fa : StaticUtils.decodeLipidNamesFromChainCombi(combi)){
        this.involvedFAs_.put(fa.getChainId(), fa);
        if (faLevels.containsKey(fa.getChainId()))
          combiLevel.addAll(faLevels.get(fa.getChainId()));
      }
      msLevels_.put(combi, combiLevel);
    }
    msLevels_.put(MSLEVEL_ALL_IDENTIFIER, allLevel);
    relativeIntensityOfCombination_ = relativeIntensityOfCombination;
    
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
    setMsnRetentionTimes(msnRetentionTimes);
  }
  
  
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
   * @param relativeIntensityOfCombination the relative share on the MS1 intensity of each chain combination - kind of percentual value
   * @param positionDefinition position definitions for fatty acid combinations - first key: id of the chain combination; second hash: lookup from the position of the fatty acid in the combination key to the assigned position
   * @param positionEvidence hash containing the fulfilled MSn evidence for a position assignment - first key: id of the chain combination; second key: position where evidence is provided; values: rules that are fulfilled by MSn evidence
   * @param numberOfPositions at how many positions the fatty acids may be assigned
   * @param basePeakValues found values for the base peak
   * @param msnRetentionTimes the retention times of the used MSn spectra; the first key is the msLevel; second key scan number; value: retention time
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  @SuppressWarnings("unchecked")
  public LipidomicsMSnSet(LipidParameterSet set, int status, float mzTolerance, Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,IntensityRuleVO> headIntensityRules,
      Hashtable<String,Hashtable<String,CgProbe>> chainFragments, Hashtable<String,Hashtable<String,IntensityChainVO>> chainIntensityRules,
      Vector<String> validChainCombinations, Hashtable<String,Double> relativeIntensityOfCombination, Hashtable<String,Hashtable<Integer,Integer>> positionDefinition,
      Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> positionEvidence, int numberOfPositions, 
      Hashtable<Integer,Float> basePeakValues, Hashtable<Integer,LinkedHashMap<Integer,Float>> msnRetentionTimes,
      HydroxyEncoding faEncoding, HydroxyEncoding lcbEncoding) throws LipidCombinameEncodingException{
    this(set, status, mzTolerance, headGroupFragments, headIntensityRules, chainFragments, chainIntensityRules, validChainCombinations, relativeIntensityOfCombination, positionDefinition,
        positionEvidence, numberOfPositions, basePeakValues, msnRetentionTimes);
    //generate the human readable lookup and the combination relative areas;
    nameLookupHumReadableToPositionInsensitve_ = new LinkedHashMap<String,String>();
    nameLookupPositionSnNomenclature_ = new LinkedHashMap<String,String>();
    nameLookupPositionInsensitve_ = new LinkedHashMap<String,String>();
    ambiguousPositionIdentifications_ = new Hashtable<String,Vector<String>>();
    chainNameLookupHumanReadable_ = new Hashtable<String,String>();
    faEncoding_ = faEncoding;
    lcbEncoding_ = lcbEncoding;
    
    Hashtable<String,Integer> faChainOccurrences = new Hashtable<String,Integer>();
    if (status_==NO_MSN_PRESENT || status_==HEAD_GROUP_DETECTED){
      nameLookupHumReadableToPositionInsensitve_.put(getNameStringWithoutRt(),getNameStringWithoutRt());
      nameLookupPositionSnNomenclature_.put(getNameStringWithoutRt(), getNameStringWithoutRt());
    } else if (status_==FRAGMENTS_DETECTED || status_==POSITION_DETECTED){
      for (String combiName : validChainCombinations_){
        //this part is for calculating the relative intensities
        Vector<FattyAcidVO> chains = StaticUtils.decodeLipidNamesFromChainCombi(combiName);
        boolean ohPresent = StaticUtils.areThereOhInCombi(chains);
        for (FattyAcidVO combiChain : chains){
          if (!chainNameLookupHumanReadable_.containsKey(combiChain.getChainId()))
            chainNameLookupHumanReadable_.put(combiChain.getChainId(), StaticUtils.getHumanReadableChainName(combiChain, faEncoding, lcbEncoding, ohPresent));
          int count = 0;
          if (faChainOccurrences.containsKey(combiChain.getChainId())) count = faChainOccurrences.get(combiChain.getChainId());
          count++;
          faChainOccurrences.put(combiChain.getChainId(),count);
        }
      }
      List<DoubleStringVO> toSort = new ArrayList<DoubleStringVO>(); 
      for (String combi : relativeIntensityOfCombination_.keySet()) {
        toSort.add(new DoubleStringVO(combi,relativeIntensityOfCombination_.get(combi)));
      }
      Collections.sort(toSort,new GeneralComparator("at.tugraz.genome.lda.vos.DoubleStringVO", "getValue", "java.lang.Double"));
      //resort valid chain combinations from highest to lowest area
      validChainCombinations_ = new Vector<String>();
      for (int i=toSort.size()-1; i!=-1; i--)
        validChainCombinations_.add(toSort.get(i).getKey());      
    }
    
    
    //now the names are sorted by area -> I can buiild my lookyp hashes
    if (status_==FRAGMENTS_DETECTED || status_==POSITION_DETECTED){
      for (String combiName : validChainCombinations_){
      	nameLookupPositionInsensitve_.put(combiName, StaticUtils.getHumanReadableCombiName(combiName,faEncoding,lcbEncoding,false));
        if (status_==POSITION_DETECTED && positionDefinition_.containsKey(combiName)){
          Vector<String> posNames = getPositionSpecificCombiNames(combiName,positionDefinition_.get(combiName),faEncoding,lcbEncoding);
          Vector<String> snPosNames = getPositionSpecificCombiName(combiName,positionDefinition_.get(combiName),faEncoding,lcbEncoding);
          if (posNames.size()==1){
        	  nameLookupHumReadableToPositionInsensitve_.put(posNames.get(0),combiName);
        	  nameLookupPositionSnNomenclature_.put(snPosNames.get(0),combiName);
          }
          else {
            ambiguousPositionIdentifications_.put(combiName, posNames);
            for (String ambigName : posNames)
              nameLookupHumReadableToPositionInsensitve_.put(ambigName,combiName);
          }
          nameLookupPositionSnNomenclature_.put(snPosNames.get(0),combiName);
          /*if (posNames.size()==1) nameLookupHumReadableToPositionInsensitve_.put(posNames.get(0),combiName);
          else {
            ambiguousPositionIdentifications_.put(combiName, posNames);
            for (String ambigName : posNames)
              nameLookupHumReadableToPositionInsensitve_.put(ambigName,combiName);
          }*/
        }else{
          nameLookupHumReadableToPositionInsensitve_.put(StaticUtils.getHumanReadableCombiName(combiName,faEncoding,lcbEncoding),combiName);
          nameLookupPositionSnNomenclature_.put(StaticUtils.getHumanReadableCombiName(combiName,faEncoding,lcbEncoding),combiName);
        }
        //this part is for calculatating the relative intensities
        Vector<FattyAcidVO> chains = StaticUtils.decodeLipidNamesFromChainCombi(combiName);
        boolean ohPresent = StaticUtils.areThereOhInCombi(chains);
        for (FattyAcidVO combiChain : chains){
          if (!chainNameLookupHumanReadable_.containsKey(combiChain.getChainId()))
            chainNameLookupHumanReadable_.put(combiChain.getChainId(), StaticUtils.getHumanReadableChainName(combiChain, faEncoding, lcbEncoding, ohPresent));
          int count = 0;
          if (faChainOccurrences.containsKey(combiChain.getChainId())) count = faChainOccurrences.get(combiChain.getChainId());
          count++;
          faChainOccurrences.put(combiChain.getChainId(),count);
        }
      }
    }    
  }
  
    
  public LipidomicsMSnSet(LipidomicsMSnSet set) throws LipidCombinameEncodingException{
    this(set,set.status_,set.mzTolerance_,set.headGroupFragments_,set.headIntensityRules_,set.chainFragments_,set.chainIntensityRules_,
        set.validChainCombinations_,set.relativeIntensityOfCombination_,set.positionDefinition_,set.positionEvidence_,set.numberOfPositions_,set.basePeakValues_,
        set.getMsnRetentionTimes());
    this.nameLookupHumReadableToPositionInsensitve_ = set.nameLookupHumReadableToPositionInsensitve_;
    this.nameLookupPositionSnNomenclature_ = set.nameLookupPositionSnNomenclature_;
    this.ambiguousPositionIdentifications_ = set.ambiguousPositionIdentifications_;
    this.nameLookupPositionInsensitve_ = set.nameLookupPositionInsensitve_;
    this.chainNameLookupHumanReadable_ = set.chainNameLookupHumanReadable_;
    this.involvedFAs_ = set.involvedFAs_;
    this.msLevels_ = set.msLevels_;    
  }
  
  /**
   * Removes a molecular species identified by the human readable name from the object.
   * @param humanReadable
   */
  public void removeMolecularSpecies(String humanReadable)
  {
  	String humanReadableWithoutDB = StaticUtils.getHumanReadableWODoubleBondPositions(humanReadable);
  	String positionInsensitive = nameLookupHumReadableToPositionInsensitve_.get(humanReadableWithoutDB);
  	
  	nameLookupHumReadableToPositionInsensitve_.remove(humanReadableWithoutDB);
  	nameLookupPositionSnNomenclature_.remove(humanReadableWithoutDB);
  	nameLookupPositionInsensitve_.remove(positionInsensitive);
  	positionDefinition_.remove(positionInsensitive);
  	positionEvidence_.remove(positionInsensitive);
  	ambiguousPositionIdentifications_.remove(positionInsensitive);
  	chainIntensityRules_.remove(positionInsensitive); //chainIntensityRules may contain single fatty acids as well as combis
  	validChainCombinations_.remove(positionInsensitive);
  	
  	//recalculating the relative area contributions
  	Double multiplier = 1.0 / (1.0 - relativeIntensityOfCombination_.get(positionInsensitive));
  	System.out.println(multiplier);
  	relativeIntensityOfCombination_.remove(positionInsensitive);
  	for (String key : relativeIntensityOfCombination_.keySet()) 
  	{
  		Double rel = relativeIntensityOfCombination_.get(key);
  		Double newRel = rel * multiplier;
  		System.out.println(newRel);
  		relativeIntensityOfCombination_.put(key, newRel);
  	}
  	
  	//adjusting the omega C=C information
  	Vector<DoubleBondPositionVO> omegaVOs = getOmegaInformation();
  	Vector<DoubleBondPositionVO> newOmegaVOs = new Vector<DoubleBondPositionVO>();
  	for (DoubleBondPositionVO omegaVO : omegaVOs)
  	{
  		String noOmegaPos = omegaVO.getEncodedDetailed(true, false);
  		if (!noOmegaPos.equals(positionInsensitive))
  		{
  			newOmegaVOs.add(omegaVO);
  		}
  	}
  	setOmegaInformation(newOmegaVOs);
  	
  	//starting here fatty acid chains are removed
  	HashSet<String> potentiallyRemovedFA = new HashSet<String>(Arrays.asList(positionInsensitive.split(LipidomicsConstants.CHAIN_COMBI_SEPARATOR)));
  	HashSet<String> remainingFA = new HashSet<String>();
  	for (String validChainCombi : validChainCombinations_) //find still valid fatty acids
  	{
  		remainingFA.addAll(Arrays.asList(validChainCombi.split(LipidomicsConstants.CHAIN_COMBI_SEPARATOR)));
  	}
  	for (String valid : remainingFA) //remove all still valid fatty acids from the ones to be removed
  	{
  		potentiallyRemovedFA.remove(valid);
  	}
  	for (String removed : potentiallyRemovedFA) //remove invalid fatty acids
  	{
  		chainFragments_.remove(removed);
  		involvedFAs_.remove(removed);
  		chainNameLookupHumanReadable_.remove(removed);
  		chainIntensityRules_.remove(removed); //chainIntensityRules may contain single fatty acids as well as combis
  	}
  	
  	//updating status if needed
  	if (chainFragments_.isEmpty())
  	{
  		status_ = headGroupFragments_.isEmpty() ? NO_MSN_PRESENT : HEAD_GROUP_DETECTED;
  	}
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
   * TODO: this can be completely replaced with the SN position one.
   * @return the lipid names for the identified lipid species based on MSn evidence - the entry may be a String, or a Vector<String> if several position specific assignments are possible
   */
  public Vector<Object> getMSnIdentificationNames() {
    Vector<Object> names = new Vector<Object>();
    Set<String> usedCombis = new HashSet<String>();
    for (String humanReadable : nameLookupHumReadableToPositionInsensitve_.keySet()) {
      String uniqueCombiId = nameLookupHumReadableToPositionInsensitve_.get(humanReadable);
      if (usedCombis.contains(uniqueCombiId))
        continue;
      usedCombis.add(uniqueCombiId);
      if (ambiguousPositionIdentifications_.containsKey(uniqueCombiId)){
        names.add(ambiguousPositionIdentifications_.get(uniqueCombiId));
      }else{
        names.add(humanReadable);
      }
    }
    return names;
  }
  
  /**
   * @return the lipid names for the identified lipid species based on MSn evidence according to new nomenclature (sn positions in brackets)
   */
  public Vector<String> getMSnIdentificationNamesWithSNPositions() {
    Vector<String> names = new Vector<String>();
    for (String humanReadable : nameLookupPositionSnNomenclature_.keySet()) {
        names.add(humanReadable);
    }
    return names;
  }
  
  /**
   * returns the position-inspecific, human readable version of a chain combination
   * @param encoded
   * @return the the position-inspecific, human readable version of a chain combination
   */
  public String getPositionInsensitiveHumanReadableCombiName(String encoded) {
    return nameLookupPositionInsensitve_.get(encoded);
  }

  
  /**
   * for position evidence, returns all the position combinations that are possible with this rules
   * @param combiName the name of the FA combination without position assignments
   * @param positions; key: the position of the chain in the combiName; value the true position on the backbone
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @return all the position combinations that are possible with this rules
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private Vector<String> getPositionSpecificCombiNames(String combiName, Hashtable<Integer,Integer> positions,
      HydroxyEncoding faEncoding, HydroxyEncoding lcbEncoding) throws LipidCombinameEncodingException{
    Vector<String> unassignedFAs = new Vector<String>();
    Vector<FattyAcidVO> chains = StaticUtils.decodeLipidNamesFromChainCombi(combiName);
    boolean hydroxylationSites = StaticUtils.areThereOhInCombi(chains);
    Hashtable<Integer,String> definedPositions = new Hashtable<Integer,String>();
    for (int i=0;i!=chains.size();i++){
      if (!positions.containsKey(i) || positions.get(i)==-1) unassignedFAs.add(StaticUtils.getHumanReadableChainName(chains.get(i), faEncoding, lcbEncoding, hydroxylationSites));
      else definedPositions.put(positions.get(i), StaticUtils.getHumanReadableChainName(chains.get(i),faEncoding, lcbEncoding, hydroxylationSites));
    }
    return getPermutedChainPositionNames(definedPositions, unassignedFAs);
  }
  
  /**
   * for position evidence, returns the name of the FA combination according to new nomenclature
   * @param combiName the name of the FA combination without position assignments
   * @param positions; key: the position of the chain in the combiName; value the true position on the backbone
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @return all the position combinations that are possible with this rules
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private Vector<String> getPositionSpecificCombiName(String combiName, Hashtable<Integer,Integer> positions,
      HydroxyEncoding faEncoding, HydroxyEncoding lcbEncoding) throws LipidCombinameEncodingException{
    Vector<String> unassignedFAs = new Vector<String>();
    Vector<FattyAcidVO> chains = StaticUtils.decodeLipidNamesFromChainCombi(combiName);
    boolean hydroxylationSites = StaticUtils.areThereOhInCombi(chains);
    Hashtable<Integer,String> definedPositions = new Hashtable<Integer,String>();
    for (int i=0;i!=chains.size();i++){
      if (!positions.containsKey(i) || positions.get(i)==-1) unassignedFAs.add(StaticUtils.getHumanReadableChainName(chains.get(i), faEncoding, lcbEncoding, hydroxylationSites));
      else definedPositions.put(positions.get(i), StaticUtils.getHumanReadableChainName(chains.get(i),faEncoding, lcbEncoding, hydroxylationSites));
    }
    return getChainPositionNames(definedPositions, unassignedFAs);
  }
  
  /**
   * returns various names if only part of the chains could be correctly assigned
   * @param definedPositions positions where a definite assignment is given
   * @param unassignedFAs fatty acids where no assignment could be made
   * @return the lipid names that are possible with these assignments
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private Vector<String> getPermutedChainPositionNames(Hashtable<Integer,String> definedPositions, Vector<String> unassignedFAs) throws LipidCombinameEncodingException{
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
      Vector<String> permutedNames = StaticUtils.getPermutedChainNames(unassigned,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
      for (String permutedName : permutedNames){
        if (!combis.containsKey(permutedName)) combis.put(permutedName, permutedName);
      }
      for (String permutedName : combis.keySet()){
        Vector<String> parts =  StaticUtils.splitChainCombiToEncodedStrings(permutedName,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
        int permutedPos = 0;
        Hashtable<Integer,String> combiHash = new Hashtable<Integer,String>();
        for (int i=0; i!=(definedPositions.size()+unassigned.size());i++){
          String fa = null;
          if (definedPositions.containsKey(i)) fa = definedPositions.get(i);
          else{
            fa = parts.get(permutedPos);
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
        if (i!=0) name+= LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS;
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
   * returns names according to new nomenclature if only part of the chains could be correctly assigned
   * @param definedPositions positions where a definite assignment is given
   * @param unassignedFAs fatty acids where no assignment could be made
   * @return the lipid names that are possible with these assignments
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private Vector<String> getChainPositionNames(Hashtable<Integer,String> definedPositions, Vector<String> unassignedFAs) throws LipidCombinameEncodingException{
    Vector<Hashtable<Integer,String>> possibleCombis = new Vector<Hashtable<Integer,String>>();
    Hashtable<Integer,String> combiHash = new Hashtable<Integer,String>();
    Vector<String> names = new Vector<String>();
    
    if (unassignedFAs.size()==0){
      possibleCombis.add(definedPositions);
      String name = "";
      for (Hashtable<Integer,String> combi : possibleCombis){
          
          for (int i=0;i!=this.numberOfPositions_;i++){
            if (i!=0) name+= LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS;
            if (combi.containsKey(i))
              name += combi.get(i);
            else
              name += LipidomicsConstants.NO_FA_LINKED;
          }
          names.add(name);
        }
    }else{
      Vector<String> unassigned = new Vector<String>(unassignedFAs);
      for (int i=(definedPositions.size()+unassignedFAs.size()); i<numberOfPositions_;i++){
        unassigned.add("-");
      }
      
      int noPositionCount = 0;
      String chainSeparator = LipidomicsConstants.CHAIN_SEPARATOR_NO_POS;
      for (int i=0; i!=(definedPositions.size()+unassigned.size());i++){
          String fa = null;
          if (definedPositions.containsKey(i)) {
        	  if(numberOfPositions_ - definedPositions.size() > 1)	// if more than 1 unknown positions (e.g. if 1 position out of 2 positions is known, the other one should also be known
        		  fa = definedPositions.get(i) + " " + LipidomicsConstants.SN_POSITION_START + (i + 1) + LipidomicsConstants.SN_POSITION_END;
        	  else {
        		  fa = definedPositions.get(i);
        		  chainSeparator = LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS;
        	  }
          }
          else {
            fa = unassigned.get(noPositionCount);
            noPositionCount++;
          }
          combiHash.put(i, fa);
      }
      possibleCombis.add(combiHash);
      
      String name = "";
      for (Hashtable<Integer,String> combi : possibleCombis){
          
          for (int i=0;i!=this.numberOfPositions_;i++){
            if (i!=0) name+= chainSeparator;
            if (combi.containsKey(i))
              name += combi.get(i);
            else
              name += "-";
          }
          names.add(name);
        }
    }


    return names;
  }
  
  /**
   * 
   * @return the relative contribution of one fatty acid combination to the total area detected in MS1 - key is the fatty acid combination id
   */
  public Hashtable<String,Double> getChainCombinationRelativeAreas(){
    return this.relativeIntensityOfCombination_;
  }
  
  /**
   * for splitting of quantified area - returns relative contribution of identification [0...1] 
   * @param fullName the name of the identified lipid species
   * @return relative contribution of this species [0...1]
   */
  public double getRelativeIntensity(String fullName){
    double relativeIntensity = 0.0;
    if (nameLookupHumReadableToPositionInsensitve_.get(fullName) != null) {
    	if (status_ < FRAGMENTS_DETECTED) return 1.0;
    	relativeIntensity = relativeIntensityOfCombination_.get(nameLookupHumReadableToPositionInsensitve_.get(fullName));
    }
    return relativeIntensity;
  }
  
  /**
   * returns the position evidence information belonging to a certain lipid (molecular) species
   * @param fullName the species name (without position information)
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @return the position evidence information belonging to a certain lipid (molecular) species
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public Hashtable<Integer,Vector<IntensityPositionVO>> getPositionEvidence(String fullName,HydroxyEncoding faEncoding,HydroxyEncoding lcbEncoding) throws LipidCombinameEncodingException{
    Hashtable<Integer,Vector<IntensityPositionVO>> posEvidence = null;
    if (!this.positionEvidence_.containsKey(fullName)){
      for (String combiName : this.validChainCombinations_){
        if (!positionDefinition_.containsKey(combiName)) continue;
       Vector<String> posCombiNames = getPositionSpecificCombiNames(combiName,positionDefinition_.get(combiName),faEncoding,lcbEncoding);
        for (String posCombiName : posCombiNames){
          if (posCombiName.equalsIgnoreCase(fullName)) return positionEvidence_.get(combiName);
        }
      }      
    }else
      posEvidence = positionEvidence_.get(fullName);
    return posEvidence;
  }
  
  
  /**
   * returns the MS-levels where detections were made for a certain (molecular) lipid species - or all levles in the case of MSLEVEL_ALL_IDENTIFIER
   * @param fullName the molecular species name; or MSLEVEL_HEAD_IDENTIFIER, or MSLEVEL_ALL_IDENTIFIER
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @return the MS-levels where detections were made for a certain (molecular) lipid species - or all levles in the case of MSLEVEL_ALL_IDENTIFIER
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public Set<Integer> getMSLevels(String fullName,HydroxyEncoding faEncoding, HydroxyEncoding lcbEncoding) throws LipidCombinameEncodingException{
    Set<Integer> msLevels = null;
    String encodedCombi = fullName;
    if (nameLookupHumReadableToPositionInsensitve_.containsKey(fullName))
      encodedCombi = nameLookupHumReadableToPositionInsensitve_.get(fullName);
    if (fullName.equalsIgnoreCase(MSLEVEL_HEAD_IDENTIFIER)||fullName.equalsIgnoreCase(MSLEVEL_ALL_IDENTIFIER) ||
        msLevels_.containsKey(encodedCombi))
      msLevels = msLevels_.get(encodedCombi);
    else {
      for (String combiName : this.validChainCombinations_){
        if (!positionDefinition_.containsKey(combiName)) continue;
        Vector<String> posCombiNames = getPositionSpecificCombiNames(combiName,positionDefinition_.get(combiName),faEncoding,lcbEncoding);
        for (String posCombiName : posCombiNames){
          if (posCombiName.equalsIgnoreCase(fullName)){
            msLevels = msLevels_.get(combiName);
            break;
          }
        }
        if (msLevels!=null)
          break;
      }            
    }
    return msLevels;
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
   * @return the retention times of the used MSn spectra; the first key is the msLevel; second key scan number; value: retention time
   */
  public Hashtable<Integer,LinkedHashMap<Integer,Float>> getMsnRetentionTimes()
  {
    return msnRetentionTimes_;
  }

  /**
   * searches for a fragment for this rule and returns the area of the base peak
   * @param rule the VO containing the rule parameters
   * @return base peak area
   */
  public float getBasePeak(IntensityRuleVO rule)
  {
    CgProbe probe = rule.getAnyNonBasepeakFragment(headGroupFragments_,chainFragments_); 
    if (probe == null)
    {
    	System.out.println("hi!");
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
    for (String biggerName : ruleVO.getBiggerNonBasePeakNames().keySet()){
      Float biggerArea = getFragmentArea(ruleVO,biggerName,headGroupFragments,chainFragments);
      if (biggerArea == null) continue;
      String identifier = biggerName;
      if (ruleVO.getBiggerPosition()>0) identifier += "["+String.valueOf(ruleVO.getBiggerPosition())+"]";
      areas.put(identifier, biggerArea);
    }
    for (String smallerName : ruleVO.getSmallerNonBasePeakNames().keySet()){
      Float smallerArea = getFragmentArea(ruleVO,smallerName,headGroupFragments,chainFragments);
      if (smallerArea == null) continue;
      String identifier = smallerName;
      if (ruleVO.getSmallerPosition()>0) identifier += "["+String.valueOf(ruleVO.getSmallerPosition())+"]";
      areas.put(identifier, smallerArea);
    }
    return areas;
  }

  /**
   * returns the area for one fragment
   * @param intRule the IntensityRuleVO holding information about the found fragments
   * @param name name of the fragment
   * @param headGroupFragments the found head group fragments with its areas; key: fragment name; value: CgProbe peak identification object
   * @param chainFragments the found chain fragments with its areas; first key: fatty acid name; second key: name of the fragment; value: CgProbe peak identification object
   * @return area of the fragment
   */
  private static Float getFragmentArea(IntensityRuleVO intRule, String name, Hashtable<String,CgProbe> headGroupFragments, Hashtable<String,Hashtable<String,CgProbe>> chainFragments){
    Float area = null;
    if (name.equalsIgnoreCase(IntensityRuleVO.BASEPEAK_NAME)) return area;
    CgProbe fragment = intRule.checkForFragmentAvailability(name, headGroupFragments, chainFragments);
    if (fragment!=null) area = fragment.Area;
    else area = -1f;
    return area;    
  }
  
  /**
   * returns vector of unique position rule information so that it can be written - each entry is one unique rule
   * @param combiName
   * @return
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public Vector<IntensityRuleVO> getFAsInSequenceAsInRule(String combiName) throws LipidCombinameEncodingException{
    Vector<IntensityRuleVO> result = new Vector<IntensityRuleVO>();
    Hashtable<Integer,Vector<IntensityPositionVO>> posRules = positionEvidence_.get(combiName);
    Hashtable<Integer,Integer> posDef = positionDefinition_.get(combiName);
    if (status_ < POSITION_DETECTED || posDef==null) return result;
    Vector<String> fas = StaticUtils.splitChainCombiToEncodedStrings(combiName, LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
    Hashtable<String,IntensityRuleVO> usedRules = new Hashtable<String,IntensityRuleVO>();
    for (int i=0;i!=fas.size();i++){
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

  /**
   * 
   * @return involved FA objects for this identification
   */
  public Hashtable<String,FattyAcidVO> getInvolvedFAs()
  {
    return involvedFAs_;
  }
  
  public String getCombiIdFromHumanReadable(String humanReadable) {
    return this.nameLookupHumReadableToPositionInsensitve_.get(humanReadable);
  }
  
  public void setMsnRetentionTimes(Hashtable<Integer,LinkedHashMap<Integer,Float>> msnRetentionTimes_) {
		this.msnRetentionTimes_ = msnRetentionTimes_;
  }
  
  
  /**
   * Compares all dynamic fields of this class with another Object. 
   * Any Object including a position insensitive annotation cannot be directly compared, as the ordering of chains is random each time this annotation is created.
   * CgProbe Objects are as of now also not taken into account as they lack the 'equals' method.
   */
  @Override
  public boolean equals(Object obj)
  {
    if (this == obj)
      return true;
    if (!super.equals(obj))
      return false;
    if (getClass() != obj.getClass())
      return false;
    LipidomicsMSnSet other = (LipidomicsMSnSet) obj;

    return Objects.equals(ambiguousPositionIdentifications_,
        other.ambiguousPositionIdentifications_)
        && Objects.equals(basePeakValues_, other.basePeakValues_)
//        && Objects.equals(chainFragments_, other.chainFragments_)
        && Objects.equals(chainIntensityRules_, other.chainIntensityRules_)
        && Objects.equals(chainNameLookupHumanReadable_,
            other.chainNameLookupHumanReadable_)
//        && Objects.equals(headGroupFragments_, other.headGroupFragments_)
        && Objects.equals(headIntensityRules_, other.headIntensityRules_)
        && Objects.equals(involvedFAs_, other.involvedFAs_)
        && Objects.equals(msLevels_, other.msLevels_)
        && Objects.equals(msnRetentionTimes_, other.msnRetentionTimes_)
        && Float.floatToIntBits(mzTolerance_) == Float
            .floatToIntBits(other.mzTolerance_)
//        && Objects.equals(nameLookupHumReadableToPositionInsensitve_,
//            other.nameLookupHumReadableToPositionInsensitve_)
//        && Objects.equals(nameLookupPositionInsensitve_,
//            other.nameLookupPositionInsensitve_)
        && numberOfPositions_ == other.numberOfPositions_
        && Objects.equals(positionDefinition_, other.positionDefinition_)
//        && Objects.equals(positionEvidence_, other.positionEvidence_)
//        && Objects.equals(relativeIntensityOfCombination_,
//            other.relativeIntensityOfCombination_)
        && status_ == other.status_ && Objects.equals(validChainCombinations_,
            other.validChainCombinations_);
  }
  
  
  
  
  
  
  
  
  
  //TODO: from here on the methods for omega
  
  public LinkedHashMap<String,String> getNameLookupHumReadableToPositionInsensitive() {
  	return nameLookupHumReadableToPositionInsensitve_;
  }
  
  /**
   * @return the lipid names for the identified lipid species based on MSn evidence as a Set of Strings
   */
  public Set<String> getHumanReadableNameSet(){
    return new HashSet<String>(getMSnIdentificationNamesWithSNPositions());  
  }
  
  public Vector<Pair<String,String>> getLabeledUnlabeledPairs()
  {
  	Vector<Pair<String,String>> labeledUnlabeledPairs = new Vector<Pair<String,String>>();
  	try
  	{
  		Set<String> humanReadableNames = getHumanReadableNameSet();
    	for (String name : humanReadableNames)
    	{
    		Vector<FattyAcidVO> chains = new Vector<FattyAcidVO>();
    		Vector<String> chainNames = new Vector<String>();
    		String[] splitName = StaticUtils.splitChainCombinationsAtChainSeparators(name);
    		for (int i=0; i<splitName.length;i++)
    		{
    			FattyAcidVO fa = StaticUtils.decodeHumanReadableChain(splitName[i],Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(),false,null);
    			chains.add(fa);
    			chainNames.add(splitName[i]);
    		}
    		boolean ohPresent = StaticUtils.areThereOhInCombi(chains);
    		
    		Vector<String> sortedChainNamesWithPrefix = StaticUtils.sortChainNames(chainNames,ohPresent,false);
    		Vector<String> sortedChainNamesWithoutPrefix = StaticUtils.sortChainNames(chainNames,ohPresent,true);
    		
    		StringBuilder nameStringWithPrefix = new StringBuilder();
    		StringBuilder nameStringWithoutPrefix = new StringBuilder();
    		for (int i=0; i<splitName.length;i++)
    		{
    			if (i>0)
    			{
    				nameStringWithPrefix.append(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    				nameStringWithoutPrefix.append(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    			}
    			nameStringWithPrefix.append(sortedChainNamesWithPrefix.get(i));
  				nameStringWithoutPrefix.append(sortedChainNamesWithoutPrefix.get(i));
    		}
    		labeledUnlabeledPairs.add(new Pair<String,String>(nameStringWithoutPrefix.toString(),nameStringWithPrefix.toString()));
    	}
  	}
  	catch (LipidCombinameEncodingException ex)
		{
			ex.printStackTrace();
		}
  	return labeledUnlabeledPairs;
  }
  
  public Set<String> getPositionInsensitiveHumanReadableNameSet(boolean excludePrefix) 
  {
  	Set<String> positionInsensitive = new HashSet<String>();
  	//INFO: nameLookupPositionInsensitive_ does not contain validChainCombinations without detected fragments, relevant for LPC etc. therefore I cannot work with that.
  	for (String combiName : validChainCombinations_)
  	{
  		try
  		{
  			Vector<FattyAcidVO> chains = StaticUtils.decodeLipidNamesFromChainCombi(combiName);
  			boolean ohPresent = StaticUtils.areThereOhInCombi(chains);
  			Vector<String> chainNames = new Vector<String>();
  			for (FattyAcidVO chain : chains)
  			{
  				chainNames.add(StaticUtils.getHumanReadableChainName(
  						chain,Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(),ohPresent));
  			}
    		Vector<String> sortedChainNames = StaticUtils.sortChainNames(chainNames,ohPresent,excludePrefix);
    		String nameString = "";
    		for (String name : sortedChainNames)
    		{
    			if (nameString.length()>0) nameString+=LipidomicsConstants.CHAIN_SEPARATOR_NO_POS;
    			nameString+=name;
    		}
    		positionInsensitive.add(nameString);
  		}
  		catch (LipidCombinameEncodingException ex)
  		{
  			ex.printStackTrace();
  		}
  	}	
  	return positionInsensitive;
  }
  
  public String getHumanReadableChainCombi(String internalRepresentation, boolean isPositionWanted) throws LipidCombinameEncodingException
  {
  	if (status_==FRAGMENTS_DETECTED || status_==POSITION_DETECTED)
  	{
  		if (isPositionWanted && status_==POSITION_DETECTED && positionDefinition_.containsKey(internalRepresentation)){
  			Vector<String> posNames = getPositionSpecificCombiNames(internalRepresentation,positionDefinition_.get(internalRepresentation),faEncoding_,lcbEncoding_);
  			if (posNames.size()==1) return posNames.get(0);
  		}
  	}
  	return StaticUtils.getHumanReadableCombiName(internalRepresentation,faEncoding_,lcbEncoding_,false);
  }
  
  public Vector<String> getValidChainCombinations()
  {
  	return this.validChainCombinations_;
  }
  
  
  
//  /**
//   * returns the human readable chain name
//   * @param encoded the LDA encoded chain id
//   * @return the human readable chain name
//   */
//  public String getHumanReadableChain(String encoded) {
//    return chainNameLookupHumanReadable_.get(encoded);
//  }
}
