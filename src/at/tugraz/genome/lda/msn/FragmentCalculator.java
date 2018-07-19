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


import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.msn.vos.FragmentVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * This class calculates the FragmentVO objects that are used for the quantitation.
 * The input consists of information about the MS1 analyte, the rest of the required
 * information is fetched from the RulesContainer and the FattyAcidsContainer
 * @author Juergen Hartler
 *
 */
public class FragmentCalculator
{
  /** directory of stored fragmentation rules */
  private String rulesDir_;
  /** the name of the lipid class */
  private String ruleName_;
  /** the name of the analyte - containing the number of C atoms and double bonds */
  private String analyteName_;
  /** the chemical formula of the analyte (precursor) */
  private String analyteFormula_;
  /** the m/z value of the precursor */
  private double precursorMass_;
  /** hash table containing the permuted chain combinations - first key: identifier for combination; values: vector containing detailed info about the involved chains */
  private Hashtable<String,Vector<FattyAcidVO>> potentialChainCombinations_;
  /** Vector storing the various alkyl-alkenyl combinations possible by the fatty acids*/ 
  private Vector<Vector<Integer>> alkylAlkenylCombinations_;
  

  /**
   * constructor requiring information about the MS1 analyte -
   * fetch of rule information is immediately started
   * @param className name of the lipid class
   * @param modName name of the adduct
   * @param analyteName the name of the analyte - containing the number of C atoms and double bonds
   * @param analyteFormula  chemical formula of the analyte (precursor)
   * @param precursorMass m/z value of the precursor
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public FragmentCalculator(String className, String modName, String analyteName, String analyteFormula, double precursorMass) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    this(null,className,modName,analyteName,analyteFormula,precursorMass);
  }

  
  /**
   * constructor requiring information about the MS1 analyte -
   * fetch of rule information is immediately started
   * @param rulesDir directory containing the fragmentation rules
   * @param className name of the lipid class
   * @param modName name of the adduct
   * @param analyteName the name of the analyte - containing the number of C atoms and double bonds
   * @param analyteFormula  chemical formula of the analyte (precursor)
   * @param precursorMass m/z value of the precursor
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public FragmentCalculator(String rulesDir, String className, String modName, String analyteName, String analyteFormula, double precursorMass) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    this.rulesDir_ = rulesDir;
    this.ruleName_ = StaticUtils.getRuleName(className, modName);
    this.analyteName_ = analyteName;
    this.analyteFormula_ = analyteFormula;
    this.precursorMass_ = precursorMass;
    initCalculator();
  }
  
  /**
   * initializes the calculator by fetching the required information (rules, fatty acids) 
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  private void initCalculator() throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    potentialChainCombinations_ = new Hashtable<String,Vector<FattyAcidVO>>();
    int amountOfChains = Integer.parseInt(RulesContainer.getAmountOfChains(ruleName_,rulesDir_));
    if (amountOfChains>0){
      String chainLib = RulesContainer.getChainlibrary(ruleName_,rulesDir_);
      try{
        Hashtable<Integer,Hashtable<Integer,Hashtable<String,FattyAcidVO>>> availableChains = checkChainsForPlausibility(FattyAcidsContainer.getFattyAcidChains(chainLib));
        int cAtoms = getIntValueFromParsingRule(RulesContainer.getCAtomsFromNamePattern(ruleName_,rulesDir_), analyteName_, ruleName_, FragRuleParser.GENERAL_CATOMS_PARSE);
        int dbs = getIntValueFromParsingRule(RulesContainer.getDoubleBondsFromNamePattern(ruleName_,rulesDir_), analyteName_, ruleName_, FragRuleParser.GENERAL_DBOND_PARSE);
        potentialChainCombinations_ = extractPotentialChainCombinations(amountOfChains,cAtoms,dbs,availableChains);
        alkylAlkenylCombinations_ = alkylAlkenylPermutations(Integer.valueOf(RulesContainer.getAmountOfChains(ruleName_, rulesDir_)),
            Integer.valueOf(RulesContainer.getAmountOfAlkylChains(ruleName_,rulesDir_)),Integer.valueOf(RulesContainer.getAmountOfAlkenylChains(ruleName_, rulesDir_)));
      } catch (NoRuleException nrx){
        throw new RulesException("Error in rule \""+ruleName_+"\"! "+nrx.getMessage());
      }
    }
  }
  
  /**
   * calculates all potential chain combinations from a given precursor (number of C atoms and double bonds must be known)
   * @param chains how many chains has the object
   * @param cs number of total C atoms
   * @param dbs number of total double bonds
   * @param availableChains available fatty acid chains
   * @return hash table containing the permuted chain combinations - first key: identifier for combination; values: vector containing detailed info about the involved chains
   */
  private Hashtable<String,Vector<FattyAcidVO>> extractPotentialChainCombinations(int chains, int cs, int dbs, Hashtable<Integer,Hashtable<Integer,Hashtable<String,FattyAcidVO>>> availableChains){
    Hashtable<String,Vector<FattyAcidVO>> combinations = new Hashtable<String,Vector<FattyAcidVO>>();
    List<Integer> availableCs = new ArrayList<Integer>(availableChains.keySet());
    Hashtable<Integer,List<Integer>> availableCHash = new Hashtable<Integer,List<Integer>>();
    for (int i=0; i!=chains; i++) availableCHash.put((i+1), availableCs);
    Vector<Vector<Integer>> combis = getCombinations(cs, chains, availableCHash, new Vector<Integer>(), 0,false);
    // this hash table is for filtering permuted double entries
    Hashtable<String,String> permutedCombinations = new Hashtable<String,String>();
    int[] isotopicLabels = new int[chains];
    int[] divisors = new int[chains];
    for (Vector<Integer> cCombis : combis){
      Hashtable<Integer,List<Integer>> availableDbsHash = new Hashtable<Integer,List<Integer>>();
      boolean doubleBonds = true;
      if (availableChains.get(cCombis.get(0)).size()==1 && availableChains.get(cCombis.get(0)).keySet().iterator().next()==-1)
        doubleBonds = false;
      Vector<Vector<Integer>> dbCombis = new Vector<Vector<Integer>>();
      if (doubleBonds){
        for (int i=0; i!=chains; i++){
          availableDbsHash.put(chains-i, new ArrayList<Integer>(availableChains.get(cCombis.get(i)).keySet()));
        }
        dbCombis = getCombinations(dbs, chains, availableDbsHash, new Vector<Integer>(), 0, true);
      }else{
        Vector<Integer> noDbsCombis = new Vector<Integer>();
        for (int i=0; i!=chains; i++) noDbsCombis.add(-1);
        dbCombis.add(noDbsCombis);
      }
      for (Vector<Integer> dbCombi : dbCombis){
        //the same CAtoms:DbsCombination may be isotopically labelled
        //to support several labels, this a little bit more complicated code was introduced
        //first, it is evaluated how many different options are present for each chain position
        for (int i=0; i!=chains; i++){
          int cAtoms = cCombis.get(i);
          int dBonds = dbCombi.get(i);
          Hashtable<String,FattyAcidVO> faVOs = availableChains.get(cAtoms).get(dBonds);
          isotopicLabels[i] = faVOs.size();
        }
        //second, it is calculated how many combinations of the different label are possible
        //the divisor is necessary to use a modulo operation to get the correct index for each
        //fatty acid position. The different labels are indexed in the following way:
        //let's assume that we have 3 fatty acid positions, where the first can contain 4 different
        //labels, the second three, and the third two. The indexing would be as follows:
        //Iteration       index1     index2    index3
        //        0            0          0         0
        //        1            0          0         1
        //        2            0          1         0
        //        3            0          1         1
        //        4            0          2         0
        //        5            0          2         1
        //        6            1          0         0
        //       ....
        //       22            3          2         0
        //       23            3          2         1
        // Thus, the index i of each fatty acid chain i is calculated as (j/divisors[i])%isotopicLabels,
        // where divisors[i] is the multiplication of number of possibilities at the higher chain position
        // E.g for 22: the first index is (22/6=3)%4 -> 3
        //             the second index is (22/2=11)%3 -> 2
        //             the third index is (22/1=22)%2 -> 0

        int iterations = 1;
        for (int i=(isotopicLabels.length-1); i!=-1; i--){
          divisors[i] = iterations;
          iterations = iterations*isotopicLabels[i];
        }
        for (int j=0; j!=iterations; j++){
          String combiName = "";
          Vector<FattyAcidVO> valueVOs = new Vector<FattyAcidVO>();
          Vector<String> singleCombiParts = new Vector<String>();
          for (int i=0; i!=chains; i++){
            int cAtoms = cCombis.get(i);
            int dBonds = dbCombi.get(i);
            Hashtable<String,FattyAcidVO> faVOs = availableChains.get(cAtoms).get(dBonds);
            int currentIndex = (j/divisors[i])%isotopicLabels[i];
            int k=0;
            for (String prefix : faVOs.keySet()){
              if (k==currentIndex){
                String partName = StaticUtils.generateLipidNameString(prefix+String.valueOf(cAtoms), dBonds);
                combiName += partName+LipidomicsConstants.FA_SEPARATOR;
                valueVOs.add(faVOs.get(prefix));
                singleCombiParts.add(partName);
                break;
              }
              k++;
            }
          }
          combiName = combiName.substring(0,combiName.length()-1);
          //this is for filtering permuted double entries
          Vector<String> permutedNames = StaticUtils.getPermutedChainNames(singleCombiParts);
          boolean isThere = false;
          for (String permutedName : permutedNames){
            if (permutedCombinations.containsKey(permutedName)){
              isThere = true;
              break;
            }
          }
          if (!isThere){
            for (String permutedName : permutedNames) permutedCombinations.put(permutedName, permutedName);
            combinations.put(combiName, valueVOs);
          }
        }        
      }
    }
    return combinations;
  }
  
  
  /**
   * recursive method to distribute the total number of C atoms/double bonds to the chains (number of chains can be defined)
   * returns all possible combinations
   * @param expectedSum total number of C atoms - the sum that has to be reached by the individual chains
   * @param currentLevel how many chains have to be added still - next recursive call decreases this value until 1 is reached
   * @param values available C/double bond fragments; first key; 
   * @param prevComb combinations that are already assigned at this current level 
   * @param currentPos current position in the values hash - prevents double iterations
   * @param currentPosDoesNotMatter for double bond calculation the double iteration of the parts is desired
   * @return vector containing the combinations; the integers are the assigned values for the parts
   */
  private Vector<Vector<Integer>> getCombinations(int expectedSum, int currentLevel, Hashtable<Integer,List<Integer>> values, Vector<Integer> prevComb, int currentPos, boolean currentPosDoesNotMatter) {
    Vector<Vector<Integer>> combinations = new Vector<Vector<Integer>>();
    int startPosition = currentPos;
    if (currentPosDoesNotMatter) startPosition = 0;
    for (int i=startPosition; i<values.get(currentLevel).size(); i++){
      Integer value = values.get(currentLevel).get(i);
      Vector<Integer> comb = new Vector<Integer>(prevComb);
      int previousSum = 0;
      for (int prev : comb) previousSum += prev;
      if ((previousSum+value)>expectedSum) continue;
      comb.add(value);
      if (currentLevel>1){
        Vector<Vector<Integer>> returnedCombis = getCombinations(expectedSum,(currentLevel-1),values,comb,i,currentPosDoesNotMatter);
        combinations.addAll(returnedCombis);
      }else{
        if (expectedSum==(previousSum+value))
          combinations.add(comb);
      }
    }
    return combinations;
  }

  /**
   * 
   * @return the FragmentVOs for the head rules - the first key: are these fragments mandatory
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Hashtable<Boolean,Vector<FragmentVO>> getHeadFragments() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Hashtable<Boolean,Vector<FragmentVO>> allHeadFragments = new Hashtable<Boolean,Vector<FragmentVO>>();
    Vector<FragmentVO> mandatoryFragments = new Vector<FragmentVO>();
    Vector<FragmentVO> addFragments = new Vector<FragmentVO>();
    Hashtable<String,FragmentRuleVO> headRules = RulesContainer.getHeadFragmentRules(ruleName_,rulesDir_);
    for (FragmentRuleVO ruleVO : headRules.values()){
      //TODO: I do not know how to handle MS3 spectra - which precursor mass and formula should I use???
      Vector<Object> formulaAndMass = ruleVO.getFormulaAndMass(analyteFormula_, precursorMass_, null, 0,ruleVO.getCharge());
      FragmentVO fragVO = new FragmentVO(ruleVO.getName(),(Double)formulaAndMass.get(1),(String)formulaAndMass.get(0),ruleVO.getCharge(),ruleVO.getMsLevel(),
          ruleVO.isMandatory());
      if (ruleVO.isMandatory()==FragmentRuleVO.MANDATORY_TRUE || ruleVO.isMandatory()==FragmentRuleVO.MANDATORY_QUANT) mandatoryFragments.add(fragVO);
      else addFragments.add(fragVO);
    }
    
    allHeadFragments.put(true, mandatoryFragments);
    allHeadFragments.put(false, addFragments);
    return allHeadFragments;
  }
  
  /**
   * 
   * @return the rules for intensity comparisons of the head group
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Vector<IntensityRuleVO> getHeadIntensityRules() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    return RulesContainer.getHeadIntensityRules(ruleName_,rulesDir_);
  }
  
 /**  
  * 
  * @return the rules for intensity comparisons of the chain group
  * @throws RulesException specifies in detail which rule has been infringed
  * @throws NoRuleException thrown if the rules are not there
  * @throws IOException exception if there is something wrong about the file
  * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
  */
  private Vector<IntensityRuleVO> getChainIntensityRules() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    return RulesContainer.getChainIntensityRules(ruleName_,rulesDir_);
  }
  
  
  /**
   * 
   * @return the rules for intensity comparisons of the chain group
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Vector<IntensityRuleVO> getChainIntensityRulesSameChain() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Vector<IntensityRuleVO> rules = getChainIntensityRules();
    Vector<IntensityRuleVO> sameRules = new Vector<IntensityRuleVO>();
    for (IntensityRuleVO rule : rules){
      if (rule.getChainType()!=IntensityRuleVO.DIFF_CHAIN_TYPES) sameRules.add(rule);
    }
    return sameRules;
  }
  
  /**
   * 
   * @return the rules for intensity comparisons of the chain group
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Vector<IntensityRuleVO> getChainIntensityRulesDiffChain() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Vector<IntensityRuleVO> rules = getChainIntensityRules();
    Vector<IntensityRuleVO> diffRules = new Vector<IntensityRuleVO>();
    for (IntensityRuleVO rule : rules){
      if (rule.getChainType()==IntensityRuleVO.DIFF_CHAIN_TYPES) diffRules.add(rule);
    }    
    return diffRules;
  }

  /**
   * 
   * @return the rules for intensity comparisons for position detection
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Vector<IntensityRuleVO> getPositionIntensityRules() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    return RulesContainer.getPositionIntensityRules(ruleName_,rulesDir_);
  }

  /**
   * 
   * @return an head or chain fragment real by its name
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public FragmentRuleVO getFragmentRuleByName(String fragmentName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    Hashtable<String,FragmentRuleVO> headRules = RulesContainer.getHeadFragmentRules(ruleName_,rulesDir_);
    if (headRules.containsKey(fragmentName)) return headRules.get(fragmentName);
    Hashtable<String,FragmentRuleVO> chainRules = RulesContainer.getChainFragmentRules(ruleName_,rulesDir_);
    if (chainRules.containsKey(fragmentName)) return chainRules.get(fragmentName);
    throw new RulesException("There exists no fragment with the name\""+fragmentName+"\" for the analyte class\""+ruleName_+"\"!");
  }
  
  /**
   * checks for which MSn levels a base peak has to be extracted - only the ones containing $BASEPEAK need to be extracted
   * @param msLevels hash containing the available MSn levels in the data
   * @return Vector of MSn levels a base peak has to be extracted - only the ones containing $BASEPEAK need to be extracted
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Vector<Integer> getBasePeakIntRuleLevels(Hashtable<Integer,Boolean> msLevels) throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    boolean hasBasePeakCutoff = (getBasePeakCutoff()>0d);
    Hashtable<Integer,Integer> levels = new Hashtable<Integer,Integer>();
    if (hasBasePeakCutoff){
      for (FragmentRuleVO ruleVO : RulesContainer.getHeadFragmentRules(ruleName_,rulesDir_).values()){
        if (msLevels.containsKey(ruleVO.getMsLevel()) && msLevels.get(ruleVO.getMsLevel()))
          levels.put(ruleVO.getMsLevel(), ruleVO.getMsLevel());
      }
      for (FragmentRuleVO ruleVO : RulesContainer.getChainFragmentRules(ruleName_,rulesDir_).values()){
        if (msLevels.containsKey(ruleVO.getMsLevel()) && msLevels.get(ruleVO.getMsLevel()))
          levels.put(ruleVO.getMsLevel(), ruleVO.getMsLevel());
      }
    }else{
      getBasePeakIntRuleLevels(levels, getHeadIntensityRules(), msLevels);
      getBasePeakIntRuleLevels(levels, getChainIntensityRules(), msLevels);
      getBasePeakIntRuleLevels(levels, getPositionIntensityRules(), msLevels);
    }
    return new Vector<Integer>(levels.keySet());
  }

  /**
   * 
   * @return cutoff ratio value relative to the base peak 
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public double getBasePeakCutoff() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    return RulesContainer.getBasePeakCutoff(ruleName_,rulesDir_);
  }

  /**
   * 
   * @return cutoff value relative to the highest chain combination found
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public double getChainCutoff() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    return RulesContainer.getChainCutoff(ruleName_,rulesDir_);
  }
  
  /**
   * 
   * @return spectrum intensity coverage that has to be fulfilled 
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public double getSpectrumCoverageMin() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    return RulesContainer.getSpectrumCoverageMin(ruleName_,rulesDir_);
  }

  
  /**
   * based on the given rules, checks for which MSn levels a base peak has to be extracted, and writes them in the "levels" hash
   * @param levels !!! is the result: hash containing the MSn levels where a base peak has to be extracted
   * @param intRules the intensity rules to be checked for $BASEPEAK
   * @param msLevels hash containing the available MSn levels in the data
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  private void getBasePeakIntRuleLevels(Hashtable<Integer,Integer> levels, Vector<IntensityRuleVO> intRules, Hashtable<Integer,Boolean> msLevels) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    for (IntensityRuleVO intRule : intRules){
      if (intRule.containsBasePeak()){      
        int msLevel = getFragmentRuleByName(intRule.getAnyNonBasePeakName()).getMsLevel();
        if (msLevels.containsKey(msLevel) && msLevels.get(msLevel)) levels.put(msLevel, msLevel);
      }
    }

  }
  
  /**
   * 
   * @return all available fatty acid chain combinations
   */
  public Vector<FattyAcidVO> getPossibleChainObjects(){
    Hashtable<String,FattyAcidVO> fas = new Hashtable<String,FattyAcidVO>();
    for (Vector<FattyAcidVO> faChain : potentialChainCombinations_.values()){
      for (FattyAcidVO fa : faChain){
        if (!fas.containsKey(fa.getName())) fas.put(fa.getName(), fa);
      }
    }
    return new Vector<FattyAcidVO>(fas.values());
  }
  
  /**
   * calculates the FragmentVOs for a given fatty acid chain - the first key: fragment type; second key are these fragments mandatory
   * @param fa fatty acid chain for which the fragments have to be calculated
   * @return the FragmentVOs for a given fatty acid chain - the first key: are these fragments mandatory
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Hashtable<Integer,Hashtable<Boolean,Vector<FragmentVO>>> getChainFragments(FattyAcidVO fa) throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Hashtable<Integer,Hashtable<Boolean,Vector<FragmentVO>>> allChainFragments = new Hashtable<Integer,Hashtable<Boolean,Vector<FragmentVO>>>();
    Hashtable<String,FragmentRuleVO> chainRules =  RulesContainer.getChainFragmentRules(ruleName_,rulesDir_);
    for (FragmentRuleVO ruleVO : chainRules.values()){
      Hashtable<Boolean,Vector<FragmentVO>> chainFragments = new Hashtable<Boolean,Vector<FragmentVO>>();
      if (allChainFragments.containsKey(ruleVO.getChainType())) chainFragments = allChainFragments.get(ruleVO.getChainType());
      Vector<FragmentVO> mandatoryFragments = new Vector<FragmentVO>();
      Vector<FragmentVO> addFragments = new Vector<FragmentVO>();
      if (chainFragments.containsKey(true)) mandatoryFragments = chainFragments.get(true);
      if (chainFragments.containsKey(false)) addFragments = chainFragments.get(false);
      
      //TODO: I do not know how to handle MS3 spectra - which precursor mass and formula should I use???
      Vector<Object> formulaAndMass = ruleVO.getFormulaAndMass(analyteFormula_, precursorMass_, fa.getFormula(), fa.getMass(), ruleVO.getCharge());
      FragmentVO fragVO = new FragmentVO(ruleVO.getName(),(Double)formulaAndMass.get(1),(String)formulaAndMass.get(0),ruleVO.getCharge(),ruleVO.getMsLevel(),
          ruleVO.isMandatory());
      if (ruleVO.isMandatory()==FragmentRuleVO.MANDATORY_TRUE || ruleVO.isMandatory()==FragmentRuleVO.MANDATORY_QUANT) mandatoryFragments.add(fragVO);
      else addFragments.add(fragVO);
      
      chainFragments.put(true, mandatoryFragments);
      chainFragments.put(false, addFragments);
      allChainFragments.put(ruleVO.getChainType(), chainFragments);
    }    
    return allChainFragments;

  }
  
  /**
   * 
   * @param foundFAs the fatty acids detected by MS evidence
   * @param chainFragments the found fragments for the fatty acid chains
   * @return only fragment combinations where MS evidence was detected
   * @throws SpectrummillParserException if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws IOException exception if there is something wrong about the file
   * @throws NoRuleException thrown if the rules are not there
   * @throws RulesException specifies in detail which rule has been infringed
   */
  public Hashtable<String,Vector<String>> getChainFragmentCombinationsWithDetectedEvidence(Collection<String> foundFAs, Hashtable<String,Hashtable<String,CgProbe>> chainFragments) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    Hashtable<String,Vector<String>> potentialChainCombinations = new Hashtable<String,Vector<String>>();
    Hashtable<String,String> fas = new Hashtable<String,String>();
    for (String fa: foundFAs){
      fas.put(fa, fa);
    }
    Hashtable<String,FragmentRuleVO> chainRules =  RulesContainer.getChainFragmentRules(ruleName_,rulesDir_);
    boolean singleChainIdentification = RulesContainer.isSingleChainIdentification(ruleName_,rulesDir_);
    Vector<String> mandatoryFrags = new Vector<String>();
    for (FragmentRuleVO ruleVO : chainRules.values()){
      if (ruleVO.isMandatory()==FragmentRuleVO.MANDATORY_TRUE || ruleVO.isMandatory()==FragmentRuleVO.MANDATORY_QUANT) mandatoryFrags.add(ruleVO.getName());
    }
    for (String key : this.potentialChainCombinations_.keySet()){
      Vector<FattyAcidVO> combFAs = this.potentialChainCombinations_.get(key);
      Vector<Vector<String>> combis = getCorrectlyLinkedCombinations(combFAs);
      for (Vector<String> combi:combis){
        boolean allFAsThere = true;
        boolean oneFAIsThere = false;
        String newKey = "";
        Hashtable<String,Boolean> mandatoryFound = new Hashtable<String,Boolean>();
        for (String fragName : mandatoryFrags){
          mandatoryFound.put(fragName, false);
        }
        for (String faName : combi){
          if (!fas.containsKey(faName)){
            allFAsThere = false;
          } else if (singleChainIdentification) {
            if (mandatoryFrags.size()>0){
              if (!chainFragments.containsKey(faName)) continue;
              for (String fragName : chainFragments.get(faName).keySet()){
                if (mandatoryFound.containsKey(fragName)) mandatoryFound.put(fragName, true);
              }
            } else oneFAIsThere = true;
          } else {
          }
          newKey += faName+LipidomicsConstants.FA_SEPARATOR;
        }
        if (mandatoryFound.size()>0){
          boolean allTrue = true;
          for (Boolean found : mandatoryFound.values()){
            if (!found) allTrue = false;
          }
          if (allTrue) oneFAIsThere = true;
        }
        if (allFAsThere || (oneFAIsThere && singleChainIdentification)){
          newKey = newKey.substring(0,newKey.length()-1);
          potentialChainCombinations.put(newKey, combi);
        }
      }
    }    
    return potentialChainCombinations;
  }
  
  
  private Vector<Vector<String>> getCorrectlyLinkedCombinations(Vector<FattyAcidVO> combFAs) {
    Vector<Vector<String>> combis = new Vector<Vector<String>>();
    for (Vector<Integer> combi : this.alkylAlkenylCombinations_){
      Vector<String> oneCombi = new Vector<String>();
      for (int i=0; i!=combFAs.size(); i++){
        FattyAcidVO fa = combFAs.get(i);
        String name = FragmentCalculator.getAcylAlkylOrAlkenylName(fa.getName(), combi.get(i));
        oneCombi.add(name);
      }
      if (!isAlreadyInAlkylAcylCombination(oneCombi, combis)){
        combis.add(oneCombi);
      }
    }
   return combis;
  }
  
  private boolean isAlreadyInAlkylAcylCombination(Vector<String> oneCombi, Vector<Vector<String>> combis){
    boolean isAlreadyThere = false;
    for (Vector<String> otherCombi: combis){
      Vector<String> other = new Vector<String>(otherCombi);
      boolean isTheSame = true;
      for (String fa: oneCombi){
        int count = 0;
        boolean found = false;
        for (count=0;count!=other.size();count++){
          String otherFA = other.get(count);
          if (fa.equalsIgnoreCase(otherFA)){
            found = true;
            break;
          }
        }
        if (found) other.remove(count);
        else{
          isTheSame = false;
          break;
        }
      }
      if (isTheSame) isAlreadyThere = true;
    }
    return isAlreadyThere;
  }
  
  /**
   * returns all possible combinations where the alkyl and/or alkenyl chain(s) may be
   * @param chains number of fatty acid chains for this molecule
   * @param alkyl number of alkyl chains for this molecule
   * @param alkenyl number of alkenyl chains for this molecule
   * @return
   */
  private Vector<Vector<Integer>> alkylAlkenylPermutations (int chains, int alkyl, int alkenyl){
    Vector<Vector<Integer>> results = new Vector<Vector<Integer>>();
    if (alkyl<1 && alkenyl<1){
      Vector<Integer> combi = new Vector<Integer>();
      for (int i=0; i!=chains;i++)combi.add(FragmentRuleVO.ACYL_CHAIN);
      results.add(combi);
    }else{
      for (int i=0;i!=alkyl;i++){
        for (int j=0; j!=chains; j++){
          Vector<Vector<Integer>> otherPermutations = alkylAlkenylPermutations(chains-1,alkyl-1,alkenyl);
          for (Vector<Integer> perms : otherPermutations){
            Vector<Integer> combi = new Vector<Integer>();
            for (int k=0;k!=chains;k++){
              if (k<j) combi.add(perms.get(k));
              else if (k==j) combi.add(FragmentRuleVO.ALKYL_CHAIN);
              else if (k>j) combi.add(perms.get(k-1));
            }
            results.add(combi);
          }
        }
      }
      for (int i=0;i!=alkenyl;i++){
        for (int j=0; j!=chains; j++){
          Vector<Vector<Integer>> otherPermutations = alkylAlkenylPermutations(chains-1,alkyl,alkenyl-1);
          for (Vector<Integer> perms : otherPermutations){
            Vector<Integer> combi = new Vector<Integer>();
            for (int k=0;k!=chains;k++){
              if (k<j) combi.add(perms.get(k));
              else if (k==j) combi.add(FragmentRuleVO.ALKENYL_CHAIN);
              else if (k>j) combi.add(perms.get(k-1));
            }
            results.add(combi);
          }
        }
      }
    }
    // this is for removal of identical combinations
    Hashtable<String,String> usedCombis = new Hashtable<String,String>();
    Vector<Vector<Integer>> cleanedResults = new Vector<Vector<Integer>>();
    for (Vector<Integer> combi:results){
      String key = "";
      for (Integer type:combi) key += String.valueOf(type)+LipidomicsConstants.FA_SEPARATOR;
      if (key.length()>0) key = key.substring(0,key.length()-1);
      if (usedCombis.containsKey(key)) continue;
      usedCombis.put(key, key);
      cleanedResults.add(combi);
    }
    return cleanedResults;
  }
  
  /**
   * 
   * @param key identifier for a chain combination
   * @return detailed FattyAcidVOs for a combination identifier 
   */
  public Vector<FattyAcidVO> getChainFragmentCombination(String key){
    return this.potentialChainCombinations_.get(key);
  }
  
  /**
   * 
   * @return the amount of chains for this class
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public int getAmountOfChains() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    return Integer.parseInt(RulesContainer.getAmountOfChains(ruleName_,rulesDir_));
  }
  
  /**
   * 
   * @return lowest and highest spectrum level required for MS2 rules
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */  
  public int[] getSpectrumLevelRange() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    return RulesContainer.getSpectrumLevelRange(ruleName_,rulesDir_);
  }
  
  /**
   * this method checks whether there are MS1 probes found for available ms levels, and whether there are fragments defined for this level
   * when both checks are OK, this msLevel is regarded valid
   * @param oldMsLevels the msLevels to check
   * @param probesWithMSnSpectra the identified probes
   * @return the checked msLevels for further rule checks
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Hashtable<Integer,Boolean> correctMsLevelsForExistingFragments(Hashtable<Integer,Boolean> oldMsLevels, Hashtable<Integer,Vector<CgProbe>> probesWithMSnSpectra) throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Hashtable<Integer,Boolean> usableLevels = new Hashtable<Integer,Boolean>();
    for (Integer msLevel : oldMsLevels.keySet()){
      if (!probesWithMSnSpectra.containsKey(msLevel)) continue;
      boolean found = false;
      for (FragmentRuleVO ruleVO : RulesContainer.getHeadFragmentRules(ruleName_,rulesDir_).values()){
        if (ruleVO.getMsLevel()==msLevel){
          found = true;
          break;
        }
      }
      for (FragmentRuleVO ruleVO : RulesContainer.getChainFragmentRules(ruleName_,rulesDir_).values()){
        if (ruleVO.getMsLevel()==msLevel){
          found = true;
          break;
        }
      }
      if (found) usableLevels.put(msLevel, true);
    }
    return usableLevels;
  }

  
  /**
   * returns the acyl/alkyl/alkenyl common name for an arbitrary fatty acid chain
   * @param faName name of fatty acid 
   * @param chainType type of chain ACYL_CHAIN/ALKYL_CHAIN/ALKENYL_CHAIN
   * @return the acyl/alkyl/alkenyl common name for an arbitrary fatty acid chain
   */
  public static String getAcylAlkylOrAlkenylName(String faName, Integer chainType){
    String chainName = new String(faName);
    if (chainType==FragmentRuleVO.ALKYL_CHAIN) chainName = LipidomicsConstants.ALKYL_PREFIX+chainName;
    else if (chainType==FragmentRuleVO.ALKENYL_CHAIN) chainName = LipidomicsConstants.ALKENYL_PREFIX+chainName;
    return chainName;
  }
  
  /**
   * 
   * @return the amount of chains for this class
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public int getAllowedChainPositions() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    return RulesContainer.getAllowedChainPositions(ruleName_,rulesDir_);
  }

  
  /**
   * This method takes a Java regular expression and a lipid species name to extract from this
   * name the number of carbon atoms or double bonds; the lClass and input key are for the error message
   * if there is something wrong with the rule
   * @param rule the Java regular expression
   * @param analyte the name of the lipid species
   * @param lClass the name of the lipid class (for error message)
   * @param inputKey the name of the rule
   * @return number of carbon atoms or double bonds
   * @throws RulesException thrown if there is somehting wrong
   */
  public static int getIntValueFromParsingRule(String rule, String analyte, String lClass, String inputKey) throws RulesException {
    Pattern cAtomsPattern =  Pattern.compile(rule);
    Matcher cAtomsMatcher = cAtomsPattern.matcher(analyte);
    if (!cAtomsMatcher.matches()) throw new RulesException("The analyte "+analyte+" does not match the "+inputKey+" pattern \""+rule+"\" of the class "+lClass+"!");
    String valueString  = cAtomsMatcher.group(1);
    try{
      int value = Integer.parseInt(valueString);
      return value;
    } catch (NumberFormatException nfx){
      throw new RulesException("The rule \""+rule+"\" of the class "+lClass+" is not correct, since it returns for the "+analyte+" the non-integer value \""+valueString+"\"!");
    }
  }
  
  /**
   * checks whether the chains might be possible for the given analyte (e.g. no deuterated chains are possible for lipids containing no deuterium in their chemical formula)
   * @param possibleChains the chains fromt the Excel file
   * @return the chains that shall be checked
   * @throws RulesException RulesException thrown if there is somehting wrong
   */
  private Hashtable<Integer,Hashtable<Integer,Hashtable<String,FattyAcidVO>>> checkChainsForPlausibility(Hashtable<Integer,Hashtable<Integer,Hashtable<String,FattyAcidVO>>> possibleChains) throws RulesException{
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,FattyAcidVO>>> chains = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,FattyAcidVO>>>();
    try{
      Hashtable<String,Integer> formulaAmounts = StaticUtils.categorizeFormula(this.analyteFormula_);
      for (Integer cAtoms : possibleChains.keySet()){
        Hashtable<Integer,Hashtable<String,FattyAcidVO>> sameCAtoms = new Hashtable<Integer,Hashtable<String,FattyAcidVO>>();
        for (Integer dbs : possibleChains.get(cAtoms).keySet()){
          Hashtable<String,FattyAcidVO> sameDbs = new Hashtable<String,FattyAcidVO>();
          for (String prefix : possibleChains.get(cAtoms).get(dbs).keySet()){
            FattyAcidVO fa = possibleChains.get(cAtoms).get(dbs).get(prefix);
            Hashtable<String,Integer> faElements = StaticUtils.categorizeFormula(fa.getFormula());
            boolean isOk = true;
            for (String element : faElements.keySet()){
              if (!formulaAmounts.containsKey(element) || formulaAmounts.get(element)<faElements.get(element)){
                isOk = false;
                break;
              }
            }
            if (isOk) sameDbs.put(prefix, fa);
          }
          if (sameDbs.size()>0) sameCAtoms.put(dbs, sameDbs);
        }
        if (sameCAtoms.size()>0) chains.put(cAtoms, sameCAtoms);
      }
    } catch (ChemicalFormulaException e) {
      throw new RulesException("The formula "+analyteFormula_+" contains fragments that have not been defined before! Fragments have to be defined in previous columns, before they can be used!");
    }
    return chains;
  }
}

