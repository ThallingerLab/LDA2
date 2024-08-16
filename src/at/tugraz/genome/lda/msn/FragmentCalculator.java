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

package at.tugraz.genome.lda.msn;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;
import at.tugraz.genome.lda.msn.utils.RulesUtils;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.msn.vos.FragmentVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.parser.ModificationParser;
import at.tugraz.genome.lda.utils.RangeInteger;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
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
  
  /** the string used to split possible distribution of OH between FA and LCB chains*/
  private final static String OH_COMBI_SEPARATOR = "-";
  /** used to split the info in one OH frequency identifier*/ 
  private final static String OH_FREQUENCY_ID_SEPRATOR = ":";
  
  /** directory of stored fragmentation rules */
  private String rulesDir_;
  /** the name of the lipid class */
  private String ruleName_;
  /** the name of the analyte - containing the number of C atoms and double bonds */
  private String analyteName_;
  /** if there is isotopic labeling used, this part contains the labeling encoding */
  private String labelInName_;
  /** the chemical formula of the analyte (precursor) */
  private String analyteFormula_;
  /** the chemical formula of the analyte (precursor) - without deductions by ionization modifications */
  private String analyteFormulaWODeducts_;
  /** the m/z value of the precursor */
  private double precursorMz_;
  /** the charge of the precursor */
  private int precursorCharge_;
  /** (total) number of hydroxylation sites present on the molecule*/
  private int ohNumber_;
  /** information about the chemical elements*/
  private ElementConfigParser elements_;
  /** hash table containing the permuted chain combinations - first key: identifier for combination; values: vector containing detailed info about the involved chains */
  private Hashtable<String,Hashtable<String,Vector<FattyAcidVO>>> potentialChainCombinations_;
  /** all the chains that are possible after the combinatorial check; key is the chain id that contains information about chain type and hydroxylation sites*/
  private Hashtable<String,FattyAcidVO> availableChains_;
  /** int[2] containing all potential hydroxylation combinations of FA and LCBs; int[0] FA hydroxylation sites; int[1] LCB hydroxylation sites*/
  private Vector<int[]> possibleOhCombinations_;
  /** this hash table stores which hydroxylations are possible/allowed for FA chains; key: the chain combination; value: Vector containing a Vector with the possibly hydroxylation combinations*/
  private Hashtable<String,Vector<Vector<Integer>>> allowedFaHydroxylationsCombinations_;
  /** this hash table stores which hydroxylations are possible/allowed for LCB chains; key: the chain combination; value: Vector containing a Vector with the possibly hydroxylation combinations*/
  private Hashtable<String,Vector<Vector<Integer>>> allowedLcbHydroxylationsCombinations_;
  /** the principally (from the chemical formula) possible FA chains*/
  private Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableFAChainsBeforeCombiCheck_;
  /** the principally (from the chemical formula) possible FA chains*/
  private Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkylChainsBeforeCombiCheck_;
  /** the principally (from the chemical formula) possible FA chains*/
  private Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkenylChainsBeforeCombiCheck_;
  /** the principally (from the chemical formula) possible FA chains*/
  private Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableLCBChainsBeforeCombiCheck_;
  /** the available isotopic labels*/
  private Hashtable<String,Integer> availableLabels_;
  /** the available isotopic labels*/
  private Hashtable<String,String> singleLabelLookup_;
  /** the single labels*/
  private Set<String> availableSingleLabels_;
  /** how many labels are allowed in the chains*/
  private Hashtable<String,Integer> allowedLabelsInChains_;
  
  private String analyteOxState_;
 
  /**
   * constructor requiring information about the MS1 analyte -
   * fetch of rule information is immediately started
   * @param rulesDir directory containing the fragmentation rules
   * @param className name of the lipid class
   * @param modName name of the adduct
   * @param analyteName the name of the analyte - containing the number of C atoms and double bonds
   * @param analyteFormula  chemical formula of the analyte (precursor)
   * @param formulaWoDeducts chemical formula of the analyte (precursor), but modifications causing element reductions are not counted
   * @param precursorMz m/z value of the precursor
   * @param precursorCharge charge of the precursor
   * @param ohNumber (total) number of hydroxylation sites present on the molecule
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   */
  public FragmentCalculator(String rulesDir, String className, String modName, String analyteName, String analyteFormula, String formulaWoDeducts,
      double precursorMz, int precursorCharge, int ohNumber, String analyteOxState) throws RulesException, NoRuleException, IOException, SpectrummillParserException, HydroxylationEncodingException, ChemicalFormulaException {
    this.rulesDir_ = rulesDir;
    this.ruleName_ = StaticUtils.getRuleName(className, modName);
    this.analyteName_ = analyteName;
    this.analyteFormula_ = analyteFormula;
    this.analyteFormulaWODeducts_ = formulaWoDeducts;
    this.precursorMz_ = precursorMz;
    this.precursorCharge_ = precursorCharge;
    this.ohNumber_ = ohNumber;
    if (ohNumber_<0)
      this.ohNumber_ = 0;
    this.analyteOxState_ = analyteOxState;
    
    initCalculator();
  }
  
  /**
   * initializes the calculator by fetching the required information (rules, fatty acids) 
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   */
  private void initCalculator() throws RulesException, NoRuleException, IOException, SpectrummillParserException, HydroxylationEncodingException, ChemicalFormulaException {
    elements_ = Settings.getElementParser();
    int amountOfChains = Integer.parseInt(RulesContainer.getAmountOfChains(ruleName_,rulesDir_));
    int[] chainAmounts = RulesUtils.getAmountOfChainsCategorized(rulesDir_, ruleName_);
    int fattyChains = chainAmounts[0];
    int lcbChains = chainAmounts[1];
    int acylChains = chainAmounts[2];
    int alkylChains = chainAmounts[3];
    int alkenylChains = chainAmounts[4];
    if (amountOfChains>0){
      try{
        String chainLib = null;
        if (fattyChains>0) chainLib = RulesContainer.getChainlibrary(ruleName_,rulesDir_);
        String lcbLib = null;
        if (lcbChains>0) lcbLib = RulesContainer.getLcbLibrary(ruleName_,rulesDir_);
        calculatePossibleFaLcbHydroxyCombinations(fattyChains,lcbChains);
        int cAtoms = getIntValueFromParsingRule(RulesContainer.getCAtomsFromNamePattern(ruleName_,rulesDir_), analyteName_, ruleName_, FragRuleParser.GENERAL_CATOMS_PARSE);
        int dbs = getIntValueFromParsingRule(RulesContainer.getDoubleBondsFromNamePattern(ruleName_,rulesDir_), analyteName_, ruleName_, FragRuleParser.GENERAL_DBOND_PARSE);
        if (LipidomicsConstants.checkChainLabelCombination()) {
          this.labelInName_ = this.analyteName_.substring(0,this.analyteName_.indexOf(String.valueOf(cAtoms)));
        }
        
        checkChainsForPlausibility(chainLib,lcbLib,cAtoms,dbs,alkylChains>0,alkenylChains>0);
        
        
//        for (String label : this.availableSingleLabels_) {
//          System.out.println("Single Label: "+label);
//        }
//        for (String key : this.availableLabels_.keySet()) {
//          System.out.println("Label: "+key+"-"+this.availableLabels_.get(key));
//        }

/*
        System.out.println(this.analyteName_);
        System.out.println("FA-chains");
        if (availableFAChainsBeforeCombiCheck_.get(0) != null)
        for (int c : this.availableFAChainsBeforeCombiCheck_.get(0).keySet()) {
          for (int db : this.availableFAChainsBeforeCombiCheck_.get(0).get(c).keySet()){
            for (String prefix: this.availableFAChainsBeforeCombiCheck_.get(0).get(c).get(db).keySet()) {
              System.out.println(prefix+c+":"+db);
            }
          }
        }
        System.out.println("LCB-chains");
        if (availableLCBChainsBeforeCombiCheck_.get(2) != null)
        for (int c : this.availableLCBChainsBeforeCombiCheck_.get(2).keySet()) {
          for (int db : this.availableLCBChainsBeforeCombiCheck_.get(2).get(c).keySet()){
            for (String prefix: this.availableLCBChainsBeforeCombiCheck_.get(2).get(c).get(db).keySet()) {
              System.out.println(prefix+c+":"+db);
            }
          }
        }
*/
        extractPotentialChainCombinations(amountOfChains,fattyChains,lcbChains,acylChains,alkylChains,alkenylChains,cAtoms,dbs);
      } catch (NoRuleException nrx){
        throw new RulesException("Error in rule \""+ruleName_+"\"! "+nrx.getMessage());
      }
    }
  }
  
  /**
   * calculates all potential chain combinations from a given precursor (number of C atoms and double bonds must be known)
   * @param chainsTotal how many chains has the object
   * @param faChains how many total acyl/alkyl/alkenyl chains has this class
   * @param lcbChains how many total LCB chains has this class
   * @param acylChains how many acyl chains has this class
   * @param alkylChains how many alkylated chains has this class
   * @param alkenylChains how many alkenylated chains has this class
   * @param cs number of total C atoms
   * @param dbs number of total double bonds
   * @param availableChains available fatty acid chains
   * @throws RulesException exception thrown when there is something wrong with the chain name
   */
  // there is a 4-step procedure to generate all of the permutations
  // 1. DONE: extract all available chain names from the provided libraries (FA and LCB)
  // 2. DONE: out of these, calculate all possible combinations - the same thing it does right now
  // 3. DONE: go over the chains extracted in checkChainsForPlausibility, and check whether these combinations are actually possible (e.g. LCB-lib does not allow for
  //    this C atoms and DBs) - remove all of the implausible hits
  // 4. DONE: extract a list of possible chain names (aggregate of FA and LCB)
  private void extractPotentialChainCombinations(int chainsTotal, int faChains, int lcbChains, int acylChains, int alkylChains,
      int alkenylChains, int cs, int dbs) throws RulesException{
    // step one: extract all available chain names from the provided libraries (FA and LCB)
	  Hashtable<Integer,Hashtable<Integer,Hashtable<String,String>>> possibleCAtomsDbsOxs = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,String>>>();
    if (availableFAChainsBeforeCombiCheck_!=null)
      this.addUniqueChainDbsCombis(possibleCAtomsDbsOxs, availableFAChainsBeforeCombiCheck_);
    if (availableLCBChainsBeforeCombiCheck_!=null)
      this.addUniqueChainDbsCombis(possibleCAtomsDbsOxs, availableLCBChainsBeforeCombiCheck_);
    // step two: calculate all possible chain combinations respecting C atoms and double bonds (OH number and isotopes are not respected)
    Vector<String> potentialCombinations = calculatePotentialChainCombinations(chainsTotal,cs,dbs,possibleCAtomsDbsOxs);
    // step three: go over all combinations of possibleOhCombinations_, check which chain combinations are actually possible by
    // respecting the FA and LCB chain constraints, and the actually present OH numbers 
    availableChains_ = new Hashtable<String,FattyAcidVO>();
    potentialChainCombinations_ = payAttentionToActuallyPresentFaAndLcbChains(chainsTotal, faChains, lcbChains, acylChains, alkylChains, alkenylChains,
        potentialCombinations,  possibleOhCombinations_,availableFAChainsBeforeCombiCheck_,availableLCBChainsBeforeCombiCheck_,
        availableAlkylChainsBeforeCombiCheck_,availableAlkenylChainsBeforeCombiCheck_,availableChains_);
  }
  
  /**
   * computes all potential permutations of C atoms and double bonds that fulfill the total C atom and double bond number of the species for a given number of chains
   * @param chains the number of chains
   * @param cs the total number of C atoms
   * @param dbs the total number of double bonds
   * @param possibleCAtomsDbsOxs the possible C atom and double bond numbers without respecting the origin (FA or LCB)
   * @return potential permutations of chains that fulfil the total C atom and double bond number of the species for a given number of chains
   */
  private Vector<String> calculatePotentialChainCombinations(int chains, int cs, int dbs, Hashtable<Integer,Hashtable<Integer,Hashtable<String,String>>> possibleCAtomsDbsOxs){
    Vector<String> combinations = new Vector<String>();
    List<Integer> availableCs = new ArrayList<Integer>(possibleCAtomsDbsOxs.keySet());
    Hashtable<Integer,List<Integer>> availableCHash = new Hashtable<Integer,List<Integer>>();
    for (int i=0; i!=chains; i++) availableCHash.put((i+1), availableCs);
    Vector<Vector<Integer>> combis = getCombinations(cs, chains, availableCHash, new Vector<Integer>(), 0,false);
    // this hash table is for filtering permuted double entries
    Hashtable<String,String> permutedCombinations = new Hashtable<String,String>();
    for (Vector<Integer> cCombis : combis){
      Hashtable<Integer,List<Integer>> availableDbsHash = new Hashtable<Integer,List<Integer>>();
      boolean doubleBonds = true;
      if (possibleCAtomsDbsOxs.get(cCombis.get(0)).size()==1 && possibleCAtomsDbsOxs.get(cCombis.get(0)).keySet().iterator().next()==-1)
        doubleBonds = false;
      Vector<Vector<Integer>> dbCombis = new Vector<Vector<Integer>>();
      if (doubleBonds){
        for (int i=0; i!=chains; i++){
          availableDbsHash.put(chains-i, new ArrayList<Integer>(possibleCAtomsDbsOxs.get(cCombis.get(i)).keySet()));
        }
        dbCombis = getCombinations(dbs, chains, availableDbsHash, new Vector<Integer>(), 0, true);
      }else{
        Vector<Integer> noDbsCombis = new Vector<Integer>();
        for (int i=0; i!=chains; i++) noDbsCombis.add(-1);
        dbCombis.add(noDbsCombis);
      }
      
      for (Vector<Integer> dbCombi : dbCombis){
          
    	    //TODO: the next few lines are required for oxidized lipids (including the for), but take an enormous amount of time - excluded until fixed
////        List<String> oxStates = new ArrayList<>(possibleCAtomsDbsOxs.get(cCombis.get(0)).get(dbCombi.get(0)).keySet());
////    	    Set<List<String>> oxCombis = getOxCombinations(oxStates,chains);
////        Set<List<String>> validOxCombis = ValidateCombinations(oxCombis); 
		  
////        for(List<String> oxCombi : validOxCombis) {
			          String combiName = "";
			          Vector<String> singleCombiParts = new Vector<String>();
			          for (int i=0; i!=chains; i++){
			            int cAtoms = cCombis.get(i);
			            int dBonds = dbCombi.get(i);
			            String oxState = null;
			            //TODO: the next line is required for oxidized lipids, but take an enormous amount of time - excluded until fixed
			            ////oxState = oxCombi.get(i);
			            String partName = StaticUtils.generateLipidNameString(String.valueOf(cAtoms), dBonds,-1,oxState);
			            combiName += partName+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS;
			            singleCombiParts.add(partName);
	          
			          }
			          combiName = combiName.substring(0,combiName.length()-1);
			          //this is for filtering permuted double entries
			          Vector<String> permutedNames = StaticUtils.getPermutedChainNames(singleCombiParts,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
			          boolean isThere = false;
			          for (String permutedName : permutedNames){
			            if (permutedCombinations.containsKey(permutedName)){
			              isThere = true;
			              break;
			            }
			          }
			          if (!isThere){
			            for (String permutedName : permutedNames) permutedCombinations.put(permutedName, permutedName);
			            combinations.add(combiName);
			          }
			        } 
    //TODO: the next line is required for oxidized lipids, but take an enormous amount of time - excluded until fixed
 ////      }
    }
    return combinations;
  }
  
  
  /**
   * this method builds the final combination hashes where the chain origin (whether it is an FA or LCB) is respected and the amount of OH numbers 
   * @param chainsTotal the total number of chains
   * @param faChains the number of FA chains
   * @param lcbChains the number of LCB chains
   * @param acylChains how many acylated chains has this class
   * @param alkylChains how many alkylated chains has this class
   * @param alkenylChains how many alkenylated chains has this class
   * @param potentialCombinations potential permutations of chains that fulfil the total C atom and double bond number of the species for a given number of chains
   * @param ohCombis the possible OH-distribution between FA and LCB chains
   * @param availableFAChains the principally (from the chemical formula) possible acyl chains 
   * @param availableLCBChains the principally (from the chemical formula) possible LCB chains
   * @param availableAlkylChains the principally (from the chemical formula) possible alkyl chains 
   * @param availableAlkenylChains the principally (from the chemical formula) possible alkenyl chains
   * @param uniqueChains all the chains that are possible after the combinatorial check; key is the chain id that contains information about chain type and hydroxylation sites
   * @return hash table containing the permuted chain combinations - first key: identifier for combination; values: vector containing detailed info about the involved chains
   * @throws RulesException exception thrown when there is something wrong with the chain name
   */
  private Hashtable<String,Hashtable<String,Vector<FattyAcidVO>>> payAttentionToActuallyPresentFaAndLcbChains(int chainsTotal, int faChains, int lcbChains,
      int acylChains, int alkylChains, int alkenylChains, Vector<String> potentialCombinations, Vector<int[]> ohCombis,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableFAChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableLCBChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkylChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkenylChains,
      Hashtable<String,FattyAcidVO> uniqueChains) throws RulesException{
      Hashtable<String,Hashtable<String,Vector<FattyAcidVO>>> potentialChainCombis = new Hashtable<String,Hashtable<String,Vector<FattyAcidVO>>>();
    //the key for these hashes is the OH number; then follows a vector of available chains
    Hashtable<Integer,Hashtable<String,String>> acylChainsInLib = null;
    Hashtable<Integer,Hashtable<String,String>> alkylChainsInLib = null;
    Hashtable<Integer,Hashtable<String,String>> alkenylChainsInLib = null;
    Hashtable<Integer,Hashtable<String,String>> lcbChainsInLib = null;

    for (int[] ohCombi : ohCombis){
      String ohCombiId = this.getOhCombiId(ohCombi);
      //System.out.println("!!!!!!!!!!!!!! "+ohCombiId);
      Vector<String> ohDistributions = calculateAllPossibleOhDistributionsAmongChainTypes(acylChains,alkylChains,alkenylChains,lcbChains,
          allowedFaHydroxylationsCombinations_.get(ohCombiId),allowedLcbHydroxylationsCombinations_.get(ohCombiId));
//      for (String combi : ohDistributions) {
//        System.out.println("The finally allowed combis: "+combi);
//      }
      Hashtable<String,Vector<FattyAcidVO>> combis = new Hashtable<String,Vector<FattyAcidVO>>();
      potentialChainCombis.put(ohCombiId, combis);
      Set<Integer> relevantFaOHs = null;
      Set<Integer> relevantLcbOHs = null;
      if (faChains>0) relevantFaOHs = getRelevantOhNumbers(allowedFaHydroxylationsCombinations_.get(ohCombiId));
      if (lcbChains>0) relevantLcbOHs = getRelevantOhNumbers(allowedLcbHydroxylationsCombinations_.get(ohCombiId));

      //first check which chains of a combination are available for which FA/LCB      
      for (String faCombi : potentialCombinations) {
        Set<String> uniqueFAs = new HashSet<String>();
        String[] fas = faCombi.split(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
        if (acylChains>0) acylChainsInLib = new Hashtable<Integer,Hashtable<String,String>>();
        if (alkylChains>0) alkylChainsInLib = new Hashtable<Integer,Hashtable<String,String>>();
        if (alkenylChains>0) alkenylChainsInLib = new Hashtable<Integer,Hashtable<String,String>>();
        if (lcbChains>0) lcbChainsInLib = new Hashtable<Integer,Hashtable<String,String>>();
        for (String fa : fas) uniqueFAs.add(fa);
        for (String fa : uniqueFAs) {
          if (acylChainsInLib!=null) checkWhetherFaIsInLib(fa,relevantFaOHs,acylChainsInLib,availableFAChains);
          if (alkylChainsInLib!=null) checkWhetherFaIsInLib(fa,relevantFaOHs,alkylChainsInLib,availableAlkylChains);
          if (alkenylChainsInLib!=null) checkWhetherFaIsInLib(fa,relevantFaOHs,alkenylChainsInLib,availableAlkenylChains);
          if (lcbChainsInLib!=null) checkWhetherFaIsInLib(fa,relevantLcbOHs,lcbChainsInLib,availableLCBChains);          
        }
        for (String ohDistri : ohDistributions) {
          permuteChainCombinationsOnOhPositionsAndCheckForValidity(faCombi,ohDistri,combis,acylChainsInLib,alkylChainsInLib,
              alkenylChainsInLib,lcbChainsInLib,availableFAChains,availableAlkylChains,availableAlkenylChains,
              availableLCBChains,uniqueChains);
        }
      }
    }
    return potentialChainCombis;
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
    if (currentLevel==0)
      return combinations;
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
   * @param ohNumber the degree of hydroxylation to check for fragments
   * @return the FragmentVOs for the head rules - the first key: are these fragments mandatory
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Hashtable<Boolean,Vector<FragmentVO>> getHeadFragments(int ohNumber) throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Hashtable<Boolean,Vector<FragmentVO>> allHeadFragments = new Hashtable<Boolean,Vector<FragmentVO>>();
    Vector<FragmentVO> mandatoryFragments = new Vector<FragmentVO>();
    Vector<FragmentVO> addFragments = new Vector<FragmentVO>();
    Hashtable<String,FragmentRuleVO> headRules = RulesContainer.getHeadFragmentRules(ruleName_,rulesDir_);
    short oh = (short)ohNumber;
    for (FragmentRuleVO ruleVO : headRules.values()){
      if (!ruleVO.hydroxylationValid(oh))
        continue;
      Vector<Object> formulaAndMass = ruleVO.getFormulaAndMz(analyteFormula_, precursorMz_*precursorCharge_, null, 0,ruleVO.getCharge());
      FragmentVO fragVO = new FragmentVO(ruleVO.getName(),(Double)formulaAndMass.get(1),(String)formulaAndMass.get(0),ruleVO.getCharge(),ruleVO.getMsLevel(),
          ruleVO.isMandatory(oh));
      if (ruleVO.isMandatory(oh)==FragmentRuleVO.MANDATORY_TRUE || ruleVO.isMandatory(oh)==FragmentRuleVO.MANDATORY_QUANT) mandatoryFragments.add(fragVO);
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
   * @return the rules for intensity comparisons of the chain group - key: chain type
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Hashtable<Short,Vector<IntensityRuleVO>> getChainIntensityRulesSameChain() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Hashtable<Short,Vector<IntensityRuleVO>> sameRules = new Hashtable<Short,Vector<IntensityRuleVO>>();    
    Vector<IntensityRuleVO> rules = getChainIntensityRules();
    for (IntensityRuleVO rule : rules){
      if (rule.getChainType()==IntensityRuleVO.DIFF_CHAIN_TYPES)
        continue;
      if (!sameRules.containsKey(rule.getChainType()))
        sameRules.put(rule.getChainType(), new Vector<IntensityRuleVO>());
      sameRules.get(rule.getChainType()).add(rule);
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
    return new Vector<FattyAcidVO>(this.availableChains_.values());
  }
  
  /**
   * calculates the FragmentVOs for a given fatty acid chain - the first key: fragment type; second key are these fragments mandatory
   * @param fa fatty acid chain for which the fragments have to be calculated
   * @return the FragmentVOs for a given fatty acid chain - key: are these fragments mandatory
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Hashtable<Boolean,Vector<FragmentVO>> getChainFragments(FattyAcidVO chain) throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Hashtable<Boolean,Vector<FragmentVO>> chainFragments = new Hashtable<Boolean,Vector<FragmentVO>>();
    Vector<FragmentVO> mandatoryFragments = new Vector<FragmentVO>();
    Vector<FragmentVO> addFragments = new Vector<FragmentVO>();
    Hashtable<String,FragmentRuleVO> chainRules =  RulesContainer.getChainFragmentRules(ruleName_,rulesDir_);
    short oh = (short)chain.getOhNumber();
    for (FragmentRuleVO ruleVO : chainRules.values()){
      if (ruleVO.getChainType()!=chain.getChainType() || !ruleVO.hydroxylationValid(oh))
        continue;
      Vector<Object> formulaAndMass = ruleVO.getFormulaAndMass(analyteFormula_, precursorMz_*precursorCharge_, chain, ruleVO.getCharge());
      FragmentVO fragVO = new FragmentVO(ruleVO.getName(),(Double)formulaAndMass.get(1),(String)formulaAndMass.get(0),ruleVO.getCharge(),ruleVO.getMsLevel(),
          ruleVO.isMandatory(oh));
      if (ruleVO.isMandatory(oh)==FragmentRuleVO.MANDATORY_TRUE || ruleVO.isMandatory(oh)==FragmentRuleVO.MANDATORY_QUANT || ruleVO.isMandatory(oh)==FragmentRuleVO.MANDATORY_CLASS) mandatoryFragments.add(fragVO);
      else addFragments.add(fragVO);
    }
    chainFragments.put(true, mandatoryFragments);
    chainFragments.put(false, addFragments);
    return chainFragments;

  }
  
  /**
   * 
   * @param foundFAs the fatty acids detected by MS evidence
   * @param chainFragments the found fragments for the fatty acid chains
   * @param forbiddenChains chains that were removed by a mandatory OR combination
   * @return only fragment combinations where MS evidence was detected
   * @throws SpectrummillParserException if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws IOException exception if there is something wrong about the file
   * @throws NoRuleException thrown if the rules are not there
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws LipidCombinameEncodingException thrown when the lipid combination cannot be decoded
   */
  public Hashtable<String,Vector<FattyAcidVO>> getChainFragmentCombinationsWithDetectedEvidence(Collection<String> foundFAs, Hashtable<String,Hashtable<String,CgProbe>> chainFragments,
      Set<String> forbiddenChains) throws RulesException, NoRuleException, IOException, SpectrummillParserException, LipidCombinameEncodingException {
    Hashtable<String,Vector<FattyAcidVO>> potentialChainCombinations = new Hashtable<String,Vector<FattyAcidVO>>();
    Hashtable<String,String> fas = new Hashtable<String,String>();
    for (String fa: foundFAs){
      fas.put(fa, fa);
    }
    Hashtable<String,FragmentRuleVO> chainRules =  RulesContainer.getChainFragmentRules(ruleName_,rulesDir_);
    boolean singleChainIdentification = RulesContainer.isSingleChainIdentification(ruleName_,rulesDir_);
    Vector<String> relevantCombinations = new Vector<String>();
    try {relevantCombinations =  getChainCombinationsReleveantForTheseChains(fas,forbiddenChains);
    }catch (LipidCombinameEncodingException e) {throw new RulesException(e);}
    short mandatory;
    Vector<String> mandatoryFrags;
    Hashtable<String,Integer> labelsOfChainsBase = new Hashtable<String,Integer>();
    Hashtable<String,Integer> labelsOfChains;
    String singleLabel;
    for (String label : this.availableSingleLabels_) {
      labelsOfChainsBase.put(label, 0);
    }
    for (String key : relevantCombinations){
      boolean allFAsThere = true;
      boolean oneFAIsThere = false;
      Vector<FattyAcidVO> chains = StaticUtils.decodeLipidNamesFromChainCombi(key);
      for (FattyAcidVO chain : chains){
        mandatoryFrags = new Vector<String>();
        for (FragmentRuleVO ruleVO : chainRules.values()){
          if (ruleVO.getChainType()!=chain.getChainType() || !ruleVO.hydroxylationValid((short)chain.getOhNumber()))
            continue;
          mandatory = ruleVO.isMandatory((short)chain.getOhNumber());
          if (mandatory==FragmentRuleVO.MANDATORY_TRUE || mandatory==FragmentRuleVO.MANDATORY_QUANT || mandatory==FragmentRuleVO.MANDATORY_CLASS)
            mandatoryFrags.add(ruleVO.getName());
        }
        if (!fas.containsKey(chain.getChainId())){
          allFAsThere = false;
        } else {
          boolean chainValid = true;
          for (String fragName : mandatoryFrags) {
            if (!chainFragments.containsKey(chain.getChainId()) || !chainFragments.get(chain.getChainId()).containsKey(fragName))
              chainValid = false;
          }
          if (chainValid)
            oneFAIsThere = true;
          else
            allFAsThere = false;
        }
      }
      if (allFAsThere || (oneFAIsThere && singleChainIdentification)) {
        boolean labelingIsOK = true;
        if (LipidomicsConstants.checkChainLabelCombination()) {
          labelsOfChains = new Hashtable<String,Integer>(labelsOfChainsBase);
          for (FattyAcidVO chain : chains){
            if (chain.getPrefix()!=null && chain.getPrefix().length()>0) {
              singleLabel = singleLabelLookup_.get(chain.getPrefix());
              labelsOfChains.put(singleLabel, labelsOfChains.get(singleLabel)+this.availableLabels_.get(chain.getPrefix()));
            }
          }
          for (String label : this.allowedLabelsInChains_.keySet()) {
            if (this.allowedLabelsInChains_.get(label) != labelsOfChains.get(label)) {
              labelingIsOK = false;
              break;
            }
          }
        }
        
        if (labelingIsOK)
          potentialChainCombinations.put(key, chains);
      }
    }
    return potentialChainCombinations;
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
    if (!cAtomsMatcher.matches()){
    	//to match also oxidation modifications e.g.: 34:3;O1
    	cAtomsPattern =  Pattern.compile(rule+"\\"+LipidomicsConstants.CHAIN_MOD_SEPARATOR + ".*");
    	cAtomsMatcher = cAtomsPattern.matcher(analyte);
    	if (!cAtomsMatcher.matches())  throw new RulesException("The analyte "+analyte+" does not match the "+inputKey+" pattern \""+rule+"\" of the class "+lClass+"!");
    }
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
   * @param chainLib the name of the FA chain library
   * @param lcbLib the name of the LCB chain library
   * @param cAtoms the number of possible carbon atoms in the chains
   * @param dbs the number of possible double bonds in the chain
   * @param alkylPresent consists this class of any alkylated chains
   * @param alkenylPresent consists this class of any alkenylated chains
   * @throws RulesException RulesException thrown if there is something wrong
   * @throws IOException exception if there is something wrong about the file
   * @throws NoRuleException thrown if the rules are not there
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws ChemicalFormulaException 
   */
  private void checkChainsForPlausibility(String chainLib, String lcbLib, int cAtoms, int dbs,
      boolean alkylPresent, boolean alkenylPresent) throws RulesException, NoRuleException, IOException, HydroxylationEncodingException, ChemicalFormulaException{
    availableLabels_ = new Hashtable<String,Integer>();
    singleLabelLookup_ = new Hashtable<String,String>();
    availableSingleLabels_ = new HashSet<String>();
    allowedLabelsInChains_ = new Hashtable<String,Integer>();
    if (chainLib==null) {
      availableFAChainsBeforeCombiCheck_ = null;
    } else {
      availableFAChainsBeforeCombiCheck_ = new Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>();
      Hashtable<String, Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>>> faHydroxies = FattyAcidsContainer.getAllFattyAcidChains(chainLib);
      addAvailableLabels(FattyAcidsContainer.getAvailableLabels(chainLib));
      for (int ohNumber : getPossibleFaHydroxylations()) {
        String encoded = HydroxyEncoding.HYDROXYLATION_ZERO;
        if (ohNumber!=0) encoded = Settings.getFaHydroxyEncoding().getEncodedPrefix((short)ohNumber);
        //System.out.println("FA: "+ohNumber);
        checkPlausibilityAndAddToHash(availableFAChainsBeforeCombiCheck_,faHydroxies.get(encoded),ohNumber,cAtoms,dbs);
      }
    }
    if (availableFAChainsBeforeCombiCheck_!=null && alkylPresent) {
      availableAlkylChainsBeforeCombiCheck_ = createAlkylAlkenylatedHash(LipidomicsConstants.CHAIN_TYPE_FA_ALKYL,availableFAChainsBeforeCombiCheck_);
    }
    if (availableFAChainsBeforeCombiCheck_!=null && alkenylPresent) {
      availableAlkenylChainsBeforeCombiCheck_ = createAlkylAlkenylatedHash(LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL,availableFAChainsBeforeCombiCheck_);
    }

    if (lcbLib==null) {
      availableLCBChainsBeforeCombiCheck_ = null;
    } else {
      availableLCBChainsBeforeCombiCheck_ = new Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>();
      Hashtable<String, Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>>> lcbHydroxies = FattyAcidsContainer.getAllLCBs(lcbLib);
      addAvailableLabels(FattyAcidsContainer.getAvailableLabels(lcbLib));
      for (int ohNumber : getPossibleLcbHydroxylations()) {
        String encoded = Settings.getLcbHydroxyEncoding().getEncodedPrefix((short)ohNumber);
        //System.out.println("LCB: "+ohNumber);
        checkPlausibilityAndAddToHash(availableLCBChainsBeforeCombiCheck_,lcbHydroxies.get(encoded),ohNumber,cAtoms,dbs);
      }
    }
    
    if (LipidomicsConstants.checkChainLabelCombination()) 
    {
      for (String label : availableSingleLabels_) {
        allowedLabelsInChains_.put(label, 0);        
      }
      Vector<String> lengthSortedLabels = new Vector<String>();
      for (String label : availableSingleLabels_) {
        int posToAdd = -1;
        for (int i=0; i!=lengthSortedLabels.size(); i++) {
          if (label.length()>lengthSortedLabels.get(i).length()) {
            posToAdd = i;
            break;
          }
        }
        if (posToAdd==-1)
          lengthSortedLabels.add(label);
        else
          lengthSortedLabels.add(posToAdd,label);
      }
      String restOfName = new String(this.labelInName_);
      for (String label : lengthSortedLabels) 
      {
        while (restOfName.indexOf(label)!=-1) 
        {
          allowedLabelsInChains_.put(label, allowedLabelsInChains_.get(label)+1);
          restOfName = restOfName.substring(0,restOfName.indexOf(label))+restOfName.substring(restOfName.indexOf(label)+label.length());
        }
        if (restOfName.length()==0)
          break;
      }
    }
  }
  
  
  /**
   * adds the possible isotopic labels stored in the chain libraries to the hashes availableLabels_ and availableSingleLabels_;
   * @param labels
   */
  private void addAvailableLabels(Set<String> labels) {
    StaticUtils.extractIsoLabelInformation(labels, availableSingleLabels_, singleLabelLookup_, availableLabels_);
  }
  
  
  /**
   * checks whether this FA/LCB chain is plausible for this particular species
   * @param hash the hash table were the possible chains shall be added
   * @param toBeChecked the chain library to be checked
   * @param ohNumber the number of oh fragments
   * @param maxCAtoms the maximum number of allowed C atoms
   * @param maxDbs the maximum number of allowed double bonds
   * @throws RulesException RulesException thrown if there is something wrong
   */
  private void checkPlausibilityAndAddToHash(Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> hash,
		  Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> toBeChecked, int ohNumber, int maxCAtoms, int maxDbs) throws RulesException {
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>> chains = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>();
    try{
      Hashtable<String,Integer> formulaAmounts = StaticUtils.categorizeFormula(this.analyteFormulaWODeducts_);
      for (Integer cAtoms : toBeChecked.keySet()){
        if (cAtoms>maxCAtoms)
          continue;
        Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>> sameCAtoms = new Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>();
        for (Integer dbs : toBeChecked.get(cAtoms).keySet()){
          if (dbs>maxDbs)
            continue;
          Hashtable<String,Hashtable<String,FattyAcidVO>> sameDbs = new Hashtable<String,Hashtable<String,FattyAcidVO>>();
          for (String prefix : toBeChecked.get(cAtoms).get(dbs).keySet())
          {
	        	Hashtable<String,FattyAcidVO> sameOxs = new Hashtable<String,FattyAcidVO>();
	        	for(String oxState : toBeChecked.get(cAtoms).get(dbs).get(prefix).keySet())
	        	{
	        	  FattyAcidVO fa = toBeChecked.get(cAtoms).get(dbs).get(prefix).get(oxState);
	        	  Hashtable<String,Integer> faElements = StaticUtils.categorizeFormula(fa.getFormula());
	        	  boolean isOk = true;
	              for (String element : faElements.keySet()){
	                if (!formulaAmounts.containsKey(element) || formulaAmounts.get(element)<faElements.get(element)){
	                  isOk = false;
	                  break;
	                }
	              }
	              if (isOk)
	            	  sameOxs.put(oxState, fa);
	        	}
	        	if (sameOxs.size()>0) sameDbs.put(prefix, sameOxs);   
          }
          if (sameDbs.size()>0) sameCAtoms.put(dbs, sameDbs);
        }
        if (sameCAtoms.size()>0) chains.put(cAtoms, sameCAtoms);
      }
    } catch (ChemicalFormulaException e) {
      throw new RulesException("The formula "+analyteFormulaWODeducts_+" contains fragments that have not been defined before! Fragments have to be defined in previous columns, before they can be used!");
    }
    hash.put(ohNumber, chains);
  }
  
  
  /**
   * calculates the potential hydroxylation site assignments to FAs and LCBs depending on the available hydroxylations, and the number of respective chains
   * @param fattyChains the number of FA chains
   * @param lcbChains the number of LCB chains
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there*/
  private void calculatePossibleFaLcbHydroxyCombinations(int faChains, int lcbChains) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    int totalPossibleFaHydroxylations = 0;
    possibleOhCombinations_ = new Vector<int[]>();
    allowedFaHydroxylationsCombinations_ = new Hashtable<String,Vector<Vector<Integer>>>();
    allowedLcbHydroxylationsCombinations_ = new Hashtable<String,Vector<Vector<Integer>>>();
    RangeInteger faRange = new RangeInteger(0,0);
    if (faChains>0 && RulesContainer.getFaHydroxyRange(ruleName_,rulesDir_)!=null) {
      faRange = RulesContainer.getFaHydroxyRange(ruleName_,rulesDir_);
      totalPossibleFaHydroxylations = faChains*faRange.getStop();
    }
    List<Integer> availableFaOHs = new ArrayList<Integer>();
    for (int i=faRange.getStart(); i!=(faRange.getStop()+1); i++) availableFaOHs.add(i);
    Hashtable<Integer,List<Integer>> availableFaHash = new Hashtable<Integer,List<Integer>>();
    for (int i=0; i!=faChains; i++) availableFaHash.put((i+1), availableFaOHs);

    int totalPossibleLcbHydroxylations = 0;
    RangeInteger lcbRange = new RangeInteger(0,0);
    if (lcbChains>0 && RulesContainer.getLcbHydroxyRange(ruleName_,rulesDir_)!=null) {
      lcbRange = RulesContainer.getLcbHydroxyRange(ruleName_,rulesDir_);
      totalPossibleLcbHydroxylations = lcbChains*lcbRange.getStop();
    }
    List<Integer> availableLcbOHs = new ArrayList<Integer>();
    for (int i=lcbRange.getStart(); i!=(lcbRange.getStop()+1); i++) availableLcbOHs.add(i);
    Hashtable<Integer,List<Integer>> availableLcbHash = new Hashtable<Integer,List<Integer>>();
    for (int i=0; i!=lcbChains; i++) availableLcbHash.put((i+1), availableLcbOHs);
    
    
    for (int i=(faRange.getStart()*faChains); i<(totalPossibleFaHydroxylations+1); i++) {
      for (int j=(lcbRange.getStart()*lcbChains); j<(totalPossibleLcbHydroxylations+1); j++) {
        if (this.ohNumber_==(i+j)) {
          this.possibleOhCombinations_.add(new int[] {i,j});
          String combiId = getOhCombiId(this.possibleOhCombinations_.lastElement());
          //now check which combinations are possible to fulfil the oh total for the FA
          allowedFaHydroxylationsCombinations_.put(combiId, getCombinations(i, faChains, availableFaHash, new Vector<Integer>(), 0,false));
          //now check which combinations are possible to fulfil the oh total for the LCB
          allowedLcbHydroxylationsCombinations_.put(combiId,getCombinations(j, lcbChains, availableLcbHash, new Vector<Integer>(), 0,false));
        }
      }
    }
//    for (int[] combi : this.possibleOhCombinations_) {
//      System.out.println(combi[0]+";"+combi[1]);
//    }
    if (this.possibleOhCombinations_.size()==0)
      throw new RulesException("The amount of "+this.ohNumber_+" hydroxylation sites, cannot be fulfilled by "+lcbChains+" LCB chains with "
          + (lcbRange.getStart()==lcbRange.getStop() ? String.valueOf(lcbRange.getStart()) : String.valueOf(lcbRange.getStart())+"-"+String.valueOf(lcbRange.getStop()))
          + " hydroxylation sites, and "+faChains+" FA chains with "
          + (faRange.getStart()==faRange.getStop() ? String.valueOf(faRange.getStart()) : String.valueOf(faRange.getStart())+"-"+String.valueOf(faRange.getStop()))
          + " hydroxylation sites!");
  }
  

  /**
   * generates the OH FA/LCB combi-id based on an int[] where int[0] = number of FA OHs, and int[1] = number of LCB OHs
   * @param combi
   * @return
   */
  private String getOhCombiId(int[] combi) {
    return String.valueOf(combi[0])+OH_COMBI_SEPARATOR+String.valueOf(combi[1]);
  }
  
  /**
   * @return the possible hydroxylation numbers of the fatty acids
   */
  private Vector<Integer> getPossibleFaHydroxylations(){
    return getPossibleHydroxylations(0);
  }

  /**
   * @return the possible hydroxylation numbers of the fatty acids
   */
  private Vector<Integer> getPossibleLcbHydroxylations(){
    return getPossibleHydroxylations(1);
  }

  
  /**
   * the possible hydroxylation sites
   * @param position should position 0 or 1 be used for the list
   * @return list of possible hydroxylation sites
   */
  private Vector<Integer> getPossibleHydroxylations(int position){
    Vector<Integer> ohIndices = new Vector<Integer>();
    Hashtable<Integer,Integer> ohs = new Hashtable<Integer,Integer>();
    int highestOh = 0;
    for (int[] ohCombi : possibleOhCombinations_) {
      String combiId = this.getOhCombiId(ohCombi);
      Vector<Vector<Integer>> allowedCombis = null;
      if (position==0) allowedCombis =  allowedFaHydroxylationsCombinations_.get(combiId);
      else if (position==1) allowedCombis =  allowedLcbHydroxylationsCombinations_.get(combiId);
      if (allowedCombis==null)continue;
      for (Vector<Integer> combi : allowedCombis) {
        for (Integer ohNumber : combi) {
          ohs.put(ohNumber, ohNumber);
          if (ohNumber>highestOh) highestOh = ohNumber;
        }
      }
    }
    for (int i=0; i<(highestOh+1); i++) {
      if (ohs.containsKey(i))
        ohIndices.add(i);
    }
    return ohIndices;
  }
  
  /**
   * adds the stored carbon atoms/dbs combinations to a unique hash - can be performed for adding several hashes to avoid duplicates, such as FA and LCB
   * @param hash the unique carbon atoms/dbs combinations hash
   * @param toBeAdded the stored hashes that have to be checked for adding
   */
  private void addUniqueChainDbsCombis(Hashtable<Integer,Hashtable<Integer,Hashtable<String,String>>> hash,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> toBeAdded) {
    for (Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> carbonHash : toBeAdded.values()) {
      for (Integer cAtoms : carbonHash.keySet()) {
        Hashtable<Integer,Hashtable<String,String>> dbHash;
        if (hash.containsKey(cAtoms)) {
          dbHash = hash.get(cAtoms);
        }else {
          dbHash = new Hashtable<Integer,Hashtable<String,String>>();
          hash.put(cAtoms, dbHash);
        }
        for (Integer dbs : carbonHash.get(cAtoms).keySet())
        {
          Hashtable<String,String> oxHash;
          if(dbHash.containsKey(dbs))
          {
        	  oxHash = hash.get(cAtoms).get(dbs);
          }
          else {
        	  oxHash = new Hashtable<String,String>();
        	  dbHash.put(dbs, oxHash);
          }
          
          for(String prefix : carbonHash.get(cAtoms).get(dbs).keySet()) 
          {
        	  for(String oxs : carbonHash.get(cAtoms).get(dbs).get(prefix).keySet())
        	  {
        		  oxHash.put(oxs,oxs);
        	  }
          }
        }
      }
    }
  }
  
  /**
   * creates the hashes when alkylated or alkenylatd chains are present
   * @param chainType the chain type to be created (alkylated or alkenylated)
   * @param faHash the acylated chain hash
   * @return the alkylated or alkenylated chain hash
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   */
  private Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>  createAlkylAlkenylatedHash(short chainType,
		  Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>  faHash) throws ChemicalFormulaException{
	Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>  newHash = new Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>();
    for (Integer oh : faHash.keySet()) {
      newHash.put(oh, new Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>());
      for (Integer c : faHash.get(oh).keySet()) {
        newHash.get(oh).put(c, new Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>());
        for (Integer dbs : faHash.get(oh).get(c).keySet()) 
        {
          newHash.get(oh).get(c).put(dbs, new Hashtable<String,Hashtable<String,FattyAcidVO>>());
          for (String prefix : faHash.get(oh).get(c).get(dbs).keySet()) 
          {
	        	newHash.get(oh).get(c).get(dbs).put(prefix, new Hashtable<String,FattyAcidVO>());  
	        	for(String oxState : faHash.get(oh).get(c).get(dbs).get(prefix).keySet())
	        	{
	        		FattyAcidVO faVO = faHash.get(oh).get(c).get(dbs).get(prefix).get(oxState);
              Hashtable<String,Integer> faElements = StaticUtils.categorizeFormula(faVO.getFormula());
              double mass = faVO.getMass();
              if (!faElements.containsKey("O"))
                continue;
              if (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ALKYL){
                mass += (2d*elements_.getElementDetails("H").getMonoMass()-elements_.getElementDetails("O").getMonoMass());
                faElements.put("H",(faElements.get("H")+2));
                faElements.put("O",(faElements.get("O")-1));
              } else if (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL){
                mass += (-1d*elements_.getElementDetails("O").getMonoMass());
                faElements.put("O",(faElements.get("O")-1));
              }
              if (faElements.get("O")<0 || mass<0d)
                continue;
              newHash.get(oh).get(c).get(dbs).get(prefix).put(oxState, new FattyAcidVO(chainType, faVO.getPrefix(), faVO.getcAtoms(),
                  faVO.getDoubleBonds(), oh, mass, StaticUtils.getFormulaInHillNotation(faElements, true),faVO.getOxState()));
	        	}
          }
        }
      }
    }
    return newHash;
  }

  
  /** returns all OH-numbers that occur in these combinations
   * @param hydroCombis the various OH combinations for an FA or LCB
   * @return a unique set of the ones that occur
   */
  private Set<Integer> getRelevantOhNumbers(Vector<Vector<Integer>> hydroCombis){
    Set<Integer> relevant = new HashSet<Integer>();
    for (Vector<Integer> ohFaCombi : hydroCombis) {
      for (Integer ohNumber : ohFaCombi) relevant.add(ohNumber);
    }
    return relevant;
  }
  
   
 /**
   * checks whether a chain name is present in the chain lib (respecting the OH number), and adds it to the hashToAdd
   * @param chain the chain name 
   * @param ohNumbers the OH numbers to check
   * @param hashToAdd the hash where the chain should be added
   * @param chainLib a library containing all available chains
   * @throws RulesException exception thrown when there is something wrong with the chain name
   */
  private void checkWhetherFaIsInLib(String chain, Set<Integer> ohNumbers, Hashtable<Integer,Hashtable<String,String>> hashToAdd,
		  Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> chainLib) throws RulesException 
  {
    String[] cAndDbAndOx = null;
    try {
      cAndDbAndOx = StaticUtils.parseCAndDbsFromChainId(chain);
    } catch (Exception e) { throw new RulesException(e.getMessage());}
    for (Integer oh : ohNumbers) 
    {
  		if (!hashToAdd.containsKey(oh)) hashToAdd.put(oh, new Hashtable<String,String>());
      if (chainLib.containsKey(oh) 
      		&& chainLib.get(oh).containsKey(Integer.parseInt(cAndDbAndOx[0])) 
      		&& chainLib.get(oh).get(Integer.parseInt(cAndDbAndOx[0])).containsKey(Integer.parseInt(cAndDbAndOx[1])) 
      		&& chainLib.get(oh).get(Integer.parseInt(cAndDbAndOx[0])).get(Integer.parseInt(cAndDbAndOx[1])).values().iterator().next().containsKey(cAndDbAndOx[2]))
      {
      	hashToAdd.get(oh).put(chain,chain);
      }
  	}   
  }
  
  /**
   * This method determines all permuted combinations how the OH groups can be spread over the chain types:
   * CHAIN_TYPE_FA_ACYL
   * CHAIN_TYPE_FA_ALKYL
   * CHAIN_TYPE_FA_ALKENYL
   * CHAIN_TYPE_LCB
   * The combinations of the chain types are separated by _
   * Each chain type that contributes to the molecule is encoded as follows:
   * $CHAIN_TYPE$:$NUMBER_OHs$:$NUMBER_CHAINS$
   * @param acylChains number of acylated chains
   * @param alkylChains number of alkylated chains
   * @param alkenylChains number of alkenylated chains
   * @param lcbChains number of LCB chains
   * @param faHydroxylationOptions which combinations of hydroxylations are allowed for the individual FA chains
   * @param lcbHydroxylationOptions which combinations of hydroxylations are allowed for the individual LCB chains
   * @return all potential hydroxylation combinations of OH groups among the chain types
   */
  private Vector<String> calculateAllPossibleOhDistributionsAmongChainTypes(int acylChains, int alkylChains, int alkenylChains, int lcbChains,
      Vector<Vector<Integer>> faHydroxylationOptions, Vector<Vector<Integer>> lcbHydroxylationOptions){
    Vector<String> encodedOptions = new Vector<String>();
    Vector<String> encodedFaMoiety = new Vector<String>();
    Vector<String> encodedLcbMoiety = new Vector<String>();
    if (faHydroxylationOptions!=null) {
      for (Vector<Integer> hydroxiesOnChain : faHydroxylationOptions) {
        //key: hydroxylation number; value: number of chains of this type
        Hashtable<Integer,Integer> numberOfChainsWithSameOh = new Hashtable<Integer,Integer>();
        for (Integer oh : hydroxiesOnChain) {
          if (numberOfChainsWithSameOh.containsKey(oh))
            numberOfChainsWithSameOh.put(oh, numberOfChainsWithSameOh.get(oh)+1);
          else
            numberOfChainsWithSameOh.put(oh, 1);
        }
//        System.out.println("chains "+acylChains+"_"+numberOfChainsWithSameOh);
        Vector<String> allPossibleAcylPermutations = getPossibleFaOHPermutations(LipidomicsConstants.CHAIN_TYPE_FA_ACYL, acylChains,numberOfChainsWithSameOh);
//        for (String combi : allPossibleAcylPermutations) {
//          System.out.println("combi Acyl: "+combi);
//        }
        Vector<String> allPossibleAlkylPermutations = getPossibleFaOHPermutations(LipidomicsConstants.CHAIN_TYPE_FA_ALKYL, alkylChains,numberOfChainsWithSameOh);
//        for (String combi : allPossibleAlkylPermutations) {
//          System.out.println("combi Alkyl: "+combi);
//        }
        Vector<String> allPossibleAlkenylPermutations = getPossibleFaOHPermutations(LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL, alkenylChains,numberOfChainsWithSameOh);
//        for (String combi : allPossibleAlkenylPermutations) {
//          System.out.println("combi Alkenyl: "+combi);
//        }
        //now calculate all possible acyl/alkyl/alkenyl permutations
        Vector<String> possibleCombis = getAcylAlkylAlkenylPermutations(allPossibleAcylPermutations,allPossibleAlkylPermutations,allPossibleAlkenylPermutations,numberOfChainsWithSameOh);
//        for (String combi : possibleCombis) {
//          System.out.println("These are OK: "+combi);
//        }
        encodedFaMoiety.addAll(possibleCombis);
      }
    }
    if (lcbHydroxylationOptions!=null) {
      int lowest;
      int highest;
      int i;
      for (Vector<Integer> hydroxiesOnChain : lcbHydroxylationOptions) {
        //System.out.println("!!!!! "+hydroxiesOnChain);
        Hashtable<Integer,Integer> numberOfChainsWithSameOh = new Hashtable<Integer,Integer>();
        lowest = Integer.MAX_VALUE;
        highest = 0;
        for (Integer oh : hydroxiesOnChain) {
          if (oh<lowest) lowest = oh;
          if (oh>highest) highest = oh;
          if (numberOfChainsWithSameOh.containsKey(oh))
            numberOfChainsWithSameOh.put(oh, numberOfChainsWithSameOh.get(oh)+1);
          else
            numberOfChainsWithSameOh.put(oh, 1);
        }
        StringBuilder encoded = new StringBuilder();
        for (i=lowest; i<(highest+1); i++) {
          if (!numberOfChainsWithSameOh.containsKey(i))
            continue;
          if (encoded.length()>0)
            encoded.append(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
          encoded.append(createTypeOHFrequencyString(LipidomicsConstants.CHAIN_TYPE_LCB,i,numberOfChainsWithSameOh.get(i)));          
        }
        encodedLcbMoiety.add(encoded.toString());
      }
    }

    //now mix the FA and LCB encoded moieties together
    //are there any FAs
    if (encodedFaMoiety.size()>0) {
      //are there FAs and LCBs
      if (encodedLcbMoiety.size()>0) {
        for (String faCombi : encodedFaMoiety) {
          for (String lcbCombi : encodedLcbMoiety) {
            encodedOptions.add(faCombi+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS+lcbCombi);
          }
        }
      }else
        encodedOptions.addAll(encodedFaMoiety);
    } else {
      //are there only LCBs
      if (encodedLcbMoiety.size()>0)
        encodedOptions.addAll(encodedLcbMoiety);
    }
    return encodedOptions;
  }
  
  /**
   * returns all possible permutations of the OH position distributed on a similar type of chain - the total number of OHs can be below the allowed ones
   * here, only unique (position-independent) permutations are returned
   * @param chainType the chain type LipidomicsConstants.CHAIN_TYPE_FA_ACYL/CHAIN_TYPE_FA_ALKYL/CHAIN_TYPE_FA_ALKENYL/CHAIN_TYPE_LCB
   * @param chains how many chains of this type are present
   * @param sameOhFrequencies the maximally allowed number of hydroxylation sites for a given hydroxylation number; key: hydroxylation number; value number of available chains
   * @return all possible permutations of the OH position distributed on a similar type of chain
   */
  private Vector<String> getPossibleFaOHPermutations(short chainType, int chains, Hashtable<Integer,Integer> sameOhFrequencies){
    Vector<String> toReturn = new Vector<String>();
    if (chains<1)
      return toReturn;
    Vector<String> permut = getOhPermutations(chains, sameOhFrequencies);
    int lowest;
    int highest;
    int i;
    Integer oh;
    String encodedStr;
    if (permut.size()>0) {
      Hashtable<String,String> unique = new Hashtable<String,String>(); 
      for (String combi : permut) {
      //the position is irrelevant -> split the String again, and allow only the unique occurrences
        lowest = Integer.MAX_VALUE;
        highest = 0;
        Hashtable<Integer,Integer> ohFrequencies = new Hashtable<Integer,Integer>();
        String[] split = combi.split(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
        for (String ohString : split) {
          oh = new Integer(ohString);
          if (oh<lowest) lowest = oh;
          if (oh>highest) highest = oh;
          if (ohFrequencies.containsKey(oh)) ohFrequencies.put(oh, ohFrequencies.get(oh)+1);
          else ohFrequencies.put(oh, 1);
        }
        StringBuilder encoded = new StringBuilder();
        for (i=lowest; i<(highest+1); i++) {
          if (!ohFrequencies.containsKey(i))
            continue;
          if (encoded.length()>0)
            encoded.append(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
          encoded.append(createTypeOHFrequencyString(chainType,i,ohFrequencies.get(i)));          
        }
        encodedStr = encoded.toString(); 
        unique.put(encodedStr, encodedStr);
      }
      for (String  encoded : unique.keySet()) {
        toReturn.add(encoded);
      }
    }    
    return toReturn;
  }
  
  
  /**
   * returns all possible permutations of the OH position distributed on a similar type of chain - the total number of OHs can be below the allowed ones
   * here, all (position-dependent) permutations are returned
   * @param chains how many chains of this type are present
   * @param sameOhFrequencies the maximally allowed number of hydroxylation sites for a given hydroxylation number; key: hydroxylation number; value number of available chains
   * @return all possible permutations of the OH position distributed on a similar type of chain
   */
  private Vector<String> getOhPermutations(int chains, Hashtable<Integer,Integer> sameOhFrequencies){
    Vector<String> toReturn = new Vector<String>();
    int amount;
    for (Integer oh : sameOhFrequencies.keySet()) {
      if (sameOhFrequencies.containsKey(oh)) {
        String combi = String.valueOf(oh);
        //a recursion is required
        if (chains>1) {
          Hashtable<Integer,Integer> newOhFrequencies = new Hashtable<Integer,Integer>(sameOhFrequencies);
          amount = newOhFrequencies.get(oh);
          amount--;
          if (amount<=0) newOhFrequencies.remove(oh);
          else newOhFrequencies.put(oh, amount);
          for (String otherPart : getOhPermutations(chains-1, newOhFrequencies))
            toReturn.add(combi+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS+otherPart);
        } else {
          toReturn.add(combi);
        }
      }
    }
    return toReturn;
  }

  
  /**
   * method for encoding an OH frequency identifier
   * @param chainType the chain type (LipidomicsConstants.CHAIN_TYPE_FA_ACYL/CHAIN_TYPE_FA_ALKYL/CHAIN_TYPE_FA_ALKENYL/CHAIN_TYPE_LCB)
   * @param oh the number of OH positions
   * @param frequency the number of chains with this number of OH positions
   * @return the encoded OH frequency identifier
   */
  private String createTypeOHFrequencyString(short chainType, int oh, int frequency) {
    return chainType+OH_FREQUENCY_ID_SEPRATOR+oh+OH_FREQUENCY_ID_SEPRATOR+frequency;
  }
  
  /**
   * method for combining the permutations from acyl, alkyl and alkenyl chains - here, the final OH calculated must be the same as in numberOfChainsWithSameOh
   * @param acyls the possible acyl permutations
   * @param alkyls the possible alkyl permutations
   * @param alkenyls the possible alkenyl permutations
   * @param numberOfChainsWithSameOh the lookup for the allowe OH chains and their amount
   * @return a combined permutation string from acyl, alkyl and alkenyl chains
   */
  private Vector<String> getAcylAlkylAlkenylPermutations(Vector<String> acyls, Vector<String> alkyls, Vector<String> alkenyls,
      Hashtable<Integer,Integer> numberOfChainsWithSameOh){
    Vector<String> faMoieties = new Vector<String>();
    String encodedTotal;
    // are there any acyl chains
    if (acyls!=null && acyls.size()>0) {
      for (String acylCombi : acyls) {
        // are there acyl and alkyl chains
        if (alkyls!=null && alkyls.size()>0) {
          for (String alkylCombi : alkyls) {
            // are there acyl-, alkyl- and alkenyl-chains
            if (alkenyls!=null && alkenyls.size()>0) {
              for (String alkenylCombi : alkenyls) {
                encodedTotal = acylCombi+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS+alkylCombi+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS+alkenylCombi;
                if (areOhNumbersFulfilled(encodedTotal, numberOfChainsWithSameOh)) faMoieties.add(encodedTotal);
              }
            }else{
              encodedTotal = acylCombi+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS+alkylCombi;
              if (areOhNumbersFulfilled(encodedTotal, numberOfChainsWithSameOh)) faMoieties.add(encodedTotal);
            }
          }
        }else{
          // are there acyl and alkenyl chains
          if (alkenyls!=null && alkenyls.size()>0) {
            for (String alkenylCombi : alkenyls) {
              encodedTotal = acylCombi+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS+alkenylCombi;
              if (areOhNumbersFulfilled(encodedTotal, numberOfChainsWithSameOh)) faMoieties.add(encodedTotal);
            }
          // there are only acyl chains 
          }else{
            encodedTotal = acylCombi;
            if (areOhNumbersFulfilled(encodedTotal, numberOfChainsWithSameOh)) faMoieties.add(encodedTotal);
          }
        }
      }
    }else{
      // are there alkyl chains
      if (alkyls!=null && alkyls.size()>0) {
        for (String alkylCombi : alkyls) {
          // are there alkyl and alkenyl chains
          if (alkenyls!=null && alkenyls.size()>0) {
            for (String alkenylCombi : alkenyls) {
              encodedTotal = alkylCombi+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS+alkenylCombi;
              if (areOhNumbersFulfilled(encodedTotal, numberOfChainsWithSameOh)) faMoieties.add(encodedTotal);
            }
          // there are only alkyl chains 
          }else{
            encodedTotal = alkylCombi;
            if (areOhNumbersFulfilled(encodedTotal, numberOfChainsWithSameOh)) faMoieties.add(encodedTotal);
          }          
        }
      }else{
        if (alkenyls!=null && alkenyls.size()>0) {
          for (String alkenylCombi : alkenyls) {
            encodedTotal = alkenylCombi;
            if (areOhNumbersFulfilled(encodedTotal, numberOfChainsWithSameOh)) faMoieties.add(encodedTotal);
          }
        // no combi is detectable 
        }else{
        }        
      }
    }
    return faMoieties;
  }
  
  /**
   * check if the allowed OH chains and their amount are the same as in an encoded String
   * @param encodedTotal the encoded combination string
   * @param numberOfChainsWithSameOh the lookup for the allowed OH chains and their amount
   * @return if they are the same and consequently, the OH distribution among the chains is valid
   */
  private boolean areOhNumbersFulfilled(String encodedTotal, Hashtable<Integer,Integer> numberOfChainsWithSameOh) {
    Hashtable<Integer,Integer> refTable = new Hashtable<Integer,Integer>();
    String[] ohNames = encodedTotal.split(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    for (String oneOhName : ohNames) {
      int[] decoded = decodeTypeOHFrequencyString(oneOhName);
      if (refTable.containsKey(decoded[1])) refTable.put(decoded[1], refTable.get(decoded[1])+decoded[2]);
      else refTable.put(decoded[1], decoded[2]);
    }
    boolean ok = true;
    if (refTable.size()==numberOfChainsWithSameOh.size()){
      for (Integer oh : numberOfChainsWithSameOh.keySet()) {
        if (!refTable.containsKey(oh) || numberOfChainsWithSameOh.get(oh)!=refTable.get(oh)) {
          ok = false;
          break;
        }
      }
    }else
      ok = false;
    
    return ok;
  }
  
  /**
   * decodes an OH frequency identifier
   * @param encoded the encoded string
   * @return int[0] chain type (LipidomicsConstants.CHAIN_TYPE_FA_ACYL/CHAIN_TYPE_FA_ALKYL/CHAIN_TYPE_FA_ALKENYL/CHAIN_TYPE_LCB), int[1] number of OHs; int[2] number of chains
   */
  private int[] decodeTypeOHFrequencyString(String encoded) {
    int[] decoded = new int[3];
    String[] parts = encoded.split(OH_FREQUENCY_ID_SEPRATOR);
    for (int i=0; i<3; i++) {
      decoded[i] = Integer.parseInt(parts[i]);
    }
    return decoded;
  }

  
  /**
   * this method permutes first the available FA chains among the defined chain type and OH combinations stored in ohDistri
   * then, it checks which combinations are possible based on the available...Chains
   * then, for such a combination, all potentially isotopically labeled "same chains" are permuted
   * finally, a unique identifier is created for storing in the toAdd hash, and the referred FattyAcidVO are added to the referenced Vector
   * @param chainCombi the combination of chains
   * @param ohDistri encoded String specifying how the OHs are distributed among the available chain types
   * @param toAdd the hash were the final results shall be added
   * @param acylChainsInLib hash for available acyl chains; key: OH number; value list of available chains defined by C atoms and dbs
   * @param alkylChainsInLib hash for available alkyl chains; key: OH number; value list of available chains defined by C atoms and dbs
   * @param alkenylChainsInLib hash for available alkenyl chains; key: OH number; value list of available chains defined by C atoms and dbs
   * @param lcbChainsInLib hash for available LCB chains; key: OH number; value list of available chains defined by C atoms and dbs
   * @param availableFAChains the principally (from the chemical formula) possible acyl chains 
   * @param availableAlkylChains the principally (from the chemical formula) possible alkyl chains 
   * @param availableAlkenylChains the principally (from the chemical formula) possible alkenyl chains
   * @param availableLCBChains the principally (from the chemical formula) possible LCB chains
   * @param uniqueChains all the chains that are possible after the combinatorial check; key is the chain id that contains information about chain type and hydroxylation sites
   * @throws RulesException thrown when there is something wrong
   */
  
  private void permuteChainCombinationsOnOhPositionsAndCheckForValidity(String chainCombi, String ohDistri, Hashtable<String,Vector<FattyAcidVO>> toAdd,
      Hashtable<Integer,Hashtable<String,String>> acylChainsInLib,  Hashtable<Integer,Hashtable<String,String>> alkylChainsInLib,
      Hashtable<Integer,Hashtable<String,String>> alkenylChainsInLib,  Hashtable<Integer,Hashtable<String,String>> lcbChainsInLib,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableFAChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkylChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkenylChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableLCBChains,
      Hashtable<String,FattyAcidVO> uniqueChains) throws RulesException {
    String[] chains = chainCombi.split(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
//    String[] ohIds = ohDistri.split(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
//    for (String ohId : ohIds) {
    boolean[] usedChains = new boolean[chains.length];
    try {
      Vector<String> chainPermutations =  getPermutedChainsSplitOnChainAndOhTypes(ohDistri,chains,usedChains,acylChainsInLib,alkylChainsInLib,
        alkenylChainsInLib,lcbChainsInLib,availableFAChains,availableAlkylChains,availableAlkenylChains,availableLCBChains);
      if (chainPermutations==null || chainPermutations.size()==0)
        return;
      //here comes the final step - permute all available isotopically labeled chains
      for (String permuteId : chainPermutations) {
        addIsotopeLabelsToMutations(permuteId,toAdd,availableFAChains,availableAlkylChains,availableAlkenylChains,availableLCBChains,
            uniqueChains);
      }
    }catch (LipidCombinameEncodingException | ChemicalFormulaException e) {
      throw new RulesException(e);
    }
  }

  //first key is the type, i.e. ohTypeFrequencyId; second key is the chain combi; third key is the indivdual chain; value: vector of FattyAcidVO (vector because might be isotopically labeled)
  private Vector<String> getPermutedChainsSplitOnChainAndOhTypes(String ohId, String[] chains, boolean[] usedChains, 
      Hashtable<Integer,Hashtable<String,String>> acylChainsInLib,  Hashtable<Integer,Hashtable<String,String>> alkylChainsInLib,
      Hashtable<Integer,Hashtable<String,String>> alkenylChainsInLib,  Hashtable<Integer,Hashtable<String,String>> lcbChainsInLib,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableFAChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkylChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkenylChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableLCBChains) throws RulesException, LipidCombinameEncodingException {
    String ohTypeFrequencyId;
    String remainingOhId = null;
    boolean lastOccurence = false;
    if (ohId.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)==-1) {
      ohTypeFrequencyId = ohId;
      lastOccurence = true;
    }else{
      // the chain composition is analysed from the end, because there are typicall less LCB
      // chains, and this reduces the combinatorial space 
      remainingOhId = ohId.substring(0, ohId.lastIndexOf(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS));
      ohTypeFrequencyId = ohId.substring(ohId.lastIndexOf(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)+1);
    }
    int[] decoded = decodeTypeOHFrequencyString(ohTypeFrequencyId);
    short type = (short)decoded[0];
    int oh = decoded[1];
    int nrOfChains = decoded[2];
    Hashtable<String,String> lookup = null;
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>> values = null;
    if (type==LipidomicsConstants.CHAIN_TYPE_FA_ACYL) {
      lookup = acylChainsInLib.get(oh);
      values = availableFAChains.get(oh);
    }else if (type==LipidomicsConstants.CHAIN_TYPE_FA_ALKYL) {
      lookup = alkylChainsInLib.get(oh);
      values = availableAlkylChains.get(oh);
    }else if (type==LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL) {
      lookup = alkenylChainsInLib.get(oh);
      values = availableAlkenylChains.get(oh);      
    }else if (type==LipidomicsConstants.CHAIN_TYPE_LCB) {
      lookup = lcbChainsInLib.get(oh);
      values = availableLCBChains.get(oh);
    }
    //System.out.println("ohTypeFrequencyId: "+type+";"+nrOfChains);

    //this vector can hold same chain configurations
    Vector<String> availableChains = new Vector<String>();
    for (int i=0; i!=chains.length; i++) {
      if (!usedChains[i] && lookup.containsKey(chains[i]))
        availableChains.add(chains[i]);
    }
    if (availableChains.size()<nrOfChains)
      return null;
    Vector<String> combis = getAllPossibleChainCombinations(nrOfChains,availableChains);
    Vector<String> permutedResults = new Vector<String>();
    Hashtable<String,String> combisEncoded = new Hashtable<String,String>();
    //encode the combinations correspondingly for the combinatorial result String
    for (String combi:combis) {
      StringBuilder encoded = new StringBuilder();
      //System.out.println("combi: "+combi);
      for (String chain : combi.split(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) {
        if (encoded.length()>0) encoded.append(LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
        String[] cAndDbs;
        try {
          cAndDbs = StaticUtils.parseCAndDbsFromChainId(chain);
          //using iterator().next() for the prefix is fine here, as no chains with prefixes end up here.
          encoded.append(values.get(Integer.parseInt(cAndDbs[0])).get(Integer.parseInt(cAndDbs[1])).values().iterator().next().get(cAndDbs[2]).getChainIdDetailed(false,true));
        }
        catch (Exception e) {
          throw new RulesException(e);
        }
      }
      combisEncoded.put(combi, encoded.toString());
    }
    
    if (lastOccurence) {
      for (String combiEncoded : combisEncoded.values())
        permutedResults.add(combiEncoded);
    }else {
      boolean[] remainingChains;
      int i;
      for (String combi : combis) {
        remainingChains = usedChains.clone();
        //mark the used chain occurrences
        for (String oneChain : combi.split(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) {
          for (i=0; i!=chains.length; i++) {
            if (!remainingChains[i] && chains[i].equalsIgnoreCase(oneChain)) {
              remainingChains[i] = true;
              break;
            }
          }
        }
        Vector<String> fromOtherParts = getPermutedChainsSplitOnChainAndOhTypes(remainingOhId, chains,
            remainingChains, acylChainsInLib, alkylChainsInLib, alkenylChainsInLib, lcbChainsInLib, availableFAChains, availableAlkylChains,
            availableAlkenylChains,availableLCBChains);
        //if the other chain types do not return a result, the combination is not valid -> go the next one
        if (fromOtherParts==null || fromOtherParts.size()==0)
          continue;
        //since this combination is valid add all results to the return vector
        for (String other : fromOtherParts) {
          permutedResults.add(other+LipidomicsConstants.CHAIN_COMBI_SEPARATOR+combisEncoded.get(combi));
        }
      }
    }    
    return permutedResults;
  }
  
  
  /**
   * recursive method for detection all potential permutations of a given number of chains
   * @param nrOfChains the total number of chains
   * @param chains the available encoded chain strings
   * @return a list of all potential permutations of a given number of chains
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private Vector<String> getAllPossibleChainCombinations(int nrOfChains,Vector<String> chains) throws LipidCombinameEncodingException{
    Vector<String> combisWithRepeats = new Vector<String>();
    String chain;
    for (int i=0; i!=chains.size(); i++) {
      chain = chains.get(i);
      if (nrOfChains>1) {
        Vector<String> remainingChains = new Vector<String>(chains);
        remainingChains.remove(i);
        for (String otherPart : getAllPossibleChainCombinations(nrOfChains-1,remainingChains)) {
          combisWithRepeats.add(chain+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS+otherPart);
        }
      }else {
        combisWithRepeats.add(chain);
      }
    }
    Hashtable<String,String> unique = new Hashtable<String,String>();
    for (String combi : combisWithRepeats) {
      String sorted = StaticUtils.sortFASequenceUnassigned(combi,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
      unique.put(sorted, sorted);
    }
    Vector<String> combis = new Vector<String>();
    for (String combi : unique.keySet()) {
      combis.add(combi);
    }
    return combis;
  }
  
  
  /**
   * this method permutes the possible isotopic labeled chains on a known combination of defined chains (chain type and number OH is given)
   * @param permuteId the encoded combinatorial name string of chains
   * @param toAdd the hash table where the results shall be added; key the new combinatorial name that includes the isotope labels (prefix in FattyAcidVO)
   * @param availableFAChains the principally (from the chemical formula) possible acyl chains 
   * @param availableAlkylChains the principally (from the chemical formula) possible alkyl chains 
   * @param availableAlkenylChains the principally (from the chemical formula) possible alkenyl chains
   * @param availableLCBChains the principally (from the chemical formula) possible LCB chains
   * @param uniqueChains all the chains that are possible after the combinatorial check; key is the chain id that contains information about chain type and hydroxylation sites
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   * @throws ChemicalFormulaException when there is something wrong with the chemical formula
   */
  private void addIsotopeLabelsToMutations(String permuteId, Hashtable<String,Vector<FattyAcidVO>> toAdd, 
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableFAChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkylChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkenylChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableLCBChains,
      Hashtable<String,FattyAcidVO> uniqueChains) throws LipidCombinameEncodingException, ChemicalFormulaException 
  {
    Vector<FattyAcidVO> chains = StaticUtils.decodeLipidNamesFromChainCombi(permuteId);
    //first, it is evaluated how many different options are present for each chain position
    int isotopicLabels[] = new int[chains.size()];
    int[] divisors = new int[chains.size()];
    FattyAcidVO chain;
    Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> available;
    Hashtable<Integer,Hashtable<String,FattyAcidVO>> chainVOsByPosition = new Hashtable<Integer,Hashtable<String,FattyAcidVO>>();
    //System.out.println("-----------------------------------");
    //System.out.println(permuteId);
    for (int i=0; i!=chains.size(); i++){
      chain = chains.get(i);
      available = getLookupHashByChainType(chain.getChainType(),availableFAChains,availableAlkylChains,availableAlkenylChains,availableLCBChains);
      Hashtable<String,FattyAcidVO> faPlusPrefix = new Hashtable<String,FattyAcidVO>();
      Set<String> prefixes = available.get(chain.getOhNumber()).get(chain.getcAtoms()).get(chain.getDoubleBonds()).keySet();
      for (String prefix : prefixes)
      {
      	faPlusPrefix.put(prefix, available.get(chain.getOhNumber()).get(chain.getcAtoms()).get(chain.getDoubleBonds()).get(prefix).get(chain.getOxState()));
      }
      chainVOsByPosition.put(i,faPlusPrefix);
      isotopicLabels[i] = chainVOsByPosition.get(i).size();
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
    Hashtable<String,FattyAcidVO> chainVOs;
    Hashtable<String,String> permutedCombinations = new Hashtable<String,String>();
    for (int j=0; j!=iterations; j++){
      String combiName = "";
      Vector<FattyAcidVO> valueVOs = new Vector<FattyAcidVO>();
      Vector<String> singleCombiParts = new Vector<String>();
      for (int i=0; i!=chains.size(); i++){
        chain = chains.get(i);
        chainVOs = chainVOsByPosition.get(i);
        int currentIndex = (j/divisors[i])%isotopicLabels[i];
        int k=0;
        for (String prefix : chainVOs.keySet()){
          if (k==currentIndex){
            String partName = chainVOs.get(prefix).getChainId();
            combiName += partName+LipidomicsConstants.CHAIN_COMBI_SEPARATOR;
            valueVOs.add(chainVOs.get(prefix));
            singleCombiParts.add(partName);
            break;
          }
          k++;
        }
      }
      combiName = combiName.substring(0,combiName.length()-LipidomicsConstants.CHAIN_COMBI_SEPARATOR.length());
      //this is for filtering permuted double entries
      Vector<String> permutedNames = StaticUtils.getPermutedChainNames(singleCombiParts,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
      boolean isThere = false;
      for (String permutedName : permutedNames){
        if (permutedCombinations.containsKey(permutedName)){
          isThere = true;
          break;
        }
      }
      if (!isThere){
        for (String permutedName : permutedNames) permutedCombinations.put(permutedName, permutedName);
        boolean ok = true;
        //this check would help to remove unlikely chain combinations from isotopically labeled probes
        //the problem is however, that the FA chains always lose and oxygen in the linkage, and I have
        //no idea how I can put this in a general rule, without possibly wreaking havoc to other classes
        //check if the total chain combination is possible at all from the chemical elements
//        Hashtable<String,Integer> combiElements = new Hashtable<String,Integer>();
//        for (FattyAcidVO chainVO : valueVOs) {
//          Hashtable<String,Integer> oneChain = StaticUtils.categorizeFormula(chainVO.getFormula());
//          for (String element : oneChain.keySet()) {
//            if (combiElements.containsKey(element))
//              combiElements.put(element, combiElements.get(element)+oneChain.get(element));
//            else
//              combiElements.put(element, oneChain.get(element));
//          }
//        }
//        boolean ok = true;
//        Hashtable<String,Integer> formulaAmounts = StaticUtils.categorizeFormula(this.analyteFormulaWODeducts_);
//        for (String element : combiElements.keySet()) {
//          if (combiElements.get(element)>0 && (!formulaAmounts.containsKey(element) || formulaAmounts.get(element)<combiElements.get(element))) {
//            ok = false;
//            break;
//          }
//        }
        if (ok) {
          toAdd.put(combiName, valueVOs);
          for (FattyAcidVO vo : valueVOs) {
            if (!uniqueChains.containsKey(vo.getChainId()))
            	uniqueChains.put(vo.getChainId(), vo);
          }
//          System.out.println(combiName);
        }
      }
    }      
  }
  
  /**
   * selects the correct lookup hash depending on the chain type
   * @param chainType the type of chain LipidomicsConstants.CHAIN_TYPE_FA_ACYL|.CHAIN_TYPE_FA_ALKYL|.CHAIN_TYPE_FA_ALKENYL|.CHAIN_TYPE_LCB
   * @param availableFAChains the principally (from the chemical formula) possible acyl chains 
   * @param availableAlkylChains the principally (from the chemical formula) possible alkyl chains 
   * @param availableAlkenylChains the principally (from the chemical formula) possible alkenyl chains
   * @param availableLCBChains the principally (from the chemical formula) possible LCB chains
   * @return correct lookup hash
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private  Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> getLookupHashByChainType(short chainType,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableFAChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkylChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableAlkenylChains,
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> availableLCBChains) throws LipidCombinameEncodingException {
    if (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ACYL)
      return availableFAChains;
    else if  (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ALKYL) 
      return availableAlkylChains;
    else if  (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL){
      return availableAlkenylChains;
    }else if  (chainType==LipidomicsConstants.CHAIN_TYPE_LCB)
      return availableLCBChains;
    else
      throw new LipidCombinameEncodingException("The chain type \""+chainType+"\" is not allowed in the combination");
  }
  
  
  /**
   * returns only combinations where fragments for chains were detected
   * @param chains the detected chains
   * @param forbiddenChains chains that were removed by a mandatory OR combination
   * @return the combinations where fragments were detected
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private Vector<String> getChainCombinationsReleveantForTheseChains(Hashtable<String,String> chains, Set<String> forbiddenChains) throws LipidCombinameEncodingException{
    Vector<String> relevantChains = new Vector<String>();
    for (Hashtable<String,Vector<FattyAcidVO>> ohCombi : this.potentialChainCombinations_.values()) {
      for (String combiId : ohCombi.keySet()) {
        boolean relevant = false;
        for (String chain : chains.keySet()) {
          if (combiId.indexOf(chain)==-1)
            continue;
          for (FattyAcidVO vo : StaticUtils.decodeLipidNamesFromChainCombi(combiId)) {
            if (vo.getChainId().equalsIgnoreCase(chain)) {
              relevant = true;
              break;
            }
          }
          if (relevant)
            break;
        }
        //check whether the combiId contains any of the forbidden chains
        if (relevant) {
          for (FattyAcidVO vo : StaticUtils.decodeLipidNamesFromChainCombi(combiId)) {
            if (forbiddenChains.contains(vo.getChainId())) {
              relevant = false;
              break;
            }
          }
        }
        if (relevant) {
          relevantChains.add(combiId);
        }
      }
    }
    return relevantChains;
  }
  
  
  /**
   * returns true when all OH combinations contain class specific fragments
   * @return true when all OH combinations contain class specific fragments
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public boolean containAllOhCombinationsClassSpecificFragments() throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    boolean foundClassFragments = false;
    boolean oneCombiHasNoClassFragments = false;
    Hashtable<String,FragmentRuleVO> chainRules =  RulesContainer.getChainFragmentRules(ruleName_,rulesDir_);
    for (int[] ohs : this.possibleOhCombinations_) {
      String combiId = getOhCombiId(ohs);
      boolean containAllFaCombisClassFragments = false;
      Vector<Vector<Integer>> faCombis = allowedFaHydroxylationsCombinations_.get(combiId);
      if (faCombis.size()>0) {
        containAllFaCombisClassFragments = true;
        for (Vector<Integer> combi : faCombis) {
          boolean containsClassFragments = false;
          for (Integer oh : combi) {
            for (FragmentRuleVO vo : chainRules.values()) {
              if (vo.getChainType()!=LipidomicsConstants.CHAIN_TYPE_FA_ACYL && vo.getChainType()!=LipidomicsConstants.CHAIN_TYPE_FA_ALKYL &&
                  vo.getChainType()!=LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL)
                continue;
              if (!vo.hydroxylationValid(oh.shortValue()))
                continue;
              if (vo.isMandatory(oh.shortValue())==FragmentRuleVO.MANDATORY_CLASS) {
                containsClassFragments = true;
                break;
              }
            }
          }
          if (!containsClassFragments) {
            containAllFaCombisClassFragments = false;
            break;
          }
        }
      }
      boolean containAllLcbCombisClassFragments = false;
      Vector<Vector<Integer>> lcbCombis = allowedLcbHydroxylationsCombinations_.get(combiId);
      if (lcbCombis.size()>0) {
        containAllLcbCombisClassFragments = true;
        for (Vector<Integer> combi : lcbCombis) {
          boolean containsClassFragments = false;
          for (Integer oh : combi) {
            for (FragmentRuleVO vo : chainRules.values()) {
              if (vo.getChainType()!=LipidomicsConstants.CHAIN_TYPE_LCB)
                continue;
              if (!vo.hydroxylationValid(oh.shortValue()))
                continue;
              if (vo.isMandatory(oh.shortValue())==FragmentRuleVO.MANDATORY_CLASS) {
                containsClassFragments = true;
                break;
              }
            }
          }
          if (!containsClassFragments) {
            containAllLcbCombisClassFragments = false;
            break;
          }
        }
      }
      if (containAllFaCombisClassFragments||containAllLcbCombisClassFragments)
        foundClassFragments = true;
      else {
        oneCombiHasNoClassFragments = true;
        break;
      }
    }
    return (foundClassFragments && !oneCombiHasNoClassFragments);
  }
  
  
  /**
   * verifies if all class specific chain fragments were found
   * @param chain the chain object
   * @param foundFragments the found fragments for this chain
   * @return true if all class specific chain fragments were found
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public boolean areAllClassFragmentsForChainFound(FattyAcidVO chain, Hashtable<String,CgProbe> foundFragments) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    boolean foundAll = true;
    Hashtable<String,FragmentRuleVO> chainRules =  RulesContainer.getChainFragmentRules(ruleName_,rulesDir_);
    for (FragmentRuleVO fragVO : chainRules.values()) {
      if (fragVO.getChainType()!=chain.getChainType() || !fragVO.hydroxylationValid((short)chain.getOhNumber()) || fragVO.isMandatory((short)chain.getOhNumber())!=FragmentRuleVO.MANDATORY_CLASS)
        continue;
      if (foundFragments==null || !foundFragments.containsKey(fragVO.getName())) {
        foundAll = false;
        break;
      }
    }
    return foundAll;
  }
  
  
  /**
   * checks if there are any combiOh requirements in these chain rules
   * @return fragment rules were there are combiOh requirements
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Hashtable<String,FragmentRuleVO> getCombiOhFragments() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Hashtable<String,FragmentRuleVO> combiOhs = new Hashtable<String,FragmentRuleVO>();
    Hashtable<String,FragmentRuleVO> allFrags = RulesContainer.getChainFragmentRules(ruleName_,rulesDir_);
    for (String name : allFrags.keySet()) {
      FragmentRuleVO rule = allFrags.get(name);
      if (rule.hasCombiOhRequirements())
        combiOhs.put(name, rule);
    }
    return combiOhs;
  }
  

	 /**
	  * creates all LipidomicsConstants.CHAIN_MOD_COLUMN_NAME combinations possible for a specified number of chains
	  * @param oxStates a list of all possible LipidomicsConstants.CHAIN_MOD_COLUMN_NAMEs on a chain 
	  * @param chains the number of chains of the analyte
	  * @return the LipidomicsConstants.CHAIN_MOD_COLUMN_NAME combinations
	  */
	public  Set<List<String>> getOxCombinations(List<String> oxStates, int chains) 
	{
	    Set<List<String>> oxStateCombinations = new HashSet<List<String>>();
	    Set<List<String>> newCombinations;

	    int index = 0;
	    for(String i: oxStates) {
	        List<String> newList = new ArrayList<String>();
	        newList.add(i);
	        oxStateCombinations.add(newList);
	    }
	    index++;
	    while(index < chains) {
	        List<String> nextList = oxStates;
	        newCombinations = new HashSet<List<String>>();
	        for(List<String> first: oxStateCombinations) {
	            for(String second: nextList) {
	                List<String> newList = new ArrayList<String>();
	                newList.addAll(first);
	                newList.add(second);
	                newCombinations.add(newList);
	            }
	        }
	        oxStateCombinations = newCombinations;

	        index++;
	    }

	    return oxStateCombinations;
	}
	
	/**
	  * validates LipidomicsConstants.CHAIN_MOD_COLUMN_NAME combinations possible on chains against the LipidomicsConstants.CHAIN_MOD_COLUMN_NAME of the analyte 
	  * @param combs the LipidomicsConstants.CHAIN_MOD_COLUMN_NAME combinations possible on chains
	  * @return only LipidomicsConstants.CHAIN_MOD_COLUMN_NAME combinations that are the same as the analyte`s LipidomicsConstants.CHAIN_MOD_COLUMN_NAME
	  */
	public  Set<List<String>> ValidateCombinations(Set<List<String>> combs) 
	{
		Set<List<String>> validCombinations = new HashSet<List<String>>();
		ModificationParser mcp; 

		for(List<String> combi : combs)
		{
			//get a string of the combination like "O,4OH"
			String combiString = "";
			for(String s : combi)
			{
				if(!combiString.equals("") && !s.equals(""))
				{
					combiString+= ","+s;
				}
				else if(combiString.equals("") && !s.equals(""))
				{
					combiString=s;
				}
			}
			combiString = combiString.replaceAll("\\s", "");
			
			//get elemental change of the chain modifications
			mcp = new ModificationParser(combiString);
			mcp.parse();
			String combiFormula = mcp.getModificationComposition();
			
			//get elemental change of the analyte modification
			mcp = new ModificationParser(analyteOxState_.replaceAll("\\s", ""));
			mcp.parse();
			String analyteOxFormula = mcp.getModificationComposition();
			
			//compare chain modifications to analyte modifications
			if(combiFormula.equals(analyteOxFormula))
	        {
				validCombinations.add(combi);
	        }
		}
		
		return validCombinations;
	} 
}

