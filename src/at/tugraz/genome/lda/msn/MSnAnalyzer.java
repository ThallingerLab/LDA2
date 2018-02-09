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
import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.dbutilities.Base64;
import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.alex123.vos.TargetlistEntry;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.FragmentVO;
import at.tugraz.genome.lda.msn.vos.IntensityChainVO;
import at.tugraz.genome.lda.msn.vos.IntensityPositionVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.msn.vos.MSnDebugVO;
import at.tugraz.genome.lda.msn.vos.SharedMS1PeakVO;
import at.tugraz.genome.lda.msn.vos.SharedPeakContributionVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.quantification.LipidomicsChromatogram;
import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.lda.utils.FloatFloatVO;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.maspectras.quantification.CgChromatogram;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.ChromatogramReader;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * Central class for the verification of MSn evidence
 * returns the composition of the analyte (precursor)
 * @author Juergen Hartler
 *
 */
public class MSnAnalyzer
{
  
  /** value to discard chain fragments that are very small - in relation to the highest found chain fragment*/ 
  private double relativeChainCutoff_ = 0.01;
  
  /** directory of stored fragmentation rules */
  private String rulesDir_;
  /** the name of the lipid class */
  private String className_;
  /** the name of the modification */
  private String modName_;
  /** MS1 data including quantitation info for the lipid to be checked */
  private LipidParameterSet set_;
  /** object that holds MS data and can quantify fragments of interest 
   * at the end of the analysis, the MSn information is transferred to this object */
  private LipidomicsAnalyzer analyzer_;
  /** if true, stores the details why a class, chain, etc. was abandoned */
  private boolean debug_;
  /** the achieved status by MSn checks */
  private int status_;
  /** class providing the potential fragments for the intact lipid */
  private FragmentCalculator fragCalc_;
  /** provides information why fragments were discarded - only available in debug mode */
  private MSnDebugVO debugVO_;
  
  /** the detected head group fragments, the key is the fragment name */
  private Hashtable<String,CgProbe> headGroupFragments_;
  /** the intensity rules that were fulfilled by the head group fragments */
  private Hashtable<String,IntensityRuleVO> fulfilledHeadIntensityRules_;
  /** the detected chain fragments, the key is the fragment name */
  private Hashtable<String,Hashtable<String,CgProbe>> chainFragments_;
  /** the intensity rules that were fulfilled by the chain fragments */
  private Hashtable<String,Hashtable<String,IntensityChainVO>> fulfilledChainIntensityRules_;
  /** stores potential fatty acid chain combinations */
  private Vector<String> validChainCombinations_;
  /** should absolute rules be ignored - whenever there is an overlap of isobars assumed*/
  private boolean ignoreAbsolute_;
  
  /** position definitions for fatty acid combinations - first key: id of the chain combination; second hash: lookup from the position of the fatty acid in the combination key to the assigned position*/
  private Hashtable<String,Hashtable<Integer,Integer>> positionDefinition_;
  /** 
   * hash containing the fulfilled MSn evidence for a position assignment
   * first key: id of the chain combination
   * second key: name of fatty acid
   * third key: position where evidence is provided
   * values: rules that are fulfilled by MSn evidence
  */
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>> positionRecommendations_;
  /** 
   * hash containing the fulfilled MSn evidence for a position assignment
   * first key: id of the chain combination
   * second key: position where evidence is provided
   * values: rules that are fulfilled by MSn evidence
  */  
  private Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> posEvidenceAccordToProposedDefinition_;

  /** base peak values at the time of interest - key is the MSn level the base peak has to be applied to */
  private Hashtable<Integer,Float> basePeakValues_;
  
  /** the found spectra where an MSn pattern matches */
  private Vector<Float> spectraFound_;
  
  /** are there are any MSn spectra found for this precursor ion m/z value*/
  private boolean msnSpectraPresent_;
  
  /** stores the probes that correspond to spectra at various MS-levels; the key is the MS-level*/
  private Hashtable<Integer,Vector<CgProbe>> probesWithMSnSpectra_;
  
  /** reports for which MS-levels spectra are present*/
  Hashtable<Integer,Boolean> msLevels_;
  
  
  
    
  /**
   * constructor requires the name for the lipid class, MS1 data including quantitation info for the lipid to be checked,
   * and object that holds MS data and can quantify fragments of interest
   * checking of rules is immediately started as soon as the constructor is called
   * @param className name of the lipid class
   * @param modName name of the adduct
   * @param set MS1 data including quantitation info for the lipid to be checked
   * @param analyzer object that holds MS data and can quantify fragments of interest
   * @param quantVO the VO containing the original parameters of the quantitation request
   * @param readMSnSpectra should the MSn spectra be read, or are they already in the cache
   * @param ignoreAbsolute should absolute rules be ignored - whenever there is an overlap of isobars assumed
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * throws CgException errors from the quantitation process
   */
  public MSnAnalyzer(String className, String modName, LipidParameterSet set, LipidomicsAnalyzer analyzer, QuantVO quantVO, boolean readMSnSpectra, boolean ignoreAbsolute) throws RulesException, IOException, SpectrummillParserException, CgException {
    this(null, className,modName,set,analyzer,quantVO,ignoreAbsolute,readMSnSpectra,false);
  }

  /**
   * constructor requires the name for the lipid class, MS1 data including quantitation info for the lipid to be checked,
   * and object that holds MS data and can quantify fragments of interest
   * checking of rules is immediately started as soon as the constructor is called
   * @param rulesDir directory containing the fragmentation rules
   * @param className name of the lipid class
   * @param modName name of the adduct
   * @param set MS1 data including quantitation info for the lipid to be checked
   * @param analyzer object that holds MS data and can quantify fragments of interest
   * @param quantVO the VO containing the original parameters of the quantitation request 
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * throws CgException errors from the quantitation process
   */
  public MSnAnalyzer(String rulesDir, String className, String modName, LipidParameterSet set, LipidomicsAnalyzer analyzer) throws RulesException, IOException, SpectrummillParserException, CgException {
    this(rulesDir, className,modName,set,analyzer,null,false,true,false);
  }

  
  /**
   * constructor requires the name for the lipid class, MS1 data including quantitation info for the lipid to be checked,
   * and object that holds MS data and can quantify fragments of interest
   * checking of rules is immediately started as soon as the constructor is called
   * @param rulesDir directory containing the fragmentation rules
   * @param className name of the lipid class
   * @param modName name of the adduct
   * @param set MS1 data including quantitation info for the lipid to be checked
   * @param analyzer object that holds MS data and can quantify fragments of interest
   * @param quantVO object that holds the quantification request
   * @param ignoreAbsolute should absolute rules be ignored - whenever there is an overlap of isobars assumed
   * @param readMSnSpectra should the MSn spectra be read, or are they already in the cache
   * @param debug should the analyzer store the details why a class, chain, etc. was abandoned?
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * throws CgException errors from the quantitation process
   */
  public MSnAnalyzer(String rulesDir, String className, String modName, LipidParameterSet set, LipidomicsAnalyzer analyzer, 
      QuantVO quantVO, boolean ignoreAbsolute, boolean readMSnSpectra, boolean debug) throws RulesException, IOException, SpectrummillParserException, CgException {
    this(className,modName,analyzer,debug,ignoreAbsolute);
    this.rulesDir_  = rulesDir;
    this.set_ = set;
    this.debug_ = debug;
    fragCalc_ = null;
    msLevels_ =  MSnAnalyzer.prepareCachedSpectra(analyzer_, set_, readMSnSpectra);
    try{
      if (Settings.useAlex() && quantVO!=null && (quantVO instanceof TargetlistEntry) && ((TargetlistEntry)quantVO).hasAlex123FragmentsForClass()){
        //This is checking if there exist any rules - used only for definition of a base peak cutoff!!!! - I am not sure if I should remove this
        try{
          fragCalc_ = new FragmentCalculator(rulesDir_,className_,modName_,set_.getNameStringWithoutRt(),set_.getChemicalFormula(),set_.Mz[0]);          
        } catch (NoRuleException nrx){
        }
        this.checkMSnByAlexFragments((TargetlistEntry)quantVO,msLevels_);
      }else{
        fragCalc_ = new FragmentCalculator(rulesDir_,className_,modName_,set_.getNameStringWithoutRt(),set_.getChemicalFormula(),set_.Mz[0]);
        this.checkMSnEvidence(msLevels_);
      }
      transferResultsToLipidParameterSet();
    } catch (NoRuleException nrx){
      // if the rule is not present, the status is NO_MSN_PRESENT
      // and there should be no further checks
      return;
    }

  }
  
  /**
   * general constructor that sets values that are required everywhere
   * @param className name of the lipid class
   * @param modName name of the adduct
   * @param analyzer object that holds MS data and can quantify fragments of interest
   * @param debug should the analyzer store the details why a class, chain, etc. was abandoned?
   * @param ignoreAbsolute should absolute rules be ignored - whenever there is an overlap of isobars assumed
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * throws CgException errors from the quantitation process
   */
  private MSnAnalyzer(String className, String modName, LipidomicsAnalyzer analyzer, boolean debug, boolean ignoreAbsolute) throws RulesException, IOException, SpectrummillParserException, CgException {
    this.rulesDir_  = null;
    this.set_ = null;
    this.className_ = className;
    this.modName_ = modName;
    this.status_ = LipidomicsMSnSet.NO_MSN_PRESENT;
    this.relativeChainCutoff_ = LipidomicsConstants.getChainCutoffValue();
    analyzer_ = analyzer;
    this.ignoreAbsolute_ = ignoreAbsolute;
  }
  
  /**
   * this is the constructor for searching single MSn spectra for identifications -
   * it returns a list of retention times where the analyte may be
   * @param className name of the lipid class
   * @param modName name of the adduct
   * @param precursorMz the mass of the precursor 
   * @param precursorTolerance the tolerance m/z value for the precursor
   * @param name name of the analyte species
   * @param doubleBonds amount of double bonds of the analyte
   * @param analyteFormula chemical formula of the analyte
   * @param modificationFormula chemical formula of modification
   * @param charge charge state of the expected hit 
   * @param analyzer object that holds MS data and can quantify fragments of interest
   * @param debug should the analyzer store the details why a class, chain, etc. was abandoned?
   * @param ignoreAbsolute should absolute rules be ignored - whenever there is an overlap of isobars assumed
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   */
  public MSnAnalyzer(String className, String modName, double precursorMz, double precursorTolerance, String name, int doubleBonds, String analyteFormula,
      String modificationFormula, int charge, LipidomicsAnalyzer analyzer, boolean debug, boolean ignoreAbsolute) throws RulesException, IOException, SpectrummillParserException, CgException {
    this(className,modName,analyzer,debug,ignoreAbsolute);
    set_ =  new LipidParameterSet((float)precursorMz, name, new Integer(doubleBonds), modName, "not appropriate", analyteFormula, modificationFormula, new Integer(charge));
    this.rulesDir_  = null;
    msnSpectraPresent_ = false;
    scanAllSpectraForCandidates(precursorMz, precursorTolerance);
  }
  
  /**
   * 
   * @param precursorMz the mass of the precursor 
   * @param precursorTolerance the tolerance m/z value for the precursor
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   */
  private void scanAllSpectraForCandidates(double precursorMz, double precursorTolerance) throws RulesException, IOException, SpectrummillParserException, CgException {
    spectraFound_ = new Vector<Float>();
    Hashtable<Integer,Boolean> msLevels =  analyzer_.prepareMSnSpectraCache((float)(precursorMz-precursorTolerance), (float)(precursorMz+precursorTolerance));
    try{
      if (msLevels.size()<1){
        this.status_ = LipidomicsMSnSet.NO_MSN_PRESENT;
        return;
      }
      fragCalc_ = new FragmentCalculator(rulesDir_,className_,modName_,set_.getNameStringWithoutRt(),set_.getChemicalFormula(),set_.Mz[0]);
      Vector<Range> ranges = analyzer_.findSingleSpectraRanges(fragCalc_.getSpectrumLevelRange());
      if (ranges.size()>0) this.msnSpectraPresent_ = false;
      for (Range range : ranges){
        Vector<CgProbe> probes = new Vector<CgProbe>();
        CgProbe probe = new CgProbe(0,1);
        probe.AreaStatus = CgAreaStatus.OK;
        probe.Area = 100f;
        probe.AreaError = 0f;
        probe.Background = 0f;
        probe.Peak = (range.getStart()+range.getStop())/2f;
        probe.LowerValley = range.getStart();
        probe.UpperValley = range.getStop();
        probe.Mz = set_.Mz[0];
        probe.LowerMzBand = LipidomicsConstants.getCoarseChromMzTolerance(set_.Mz[0]);
        probe.UpperMzBand = LipidomicsConstants.getCoarseChromMzTolerance(set_.Mz[0]);
        probe.isotopeNumber = 0;
        probes.add(probe);
        set_.setProbes(probes);
        this.status_ = LipidomicsMSnSet.NO_MSN_PRESENT;
        this.checkMSnEvidence(msLevels);
        if (this.status_>LipidomicsMSnSet.DISCARD_HIT){
          spectraFound_.add(probe.Peak);
        }
      }
    } catch (NoRuleException nrx){
      // if the rule is not present, the status is NO_MSN_PRESENT
      // and there should be no further checks
      return;
    }
  }
  
  
  
  /**
   * verifies the evidence for the lipid class, the chain composition, and the position at the backbone
   * according to the defined fragmentation rules
   * @param msLevelsFromSpectralData the MS-levels that might be possible due to the available spectra
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   * @throws NoRuleException thrown if the rules are not there
   */
  private void checkMSnEvidence(Hashtable<Integer,Boolean> msLevelsFromSpectralData) throws RulesException, IOException, SpectrummillParserException, CgException, NoRuleException {
    probesWithMSnSpectra_ = performStandardInitialProcesses(msLevelsFromSpectralData);
    Hashtable<Integer,Boolean> msLevels = fragCalc_.correctMsLevelsForExistingFragments(msLevelsFromSpectralData,probesWithMSnSpectra_);
    if (probesWithMSnSpectra_==null || probesWithMSnSpectra_.size()<1 || msLevels.size()<1){
      this.status_ = LipidomicsMSnSet.NO_MSN_PRESENT;
      return;
    }
    if (fragCalc_.getChainCutoff()>=0) relativeChainCutoff_ = fragCalc_.getChainCutoff();
    basePeakValues_ = calculateBasePeakValuesIfRequired(msLevels);
    checkHeadGroupFragments(probesWithMSnSpectra_);
    if (status_== LipidomicsMSnSet.DISCARD_HIT && !debug_) return;
    checkChainFragments(probesWithMSnSpectra_);
    if (debug_) debugVO_.setSpectrumCoverageFulfilled(true);
    checkSpectrumCoverage(msLevels);
    if (status_!= LipidomicsMSnSet.FRAGMENTS_DETECTED && !debug_) return;
    if (fragCalc_.getAllowedChainPositions()>1)checkPositions();
    else if (status_!=LipidomicsMSnSet.DISCARD_HIT)status_ = LipidomicsMSnSet.POSITION_DETECTED;
  }
  

  /**
   * initializes the hashtables, checks if MSn spectra are present, and if they are within the region of the quantified peak
   * @param msLevels the msLevels to be checked
   * @return only probes where MSn spectra are present; key is the msLevel
   * @throws CgException
   */
  private Hashtable<Integer,Vector<CgProbe>> performStandardInitialProcesses(Hashtable<Integer,Boolean> msLevels) throws CgException{
    headGroupFragments_ = new Hashtable<String,CgProbe>();
    fulfilledHeadIntensityRules_ = new Hashtable<String,IntensityRuleVO>();
    chainFragments_ = new Hashtable<String,Hashtable<String,CgProbe>>();
    fulfilledChainIntensityRules_ = new Hashtable<String,Hashtable<String,IntensityChainVO>>();
    validChainCombinations_ = new Vector<String>();
    basePeakValues_ = new Hashtable<Integer,Float>();
    positionDefinition_ = new Hashtable<String,Hashtable<Integer,Integer>>();
    positionRecommendations_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>>();
    posEvidenceAccordToProposedDefinition_ = new Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>();
    if (debug_) debugVO_ = new MSnDebugVO();
    if (msLevels.size()<1){
      this.status_ = LipidomicsMSnSet.NO_MSN_PRESENT;
      return null;
    }
    Hashtable<Integer,Vector<CgProbe>> probesWithMSnSpectra = new Hashtable<Integer,Vector<CgProbe>>();
    for (int msLevel : msLevels.keySet()){
      Vector<CgProbe> probes = new Vector<CgProbe>();
      for (CgProbe probe : set_.getIsotopicProbes().get(0)){
        if (analyzer_.areMSnSpectraInThisRegion(probe.LowerValley,probe.UpperValley,msLevel)){
          probes.add(probe);
        }
        if (probes.size()>0)probesWithMSnSpectra.put(msLevel, probes);
      }
    }
    return probesWithMSnSpectra;
  }
  
  /**
   * checks for fragments by Alex123 target lists
   * @param quantVO object that holds the quantification request
   * @param msLevels the msLevels to be checked
   * @throws CgException errors from the quantitation process
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  private void checkMSnByAlexFragments(TargetlistEntry quantVO, Hashtable<Integer,Boolean> msLevels) throws CgException,
    RulesException, NoRuleException, IOException, SpectrummillParserException{
    if (quantVO.getMsnFragments()==null || quantVO.getMsnFragments().size()==0){
      this.status_ = LipidomicsMSnSet.NO_MSN_PRESENT;
      return;      
    }
    probesWithMSnSpectra_ = performStandardInitialProcesses(msLevels);
    if (probesWithMSnSpectra_==null || probesWithMSnSpectra_.size()<1){
      this.status_ = LipidomicsMSnSet.NO_MSN_PRESENT;
      return;
    }
    if (fragCalc_!=null)
      basePeakValues_ = calculateBasePeakValuesIfRequired(msLevels);
//    for (Integer level : basePeakValues_.keySet()){
//      System.out.println("Class: "+quantVO.getAnalyteClass());
//      System.out.println("Basepeakvalue "+level+": "+basePeakValues_.get(level));
//    }
    
    Hashtable<String,String> verifiedMolecularSpecies = new Hashtable<String,String>();
    Hashtable<Integer,Hashtable<String,Hashtable<String,TargetlistEntry>>> allFragments = quantVO.getMsnFragments();
    boolean foundAnyFragments = false;
    boolean foundHeadFragments = false;
    boolean foundChainFragments = false;
    for (Integer msLevel : allFragments.keySet()){
      if (!msLevels.get(msLevel)) continue;
      Hashtable<String,Hashtable<String,TargetlistEntry>> fragments = allFragments.get(msLevel);
      for (String fragmentName : fragments.keySet()){
        //System.out.println("Fragment: "+fragmentName);
        Hashtable<String,TargetlistEntry> molSpecies = fragments.get(fragmentName);
        TargetlistEntry oneFragment = molSpecies.values().iterator().next();
        CgProbe probe = analyzer_.calculateMs2Area(oneFragment.getAnalyteMass(), oneFragment.getFragmentFormula(),
            oneFragment.getMsLevel(), oneFragment.getCharge(), false, probesWithMSnSpectra_.get(oneFragment.getMsLevel()));
        if (probe.AreaStatus == CgAreaStatus.OK && (fragCalc_==null || checkBasePeakCutoff(probe,oneFragment.getMsLevel())) ){
          foundAnyFragments = true;
          //if (oneFragment.getStructure().equalsIgnoreCase(quantVO.getAnalyteClass()) || quantVO.getAnalyteClass().equalsIgnoreCase("IS "+oneFragment.getStructure())){
          if (!oneFragment.getStructure().startsWith("FA ") && !oneFragment.getStructure().startsWith("LCB ")){
            headGroupFragments_.put(fragmentName,probe);
            foundHeadFragments = true;
          }else{
            String chainName = oneFragment.getStructure();////oneFragment.getFragmentSumComposition();
            Hashtable<String,CgProbe> fragmentsOfChain = new Hashtable<String,CgProbe>();
            if (chainFragments_.containsKey(chainName)) fragmentsOfChain = chainFragments_.get(chainName);
            fragmentsOfChain.put(fragmentName, probe);
            chainFragments_.put(chainName, fragmentsOfChain);
            for (String molSpec : molSpecies.keySet()){
              foundChainFragments = true;
              verifiedMolecularSpecies.put(molSpec, molSpec);            
            }
          }
        }
      }
    }
    if (!foundAnyFragments) this.status_ = LipidomicsMSnSet.DISCARD_HIT;
    if (foundHeadFragments && status_!=LipidomicsMSnSet.DISCARD_HIT) this.status_ = LipidomicsMSnSet.HEAD_GROUP_DETECTED;
    if (foundChainFragments && status_!=LipidomicsMSnSet.DISCARD_HIT) this.status_ = LipidomicsMSnSet.FRAGMENTS_DETECTED;
    for (String molSpecies : verifiedMolecularSpecies.keySet()){
      validChainCombinations_.add(molSpecies);
    }
  }
  
  /**
   * verifies if the head group fragment rules are fulfilled
   * if mandatory rules are infringed, the status changes to DISCARD_HIT
   * @param probesWithMSnSpectra MS1 identifications where MSn spectra can be extracted
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   * @throws NoRuleException thrown if the rules are not there
   */
  private void checkHeadGroupFragments(Hashtable<Integer,Vector<CgProbe>> probesWithMSnSpectra) throws RulesException, IOException, SpectrummillParserException, CgException, NoRuleException {
    //calculate the areas of the head groups
    Hashtable<Boolean,Vector<FragmentVO>> headFragments = fragCalc_.getHeadFragments();
    Vector<FragmentVO> mandatoryHeadFragments = headFragments.get(true);
    Vector<FragmentVO> addHeadFragments = headFragments.get(false);
    boolean foundHeadFragments = false;
    for (FragmentVO fragment : mandatoryHeadFragments){
      if (!probesWithMSnSpectra.containsKey(fragment.getMsLevel())) continue;
      CgProbe probe = analyzer_.calculateMs2Area(fragment.getMass(), fragment.getFormula(), fragment.getMsLevel(), fragment.getCharge(), fragment.isFromOtherSpecies(), probesWithMSnSpectra.get(fragment.getMsLevel()));
//      System.out.println("Mand: "+fragment.getName()+";"+fragment.getMass()+";"+probe.Area+";"+probe.AreaStatus+";"+checkBasePeakCutoff(probe,fragment.getMsLevel()));
      if (probe.AreaStatus == CgAreaStatus.OK && checkBasePeakCutoff(probe,fragment.getMsLevel())){
        foundHeadFragments = true;
        headGroupFragments_.put(fragment.getName(),probe);
      }else{
        status_ = LipidomicsMSnSet.DISCARD_HIT;
        if (debug_){
          debugVO_.addDiscardedHeadGroupFragment(fragment.getName(),findAreaDiscardReason(probe, fragment));
        }else{
          return;
        }
      }
    }
    for (FragmentVO fragment : addHeadFragments){
      if (!probesWithMSnSpectra.containsKey(fragment.getMsLevel())) continue;
      CgProbe probe = analyzer_.calculateMs2Area(fragment.getMass(), fragment.getFormula(), fragment.getMsLevel(), fragment.getCharge(), fragment.isFromOtherSpecies(), probesWithMSnSpectra.get(fragment.getMsLevel()));
//      System.out.println("Add: "+fragment.getName()+";"+fragment.getMass()+";"+probe.Area+";"+probe.AreaStatus+";"+checkBasePeakCutoff(probe,fragment.getMsLevel()));
      if (probe.AreaStatus == CgAreaStatus.OK && checkBasePeakCutoff(probe,fragment.getMsLevel())){
        foundHeadFragments = true;
        headGroupFragments_.put(fragment.getName(),probe);
      } else if (debug_){
        debugVO_.addDiscardedHeadGroupFragment(fragment.getName(),findAreaDiscardReason(probe, fragment));
      }
    }
    if (foundHeadFragments && status_!=LipidomicsMSnSet.DISCARD_HIT) this.status_ = LipidomicsMSnSet.HEAD_GROUP_DETECTED;
    
    //check if intensity rules are OK
    Vector<IntensityRuleVO> intRules = fragCalc_.getHeadIntensityRules();
    for (IntensityRuleVO intRule : intRules){
      if (intRule.enoughFragmentsForRuleEvaluationFound(headGroupFragments_)){
        if (this.ignoreAbsolute_ && intRule.isAbsoluteComparison()){
          if (intRule.isRuleFulfilled(headGroupFragments_,getBasepeakIfRequired(intRule)))
            fulfilledHeadIntensityRules_.put(intRule.getRuleIdentifier(), intRule);
          continue;
        }
        if (intRule.isRuleFulfilled(headGroupFragments_,getBasepeakIfRequired(intRule))){
          fulfilledHeadIntensityRules_.put(intRule.getRuleIdentifier(), intRule);
        }else{
          if (debug_){
            this.debugVO_.addViolatedHeadRule(intRule);
          }
          if (intRule.isMandatory()){
            status_ = LipidomicsMSnSet.DISCARD_HIT;
            return;
          }
        }
      }
    }
  }
  
  /**
   * returns the area of the base peak, if this intensity rule contains a base peak place holder
   * @param ruleVO the IntensityRuleVO where we want to get the base peak
   * @return base peak area
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  private Float getBasepeakIfRequired(IntensityRuleVO ruleVO) throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Float basePeak = null;
    if (ruleVO.containsBasePeak()){
      int msLevel = fragCalc_.getFragmentRuleByName(ruleVO.getAnyNonBasePeakName()).getMsLevel();
      if (basePeakValues_.containsKey(msLevel)) basePeak = basePeakValues_.get(msLevel);
    }
    return basePeak;
  }
  
  /**
   * verifies if the chain fragment rules are fulfilled
   * stores all verified chain combinations, and discards hits below the cutoff
   * @param probesWithMSnSpectra MS1 identifications where MSn spectra can be extracted
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   */
  private void checkChainFragments(Hashtable<Integer,Vector<CgProbe>> probesWithMSnSpectra) throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException {
    Vector<FattyAcidVO> fas = fragCalc_.getPossibleChainObjects();
    Vector<IntensityRuleVO> intRules = fragCalc_.getChainIntensityRulesSameChain();
    for (FattyAcidVO fa : fas){
      Hashtable<Integer,Hashtable<Boolean,Vector<FragmentVO>>> allChainFragments = fragCalc_.getChainFragments(fa);
      for (Integer chainType : allChainFragments.keySet()){
        Hashtable<Boolean,Vector<FragmentVO>> chainFragments = allChainFragments.get(chainType);
        Vector<FragmentVO> mandatoryChainFragments = chainFragments.get(true);
        Vector<FragmentVO> addChainFragments = chainFragments.get(false);
        boolean discardChain = false;
        boolean foundChainFragments = false;
        Hashtable<String,CgProbe> foundFragments = new Hashtable<String,CgProbe>();
//      System.out.println(fa.getName()+";"+mandatoryChainFragments.size());
        for (FragmentVO fragment : mandatoryChainFragments){
        //System.out.println(fa+";"+fragment.getName()+" ; "+fragment.getMass());
          if (!probesWithMSnSpectra.containsKey(fragment.getMsLevel())) continue;
          CgProbe probe = analyzer_.calculateMs2Area(fragment.getMass(), fragment.getFormula(), fragment.getMsLevel(), fragment.getCharge(), fragment.isFromOtherSpecies(), probesWithMSnSpectra.get(fragment.getMsLevel()));
//        if (fa.getName().equalsIgnoreCase("16:0")) System.out.println("!!! "+fa.getName()+";"+fragment.getMass()+";"+fragment.getName()+";"+probe.Area);
          if (probe.AreaStatus == CgAreaStatus.OK && checkBasePeakCutoff(probe,fragment.getMsLevel())){
            foundChainFragments = true;
            foundFragments.put(fragment.getName(),probe);
//          System.out.println(fa.getName()+";"+fragment.getName()+";"+fragment.getMass()+";"+probe.Area+"; Mand");
          }else{
            if (debug_) debugVO_.addViolatedChainFragment(fa.getName(), fragment.getName(), findAreaDiscardReason(probe, fragment));
            discardChain = true;
            break;
          }
        }
        if (discardChain) continue;
        for (FragmentVO fragment : addChainFragments){
//          if (chainType==FragmentRuleVO.ALKYL_CHAIN && fa.getName().equalsIgnoreCase("18:1")) System.out.println("!!! 18:1: "+fragment.getMass());
          if (!probesWithMSnSpectra.containsKey(fragment.getMsLevel())) continue;
          CgProbe probe = analyzer_.calculateMs2Area(fragment.getMass(), fragment.getFormula(), fragment.getMsLevel(), fragment.getCharge(), fragment.isFromOtherSpecies(), probesWithMSnSpectra.get(fragment.getMsLevel()));
          if (probe.AreaStatus == CgAreaStatus.OK && checkBasePeakCutoff(probe,fragment.getMsLevel())){
//          System.out.println(fa.getName()+";"+fragment.getName()+";"+fragment.getMass()+";"+probe.Area+";"+probe.AreaStatus+"; Add");
            foundChainFragments = true;
            foundFragments.put(fragment.getName(),probe);
          } else {
            if (debug_) debugVO_.addViolatedChainFragment(fa.getName(), fragment.getName(), findAreaDiscardReason(probe, fragment));          
          }
        }
        if (foundChainFragments){
          Hashtable<String,IntensityChainVO> fulfilledChaindIntensityRules = new Hashtable<String,IntensityChainVO>();
          for (IntensityRuleVO intRule : intRules){
            if (intRule.getChainType()!=chainType) continue;
            Hashtable<String,CgProbe> allFragments = new Hashtable<String,CgProbe>(foundFragments);
            allFragments.putAll(headGroupFragments_);
            if (intRule.enoughFragmentsForRuleEvaluationFound(allFragments)){
              if (this.ignoreAbsolute_ && intRule.isAbsoluteComparison()){
                if (intRule.isRuleFulfilled(headGroupFragments_,getBasepeakIfRequired(intRule)))
                  fulfilledChaindIntensityRules.put(intRule.getRuleIdentifier(), new IntensityChainVO(intRule,fa.getName()));
                continue;
              }
              if (intRule.isRuleFulfilled(allFragments,getBasepeakIfRequired(intRule))){
                fulfilledChaindIntensityRules.put(intRule.getRuleIdentifier(), new IntensityChainVO(intRule,fa.getName()));
              }else{
                if (debug_) debugVO_.addViolatedChainRule(fa.getName(), intRule); 
                if (intRule.isMandatory()){
                  discardChain = true;
                  break;
                }
              }
            }
          }
//        System.out.println("!! "+fa.getName()+" ; "+discardChain);
          if (discardChain) continue;
          String chainName = FragmentCalculator.getAcylAlkylOrAlkenylName(fa.getName(), chainType);
          chainFragments_.put(chainName, foundFragments);
          if (fulfilledChaindIntensityRules.size()>0)
            fulfilledChainIntensityRules_.put(chainName, fulfilledChaindIntensityRules);
        }
      }
    }
    if (chainFragments_.size()==0) return;
    
    // this checks if intensity relationships of different chains are fulfilled
    Hashtable<String,Vector<String>> combis = fragCalc_.getChainFragmentCombinationsWithDetectedEvidence(chainFragments_.keySet(),chainFragments_);
    Vector<IntensityRuleVO> diffRules = fragCalc_.getChainIntensityRulesDiffChain();
    if (combis.size()>0 && combis.values().iterator().next().size()>1 && diffRules.size()>0){
      Hashtable<String,Vector<String>> intRulesFulFilled = new Hashtable<String,Vector<String>>();
      for (String combiKey: combis.keySet()){
        Vector<String> combi = combis.get(combiKey);
        if (combi.size()<2) continue;
        boolean discardCombi = false;
        Hashtable<String,IntensityChainVO> fulfilledChaindIntensityRules = new Hashtable<String,IntensityChainVO>();
        for (int i=0; i<combi.size(); i++){
          for (int j=i+1; j<combi.size(); j++){
            String fa1 = combi.get(i);
            String fa2 = combi.get(j);
            Hashtable<String,CgProbe> allFragments1 = new Hashtable<String,CgProbe>();
            if (chainFragments_.containsKey(fa1)) allFragments1 = new Hashtable<String,CgProbe>(chainFragments_.get(fa1));
            allFragments1.putAll(headGroupFragments_);
            Hashtable<String,CgProbe> allFragments2 = new Hashtable<String,CgProbe>();
            if (chainFragments_.containsKey(fa2)) allFragments2 = new Hashtable<String,CgProbe>(chainFragments_.get(fa2));
            allFragments2.putAll(headGroupFragments_);
            for (IntensityRuleVO intRule : diffRules){
              boolean[] enoughAndFirst = intRule.enoughFragmentsForRuleEvaluationFound(allFragments1,allFragments2);
              if (enoughAndFirst[0]){
                if (this.ignoreAbsolute_ && intRule.isAbsoluteComparison()){
                  if (intRule.isRuleFulfilled(headGroupFragments_,getBasepeakIfRequired(intRule))){
                    String biggerFA = fa2;
                    String smallerFA = fa1;
                    if (enoughAndFirst[1]){
                      biggerFA = fa1;
                      smallerFA = fa2;
                    }
                    fulfilledChaindIntensityRules.put(intRule.getRuleIdentifier(), new IntensityChainVO(intRule,biggerFA,smallerFA));                    
                  }
                  continue;
                }
                if (intRule.isRuleFulfilled(allFragments1,allFragments2,getBasepeakIfRequired(intRule))){
                  String biggerFA = fa2;
                  String smallerFA = fa1;
                  if (enoughAndFirst[1]){
                    biggerFA = fa1;
                    smallerFA = fa2;
                  }
                  fulfilledChaindIntensityRules.put(intRule.getRuleIdentifier(), new IntensityChainVO(intRule,biggerFA,smallerFA));
                }else{
                  if (debug_) debugVO_.addViolatedChainRule(combiKey, intRule); 
                  if (intRule.isMandatory()){
                    discardCombi = true;
                    break;
                  }
                }
              }
            }
          }
        }
        if (!discardCombi){
          if (fulfilledChaindIntensityRules.size()>0)
            fulfilledChainIntensityRules_.put(combiKey, fulfilledChaindIntensityRules);
          intRulesFulFilled.put(combiKey, combi);
        }
      }
      combis = intRulesFulFilled;
    }
    
    // remove hits that cannot occur in any combination
    // and count if the FA occurs in more than one combination
    int removedKeys = 1;
    // the while is necessary for classes with more than 2 fatty acyls where on chain fragment might be shared by several combinations
    while (removedKeys>0){
      validChainCombinations_ = new Vector<String>();
      removedKeys = 0;
      Hashtable<String,String> allowedFAs = new Hashtable<String,String>();
      Hashtable<String,Integer> faChainOccurrences = new Hashtable<String,Integer>();
      for (Vector<String> combiFAs : combis.values()){
        for (String combiFA : combiFAs){
          allowedFAs.put(combiFA, combiFA);
          int count = 0;
          if (faChainOccurrences.containsKey(combiFA)) count = faChainOccurrences.get(combiFA);
          count++;
          faChainOccurrences.put(combiFA,count);
        }
      }
      this.removeNotNecessaryFragments(allowedFAs,true);
      this.removeNotNecessaryDiffIntensityRules(combis.keySet());
      if (combis.size()==0) return;
      // area values for each fragment
      Hashtable<String,Double> areas = new Hashtable<String,Double>();
      // calculate a total area for each chain fragment
      for (String key : chainFragments_.keySet()){
        double area = 0f;
        for (CgProbe probe : chainFragments_.get(key).values()){
          double oneArea = (double)probe.Area;
          area += oneArea;
        }
        areas.put(key, area);
      }
    
      Hashtable<String,Double> combiAreas = new Hashtable<String,Double>();
      double highestIntensity = 0d;
      for (String key : combis.keySet()){
        double relative = 0d;
        for (String combiFA : combis.get(key)){
          if (areas.containsKey(combiFA))
            relative += areas.get(combiFA)/((double)faChainOccurrences.get(combiFA));
        }
        if (relative>highestIntensity) highestIntensity = relative;
        combiAreas.put(key, relative);
      }
//    double totalArea = 0d;
      for (String key: new Vector<String>(combiAreas.keySet())){
        double intensity = combiAreas.get(key);
        if (intensity>=(relativeChainCutoff_*highestIntensity)){
//        totalArea+=intensity;
        }else{
          if (debug_) debugVO_.addViolatedCombinations(key, MSnDebugVO.COMBINATION_LOWER_CHAIN_CUTOFF);
          removedKeys++;
          combiAreas.remove(key);
          combis.remove(key);
        }
      }
      allowedFAs = new Hashtable<String,String>();
      for (String key: new Vector<String>(combiAreas.keySet())){
        validChainCombinations_.add(key);
//      relativeIntensityOfCombination_.put(key, combiAreas.get(key)/totalArea);
        for (String combiFA : combis.get(key)){
          allowedFAs.put(combiFA, combiFA);
        }      
      }
      removeNotNecessaryFragments(allowedFAs,false);
      this.removeNotNecessaryDiffIntensityRules(combiAreas.keySet());
    }
    if (status_!=LipidomicsMSnSet.DISCARD_HIT) this.status_ = LipidomicsMSnSet.FRAGMENTS_DETECTED;
  }
  
  /**
   * removes fatty acid chains from the result that are not possible
   * @param allowedFAs fatty acids where the evidence is OK
   * @param addToDiscard for debugging purposes only - chains that are not allowed anymore are stored in the debugVO_
   */
  private void removeNotNecessaryFragments(Hashtable<String,String> allowedFAs, boolean addToDiscard){
    Vector<String> chains = new Vector<String>(chainFragments_.keySet());
    for (String chain: chains){
      if (!allowedFAs.containsKey(chain)){
        if (debug_ && addToDiscard) debugVO_.addViolatedChainFragment(chain, MSnDebugVO.ID_COMBI_REMOVE, MSnDebugVO.NO_CHAIN_COMINATION_POSSIBLE);
        chainFragments_.remove(chain);
//        String chainWithoutLinkInfo = new String(chain);
//        if (chain.startsWith(FragmentCalculator.ALKYL_PREFIX)) chainWithoutLinkInfo = chainWithoutLinkInfo.substring(FragmentCalculator.ALKYL_PREFIX.length());
//        if (chain.startsWith(FragmentCalculator.ALKENYL_PREFIX)) chainWithoutLinkInfo = chainWithoutLinkInfo.substring(FragmentCalculator.ALKENYL_PREFIX.length());
        fulfilledChainIntensityRules_.remove(chain);
      }
    }
  }
  
  /**
   * removes rules that have belonged to a discarded chain combination
   * @param allowedCombis the allowed chain combinations
   */
  private void removeNotNecessaryDiffIntensityRules(Set<String> allowedCombis){
    Hashtable<String,String> allowed = new Hashtable<String,String>();
    for (String combi : allowedCombis) allowed.put(combi, combi);
    Vector<String> combisToRemove = new Vector<String>();
    for (String key : fulfilledChainIntensityRules_.keySet()){
      if (key.indexOf(LipidomicsConstants.FA_SEPARATOR)==-1) continue;
      if (!allowed.containsKey(key)) combisToRemove.add(key);
    }
    
    for (String combi : combisToRemove) fulfilledChainIntensityRules_.remove(combi);
  }
  
  /**
   * 
   * @return the status if the verification process - corresponds to which level of evidence could have been extracted
   */
  public int checkStatus(){
    return this.status_;
  }
  
  /**
   * 
   * @return the found and quantified head group fragments - the key is the fragment name
   */
  public Hashtable<String,CgProbe> getHeadGroupFragments(){
    return this.headGroupFragments_;
  }
  
  /**
   * 
   * @return intensity rules that were fulfilled by the head group fragments
   */
  public Hashtable<String,IntensityRuleVO> getFulfilledHeadIntensityRules(){
    return this.fulfilledHeadIntensityRules_;
  }
  
  /**
   * 
   * @return the found and quantified chain fragments - the key is the fragment name
   */
  public Hashtable<String,Hashtable<String,CgProbe>> getChainFragments(){
    return this.chainFragments_;
  } 
  
  /**
   * 
   * @return intensity rules that were fulfilled by the chain fragments
   */
  public Hashtable<String,Hashtable<String,IntensityChainVO>> getFulfilledChainIntensityRules(){
  return this.fulfilledChainIntensityRules_;
  }
  
  /**
   * 
   * @return position recommendations
   */
  public Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>> getPositionRecommendations(){
    return this.positionRecommendations_;
  }
  
  
  /**
   * ATTENTION: prepareMSnSpectraCache has to be called before
   * calculates base peaks if they are required
   * this is the case if in the rules the BasePeakCutoff is greater than 0 or
   * if any of the rules contains the $BASEPEAK fragment
   * @param msLevels for which MSn levels should the base peak values be extracted
   * @return intensities for the base peaks - key is the MSn level where the peak was detected
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   */
  private Hashtable<Integer,Float> calculateBasePeakValuesIfRequired(Hashtable<Integer,Boolean> msLevels) throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException {
    Hashtable<Integer,Float> basePeakValues = new Hashtable<Integer,Float>();
    Vector<Integer> levels = fragCalc_.getBasePeakIntRuleLevels(msLevels);
    if (levels.size()>0){
      basePeakValues = analyzer_.extractBasePeakValues(levels,set_.getIsotopicProbes().get(0));
    }
    return basePeakValues;
  }
  
  /**
   * checks if a fragment is strong enough in relation the base peak
   * the relative threshold is defined by BasePeakCutoff in the rules file
   * @param probe identified fragment
   * @param msLevel MSn level of the fragment
   * @return true if the intensity of the fragment is sufficiently strong
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   */
  private boolean checkBasePeakCutoff(CgProbe probe, int msLevel)throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException {
    double cutoff = fragCalc_.getBasePeakCutoff();
    if (cutoff>0){
      if (probe.Area > basePeakValues_.get(msLevel)*((float)cutoff)) return true;
      else return false;
    }else return true;
  }
    
  /**
   * verifies if the position rules are fulfilled
   * all mandatory rules must be fulfilled (if the affected fragments are present)
   * if rules are not mandatory, the position for the combinations where more rules are fulfilled
   * if there is equality in contradicting assignments - no position is assigned
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   */
  private void checkPositions() throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException {
    Vector<IntensityRuleVO> intRules = fragCalc_.getPositionIntensityRules();
    Hashtable<String,Hashtable<Boolean,Vector<IntensityRuleVO>>> positionSpecificRules = new Hashtable<String,Hashtable<Boolean,Vector<IntensityRuleVO>>>();
    Hashtable<String,Hashtable<Boolean,Vector<IntensityRuleVO>>> positionSpecificSingleRules = new Hashtable<String,Hashtable<Boolean,Vector<IntensityRuleVO>>>();
    for (IntensityRuleVO rule : intRules){
      List<Integer> affectedPositions = new ArrayList<Integer>();
      affectedPositions.add(rule.getBiggerPosition());
      if (rule.getBiggerPosition() != rule.getSmallerPosition())
        affectedPositions.add(rule.getSmallerPosition());
      String identifier = getIdentifierForPositions(affectedPositions);
      if (rule.getBiggerPosition() == rule.getSmallerPosition())
        addToPositionSpecificRuleHash(positionSpecificSingleRules,identifier,rule);
      else
        addToPositionSpecificRuleHash(positionSpecificRules,identifier,rule);
    }
    boolean foundOnePosition = false;
    for (String combiKey : validChainCombinations_){
      String[] fas = StaticUtils.getFAsFromCombiName(combiKey);
      Hashtable<String,Integer> chainOccurenceInCombi = new Hashtable<String,Integer>();
      for (String fa: fas){
        int count = 0;
        if (chainOccurenceInCombi.containsKey(fa)) count = chainOccurenceInCombi.get(fa);
        count++;
        chainOccurenceInCombi.put(fa, count);
      }
      // in this case, the chain contains the same fragments only
      if (chainOccurenceInCombi.size()==1 && fragCalc_.getAllowedChainPositions()==fragCalc_.getAmountOfChains()){
        foundOnePosition = true;
        Hashtable<Integer,Integer> definitions = new Hashtable<Integer,Integer>();
        for (int i=0; i!=fas.length;i++){
          definitions.put(i, i);
        }
        positionDefinition_.put(combiKey, definitions);
        posEvidenceAccordToProposedDefinition_.put(combiKey, new Hashtable<Integer,Vector<IntensityPositionVO>>());
        continue;
      }
      // this is the normal case, where the rules have to be applied
      int[] positions = detectPostions(combiKey,chainOccurenceInCombi,positionSpecificRules,positionSpecificSingleRules);
      if (positions!=null){
        boolean foundAtLeastOnePosition = false;
        for (int position : positions){
          if (position>-1) foundAtLeastOnePosition = true;
        }
        if (foundAtLeastOnePosition){
          foundOnePosition = true;
          Hashtable<Integer,Integer> definitions = new Hashtable<Integer,Integer>();
          for (int i=0; i!=fas.length;i++){
            if (positions[i]>-1)
              definitions.put(i, (positions[i]-1));
          }
          positionDefinition_.put(combiKey, definitions); 
        }
      }
      
    }
    if (foundOnePosition && status_!=LipidomicsMSnSet.DISCARD_HIT) status_ =  LipidomicsMSnSet.POSITION_DETECTED;
    
//    for (String combi : posEvidenceAccordToProposedDefinition_.keySet()){
//      System.out.println("1. "+combi);
//      Hashtable<Integer,Vector<IntensityRuleVO>> positions = posEvidenceAccordToProposedDefinition_.get(combi);
//      String[] fas = combi.split("_");
//      for (int i=0;i!=fas.length;i++){
//        if (positions.containsKey(i+1)){
//          System.out.println("  2.  "+(i+1)+"="+fas[i]);
//          Vector<IntensityRuleVO> rules = positions.get(i+1);
//          for (int j=0;j!=rules.size();j++){
//            System.out.println("   3. "+rules.get(j));
//          }
//        }
//      }
//    }
  }
  
  /**
   * 
   * @param positions list of possible positions
   * @return String identifier for possible position combinations - String is separated by ","
   */
  private String getIdentifierForPositions(List<Integer> positions){
    Collections.sort(positions);
    String id = "";
    for (Integer position: positions) id += position.toString()+",";
    if (id.length()>0) id = id.substring(0,id.length()-1);
    return id;
  }
  
  /**
   * verifies if the position rules are fulfilled for one specific fatty acid combination
   * decides which position assignment is the most likely one
   * all mandatory rules must be fulfilled (if the affected fragments are present)
   * if rules are not mandatory, the position for the combinations where more rules are fulfilled
   * if there is equality in contradicting assignments - no position is assigned
   * if there is no assignment possible - an empty array is returned
   * @param combiName id for the chain combination
   * @param chainOccurenceInCombi how often occurs the same chain in this possible combination
   * @param positionSpecificRules the rules for position assignment
   * @return lookup in int[] - the sequence of the occurence of fatty acids in the combination defines the position in the array, the entry is the assigned position
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   */
  private int[] detectPostions(String combiName, Hashtable<String,Integer> chainOccurenceInCombi, Hashtable<String,Hashtable<Boolean,Vector<IntensityRuleVO>>> positionSpecificRules,
      Hashtable<String,Hashtable<Boolean,Vector<IntensityRuleVO>>> positionSpecificSingleRules) throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException {
    int[] positions = null;
    List<String> faList = new ArrayList<String>(chainOccurenceInCombi.keySet());
    // the first String is the FA, the second Integer is the position, the Vector is the amount of rules that is fulfilled
    Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> positionRecommendations = new Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>();
    Hashtable<String,Hashtable<Integer,IntensityPositionVO>> impossiblePositions = new Hashtable<String,Hashtable<Integer,IntensityPositionVO>>();
    // check first for same chain position identifications
    for (int i=0;i!=faList.size(); i++){
      String fa = faList.get(i);
      Integer occurrence = chainOccurenceInCombi.get(fa);
      Hashtable<String,CgProbe> fragments = chainFragments_.get(fa);
      for (String position : positionSpecificSingleRules.keySet()){
        Hashtable<Boolean,Vector<IntensityRuleVO>> mandAndAddRules = positionSpecificSingleRules.get(position);
        Vector<IntensityRuleVO> mandRules = mandAndAddRules.get(true);
        Vector<IntensityRuleVO> addRules = mandAndAddRules.get(false);
        int proposedPosition = -1;
        Vector<IntensityPositionVO> fulfilledRules = new Vector<IntensityPositionVO>();
        int notEnoughFragmentsForEvaluation = 0;
        for (IntensityRuleVO rule : mandRules){
          if (rule.enoughFragmentsForRuleEvaluationFound(fragments)){
            int[] poss = getPositionRecommendation(fa,occurrence,fragments,fa,occurrence,fragments,rule);
            if (poss==null || (proposedPosition>-1 && poss[0]!=proposedPosition)){
              if (debug_) debugVO_.addUnfulfilledPositionRule(combiName, rule);
              proposedPosition = -1;
              fulfilledRules = new Vector<IntensityPositionVO>();
              Hashtable<Integer,IntensityPositionVO> imp = new Hashtable<Integer,IntensityPositionVO>();
              if (impossiblePositions.containsKey(fa)) imp = impossiblePositions.get(fa);
              imp.put(rule.getBiggerPosition(), createIntensityPositionVO(rule,new int[]{rule.getBiggerPosition(),rule.getBiggerPosition()},fa,fa));
              impossiblePositions.put(fa, imp);
              break;
            }else{
              proposedPosition = poss[0];
              fulfilledRules.add(createIntensityPositionVO(rule,poss,fa,fa));
            }               
          } else notEnoughFragmentsForEvaluation++;
        }
        // if there are mandatory rules, there must be a proposed position - otherwise no detection possible
        if (mandRules.size()>0 && notEnoughFragmentsForEvaluation<mandRules.size() && proposedPosition==-1) continue;
        int fulfilled = 0;
        int contradicting = 0;
        Vector<IntensityPositionVO> addFulfilled = new Vector<IntensityPositionVO>();
        for (IntensityRuleVO rule : addRules){
          if (rule.enoughFragmentsForRuleEvaluationFound(fragments)){
            int[] poss = getPositionRecommendation(fa,occurrence,fragments,fa,occurrence,fragments,rule);
            if (poss!=null){
              // if there is a mandatory rule, the recommendation of the additional rule must fulfil it - otherwise the add rule is wrong
              if (mandRules.size()>0 && notEnoughFragmentsForEvaluation<mandRules.size()){              
                if (proposedPosition==poss[0]) fulfilledRules.add(createIntensityPositionVO(rule,poss,fa,fa));
                else if (debug_) debugVO_.addContradictingPositionRules(combiName, fulfilledRules.lastElement(),createIntensityPositionVO(rule,poss,fa,fa));
              // if there is no mandatory rule - the fulfilled add rules are stored position specific -
              // the positions with more fulfilled rules are said to be correct - if there is equality - no position specific
              // prediction is possible
              } else {
                addFulfilled.add(createIntensityPositionVO(rule,poss,fa,fa));
                fulfilled++;
              }
            }else { 
              contradicting++;
              if (debug_) debugVO_.addUnfulfilledPositionRule(combiName, new IntensityPositionVO(rule,fa,fa));
            }
          }
        }
        if ((mandRules.size()==0 || mandRules.size()==notEnoughFragmentsForEvaluation) && fulfilled>contradicting){
          proposedPosition = new Integer(position);
          fulfilledRules.addAll(addFulfilled);
        }
        if (proposedPosition>-1){
          Hashtable<Integer,Vector<IntensityPositionVO>> positionRuleEvidence = new Hashtable<Integer,Vector<IntensityPositionVO>>();
          if (positionRecommendations.containsKey(fa)) positionRuleEvidence = positionRecommendations.get(fa);
          Vector<IntensityPositionVO> rules = new Vector<IntensityPositionVO>();
          if (positionRuleEvidence.containsKey(proposedPosition)) rules = positionRuleEvidence.get(proposedPosition);
          rules.addAll(fulfilledRules);
          positionRuleEvidence.put(proposedPosition, rules);
          positionRecommendations.put(fa, positionRuleEvidence);
        }
      }
    }
//    for (String fa: impossiblePositions.keySet()){
//      System.out.println("222222222222222222 "+fa+": "+impossiblePositions.get(fa).keySet());
//    }
    for (int i=0;i!=faList.size(); i++){
      String firstFA = faList.get(i);
      Integer firstOccurrence = chainOccurenceInCombi.get(firstFA);
      Hashtable<String,CgProbe> firstFragments = chainFragments_.get(firstFA);
      for (int j=i+1; j!=faList.size(); j++){
        String secondFA = faList.get(j);
        Integer secondOccurrence = chainOccurenceInCombi.get(secondFA);
        Hashtable<String,CgProbe> secondFragments = chainFragments_.get(secondFA);
        for (String positionId : positionSpecificRules.keySet()){
          Hashtable<Boolean,Vector<IntensityRuleVO>> mandAndAddRules = positionSpecificRules.get(positionId);
          Vector<IntensityRuleVO> mandRules = mandAndAddRules.get(true);
          Vector<IntensityRuleVO> addRules = mandAndAddRules.get(false);
          int[] proposedPositions = null;
          Vector<IntensityPositionVO> fulfilledRules = new Vector<IntensityPositionVO>();
          int notEnoughFragmentsForEvaluation = 0;
          Vector<IntensityRuleVO> inversibleRules = new Vector<IntensityRuleVO>();
          for (IntensityRuleVO rule : mandRules){
            if (rule.enoughFragmentsForRuleEvaluationFound(firstFragments) && rule.enoughFragmentsForRuleEvaluationFound(secondFragments)){
              int[] poss = getPositionRecommendation(firstFA,firstOccurrence,firstFragments,secondFA,secondOccurrence,secondFragments,rule);
              int[] poss2 = getPositionRecommendation(secondFA,secondOccurrence,secondFragments,firstFA,firstOccurrence,firstFragments,rule);
              boolean isMandRuleOK = true;
              if (poss!=null){
                if (poss[0]==poss2[0] && poss[1]==poss2[1]){
                  boolean firstOK = true;
                  if (impossiblePositions.containsKey(firstFA) && impossiblePositions.get(firstFA).containsKey(poss[0])||
                      impossiblePositions.containsKey(secondFA) && impossiblePositions.get(secondFA).containsKey(poss[1]))
                    firstOK = false;
                  boolean secondOK = true;
                  if (impossiblePositions.containsKey(firstFA) && impossiblePositions.get(firstFA).containsKey(poss[1])||
                      impossiblePositions.containsKey(secondFA) && impossiblePositions.get(secondFA).containsKey(poss[0]))
                    secondOK = false;
                  if (firstOK && secondOK){
                    inversibleRules.add(rule);
                    continue;
                  } else if (firstOK){
                  } else if (secondOK){
                    poss = poss2;
                  } else isMandRuleOK = false;
                }
                if (impossiblePositions.containsKey(firstFA) && impossiblePositions.get(firstFA).containsKey(poss[0])||
                    impossiblePositions.containsKey(secondFA) && impossiblePositions.get(secondFA).containsKey(poss[1]))
                  isMandRuleOK = false;
                if (proposedPositions!=null){
                  for (int k=0;k!=poss.length;k++){
                    if (poss[k] != proposedPositions[k]) isMandRuleOK = false;
                  }
                  if (!isMandRuleOK && debug_) debugVO_.addContradictingPositionRules(combiName, fulfilledRules.lastElement(),createIntensityPositionVO(rule,poss,firstFA,secondFA));
                }
              }else{
                isMandRuleOK = false;
                if (debug_) debugVO_.addUnfulfilledPositionRule(combiName, rule);
              }
              if (isMandRuleOK){
                proposedPositions = poss;
                fulfilledRules.add(createIntensityPositionVO(rule,poss,firstFA,secondFA));
              }else{
                proposedPositions = null;
                fulfilledRules = new Vector<IntensityPositionVO>();
                break;
              } 
            } else notEnoughFragmentsForEvaluation++;
          }
          // if there are inversible rules, add them to the evidence
          if (fulfilledRules.size()>0 && inversibleRules.size()>0){
            for (IntensityRuleVO rule : inversibleRules){
              fulfilledRules.add(createIntensityPositionVO(rule,proposedPositions,firstFA,secondFA));
            }
            inversibleRules = new Vector<IntensityRuleVO>();
          }
          // if there are mandatory rules, there must be a proposed position - otherwise no detection possible
          if (mandRules.size()>0 && notEnoughFragmentsForEvaluation<mandRules.size() && proposedPositions==null) continue;
          Hashtable<String,Vector<IntensityPositionVO>> addRulesHash = new Hashtable<String,Vector<IntensityPositionVO>>();
          for (IntensityRuleVO rule : addRules){
            if (rule.enoughFragmentsForRuleEvaluationFound(firstFragments) && rule.enoughFragmentsForRuleEvaluationFound(secondFragments)){
              int[] poss = getPositionRecommendation(firstFA,firstOccurrence,firstFragments,secondFA,secondOccurrence,secondFragments,rule);
              int[] poss2 = getPositionRecommendation(secondFA,secondOccurrence,secondFragments,firstFA,firstOccurrence,firstFragments,rule);
              if (poss!=null){
                boolean discard = false;
                if (poss[0]==poss2[0] && poss[1]==poss2[1]){
                  boolean firstOK = true;
                  if (impossiblePositions.containsKey(firstFA) && impossiblePositions.get(firstFA).containsKey(poss[0])||
                      impossiblePositions.containsKey(secondFA) && impossiblePositions.get(secondFA).containsKey(poss[1]))
                    firstOK = false;
                  boolean secondOK = true;
                  if (impossiblePositions.containsKey(firstFA) && impossiblePositions.get(firstFA).containsKey(poss[1])||
                      impossiblePositions.containsKey(secondFA) && impossiblePositions.get(secondFA).containsKey(poss[0]))
                    secondOK = false;
                  if (firstOK && secondOK){
                    inversibleRules.add(rule);
                    continue;
                  } else if (firstOK){
                  } else if (secondOK){
                    poss = poss2;
                  } else discard = true;
                }
                if (impossiblePositions.containsKey(firstFA) && impossiblePositions.get(firstFA).containsKey(poss[0])||
                    impossiblePositions.containsKey(secondFA) && impossiblePositions.get(secondFA).containsKey(poss[1]))
                  discard = true;
                if (discard){
                  if (debug_) debugVO_.addUnfulfilledPositionRule(combiName, new IntensityPositionVO(rule,firstFA,secondFA));
                  continue;
                }
                // if there is a mandatory rule, the recommendation of the additional rule must fulfil it - otherwise the add rule is wrong
                if (mandRules.size()>0 && notEnoughFragmentsForEvaluation<mandRules.size()){
                
                  boolean theSame = true;
                  for (int k=0;k!=poss.length;k++){
                    if (poss[k] != proposedPositions[k]) theSame = false;
                  }
                  if (theSame) fulfilledRules.add(createIntensityPositionVO(rule,poss,firstFA,secondFA));
                  else if (debug_) debugVO_.addContradictingPositionRules(combiName, fulfilledRules.lastElement(),createIntensityPositionVO(rule,poss,firstFA,secondFA));
                // if there is no mandatory rule - the fulfilled add rules are stored position specific -
                // the positions with more fulfilled rules are said to be correct - if there is equality - no position specific
                // prediction is possible
                } else {
                  String id = String.valueOf(poss[0])+","+String.valueOf(poss[1]);
                  Vector<IntensityPositionVO> rules = new Vector<IntensityPositionVO>();
                  if (addRulesHash.containsKey(id)) rules = addRulesHash.get(id);
                  rules.add(createIntensityPositionVO(rule,poss,firstFA,secondFA));
                  addRulesHash.put(id, rules);
                }
              }else if (debug_){
                debugVO_.addUnfulfilledPositionRule(combiName, new IntensityPositionVO(rule,firstFA,secondFA));
              }
            }
          }
          // if there is no mandatory rule - the fulfilled add rules are stored position specific -
          // the positions with more fulfilled rules are said to be correct - if there is equality - no position specific
          // prediction is possible
          if ((mandRules.size()==0 || mandRules.size()==notEnoughFragmentsForEvaluation) && addRulesHash.size()>0){
            int maxFoundRules = 0;
            for (String id : addRulesHash.keySet()){
              Vector<IntensityPositionVO> rules = addRulesHash.get(id);
              // there is equality -> no decision
              if (rules.size() == maxFoundRules){
                proposedPositions = null;
                fulfilledRules = new Vector<IntensityPositionVO>();
              } else if (rules.size()>maxFoundRules){
                maxFoundRules = rules.size();
                fulfilledRules = rules;
                proposedPositions = new int[2];
                String[] poss = id.split(",");
                proposedPositions[0] = Integer.parseInt(poss[0]);
                proposedPositions[1] = Integer.parseInt(poss[1]);
              }
            }
            if (debug_ && addRulesHash.size()>1){
              Iterator<String> keys = addRulesHash.keySet().iterator();
              String key1 = keys.next();
              String key2 = keys.next();
              IntensityPositionVO rule1 = addRulesHash.get(key1).lastElement();
              IntensityPositionVO rule2 = addRulesHash.get(key2).lastElement();
              debugVO_.addContradictingPositionRules(combiName, rule1, rule2);
            }
          }
          if (proposedPositions!=null){
            for (IntensityRuleVO rule : inversibleRules) fulfilledRules.add(createIntensityPositionVO(rule,proposedPositions,firstFA,secondFA));
            Hashtable<Integer,Vector<IntensityPositionVO>> positionRuleEvidence = new Hashtable<Integer,Vector<IntensityPositionVO>>();
            if (positionRecommendations.containsKey(firstFA)) positionRuleEvidence = positionRecommendations.get(firstFA);
            Vector<IntensityPositionVO> rules = new Vector<IntensityPositionVO>();
            if (positionRuleEvidence.containsKey(proposedPositions[0])) rules = positionRuleEvidence.get(proposedPositions[0]);
            rules.addAll(fulfilledRules);
            positionRuleEvidence.put(proposedPositions[0], rules);
            positionRecommendations.put(firstFA, positionRuleEvidence);
            positionRuleEvidence = new Hashtable<Integer,Vector<IntensityPositionVO>>();
            if (positionRecommendations.containsKey(secondFA)) positionRuleEvidence = positionRecommendations.get(secondFA);
            rules = new Vector<IntensityPositionVO>();
            if (positionRuleEvidence.containsKey(proposedPositions[1])) rules = positionRuleEvidence.get(proposedPositions[1]);
            rules.addAll(fulfilledRules);
            positionRuleEvidence.put(proposedPositions[1], rules);
            positionRecommendations.put(secondFA, positionRuleEvidence);

          }

        }
      }
    }
    if (positionRecommendations.size()>0){
      for (int i=0;i!=faList.size(); i++){
        String fa = faList.get(i);
        if (!positionRecommendations.containsKey(fa)) continue;
//        Hashtable<Integer,Vector<IntensityPositionVO>> recomms = positionRecommendations.get(fa);
//        System.out.println(fa+": "+recomms);
      }  
      
      positions = getMostLikelyPositionsFromRecommendations(combiName, positionRecommendations, impossiblePositions, fragCalc_.getAllowedChainPositions());
      positionRecommendations_.put(combiName, positionRecommendations);
    }
//    String[] fas = LipidomicsMSnSet.getFAsFromCombiName(combiName);
//    for (int i=0;i!=fas.length; i++){
//      String fa = fas[i];
//      System.out.println(fa+": "+positions[i]);
//    }  
    return positions;
  }
  
  /**
   * the rules return recommendations for each position
   * this algorithm defines the most likely position from the rules
   * if a definite recommendation for one position cannot be made, -1 is entered at this entry of the int array
   * the algorithm removes iteratively the unique assignments
   * then it assigns the positions to the most fatty acids where more rules match
   * if there is no "winner" -1 is assigned
   * @param combiName String containing the fatty acids separated by "_" - the returned int array with the positions corresponds to this fatty acid sequence
   * @param positionRecommendations the found rules for the assignment - first key: fatty acid name, second key: possible position, thrid key: rules corroborating this assignment
   * @param impossiblePositions positions that were excluded by tests on fragments of the same FA
   * @param the amount of positions that are possible
   * @return the assigned positions for the sequence in the fatty acid string array
   */
  private int[] getMostLikelyPositionsFromRecommendations(String combiName, Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> positionRecommendations, Hashtable<String,Hashtable<Integer,IntensityPositionVO>> impossiblePositions,
      int allowedChainPositions){
    String[] fas = LipidomicsMSnSet.getFAsFromCombiName(combiName);
    Hashtable<Integer,Vector<IntensityPositionVO>> positionEvidence = new Hashtable<Integer,Vector<IntensityPositionVO>>();
    int[] positions = new int[fas.length];
    ////Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> posRec = new Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>(positionRecommendations);
    Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> posRec = new Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>();
    for (int i=0; i!= fas.length; i++){
      String fa = fas[i];
      if (positionRecommendations.containsKey(fa))
        posRec.put(fa+"_"+String.valueOf(i), new Hashtable<Integer,Vector<IntensityPositionVO>>(positionRecommendations.get(fa)));
    }
    Hashtable<Integer,Integer> unassignedPositions = new Hashtable<Integer,Integer>();
    for (int i=0;i!=fas.length;i++) positions[i] = -1;
    for (int i=0;i!=allowedChainPositions;i++){
      int position = i+1;
      unassignedPositions.put(position, position);
    }
    while (unassignedPositions.size()!=allowedChainPositions-fas.length){
      Hashtable<Integer,String> singles = getSinglePositionAssignments(posRec);
      // single identifications are thought to be valid
      // they are assigned and removed from the remaining positions
      if (singles.size()>0){
        for (Integer position : singles.keySet()){
          String fa = singles.get(position);
          for (int i=0;i!=fas.length;i++){
            if (fas[i].equalsIgnoreCase(fa)&&positions[i]==-1){
              positions[i] = position;
              positionEvidence.put(position, new Vector<IntensityPositionVO>(positionRecommendations.get(fa).get(position)));
              break;
            }
          }
          String keyToRemove = "";
          for (String key : posRec.keySet()){
            if (key.startsWith(fa+"_")){
              keyToRemove = key;
              break;
            }
          }
          posRec.remove(keyToRemove);
////          if (posRec.get(fa).size()<2)
////            posRec.remove(fa);
          for (String otherFA : posRec.keySet()){
            posRec.get(otherFA).remove(position);
          }
          removeWrongEvidence(fa,position,posRec);
          unassignedPositions.remove(position);
        }
      // if there are no singles -> more than one FA can have the position -> find one with more evidence if not -> stop
      }else{
        Hashtable<Integer,String> multiples = getMultiplePositionsWithHigherEvidence(posRec);
        // do the same as for singles
        if (multiples.size()>0){
          for (Integer position : multiples.keySet()){
            String fa = multiples.get(position);
            for (int i=0;i!=fas.length;i++){
              if (fas[i].equalsIgnoreCase(fa)){
                positions[i] = position;
                String key = fa+"_"+String.valueOf(i);
                positionEvidence.put(position, new Vector<IntensityPositionVO>(posRec.get(key).get(position)));
                posRec.remove(key);
                break;
              }
            }
            for (String otherFA : posRec.keySet()){
              posRec.get(otherFA).remove(position);
            }
            unassignedPositions.remove(position);
            removeWrongEvidence(fa,position,posRec);
          }
        // stop the while - further positions cannot be assigned
        }else{
          break;
        }
      }
    }
    // this if is to check if there is only one position possible for the assignment,
    // since the others were removed by tests of chains of the same fragment
    if (unassignedPositions.size()==(allowedChainPositions-fas.length+1)){
      int faPos = -1;
      for (int i=0; i!= fas.length; i++){
        if (positions[i]==-1){
          faPos = i;
          break;
        }
      }
      String fa = fas[faPos];
      for (int i=0;i!=fas.length;i++){
        if (positions[i]>0) unassignedPositions.remove(positions[i]);
      }
      if (impossiblePositions.containsKey(fa)){
        for (Integer pos : impossiblePositions.get(fa).keySet()){
          if (unassignedPositions.containsKey(pos)) unassignedPositions.remove(pos);
        }
      }
      if (unassignedPositions.size()==1){
        positions[faPos] = unassignedPositions.keySet().iterator().next();
        String key = fa+"_"+String.valueOf(faPos);
        Vector<IntensityPositionVO> evidence = new Vector<IntensityPositionVO>();
        if (impossiblePositions.containsKey(fa)){
          for (IntensityPositionVO posVO : impossiblePositions.get(fa).values()){
            evidence.add(IntensityPositionVO.createNegatedVO(posVO,positions[faPos]));
          }
        }
        positionEvidence.put(positions[faPos], evidence);
        posRec.remove(key);
        unassignedPositions.remove(positions[faPos]);
        removeWrongEvidence(fa,positions[faPos],posRec);
      }
    }
    positionEvidence = cleanWrongEvidence(fas, positions, positionEvidence);
    posEvidenceAccordToProposedDefinition_.put(combiName, positionEvidence);
    return positions;
  }
  
  /**
   * removes evidence that is based on a position, that was assigned later to something else
   * @param fas the fatty acids in the same sequence as in the combi name
   * @param positions the assigned positions
   * @param positionEvidence the current evidence
   * @return the cleared position evidence hash
   */
  private Hashtable<Integer,Vector<IntensityPositionVO>> cleanWrongEvidence(String[] fas, int[] positions, Hashtable<Integer,Vector<IntensityPositionVO>> positionEvidence){
    Hashtable<Integer,Vector<IntensityPositionVO>> cleaned = new Hashtable<Integer,Vector<IntensityPositionVO>>();
    Hashtable<String,Hashtable<Integer,Integer>> faPosOK = new Hashtable<String,Hashtable<Integer,Integer>>();
    Hashtable<String,String> unassigned = new Hashtable<String,String>();
    for (int i=0; i!=fas.length; i++){
      String fa = fas[i];
      if (positions[i]>-1){
        Hashtable<Integer,Integer> poss = new Hashtable<Integer,Integer>();
        if (faPosOK.containsKey(fa)) poss = faPosOK.get(fa);
        poss.put(positions[i], positions[i]);
        faPosOK.put(fa, poss);
      } else unassigned.put(fa, fa);
    }
    for (Integer pos : positionEvidence.keySet()){
      Vector<IntensityPositionVO> ok = new Vector<IntensityPositionVO>();
      for (IntensityPositionVO posVO : positionEvidence.get(pos)){
        String biggerFA = posVO.getBiggerFA();
        int biggerPos = posVO.getBiggerPosition();
        String smallerFA = posVO.getSmallerFA();
        int smallerPos = posVO.getSmallerPosition();
        if (((unassigned.containsKey(biggerFA)||(faPosOK.containsKey(biggerFA)&&faPosOK.get(biggerFA).containsKey(biggerPos)))&&
            (unassigned.containsKey(smallerFA)||(faPosOK.containsKey(smallerFA)&&faPosOK.get(smallerFA).containsKey(smallerPos))))||
            posVO.isNegated()){
          ok.add(posVO);
        }
      }
      cleaned.put(pos, ok);
    }
    return cleaned;
  }
  
  private void removeWrongEvidence(String fa, int position, Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> posRec){
    for (String faWithPos :  posRec.keySet()){
      Hashtable<Integer,Vector<IntensityPositionVO>> recomm = posRec.get(faWithPos);
      for (Integer pos : recomm.keySet()){
        Vector<IntensityPositionVO> old = recomm.get(pos);
        Vector<IntensityPositionVO> newEv = new Vector<IntensityPositionVO>();
        for (IntensityPositionVO ev : old){
          if (ev.getBiggerPosition()!=position && ev.getSmallerPosition()!=position)
            newEv.add(ev);
          else {
            if ((ev.getBiggerPosition()==position && ev.getBiggerFA().equalsIgnoreCase(fa)) || (ev.getSmallerPosition()==position && ev.getSmallerFA().equalsIgnoreCase(fa)))
              newEv.add(ev);
          }
        }
        if (newEv.size()>0) recomm.put(pos, newEv);
      }
      posRec.put(faWithPos, recomm);
    }
  }
  
  /**
   * assigns fatty acids to positions where evidence for the same position is returned from more than one fatty acid
   * @param posRec the found rules for the assignment - first key: fatty acid name, second key: possible position, thrid key: rules corroborating this assignment
   * @return proposed assignments - key: assigned position, value: assigned fatty acid - if nothing is assigned -> the hash is empty
   */
  private Hashtable<Integer,String> getMultiplePositionsWithHigherEvidence(Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> posRec){
    Hashtable<Integer,String> moreLikelyMultiples = new Hashtable<Integer,String>();
    Hashtable<Integer,Integer> occurences = new Hashtable<Integer,Integer>();
    Hashtable<Integer,Vector<String>> fasForPosition = buildFasForPositionHash(posRec);
    int lowestOccurence = Integer.MAX_VALUE;

    for (Vector<String> fas : fasForPosition.values()){
      if (fas.size()<lowestOccurence) lowestOccurence = fas.size();
    }
    int rulesForPosition = 0;
    Hashtable<String,Integer> allowedOccurences = getAllowedOccurences(posRec);
    for (Integer position: fasForPosition.keySet()){
      Vector<String> fas = fasForPosition.get(position);
      Hashtable<String,String> usedFAs = new Hashtable<String,String>();
      for (String fa : fas){
        if (usedFAs.containsKey(fa)) continue;
        usedFAs.put(fa, fa);
        String key = "";
        for (String faWithPos : posRec.keySet()){
          if (faWithPos.startsWith(fa+"_")){
            key = faWithPos;
            break;
          }
        }
        Vector<IntensityPositionVO> rules = posRec.get(key).get(position);
        // there is equality -> remove
        if (rules.size()==rulesForPosition){
          moreLikelyMultiples.remove(position);
          occurences.remove(position);
        }else if (rules.size()>rulesForPosition){
          int foundElse = 0;
          for (Integer pos2 : fasForPosition.keySet()){
            if (position==pos2 || !posRec.get(key).containsKey(pos2) || rules.size()!=posRec.get(key).get(pos2).size()) continue;
            foundElse++;
          }
          if (foundElse<allowedOccurences.get(fa)){
            rulesForPosition = rules.size();
            moreLikelyMultiples.put(position, fa);
            occurences.put(position, rulesForPosition);
          }
        }
      }
    }
    int highestOccurence = 0;
    for (Integer pos : moreLikelyMultiples.keySet()){
      if (occurences.get(pos)>highestOccurence) highestOccurence = occurences.get(pos);
    }
    for (Integer pos : occurences.keySet()){
      if (occurences.get(pos)!=highestOccurence) moreLikelyMultiples.remove(pos);
    }
    return moreLikelyMultiples;
  }
  
  /**
   * assigns fatty acids to positions where evidence for the same position is returned by only one fatty acid
   * @param posRec the found rules for the assignment - first key: fatty acid name, second key: possible position, thrid: rules corroborating this assignment
   * @return proposed assignments - key: assigned position, value: assigned fatty acid - if nothing is assigned -> the hash is empty
   */
  private Hashtable<Integer,String> getSinglePositionAssignments(Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> posRec){
    Hashtable<Integer,String> singles = new Hashtable<Integer,String>();
    Hashtable<Integer,Vector<String>> fasForPosition = buildFasForPositionHash(posRec);
    Hashtable<String,Integer> allowedOccurences = getAllowedOccurences(posRec);
    for (Integer position : fasForPosition.keySet()){
      if (fasForPosition.get(position).size()==1){
        String fa = fasForPosition.get(position).get(0);
        int foundElse = 0;
        for (Integer pos2 : fasForPosition.keySet()){
          if (position==pos2 || fasForPosition.get(position).size()>1) continue;
          for (String other : fasForPosition.get(pos2)){
            if (fa.equalsIgnoreCase(other)){              
              foundElse++;
              break;
            }
          }
        }
        if (foundElse<allowedOccurences.get(fa))
          singles.put(position, fa);
      }
    }
    return singles;
  }
  
  private Hashtable<String,Integer> getAllowedOccurences (Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> posRec){
    Hashtable<String,Integer> allowedOccurences = new Hashtable<String,Integer>();
    for (String faWithPos : posRec.keySet()){
      String fa = faWithPos.substring(0,faWithPos.lastIndexOf("_"));
      int occur = 0;
      if (allowedOccurences.containsKey(fa)) occur = allowedOccurences.get(fa);
      occur++;
      allowedOccurences.put(fa, occur);
    }
    return allowedOccurences;
  }
  
  /**
   * builds a hash table: key: the position; value: the fatty acids that are possible for this position
   * @param posRec the found rules for the assignment - first key: fatty acid name, second key: possible position, thrid: rules corroborating this assignment
   * @return hash table: key: the position; value: the fatty acids that are possible for this position
   */
  private Hashtable<Integer,Vector<String>> buildFasForPositionHash(Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> posRec){
    Hashtable<Integer,Vector<String>> fasForPosition = new Hashtable<Integer,Vector<String>>();
    for (String faWithPos : posRec.keySet()){
      String fa = faWithPos.substring(0,faWithPos.lastIndexOf("_"));
      Hashtable<Integer,Vector<IntensityPositionVO>> intToRules = posRec.get(faWithPos);
      for (Integer position: intToRules.keySet()){
        Vector<String> fas = new Vector<String>();
        if (fasForPosition.containsKey(position)) fas = fasForPosition.get(position);
        fas.add(fa);
        fasForPosition.put(position, fas);
      }
    }
    return fasForPosition;
  }
  
  
  /**
   * verifies if one position rule for two affected fatty acids is fulfilled
   * @param firstFA first fatty acid of the two
   * @param firstOccurence how often occurs this fatty acid on the proposed combination
   * @param firstFragments the found fragments of the first fatty acid
   * @param secondFA second fatty acid of the two
   * @param secondOccurence how often occurs this fatty acid on the proposed combination
   * @param secondFragments the found fragments of the first fatty acid
   * @param rule intensity rule to be verified
   * @return the position recommendations according to this rule - if no recommendation is made, an empty int array is returned
   * @throws SpectrummillParserException 
   * @throws IOException 
   * @throws NoRuleException 
   * @throws RulesException 
   */
  private int[] getPositionRecommendation(String firstFA, int firstOccurence, Hashtable<String,CgProbe> firstFragments, String secondFA, int secondOccurence, Hashtable<String,CgProbe> secondFragments, IntensityRuleVO rule) throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    int[] recommendation = null;
    // check if the firstFA is the biggerArea
    Float basepeak = getBasepeakIfRequired(rule);
    double biggerArea = rule.getBiggerArea(firstFragments, basepeak)/((float)firstOccurence);
    double smallerArea = rule.getSmallerArea(secondFragments, basepeak)/((float)secondOccurence);
    if (biggerArea > smallerArea){
      recommendation = new int[2];
      recommendation[0] = rule.getBiggerPosition();
      recommendation[1] = rule.getSmallerPosition();     
    }
    if (recommendation!=null) return recommendation;
    biggerArea = rule.getBiggerArea(secondFragments, basepeak)/((float)secondOccurence);
    smallerArea = rule.getSmallerArea(firstFragments, basepeak)/((float)firstOccurence);
    if (biggerArea > smallerArea){
      recommendation = new int[2];
      recommendation[1] = rule.getBiggerPosition();
      recommendation[0] = rule.getSmallerPosition();     
    }   
    return recommendation;
  }  
  
  
  private void transferResultsToLipidParameterSet() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    if (status_ > LipidomicsMSnSet.DISCARD_HIT || debug_){
      Hashtable<Integer,Vector<Float>> msnRetentionTimes = new Hashtable<Integer,Vector<Float>>();
      for (int msLevel : msLevels_.keySet()){
        if (!msLevels_.get(msLevel) || !this.probesWithMSnSpectra_.containsKey(msLevel)) continue;
        msnRetentionTimes.put(msLevel, analyzer_.getMSnSpectraRetentionTimes(msLevel,probesWithMSnSpectra_.get(msLevel)));
      }      
      set_ = new LipidomicsMSnSet(set_,status_,LipidomicsConstants.getMs2MzTolerance(),headGroupFragments_,fulfilledHeadIntensityRules_,
          chainFragments_, fulfilledChainIntensityRules_, validChainCombinations_,positionDefinition_, posEvidenceAccordToProposedDefinition_,
          fragCalc_!=null ? fragCalc_.getAllowedChainPositions() : 1, basePeakValues_, msnRetentionTimes);
    }
  }
  
  public LipidParameterSet getResult(){
    return this.set_;
  }
  
  
  /**
   * creates an IntensityPositionVO from the IntensityRuleVO, the assigned positions, and the name of the two fatty acids
   * @param rule the rule that is fulfilled
   * @param poss the assigned positions - pos[0] is the position of the firstFA
   * @param firstFA the first assigned fatty acid
   * @param secondFA the second fatty acid
   * @return created IntensityPositionVO
   */
  private IntensityPositionVO createIntensityPositionVO(IntensityRuleVO rule, int[] poss, String firstFA, String secondFA){
    String biggerFA = null;
    String smallerFA = null;
    if (poss[0]==rule.getBiggerPosition()){
      biggerFA = firstFA;
      smallerFA = secondFA;
    }else{
      biggerFA = secondFA;
      smallerFA = firstFA;
    }
    return new IntensityPositionVO(rule,biggerFA,smallerFA);
  }
    
  /**
   * 
   * @param msLevels levels where the coverage has to be checked
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   */
  private void checkSpectrumCoverage(Hashtable<Integer,Boolean> msLevels) throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException{
    Vector<CgProbe> foundHits = new Vector<CgProbe>();
    Vector<CgProbe> notCounted = new Vector<CgProbe>();
    for (CgProbe hit : this.headGroupFragments_.values()){
      if (hit.isFromOtherSpecies()) notCounted.add(hit);
      else foundHits.add(hit);
    }
    for (Hashtable<String,CgProbe> frags : this.chainFragments_.values()){
      for (CgProbe hit : frags.values()){
        if (hit.isFromOtherSpecies()) notCounted.add(hit);
        else foundHits.add(hit);
      }
    }
    if (!isSpectrumCovered(this.analyzer_,this.set_,msLevels,foundHits,new Vector<CgProbe>(),(float)RulesContainer.getBasePeakCutoff(StaticUtils.getRuleName(this.className_, this.modName_),this.rulesDir_),
        (float)RulesContainer.getSpectrumCoverageMin(StaticUtils.getRuleName(this.className_, this.modName_),this.rulesDir_),this.ignoreAbsolute_,this.debug_,this.debugVO_))
      this.status_ = LipidomicsMSnSet.DISCARD_HIT;
  }


  /**
   * checks if the spectrum is sufficiently covered by found peaks
   * @param analyzer object that holds MS data and can quantify fragments of interest
   * @param set MS1 data including quantitation info for the lipid to be checked
   * @param msLevels levels where the coverage has to be checked
   * @param foundHits identified m/z values
   * @param notCountedHits hits that originate from another species; these hits should be excluded from the spectrum coverage
   * @param bpCutoff the base peak cutoff that shall be used for the spectrum coverage calculation
   * @param coverageMin the spectrum coverage that must be fulfilled
   * @param ignoreAbsolute should absolute rules be ignored - whenever there is an overlap of isobars assumed
   * @param debug should debug messages be written
   * @param debugVO provides information why fragments were discarded - only available in debug mode
   * @return true if the spectrum is sufficiently covered, false otherwise
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rule does not exist
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * throws CgException errors from the quantitation process
   */
  private static boolean isSpectrumCovered(LipidomicsAnalyzer analyzer, LipidParameterSet set, Hashtable<Integer,Boolean> msLevels, Vector<CgProbe> foundHits, Vector<CgProbe> notCountedHits, float bpCutoff, float coverageMin, boolean ignoreAbsolute, boolean debug, MSnDebugVO debugVO) throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException {
    Vector<Integer> levels = new Vector<Integer>();
    for (int level : msLevels.keySet()){
      //TODO: spectrum coverage is currently only allowed for MS-level: 2
      //to extend this, spectrum coverage values have to be added to be added for each MS-Level in the [GENERAL] section of the fragmentation rules
      if (msLevels.get(level) && level==2) levels.add(level);
    }
    Hashtable<Integer,Float> coverages = analyzer.calculateSpectrumCoverage(set.getIsotopicProbes().get(0), foundHits, notCountedHits, levels,
        bpCutoff);
    boolean covered = true;
    for (int key : coverages.keySet()){
      //TODO: if in any msLevel the hit is lower than the required coverage -> the hit is discarded; 
      //      to implement: discard only fragments of the corresponding level and check if mandatory rules are violated
      //System.out.println(key+": "+coverages.get(key));
      if (ignoreAbsolute) coverageMin = 0f;
//      if (notCountedHits.size()>0) System.out.println(className+"/"+modName+"-"+set.getNameString()+";"+set.getRt()+"    "+coverages.get(key));      
      if (!(coverages.get(key) > coverageMin || (!ignoreAbsolute && coverageMin==0f))){
        if (debug) debugVO.setSpectrumCoverageFulfilled(false);
        covered = false;
        break;
      }
    }
    return covered;
  }
    
  /**
   * 
   * @param analyzer object that holds MS data and can quantify fragments of interest
   * @param set MS1 data including quantitation info for the lipid to be checked
   * @param msLevels levels where the coverage has to be checked
   * @param foundHits identified m/z values
   * @param notCountedHits hits that originate from another species; these hits should be excluded from the spectrum coverage
   * @param bpCutoff the base peak cutoff that shall be used for this calculation
   * @param coverageMin spectrum coverage that need to be fulfilled
   * @return true if the spectrum is sufficiently covered, false otherwise
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rule does not exist
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * throws CgException errors from the quantitation process
   */
  public static boolean isSpectrumCovered(LipidomicsAnalyzer analyzer, LipidParameterSet set, Hashtable<Integer,Boolean> msLevels, Vector<CgProbe> foundHits, Vector<CgProbe> notCountedHits, float bpCutoff, float coverageMin) throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException {
    return isSpectrumCovered(analyzer, set, msLevels, foundHits, notCountedHits, bpCutoff, coverageMin, false, false, null);
  }


  /**
   * returns information why fatty acid chains were discarded - can be called only if "debug" in the constructor is true
   * @return the debg information object
   * @throws RulesException
   */
  public MSnDebugVO getDebugInfo() throws RulesException
  {
    if (debug_){
      debugVO_.setStatus(status_);
      return debugVO_;
    }else throw new RulesException("The analyzer is not in debug mode -> there is no debug information available!");
  }
  
  /**
   * @param probe VO storing information about the peak
   * @param fragment the fragment that should be detected
   * @return the status why the probe was not detected
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   */
  private int findAreaDiscardReason(CgProbe probe, FragmentVO fragment) throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException{
    int discardStatus = MSnDebugVO.NO_PEAK_THERE;
    if (probe.AreaStatus==CgAreaStatus.NothingThere) discardStatus = MSnDebugVO.NO_PEAK_THERE;
    else if (!checkBasePeakCutoff(probe,fragment.getMsLevel())) discardStatus = MSnDebugVO.BELOW_BASE_PEAK_CUTOFF;
    return discardStatus;
  }
  
  /**
   * 
   * @return retention times where MSn spectra matched (only useful if scanAllSpectraForCandidates was called by a certain constructor)
   */
  public Vector<Float> getFoundMatchingSpectraTimes(){
    return spectraFound_;
  }
  
  /**
   * 
   * @return whether there are any MSn spectra found for this precursor ion m/z value
   */
  public boolean areMSnSpectraPresent(){
    return msnSpectraPresent_;
  }
  
  public Hashtable<String,Hashtable<Integer,Integer>> getPositionDefinition()
  {
    return this.positionDefinition_;
  }
  
  /**
   * Returns the rules direction
   * @return
   */
  public String getRulesDir()
  {
    return this.rulesDir_;
  }  
  
  /**
   * adds an IntensityRuleVO to an hash that stores a key for the position comparision, and a key if the rule is mandatory
   * @param positionSpecificRules the hash where the rule shall be added
   * @param identifier the position key
   * @param rule the rule to be added
   */
  private void addToPositionSpecificRuleHash (Hashtable<String,Hashtable<Boolean,Vector<IntensityRuleVO>>> positionSpecificRules, String identifier, IntensityRuleVO rule){
    Hashtable<Boolean,Vector<IntensityRuleVO>> mandAndAddRules = new Hashtable<Boolean,Vector<IntensityRuleVO>>();
    Vector<IntensityRuleVO> mandRules = new Vector<IntensityRuleVO>();
    Vector<IntensityRuleVO> addRules = new Vector<IntensityRuleVO>();
    if (positionSpecificRules.containsKey(identifier)){
      mandAndAddRules = positionSpecificRules.get(identifier);
      mandRules = mandAndAddRules.get(true);
      addRules = mandAndAddRules.get(false);
    }      
    if (rule.isMandatory()) mandRules.add(rule);
    else addRules.add(rule);
    mandAndAddRules.put(true,mandRules);
    mandAndAddRules.put(false,addRules);
    positionSpecificRules.put(identifier, mandAndAddRules);
  }
  
  /**
   * if there are not more than two isobaric species under one MS1 peak, the peaks are split according to their fragments
   * @param analyzer the access to the chrom files
   * @param start start time for checking MSn spectra
   * @param stop stop time for checking MSn spectra
   * @param startRelative start time for the relative check of MSns spectra (only spectra at higher intensity reasons should be taken for this, otherwise false positives are possible)
   * @param stopRelative stop time for the relative check of MSns spectra (only spectra at higher intensity reasons should be taken for this, otherwise false positives are possible)
   * @param shared VO containing information about the species that are in the peak and their distinct fragments
   * @param upperMsLevel the level of the peak to be split
   * @return the split peak intensities if possible
   * @throws CgException 
   */
  public static Vector<SharedPeakContributionVO> splitTwoIsobaricPeaks(LipidomicsAnalyzer analyzer, float start, float stop, float startRelative, float stopRelative, SharedMS1PeakVO shared, int upperMsLevel) throws CgException{
    // this is for the preparation of the data hashes for processing
    Vector<SharedPeakContributionVO> splitted = new Vector<SharedPeakContributionVO>();
    Hashtable<Integer,Hashtable<Integer,String>> relevantSpectra = new Hashtable<Integer,Hashtable<Integer,String>>();
    Hashtable<Integer,Hashtable<Integer,String>> allSpectra = analyzer.getMSnSpectraCache();
    Hashtable<Integer,Hashtable<Integer,Float>> allNoise = analyzer.getMSnSpectraNoise();
    Hashtable<Integer,Hashtable<Integer,Float>> retTimeLookup = new Hashtable<Integer,Hashtable<Integer,Float>>();
    Hashtable<Integer,List<Integer>> scanNumbersSorted = new Hashtable<Integer,List<Integer>>();
    for (Integer msLevel : allSpectra.keySet()){
      Hashtable<Integer,Float> retTimes = analyzer.getRetentionTimes(msLevel);
      Hashtable<Integer,String> spectra = allSpectra.get(msLevel);
      Hashtable<Integer,Float> relRetTimes = new Hashtable<Integer,Float>();
      Hashtable<Integer,String> relSpectra = new Hashtable<Integer,String>();
      List<Integer> scansSorted = new  ArrayList<Integer>();
      for (Integer consScanNumber : spectra.keySet()){
        float rt = retTimes.get(consScanNumber);
        if (start<=rt && rt<=stop){
          relRetTimes.put(consScanNumber, retTimes.get(consScanNumber));
          relSpectra.put(consScanNumber, spectra.get(consScanNumber));
          scansSorted.add(consScanNumber);
        }
      }
      Collections.sort(scansSorted);
      
      retTimeLookup.put(msLevel, relRetTimes);
      relevantSpectra.put(msLevel, relSpectra);
      scanNumbersSorted.put(msLevel, scansSorted);      
    }
    
    // now figure out which m/z values are relevant for the extraction
    Hashtable<Integer,Hashtable<QuantVO,Hashtable<String,CgProbe>>> relevantMzHash =  new Hashtable<Integer,Hashtable<QuantVO,Hashtable<String,CgProbe>>>();
    for (int i=0; i!=shared.getPartners().size();i++){
      SharedPeakContributionVO partner = shared.getPartners().get(i);
      QuantVO quant = partner.getQuantVO();
      for (String key : partner.getDistinctFragments().keySet()){
        CgProbe probe = partner.getDistinctFragments().get(key);
        Hashtable<QuantVO,Hashtable<String,CgProbe>> allOfMsLevel = new Hashtable<QuantVO,Hashtable<String,CgProbe>>();
        if (relevantMzHash.containsKey(probe.getMsLevel())) allOfMsLevel = relevantMzHash.get(probe.getMsLevel());
        Hashtable<String,CgProbe> allOfQuant = new Hashtable<String,CgProbe>();
        if (allOfMsLevel.containsKey(quant)) allOfQuant = allOfMsLevel.get(quant);
        allOfQuant.put(key, probe);
        allOfMsLevel.put(quant, allOfQuant);
        relevantMzHash.put(probe.getMsLevel(),allOfMsLevel);
      }
    }
    // now extract the chromatograms for each distinct mz Value
    Hashtable<QuantVO,Hashtable<String,Vector<LipidomicsChromatogram>>> chromsForMzs = new Hashtable<QuantVO,Hashtable<String,Vector<LipidomicsChromatogram>>>();
    for (Integer msLevel : relevantMzHash.keySet()){
      Hashtable<QuantVO,Hashtable<String,CgProbe>> mzsForChroms = relevantMzHash.get(msLevel);
      List<Integer> scansSorted = scanNumbersSorted.get(msLevel);
      Hashtable<Integer,Float> rts = retTimeLookup.get(msLevel);
      Hashtable<Integer,Float> noiseLevels = allNoise.get(msLevel);
      Hashtable<Integer,String> spectra = relevantSpectra.get(msLevel);
      //init the chromatograms that should be filled with the values
      for (QuantVO quant : mzsForChroms.keySet()){
        Hashtable<String,CgProbe> allOfQuant = mzsForChroms.get(quant);
        Hashtable<String,Vector<LipidomicsChromatogram>> chromsOfPartner = new Hashtable<String,Vector<LipidomicsChromatogram>>();
        if (chromsForMzs.containsKey(chromsOfPartner)) chromsOfPartner = chromsForMzs.get(chromsOfPartner);
        for (String key: allOfQuant.keySet()){
          LipidomicsChromatogram chrom = new LipidomicsChromatogram(new CgChromatogram(scansSorted.size()));
          chrom.LowerMzBand = analyzer.getMsnMzTolerance();
          chrom.UpperMzBand = analyzer.getMsnMzTolerance();
          LipidomicsChromatogram chromRelative = new LipidomicsChromatogram(new CgChromatogram(scansSorted.size()));
          for (int i=0; i!=scansSorted.size(); i++){
            chrom.Value[i][0] = rts.get(scansSorted.get(i));
            chrom.Value[i][1] = 0f;
            chrom.Value[i][2] = 0f;
            chrom.Value[i][3] = 0f;
            chromRelative.Value[i][0] = rts.get(scansSorted.get(i));
            chromRelative.Value[i][1] = 0f;
            chromRelative.Value[i][2] = 0f;
            chromRelative.Value[i][3] = 0f;
          }
          //LipidomicsChromatogram chromRelative = new LipidomicsChromatogram(chrom);
          Vector<LipidomicsChromatogram> chroms = new Vector<LipidomicsChromatogram>();
          chroms.add(chrom);
          chroms.add(chromRelative);
          chromsOfPartner.put(key, chroms);
        }
        chromsForMzs.put(quant, chromsOfPartner);
      }
      // now write the values into the files
      for (int i=0; i!=scansSorted.size();i ++){
        int consScanNumber = scansSorted.get(i);
        float noise = noiseLevels.get(consScanNumber);
        float noiseThreshold = noise*ChromatogramReader.NOISE_CUTOFF_MULTIPLICATOR;
        String spectrum = spectra.get(consScanNumber);
        spectrum = spectrum.substring(spectrum.indexOf(" ")+1);
        FloatBuffer buffer = ByteBuffer.wrap(Base64.decode(spectrum)).asFloatBuffer();
        int limit = buffer.limit();
        float mz = 0f;
        float intensity = 0f;
        float highestInt = 0f;
        for(int iItem = 0; iItem < limit; iItem++){
          if (iItem%2==0) mz = buffer.get();
          else{
            intensity = buffer.get();
            if (intensity>noiseThreshold){
              if (intensity>highestInt) highestInt = intensity;
              //this is for the absolute intensites
              for (QuantVO quant : mzsForChroms.keySet()){
                Hashtable<String,CgProbe> mzs = mzsForChroms.get(quant);
                for (String key : mzs.keySet()){
                  CgProbe probe = mzs.get(key);
                  if ((probe.Mz-probe.LowerMzBand)<mz && mz<(probe.Mz+probe.UpperMzBand)){
                    chromsForMzs.get(quant).get(key).get(0).Value[i][1] = intensity;
                  }
                }
              }
            }
          }
        }
        if (highestInt<=0f) continue;
        //this is for the relative intensity values;
        for (QuantVO quant : chromsForMzs.keySet()){
          Hashtable<String,Vector<LipidomicsChromatogram>> chromsForQuant = chromsForMzs.get(quant);
          for (String key : chromsForQuant.keySet()){
            Vector<LipidomicsChromatogram> chroms = chromsForQuant.get(key);
            if (startRelative<=chroms.get(1).Value[i][0] && chroms.get(1).Value[i][0]<=stopRelative)
              chroms.get(1).Value[i][1] = chroms.get(0).Value[i][1]/highestInt;
          }
        }
      }
    }
    
    //now calculate for each QuantVO where are the highest intensites reached for each fragment (absoluteValues)
    Hashtable<QuantVO,Hashtable<String,FloatFloatVO>> rtPeaksAbsolute = new Hashtable<QuantVO,Hashtable<String,FloatFloatVO>>();
    Hashtable<QuantVO,Float> proposedRts = new Hashtable<QuantVO,Float>();
    for (QuantVO quant : chromsForMzs.keySet()){
      Hashtable<String,FloatFloatVO> rts = new Hashtable<String,FloatFloatVO>();
      //String name = quant.getAnalyteClass()+quant.getAnalyteName()+"_"+quant.getDbs()+"_"+quant.getModName();
      float sumWeightings = 0f;
      float weightedRt = 0f;
      for (String key: chromsForMzs.get(quant).keySet()){
        LipidomicsChromatogram chrom = chromsForMzs.get(quant).get(key).get(0);
        int index = chrom.findIndexOfHighestIntensity();
        if (index>=0){
          rts.put(key, new FloatFloatVO(chrom.Value[index][0],chrom.Value[index][1]));
          weightedRt += chrom.Value[index][0]*chrom.Value[index][1];
          sumWeightings += chrom.Value[index][1];
//          System.out.println(name+": "+key+" "+chrom.Value[index][0]/60f);
        }
      }
      float proposedRt = weightedRt/sumWeightings;
      proposedRts.put(quant, proposedRt);
      rtPeaksAbsolute.put(quant, rts);
    }
    if (proposedRts.size()<2){
      String classes = "";
      String rts = "";
      for (SharedPeakContributionVO partner : shared.getPartners()){
        classes += partner.getQuantVO().getAnalyteClass()+"_"+partner.getSet().getNameString()+"-";
        rts += partner.getSet().getRt()+"/";
      }
      classes = classes.substring(0,classes.length()-1);
      rts = rts.substring(0,rts.length()-1);
//      System.out.println("!!!!!!!!!!!!!!! "+classes+" "+rts);
    }
    Iterator<QuantVO> it = proposedRts.keySet().iterator();
    QuantVO quant1 = it.next();
    QuantVO quant2 = it.next();
    checkProposedRtsInsidePeakBorders(proposedRts,shared);
    float threshold = checkForThresholdRt(proposedRts.get(quant1),proposedRts.get(quant2),scanNumbersSorted, retTimeLookup);
    if (threshold<0){
      Hashtable<QuantVO,Float> proposedRelativeRts = new Hashtable<QuantVO,Float>();
      for (QuantVO quant : chromsForMzs.keySet()){
//        Hashtable<String,FloatFloatVO> rts = new Hashtable<String,FloatFloatVO>();
//        String name = quant.getAnalyteClass()+quant.getAnalyteName()+"_"+quant.getDbs()+"_"+quant.getModName();
        float sumWeightings = 0f;
        float weightedRt = 0f;
        for (String key: chromsForMzs.get(quant).keySet()){
          LipidomicsChromatogram absChrom = chromsForMzs.get(quant).get(key).get(0);
          LipidomicsChromatogram relChrom = chromsForMzs.get(quant).get(key).get(1);
          Vector<Double> basePeaks = new Vector<Double>();
          for (int i=0; i!=relChrom.Value.length;i++){
            if (relChrom.Value[i][1]>0.999999f) basePeaks.add((double)relChrom.Value[i][0]);
          }
          int index = -1;
          float rt = 0f;
          if (basePeaks.size()>0){
            if ((basePeaks.size()*5)/4>=relChrom.Value.length){
              index = absChrom.findIndexOfHighestIntensity();
              rt = absChrom.Value[index][0];
            } else {
              rt = (float)Calculator.mean(basePeaks);
              index = absChrom.findIndexOfHighestIntensity();
            }
          } else {
            index = relChrom.findIndexOfHighestIntensity();
            if (index>0)
              rt = relChrom.Value[index][0];
          }
          if (index>=0){
//            rts.put(key, new FloatFloatVO(rt,relChrom.Value[index][1]));
            weightedRt += rt*relChrom.Value[index][1];
            sumWeightings += relChrom.Value[index][1];
//            System.out.println(name+": "+key+" "+relChrom.Value[index][0]/60f+" "+relChrom.Value[index][1]);
          }
        }
        float proposedRt = weightedRt/sumWeightings;
//        System.out.println(name+": "+proposedRt/60f);
        proposedRelativeRts.put(quant, proposedRt);
//        rtPeaksAbsolute.put(quant, rts);
      }
      if ((proposedRts.get(quant1)*1.02f>proposedRts.get(quant2) && proposedRts.get(quant2)*1.02f>proposedRts.get(quant1)) ||
          (proposedRelativeRts.get(quant1)<proposedRelativeRts.get(quant2) && proposedRts.get(quant1)<proposedRts.get(quant2)) ||
          (proposedRelativeRts.get(quant1)>proposedRelativeRts.get(quant2) && proposedRts.get(quant1)>proposedRts.get(quant2))){
        //System.out.println("111111111111111111111111111");
        checkProposedRtsInsidePeakBorders(proposedRelativeRts,shared);
        threshold = checkForThresholdRt(proposedRelativeRts.get(quant1),proposedRelativeRts.get(quant2),scanNumbersSorted, retTimeLookup);
        if (threshold>=0) proposedRts = proposedRelativeRts;
      }
    }
    if (threshold>=0){
      if (!areIntensitesOKAccordingToSeparation(shared,proposedRts,chromsForMzs)){
        threshold = -1;
      }
    }
    // now comes the physical split of the peak by the calculated threshold
    if (threshold>-1){
      SharedPeakContributionVO[] lowerHigher = getLowerAndHigherContribution(shared, proposedRts);
      SharedPeakContributionVO lowerContr = lowerHigher[0];
      lowerContr.makeCopyOfOldIdentification();
      SharedPeakContributionVO upperContr = lowerHigher[1];
      upperContr.makeCopyOfOldIdentification();
      lowerContr.getSet().setUpperRtHardLimit(threshold);
      analyzer.recalculatePeaksAccordingToHardLimits(lowerContr.getSet(),upperMsLevel,proposedRts.get(lowerContr.getQuantVO()));
      try {
        recalculateMS2Identification(lowerContr,analyzer);
      }
      catch (RulesException | IOException | SpectrummillParserException e) {
        e.printStackTrace();
      }
      upperContr.getSet().setLowerRtHardLimit(threshold);
      analyzer.recalculatePeaksAccordingToHardLimits(upperContr.getSet(),upperMsLevel,proposedRts.get(upperContr.getQuantVO()));
      try {
        recalculateMS2Identification(upperContr,analyzer);
      }
      catch (RulesException | IOException | SpectrummillParserException e) {
        e.printStackTrace();
      }
      splitted.add(lowerContr);
      splitted.add(upperContr);
    }
//    System.out.println("Threshold detected: "+threshold/60f);
    
//    System.out.println("!!!!!!!!!! Size: "+scanNumbersSorted.get(2).size());
//    for (Integer consScanNumber : scanNumbersSorted.get(2)){
//      System.out.println(retTimeLookup.get(2).get(consScanNumber)/60f);
//    }
    return splitted;
  }
  
  
  private static void checkProposedRtsInsidePeakBorders(Hashtable<QuantVO,Float> proposedRts, SharedMS1PeakVO shared){
    for (SharedPeakContributionVO partner : shared.getPartners()){
      if (!proposedRts.containsKey(partner.getQuantVO())) continue;
      float lowest = Float.MAX_VALUE;
      float highest = 0f;
      for (CgProbe probe : partner.getSet().getIsotopicProbes().get(0)){
        if (probe.LowerValley<lowest) lowest = probe.LowerValley;
        if (probe.UpperValley>highest) highest = probe.UpperValley;
      }
      float proposedRt = proposedRts.get(partner.getQuantVO());
      if (proposedRt<lowest) proposedRts.put(partner.getQuantVO(), lowest);
      if (proposedRt>highest) proposedRts.put(partner.getQuantVO(), highest);
    }
    
  }
  
  /**
   * rechecks the MS2 identfication based on absolute values
   * @param contr a species contributing to the shared peak
   * @param analyzer the access to the chrom files
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * throws CgException errors from the quantitation process
   */
  private static void recalculateMS2Identification(SharedPeakContributionVO contr, LipidomicsAnalyzer analyzer) throws RulesException, IOException, SpectrummillParserException, CgException{
    //TODO: the parameter before the last one is set to true in the meantime - maybe play around with caching in the future to improve calculation time
    MSnAnalyzer msnAnalyzer = new MSnAnalyzer(contr.getQuantVO().getAnalyteClass(),contr.getQuantVO().getModName(),contr.getSet(),analyzer,contr.getQuantVO(),true,true);  
    if (msnAnalyzer.checkStatus()==LipidomicsMSnSet.DISCARD_HIT) contr.setSet(null);
    else contr.setSet(msnAnalyzer.getResult());

  }
  

  

  /**
   * checks if it is possible to calculate a retention time, where the two peaks should be separated
   * @param rt1 retention time of peak summit of first species
   * @param rt2 retention time of peak summit of second species
   * @param scans a sorted list of consecutive scan numbers - key is the MS level
   * @param retTimeLookup the retention time lookup table for the found scans - first key is the MS level; second key is the consecutive scan number
   * @return the retention time threshold - if the value is negative: no threshold is detectable
   */
  private static float checkForThresholdRt(float rt1, float rt2, Hashtable<Integer,List<Integer>> scans,  Hashtable<Integer,Hashtable<Integer,Float>> retTimeLookup){
    float threshold = -1f;
    float lowerRt = rt1;
    float higherRt = rt2;
    if (lowerRt>higherRt){
      lowerRt = rt2;
      higherRt = rt1;
    }
    if (scanInBetween(lowerRt,higherRt,scans,retTimeLookup) || rtsCloseToNextScan(lowerRt,higherRt,scans,retTimeLookup)){
      threshold = (rt1+rt2)/2f;
    }
    return threshold;
  }

  /**
   * checks if there is at least one scan between tow retention time points  
   * @param lowerRt retention time of peak summit of the species that comes chromatographically first
   * @param higherRt retention time of peak summit of the species that comes chromatographically later
   * @param scans a sorted list of consecutive scan numbers - key is the MS level
   * @param retTimeLookup the retention time lookup table for the found scans - first key is the MS level; second key is the consecutive scan number
   * @return true if there is a scan in between
   */
  private static boolean scanInBetween(float lowerRt, float higherRt, Hashtable<Integer,List<Integer>> scans,  Hashtable<Integer,Hashtable<Integer,Float>> retTimeLookup){
    boolean scanInBetween = false;
    for (Integer msLevel : scans.keySet()){
      for (Integer scanNumber : scans.get(msLevel)){
        float rt = retTimeLookup.get(msLevel).get(scanNumber);
        if (lowerRt<rt && rt<higherRt){
          scanInBetween = true;
          return scanInBetween;
        }
      }
    }
    return scanInBetween;
  }
  
  /**
   * checks if the retention times of the species summits are sufficiently close to their next scan (within 20% time distance)
   * @param lowerRt retention time of peak summit of the species that comes chromatographically first
   * @param higherRt retention time of peak summit of the species that comes chromatographically later
   * @param scans a sorted list of consecutive scan numbers - key is the MS level
   * @param retTimeLookup the retention time lookup table for the found scans - first key is the MS level; second key is the consecutive scan number
   * @return true if the scans are close enough
   */
  private static boolean rtsCloseToNextScan(float lowerRt, float higherRt, Hashtable<Integer,List<Integer>> scans,  Hashtable<Integer,Hashtable<Integer,Float>> retTimeLookup){
    boolean closeEnough = false;
    float nextLowerScan = -1f;
    float nextUpperScan = -1f;
    float timeDistanceLower = Float.MAX_VALUE;
    float timeDistanceUpper = Float.MAX_VALUE;
    for (Integer msLevel : scans.keySet()){
      for (Integer scanNumber : scans.get(msLevel)){
        float rt = retTimeLookup.get(msLevel).get(scanNumber);
        if (rt<=lowerRt && (lowerRt-rt)<timeDistanceLower){
          nextLowerScan = rt;
          timeDistanceLower = lowerRt-rt;
        }
        if (rt>=higherRt && (rt-higherRt)<timeDistanceUpper){
          nextUpperScan = rt;
          timeDistanceUpper = rt-higherRt;
        }
      }
    }
    float scanDistance = nextUpperScan-nextLowerScan;
    float tolerance = scanDistance/5f;
    if (lowerRt<=(nextLowerScan+tolerance) && higherRt>=(nextUpperScan-tolerance))
      closeEnough = true;
    return closeEnough;
  }
  
  /**
   * this method checks if fragments of the species that comes later are comparably higher
   * in the region after the MSn peak summit of the later species to the same intensities
   * in the region before the MSn peak summit of the earlier species
   * @param shared VO containing information about the species that are in the peak and their distinct fragments
   * @param proposedRts the retention times of each analyte
   * @param chromsForMzs the chromatograms of the distinct fragments of each analyte
   * @return
   */
  private static boolean areIntensitesOKAccordingToSeparation(SharedMS1PeakVO shared, Hashtable<QuantVO,Float> proposedRts, 
      Hashtable<QuantVO,Hashtable<String,Vector<LipidomicsChromatogram>>> chromsForMzs){
    boolean ok = false;
    float lowerAreaOfLowerContr = 0f;
    float lowerAreaOfUpperContr = 0f;
    float upperAreaOfLowerContr = 0f;
    float upperAreaOfUpperContr = 0f;
    SharedPeakContributionVO[] lowerHigher = getLowerAndHigherContribution(shared, proposedRts);
    SharedPeakContributionVO lowerContr = lowerHigher[0];
    SharedPeakContributionVO upperContr = lowerHigher[1];
    float lowerRt = proposedRts.get(lowerContr.getQuantVO());
    float upperRt = proposedRts.get(upperContr.getQuantVO());

    for (Vector<LipidomicsChromatogram> chroms : chromsForMzs.get(lowerContr.getQuantVO()).values()){
      LipidomicsChromatogram chrom = chroms.get(0);
      for (int i=0; i!=chrom.Value.length; i++){
        if (chrom.Value[i][0]<=lowerRt) lowerAreaOfLowerContr += chrom.Value[i][1];
        if (chrom.Value[i][0]>=upperRt) upperAreaOfLowerContr += chrom.Value[i][1];
      }
    }
    for (Vector<LipidomicsChromatogram> chroms : chromsForMzs.get(upperContr.getQuantVO()).values()){
      LipidomicsChromatogram chrom = chroms.get(0);
      for (int i=0; i!=chrom.Value.length; i++){
        if (chrom.Value[i][0]<=lowerRt) lowerAreaOfUpperContr += chrom.Value[i][1];
        if (chrom.Value[i][0]>=upperRt) upperAreaOfUpperContr += chrom.Value[i][1];
      }
    }
    float lowerRatio = lowerAreaOfUpperContr/lowerAreaOfLowerContr;
    float upperRatio = upperAreaOfUpperContr/upperAreaOfLowerContr;
    if (upperRatio/lowerRatio>1.5f){
      ok = true;
    }
    return ok;
  }
  
  /**
   * from two SharedPeakContributionVO (species contributing to one MS1 peak), it is found out which one is the one that comes earlier and the one that comes later
   * ATTENTION: this is only applicable for two shared objects
   * @param shared VO containing information about the species that are in the peak and their distinct fragments
   * @param proposedRts the retention times of each analyte
   * @return
   */
  private static SharedPeakContributionVO[] getLowerAndHigherContribution(SharedMS1PeakVO shared, Hashtable<QuantVO,Float> proposedRts){
    Iterator<SharedPeakContributionVO> it = shared.getPartners().iterator();
    SharedPeakContributionVO contr1 = it.next();
    SharedPeakContributionVO contr2 = it.next();
    SharedPeakContributionVO lowerContr = contr1;
    SharedPeakContributionVO upperContr = contr2;
    float lowerRt = proposedRts.get(lowerContr.getQuantVO());
    float upperRt = proposedRts.get(upperContr.getQuantVO());
    if (lowerRt>upperRt){
      lowerContr = contr2;
      upperContr = contr1;
      lowerRt = proposedRts.get(lowerContr.getQuantVO());
      upperRt = proposedRts.get(upperContr.getQuantVO());
    }
    SharedPeakContributionVO[] lowerHigher = new SharedPeakContributionVO[2];
    lowerHigher[0] = lowerContr;
    lowerHigher[1] = upperContr;
    return lowerHigher;
  }
  
  public static Hashtable<Integer,Boolean> prepareCachedSpectra(LipidomicsAnalyzer analyzer, LipidParameterSet set, boolean readMSnSpectra) throws CgException{
    if (readMSnSpectra)
      return analyzer.prepareMSnSpectraCache(set.Mz[0]-LipidomicsConstants.getMs2PrecursorTolerance(), set.Mz[0]+LipidomicsConstants.getMs2PrecursorTolerance());
    else{
      return analyzer.checkMSnLevels();
    }
  }
}

