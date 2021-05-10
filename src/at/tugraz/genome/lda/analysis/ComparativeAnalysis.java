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


package at.tugraz.genome.lda.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.Vector;



import at.tugraz.genome.lda.LDAResultReader;
import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.DoubleCalculator;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.AbsoluteSettingsVO;
import at.tugraz.genome.lda.vos.DoubleStringVO;
import at.tugraz.genome.lda.vos.InternalStandardStatistics;
import at.tugraz.genome.lda.vos.IsotopicLabelVO;
import at.tugraz.genome.lda.vos.LipidClassSettingVO;
import at.tugraz.genome.lda.vos.ProbeVolConcVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.lda.vos.ResultCompGroupVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.VolumeConcVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.voutils.GeneralComparator;


/**
 * 
 * @author Juergen Hartler
 *
 */
public class ComparativeAnalysis extends ComparativeNameExtractor implements ComparativeResultsLookup
{

  private Hashtable<String,Hashtable<String,Vector<Hashtable<String,ResultAreaVO>>>> unprocessedResults_;
  private Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>> isNullResult_;

  
  // the results-hash is built in the following manner; the first key
  // is the sheet-name, since the sheet is one molecule group and represents
  // one logical entity for the analysis; the second one is the name of the file, 
  // and the contents are the found molecules (ResultAreaVO) and their isotopic areas
  private Hashtable<String,Hashtable<String,Vector<ResultAreaVO>>> allResults_;
  private Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> allResultsHash_;
  // the results-hash is built in the following manner; the first key
  // is the sheet-name, since the sheet is one molecule group and represents
  // one logical entity for the analysis; the second one is the name of the internal 
  // standard in this group, then comes the name of the internal standard and then
  // results of single searches (ResultAreaVO)
  private Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> isResults_;
  private Hashtable<String,Hashtable<String,Hashtable<String,Double>>> isCorrectionFactors_ ;

  private Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> esResults_;
  private Hashtable<String,Hashtable<String,Hashtable<String,Double>>> esCorrectionFactors_ ;
  
  private Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactorsToBestIS_ ;
  
  private Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactorsToBestES_ ;
  private Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactorsToBestESISCorr_ ;
  private Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactorsToBestESMedianCorr_ ;
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,Double>>>> correctionFactorsToBestESSingleCorr_;
  // contains statistics about the internal standards used; the first key
  // is the sheet-name, since the sheet is one molecule group and represents
  // one logical entity for the analysis; the second one is the name of the internal 
  // standard in this group, then the statistics about this standard
  private Hashtable<String,Hashtable<String,InternalStandardStatistics>> intStandStatistics_;
  
  private Hashtable<String,Hashtable<String,InternalStandardStatistics>> extStandStatistics_;
  private Hashtable<String,Hashtable<String,InternalStandardStatistics>> extStandStatisticsISCorrected_;
  private Hashtable<String,Hashtable<String,InternalStandardStatistics>> extStandStatisticsISMedianCorrected_;
  private Hashtable<String,Hashtable<String,Hashtable<String,InternalStandardStatistics>>> extStandStatisticsISSingleCorrected_;
  
  // standards in preferable sequence
  private Hashtable<String,Vector<String>> standardsOrderedConcerningReliability_;
  
  private Hashtable<String,Vector<String>> extstandsOrderedConcerningReliability_;
  private Hashtable<String,Vector<String>> extstandsISCorrOrderedConcerningReliability_;
  private Hashtable<String,Vector<String>> extstandsMedianCorrOrderedConcerningReliability_;
  private Hashtable<String,Hashtable<String,Vector<String>>> extstandsSingleCorrOrderedConcerningReliability_;

  private Hashtable<String,Hashtable<String,String>> modifications_;
//  Hashtable<String,Hashtable<String,Integer>> comparativeMaxIsotopes_;
  
  // this hash stores if a standard in the probe should be used for the comparison or
  // if it should be regarded as outlier; the first string is the sheet-name; the second
  // the experiment name; and the third the standard;
  private Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>> applicableStandards_;
  private Hashtable<String,Hashtable<String,String>> bestExpForStandard_;
  private Hashtable<String,Hashtable<String,String>> inProperForBestComparableProbe_;
  
  private Hashtable<String,Hashtable<String,String>> inProperForBestComparableProbeES_;
  private Hashtable<String,Hashtable<String,String>> inProperForBestComparableProbeESISCorr_;
  private Hashtable<String,Hashtable<String,String>> inProperForBestComparableProbeESMedianCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<String,String>>> inProperForBestComparableProbeESSingleCorr_;
  
  private Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>> applicableStandardsES_;
  private Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>> applicableStandardsESISCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>> applicableStandardsESMedianCorr_;
  
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>>> applicableStandardsESSingleCorr_;
  
  
  private Hashtable<String,Hashtable<String,String>> bestExpForExtStandard_;
  private Hashtable<String,Hashtable<String,String>> bestExpForExtStandardISCorr_;
  private Hashtable<String,Hashtable<String,String>> bestExpForExtStandardMedianCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<String,String>>> bestExpForExtStandardSingleCorr_;
  
  // values to which all of the other molecules should be compared; to relatively compare different probes
  // first key is the sheet name; second one is the experiment name
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> referenceValues_;
  
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> referenceValuesES_;
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> referenceValuesESISCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> referenceValuesESMedianCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>> referenceValuesESSingleCorr_;
  
  private Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> comparativeRatios_;
  private Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> comparativeRatiosGroups_;
  
  private Hashtable<String,Vector<String>> allMoleculeNames_;
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> medianOfRatios_;
  
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> medianOfESRatios_;
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> medianOfESRatiosISCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> medianOfESRatiosMedianCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>> medianOfESRatiosSingleCorr_;
  
  private Hashtable<String,Vector<String>> expNamesOfGroup_;
  
  private Hashtable<String,Integer> maxIsotopesOfGroup_;
  
  private ElementConfigParser elementParser_;
  
  private AbsoluteSettingsVO absSetting_;
  private Hashtable<String,Double> classCutoffs_;
  private int maxCutoffIsotope_;
  
  private Hashtable<String,Hashtable<Integer,Vector<Double>>> isSingleRefAreas_;
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<Double>>>> isSingleCorrectiveFactors_;
  private Hashtable<String,Hashtable<String,Integer>> correctionTypeISLookup_;

  private Hashtable<String,Hashtable<Integer,Vector<Double>>> esSingleRefAreasNoCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<Double>>>> esSingleCorrectiveFactorsNoCorr_;
  private Hashtable<String,Hashtable<String,Integer>> correctionTypeESLookup_;
//  Hashtable<String,Hashtable<String,Integer>> correctionTypeESLookupNoCorr_;
  private Hashtable<String,Hashtable<Integer,Vector<Double>>> esSingleRefAreasIntCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<Double>>>> esSingleCorrectiveFactorsIntCorr_;
//  Hashtable<String,Hashtable<String,Integer>> correctionTypeESLookupIntCorr_;
  private Hashtable<String,Hashtable<Integer,Vector<Double>>> esSingleRefAreasMedCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<Double>>>> esSingleCorrectiveFactorsMedCorr_;
//  Hashtable<String,Hashtable<String,Integer>> correctionTypeESLookupMedCorr_;
  
  private Hashtable<String,Hashtable<Integer,Hashtable<Integer,Vector<Double>>>> esSingleRefAreasSingleCorr_;
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Hashtable<Integer,Vector<Double>>>>> esSingleCorrectiveFactorsSingleCorr_;
  
  private Hashtable<String,Hashtable<Integer,VolumeConcVO>> isAmountLookup_;
  private Hashtable<String,Hashtable<Integer,VolumeConcVO>> esAmountLookup_;
  
  /** the sequence of the classes as in the quantification file*/
  private LinkedHashMap<String,Integer> classSequence_;
  /** the sequence of the analytes as in the quantification file*/
  private Hashtable<String,Vector<String>> correctAnalyteSequence_;
  /** the original quantification file objects*/
  private Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects_;
  
  private double expRtGroupingTime_;
  
  /** list of potential isotopic labels*/
  private List<IsotopicLabelVO> isoLabels_;
  
  /** the character encoding of the number of hydroxylation sites for the FA*/
  private HydroxyEncoding faHydroxyEncoding_;
  /** the character encoding of the number of hydroxylation sites for the LCB*/
  private HydroxyEncoding lcbHydroxyEncoding_;
  /** how many attached chains contains this class*/
  private Hashtable<String,Integer> chainsOfClass_; 
  
  
  public ComparativeAnalysis(Vector<File> resultFiles, String isSelectionPrefix, String esSelectionPrefix, AbsoluteSettingsVO absSetting, Hashtable<String,Double> classCutoffs, int maxCutoffIsotope, 
      LinkedHashMap<String,Integer> classSequence, Hashtable<String,Vector<String>> correctAnalyteSequence, Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects,
      double expRtGroupingTime){
    this(resultFiles, isSelectionPrefix, esSelectionPrefix,absSetting,classCutoffs,maxCutoffIsotope,null,null,classSequence,correctAnalyteSequence,
        quantObjects,expRtGroupingTime);
  }
  
  public ComparativeAnalysis(Vector<File> resultFiles, String isSelectionPrefix, String esSelectionPrefix,
      AbsoluteSettingsVO absSetting, Hashtable<String,Double> classCutoffs, int maxCutoffIsotope, Vector<String> groups, Hashtable<String,Vector<File>> filesOfGroup, 
      LinkedHashMap<String,Integer> classSequence, Hashtable<String,Vector<String>> correctAnalyteSequence,
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects, double expRtGroupingTime){
    super(resultFiles, isSelectionPrefix, esSelectionPrefix,groups,filesOfGroup);
    absSetting_ = absSetting;
    this.classSequence_ = classSequence;
    correctAnalyteSequence_ = correctAnalyteSequence;
    this.quantObjects_ = quantObjects;
    expRtGroupingTime_ = expRtGroupingTime;
    classCutoffs_ = classCutoffs;
    maxCutoffIsotope_ = maxCutoffIsotope;
    faHydroxyEncoding_ = null;
    lcbHydroxyEncoding_ = null;
    chainsOfClass_ = new Hashtable<String,Integer>();
  }
  
  public void parseInput() throws ExcelInputFileException, LipidCombinameEncodingException{
    elementParser_ = Settings.getElementParser();
    unprocessedResults_ = new Hashtable<String,Hashtable<String,Vector<Hashtable<String,ResultAreaVO>>>>();
    isNullResult_ = new Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>>();
    allResults_ = new Hashtable<String,Hashtable<String,Vector<ResultAreaVO>>>();
    allResultsHash_ = new Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>();
    modifications_ = new Hashtable<String,Hashtable<String,String>>();
    extractInformation();
    // this is to find out the experiments for the files of the group
    if (groups_!=null){
      expNamesOfGroup_ = new Hashtable<String,Vector<String>>();
      for (String group:groups_){
        Vector<String> expNames = new Vector<String>();
        for (File file :filesOfGroup_.get(group)){
          for (String expName: expNameToFile_.keySet()){
            if (file.getAbsolutePath().equalsIgnoreCase(expNameToFile_.get(expName).getAbsolutePath())){
              expNames.add(expName);
              break;
            }
          }
        }
        expNamesOfGroup_.put(group, expNames);
      }
    }
  }
  
  public void calculateStatistics() throws LipidCombinameEncodingException  {
    this.extractInternalStandardStatistics();
    this.extractExternalStandardStatistics();
    this.calculateRelativeComparisonValues();
  }
  
  
  private void extractExternalStandardStatistics(){
    esResults_ = new Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>();
    allESNames_ = new Hashtable<String,Hashtable<String,String>>();
    extractStandardVOsFromResults(esResults_, allESNames_, esSelectionPrefix_);
    calculteESRelativeCorrectionFactors();
    calculateStatisticsForESStandards();
    selectMostReliableESStandard();
    getApplicableESStandardsForExperiment();
    selectBestComparableESProbeForStandard();
    calculateRelativeCorrectiveValuesComparedToBestES();
    calculateESMedians();
    extractESAmountValues();
    calculateESReferenceValueForExperiments();
  }
  
  private void extractInternalStandardStatistics() /*throws NoISForGroupFoundException,NoISFoundException*/ {
    isResults_ = new Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>();
    allISNames_ = new Hashtable<String,Hashtable<String,String>>();
    
    maxIsotopesOfGroup_ = new Hashtable<String,Integer>();
    this.extractStandardVOsFromResults(isResults_, allISNames_, isSelectionPrefix_);

    calculteISRelativeCorrectionFactors();
    this.calculateStatisticsForISStandards();
    this.selectMostReliableISStandard();
    this.getApplicableISStandardsForExperiment();
    this.selectBestComparableISProbeForStandard();
    this.calculateRelativeCorrectiveValuesComparedToBestIS();
    this.calculateISMedians();
    this.extractISAmountValues();
    this.calculateISReferenceValueForExperiments();
  }
  
  private void extractStandardVOsFromResults(Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> standResults,
      Hashtable<String,Hashtable<String,String>> standNames, String selectionPrefix){
    for (String molGroupName : allResults_.keySet()){
      Hashtable<String,Vector<ResultAreaVO>> resultsMoleculeGroup = allResults_.get(molGroupName);
      Hashtable<String,Hashtable<String,ResultAreaVO>> isResultOfOneGroup = new Hashtable<String,Hashtable<String,ResultAreaVO>>();
      Hashtable<String,String> isOfOneMolGroup = new Hashtable<String,String>();
      int highestIsoValue = 0;
      for (String experimentName : resultsMoleculeGroup.keySet()){
        for (ResultAreaVO molecule : resultsMoleculeGroup.get(experimentName)){
          if (molecule.getMaxIsotope()>highestIsoValue)
            highestIsoValue = molecule.getMaxIsotope();
          if (selectionPrefix!=null && molecule.getMoleculeName().startsWith(selectionPrefix)){
            String isName = molecule.getMoleculeName();
            Hashtable<String,ResultAreaVO> oneIsResults = new Hashtable<String,ResultAreaVO>();
            if (isResultOfOneGroup.containsKey(isName)) oneIsResults = isResultOfOneGroup.get(isName);
            oneIsResults.put(experimentName, molecule);
            isResultOfOneGroup.put(isName,oneIsResults);
            isOfOneMolGroup.put(isName, isName);
          }
        }
      }
      if (isResultOfOneGroup.size()>0){
        standResults.put(molGroupName, isResultOfOneGroup);
      }
      standNames.put(molGroupName, isOfOneMolGroup);
      maxIsotopesOfGroup_.put(molGroupName, highestIsoValue);
    }
   
  }
  
  private boolean isISAvailable(){
    boolean isFound = false;
    for (String molGroupName : isResults_.keySet()){
      if (this.allISNames_.get(molGroupName).size()>0) isFound = true;
    }
    return isFound;
  }
 
  private void calculateStatisticsForESStandards(){
    extStandStatistics_ = new Hashtable<String,Hashtable<String,InternalStandardStatistics>>();
    this.calculateStatisticsForStandard(esResults_, allESNames_, esCorrectionFactors_,extStandStatistics_, true, ResultCompVO.NO_STANDARD_CORRECTION);
    boolean isFound = isISAvailable();
    if (isFound){
      extStandStatisticsISCorrected_ = new Hashtable<String,Hashtable<String,InternalStandardStatistics>>();
      this.calculateStatisticsForStandard(esResults_, allESNames_, esCorrectionFactors_,extStandStatisticsISCorrected_,true,ResultCompVO.STANDARD_CORRECTION_INTERNAL);
      extStandStatisticsISMedianCorrected_ = new Hashtable<String,Hashtable<String,InternalStandardStatistics>>();
      this.calculateStatisticsForStandard(esResults_, allESNames_, esCorrectionFactors_,extStandStatisticsISMedianCorrected_,true,ResultCompVO.STANDARD_CORRECTION_MEDIAN);
      extStandStatisticsISSingleCorrected_ = new Hashtable<String,Hashtable<String,Hashtable<String,InternalStandardStatistics>>>();
      for (String molGroupName : allResults_.keySet()){
        Hashtable<String,Hashtable<String,InternalStandardStatistics>> extStandStatisticsISSingleCorrectedGroup = new Hashtable<String,Hashtable<String,InternalStandardStatistics>>();
        for (String isName : allISNames_.get(molGroupName).keySet()){
          Hashtable<String,InternalStandardStatistics> molGroupStatistics = new Hashtable<String,InternalStandardStatistics>();
          if (allESNames_.get(molGroupName).size()>0)
            calculateStatisticsForStandardGroup(molGroupName,esResults_.get(molGroupName), esCorrectionFactors_.get(molGroupName), molGroupStatistics, true, correctionTypeISLookup_.get(molGroupName).get(isName).intValue(),isName);
          extStandStatisticsISSingleCorrectedGroup.put(isName, molGroupStatistics);
        }
        extStandStatisticsISSingleCorrected_.put(molGroupName, extStandStatisticsISSingleCorrectedGroup);
      }
    }  
  }
  
  private void calculateStatisticsForISStandards(){
    intStandStatistics_ = new Hashtable<String,Hashtable<String,InternalStandardStatistics>>();
    this.calculateStatisticsForStandard(isResults_, allISNames_,isCorrectionFactors_,intStandStatistics_, false, ResultCompVO.NO_STANDARD_CORRECTION);
    
  }
  
  private void calculateStatisticsForStandard(Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> standResults,
      Hashtable<String,Hashtable<String,String>> standNames, Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactors,
      Hashtable<String,Hashtable<String,InternalStandardStatistics>> standardsStatistics, boolean respectDilution, int standMethod){
    for (String molGroupName : standResults.keySet()){
      if (standNames.get(molGroupName).size()>0){
        Hashtable<String,Hashtable<String,ResultAreaVO>> standards = standResults.get(molGroupName);
        Hashtable<String,InternalStandardStatistics> molGroupISStatistics = new Hashtable<String,InternalStandardStatistics>();      
        calculateStatisticsForStandardGroup(molGroupName,standards, correctionFactors.get(molGroupName), molGroupISStatistics, respectDilution, standMethod,null);
        standardsStatistics.put(molGroupName, molGroupISStatistics);
      }
    }
  }
  
  private void calculateStatisticsForStandardGroup(String molGroupName,Hashtable<String,Hashtable<String,ResultAreaVO>> standards,
      Hashtable<String,Hashtable<String,Double>> correctionFactors,
      Hashtable<String,InternalStandardStatistics> molGroupISStatistics, boolean respectDilution, int standMethod, String correctionStandard){
        
    for (String standardName : standards.keySet()){
      Hashtable<String,ResultAreaVO> standardElements = standards.get(standardName);
      int howOftenFound = standardElements.size();
      float mean = 0;
      float median = 0;
      float stdev = Float.NaN;
      float coefficentOfVariation = Float.NaN;
      float lowerThreshold = Float.NaN;
      float upperThreshold = Float.NaN;
      int amountIsotopesFoundInAll = 0;
      if (howOftenFound>0){
        amountIsotopesFoundInAll = standardElements.values().iterator().next().getMaxIsotope();
      }
      Hashtable<String,Double> correctiveFactors = correctionFactors.get(standardName);
      howOftenFound = 0;
      for (String expName : standardElements.keySet()){
        ResultAreaVO resultVO = standardElements.get(expName);
        double area = correctAreaCorrespondingly(resultVO.getTotalArea(amountIsotopesFoundInAll),molGroupName,expName,amountIsotopesFoundInAll, correctiveFactors, respectDilution, standMethod, correctionStandard);
        if (area>0 && !Double.isNaN(area) && !Double.isInfinite(area))
          howOftenFound++;    
      }
      if (howOftenFound>1){
        for (ResultAreaVO resultVO : standardElements.values()){
          if (resultVO.getMaxIsotope()<amountIsotopesFoundInAll)
            amountIsotopesFoundInAll = resultVO.getMaxIsotope();
        }
      }
      if (howOftenFound>0){
        Vector<Float> valuesVector = new Vector<Float>();
        for (String expName : standardElements.keySet()){
          ResultAreaVO resultVO = standardElements.get(expName);
          float value = (float)this.correctAreaCorrespondingly(resultVO.getTotalArea(amountIsotopesFoundInAll),molGroupName,expName,amountIsotopesFoundInAll, correctiveFactors, respectDilution, standMethod, correctionStandard);
          if (value>0 && !Float.isNaN(value) && !Float.isInfinite(value))
            valuesVector.add(value);
        }
        Float[] medians = ComparativeAnalysis.medianPlusUpperLowerValues(valuesVector);
        float[] lowerUpperThreshold = ComparativeAnalysis.getLowerUpperOutlierThresholds(medians,2f);
        lowerThreshold = lowerUpperThreshold[0];
        upperThreshold = lowerUpperThreshold[1];
        Vector<Float> valuesVector2 = new Vector<Float>();
        // this removes the outlier from the statistics
        for (Float value : valuesVector){
          if (lowerThreshold<=value&&value<=upperThreshold)
            valuesVector2.add(value);
        }
        float[] values = new float[howOftenFound];
        howOftenFound = valuesVector2.size();
        for (int i=0;i!=howOftenFound;i++) values[i] = valuesVector2.get(i);
        mean = Calculator.mean(values);
        median = Calculator.median(valuesVector2);
        stdev = Calculator.stddeviation(values);
        coefficentOfVariation = stdev/mean;
      }
      InternalStandardStatistics statVO = new InternalStandardStatistics(standardName, howOftenFound, amountIsotopesFoundInAll, 
          mean,median, stdev, coefficentOfVariation,lowerThreshold,upperThreshold);
            molGroupISStatistics.put(standardName, statVO);
    }
  }
  
  private double correctAreaCorrespondingly(double area, String molGroupName, String expName, int amountIsotopesFoundInAll,
      Hashtable<String,Double> correctiveFactors, boolean respectDilution, int standMethod, String preferredStandard){
    double value = area * correctiveFactors.get(expName);
    if (this.absSetting_ != null && respectDilution){
      value = value * absSetting_.getClassSettings().get(molGroupName).getDilutionFactors().get(expName).floatValue();
    }
    if (standMethod == ResultCompVO.STANDARD_CORRECTION_INTERNAL || standMethod == ResultCompVO.STANDARD_CORRECTION_MEDIAN ||
        standMethod>ResultCompVO.STANDARD_CORRECTION_MEDIAN){
      String mostReliableStandard = preferredStandard;
      if ( (mostReliableStandard==null || mostReliableStandard.length()<1)&&standardsOrderedConcerningReliability_.containsKey(molGroupName)&& standardsOrderedConcerningReliability_.get(molGroupName).size()>0)
        //this is valid because we can correct the values of the ES only against the ones of IS
        mostReliableStandard = standardsOrderedConcerningReliability_.get(molGroupName).get(0);
      if (mostReliableStandard!=null&&mostReliableStandard.length()>0){
        String bestExp = bestExpForStandard_.get(molGroupName).get(mostReliableStandard);
        double correctionFactor = 1d;
        if (standMethod == ResultCompVO.STANDARD_CORRECTION_INTERNAL){
          correctionFactor = isResults_.get(molGroupName).get(mostReliableStandard).get(bestExp).getTotalArea(amountIsotopesFoundInAll)/referenceValues_.get(molGroupName).get(expName).get(amountIsotopesFoundInAll-1);
        }else if (standMethod == ResultCompVO.STANDARD_CORRECTION_MEDIAN){
          correctionFactor = medianOfRatios_.get(molGroupName).get(bestExp).get(amountIsotopesFoundInAll-1)/medianOfRatios_.get(molGroupName).get(expName).get(amountIsotopesFoundInAll-1);
        }else{
          if (isSingleCorrectiveFactors_.get(molGroupName).get(expName)!=null && isSingleCorrectiveFactors_.get(molGroupName).get(expName).get(standMethod)!=null)
            correctionFactor = isSingleCorrectiveFactors_.get(molGroupName).get(expName).get(standMethod).get(amountIsotopesFoundInAll-1);
          else return 0;
        }
      if (correctionFactor>0)
        value = value * correctionFactor;
      else
        value = 0;
      }
    }
    
    return value;
  }
  
  private void selectMostReliableESStandard(){
    extstandsOrderedConcerningReliability_ = new Hashtable<String,Vector<String>>();
    this.selectMostReliableStandard(extStandStatistics_,allESNames_,extstandsOrderedConcerningReliability_);
    boolean isFound = isISAvailable();
    if (isFound){
      extstandsISCorrOrderedConcerningReliability_ = new Hashtable<String,Vector<String>>();
      this.selectMostReliableStandard(extStandStatisticsISCorrected_,allESNames_,extstandsISCorrOrderedConcerningReliability_);
      extstandsMedianCorrOrderedConcerningReliability_ = new Hashtable<String,Vector<String>>();
      this.selectMostReliableStandard(extStandStatisticsISMedianCorrected_,allESNames_,extstandsMedianCorrOrderedConcerningReliability_);
      extstandsSingleCorrOrderedConcerningReliability_ = new Hashtable<String,Hashtable<String,Vector<String>>>();
      for (String molGroupName : allResults_.keySet()){
        Hashtable<String,Vector<String>> standardsOrderedCorrSingle = new Hashtable<String,Vector<String>>();
        for (String isName : allISNames_.get(molGroupName).keySet()){
          Hashtable<String,Vector<String>> standardsGroupCorrSingle = new Hashtable<String,Vector<String>>();
          if (allESNames_.get(molGroupName).size()>0){
            selectMostReliableStandardGroup(molGroupName,extStandStatisticsISSingleCorrected_.get(molGroupName).get(isName),
              allESNames_.get(molGroupName),standardsGroupCorrSingle);
            standardsOrderedCorrSingle.put(isName, standardsGroupCorrSingle.get(molGroupName));
          }  
        }
        extstandsSingleCorrOrderedConcerningReliability_.put(molGroupName, standardsOrderedCorrSingle);
      }
    }
  }
  
  private void selectMostReliableISStandard(){
    standardsOrderedConcerningReliability_ = new Hashtable<String,Vector<String>>();
    this.selectMostReliableStandard(intStandStatistics_,allISNames_,standardsOrderedConcerningReliability_);
  }
  
  // this method selects the most reliable of the standards, to which
  // the comparison of the other molecules should be made;
  private void selectMostReliableStandard(Hashtable<String,Hashtable<String,InternalStandardStatistics>> standStatistics,
      Hashtable<String,Hashtable<String,String>> standNames,Hashtable<String,Vector<String>> standardsOrdered){
    for (String molGroupName : standStatistics.keySet()){
      this.selectMostReliableStandardGroup(molGroupName,standStatistics.get(molGroupName), standNames.get(molGroupName), standardsOrdered);
    }
  }
  
  @SuppressWarnings("unchecked")
  private void selectMostReliableStandardGroup(String molGroupName,Hashtable<String,InternalStandardStatistics> standStatistics,
      Hashtable<String,String> standNames,Hashtable<String,Vector<String>> standardsOrdered){    
    if (standNames.size()>0){
      Vector<String> orderedConcerningReliability = new Vector<String>();
      Hashtable<String,InternalStandardStatistics> molGroupISStatistics = standStatistics;
      // first we have to figure out the standard that is most often found;
      // split then the standards in 2 groups first the the ones who have been found 90-100%
      // if the most often found standards and the rest
      Vector<InternalStandardStatistics> oftenFound = new Vector<InternalStandardStatistics>();
      List<InternalStandardStatistics> rarelyFound = new ArrayList<InternalStandardStatistics>();
      int mostOftenFound = 0;
      for (InternalStandardStatistics stat : molGroupISStatistics.values()){
        if (stat.getHowOftenFound()>mostOftenFound)
          mostOftenFound = stat.getHowOftenFound();
      }
      int foundThreshold = (int)Calculator.roundFloat(((float)mostOftenFound*0.9f),0);
      for (InternalStandardStatistics stat : molGroupISStatistics.values()){
        if (stat.getHowOftenFound()>=foundThreshold) oftenFound.add(stat);
        else rarelyFound.add(stat);
      }
      // now find the one with the smallest coefficient of variation
      // if the CV differs just less than 5% between neighbouring IS the one with the higher mean is taken 
      Collections.sort(oftenFound, new GeneralComparator("at.tugraz.genome.lda.vos.InternalStandardStatistics", "getCoefficentOfVariation", "java.lang.Float"));
      for (int i=0;i!=oftenFound.size();i++){
        Vector<InternalStandardStatistics> sortByArea = new Vector<InternalStandardStatistics>();
        sortByArea.add(oftenFound.get(i));
        for (int j=(i+1); j<oftenFound.size();j++){
          if (oftenFound.get(j).getCoefficentOfVariation()<(oftenFound.get(i).getCoefficentOfVariation()+0.05))
            sortByArea.add(oftenFound.get(j));
          else
            break;
        }
        Collections.sort(sortByArea, new GeneralComparator("at.tugraz.genome.lda.vos.InternalStandardStatistics", "getMedian", "java.lang.Float"));
        for (int j=(sortByArea.size()-1);j!=-1;j--){
          orderedConcerningReliability.add(sortByArea.get(j).getStandardName());
        }
        i += (sortByArea.size()-1);
      }
      if (rarelyFound.size()>0){
        Collections.sort(rarelyFound, new GeneralComparator("at.tugraz.genome.lda.vos.InternalStandardStatistics", "getHowOftenFound", "java.lang.Integer"));
        for (int i=(rarelyFound.size()-1);i!=-1;i--){
          orderedConcerningReliability.add(rarelyFound.get(i).getStandardName());
        }
      }
      standardsOrdered.put(molGroupName,orderedConcerningReliability );
    }
  }
  
  // selects against which probe the standard should be compared to
  private void selectBestComparableESProbeForStandard(){
    bestExpForExtStandard_ = new Hashtable<String,Hashtable<String,String>>();
    selectBestComparableISProbeForStandard(esResults_,extStandStatistics_,allESNames_,esCorrectionFactors_,extstandsOrderedConcerningReliability_,applicableStandardsES_,
        inProperForBestComparableProbeES_,bestExpForExtStandard_,true,ResultCompVO.NO_STANDARD_CORRECTION);
    boolean isFound = this.isISAvailable();
    if (isFound){
      bestExpForExtStandardISCorr_ = new Hashtable<String,Hashtable<String,String>>();
      selectBestComparableISProbeForStandard(esResults_,extStandStatisticsISCorrected_,allESNames_,esCorrectionFactors_,extstandsISCorrOrderedConcerningReliability_,applicableStandardsESISCorr_,
          inProperForBestComparableProbeESISCorr_,bestExpForExtStandardISCorr_,true,ResultCompVO.STANDARD_CORRECTION_INTERNAL);
      bestExpForExtStandardMedianCorr_ = new Hashtable<String,Hashtable<String,String>>();
      selectBestComparableISProbeForStandard(esResults_,extStandStatisticsISMedianCorrected_,allESNames_,esCorrectionFactors_,extstandsMedianCorrOrderedConcerningReliability_,applicableStandardsESMedianCorr_,
          inProperForBestComparableProbeESMedianCorr_,bestExpForExtStandardMedianCorr_,true,ResultCompVO.STANDARD_CORRECTION_MEDIAN);
      bestExpForExtStandardSingleCorr_ = new Hashtable<String,Hashtable<String,Hashtable<String,String>>>();
      for (String molGroupName : allResults_.keySet()){
        Hashtable<String,Hashtable<String,String>> bestExpForExtStandardGroup = new Hashtable<String,Hashtable<String,String>>();
        for (String isName : allISNames_.get(molGroupName).keySet()){
          Hashtable<String,String> bestExpForExtStandard = new Hashtable<String,String>();
          if (allESNames_.get(molGroupName).size()>0)
            selectBestComparableISProbeForStandardGroup(molGroupName, esResults_.get(molGroupName), extStandStatisticsISSingleCorrected_.get(molGroupName).get(isName), 
                esCorrectionFactors_.get(molGroupName), extstandsSingleCorrOrderedConcerningReliability_.get(molGroupName).get(isName), applicableStandardsESSingleCorr_.get(molGroupName).get(isName), 
                inProperForBestComparableProbeESSingleCorr_.get(molGroupName).get(isName), bestExpForExtStandard, true, correctionTypeISLookup_.get(molGroupName).get(isName).intValue(),isName);
          if (bestExpForExtStandard.size()>0)
            bestExpForExtStandardGroup.put(isName, bestExpForExtStandard);
        }
        if (bestExpForExtStandardGroup.size()>0)
          bestExpForExtStandardSingleCorr_.put(molGroupName, bestExpForExtStandardGroup);
      }      
    }
  }
  
  private void selectBestComparableISProbeForStandard(){
    bestExpForStandard_ = new Hashtable<String,Hashtable<String,String>>();
    selectBestComparableISProbeForStandard(isResults_,intStandStatistics_,allISNames_,isCorrectionFactors_,standardsOrderedConcerningReliability_,applicableStandards_,
        inProperForBestComparableProbe_,bestExpForStandard_,false,ResultCompVO.NO_STANDARD_CORRECTION);
  }
  
  private void selectBestComparableISProbeForStandard(Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> standResults,
      Hashtable<String,Hashtable<String,InternalStandardStatistics>> standStatistics, Hashtable<String,Hashtable<String,String>> standNames, 
      Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactors, Hashtable<String,Vector<String>> standardsOrderedConcerningReliability, 
      Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>> applicableStandards,Hashtable<String,Hashtable<String,String>> inProperForBestComparableProbe,
      Hashtable<String,Hashtable<String,String>> bestExpForStandard, boolean respectDilution, int standMethod){
  
  // selects against which probe the standard should be compared to
    for (String molGroupName : standardsOrderedConcerningReliability.keySet()){
      if (standNames.get(molGroupName).size()>0){
        Hashtable<String,String> expForStandard = new Hashtable<String,String>();
        this.selectBestComparableISProbeForStandardGroup(molGroupName,standResults.get(molGroupName), standStatistics.get(molGroupName),
            correctionFactors.get(molGroupName),standardsOrderedConcerningReliability.get(molGroupName),
            applicableStandards.get(molGroupName),inProperForBestComparableProbe.get(molGroupName),expForStandard,
            respectDilution, standMethod,null);
        if (expForStandard.size()>0)
          bestExpForStandard.put(molGroupName, expForStandard);
      }
    }
  }  
      
  private void selectBestComparableISProbeForStandardGroup(String molGroupName, Hashtable<String,Hashtable<String,ResultAreaVO>> standResults,
    Hashtable<String,InternalStandardStatistics> standStatistics, Hashtable<String,Hashtable<String,Double>> correctionFactors,
    Vector<String> orderedConcerningReliability, Hashtable<String,Hashtable<String,Boolean>> applicableStandards,
    Hashtable<String,String> inProperForBestComparableProbe, Hashtable<String,String> expForStandard,
    boolean respectDilution, int standMethod, String corrStandard){    
    for (String isName : orderedConcerningReliability){
      // first we have to find experiments that have found most of the standards
      Hashtable<String,Integer> experimentsThatFoundStandard = new Hashtable<String,Integer>();
      int expFoundHowManyStandardsMax = 0;
      for (String expName : applicableStandards.keySet()){
        if (applicableStandards.get(expName).get(isName)&&!inProperForBestComparableProbe.containsKey(expName)){
          experimentsThatFoundStandard.put(expName,applicableStandards.get(expName).size());
          if (experimentsThatFoundStandard.get(expName)>expFoundHowManyStandardsMax)
            expFoundHowManyStandardsMax = experimentsThatFoundStandard.get(expName);
        }
      }
      int lowerThreshold = (int)Calculator.roundFloat((float)expFoundHowManyStandardsMax*0.8f,0);
      Hashtable<String,String> expThatAreFine = new Hashtable<String,String>();
      for (String expName : experimentsThatFoundStandard.keySet()){
        if (experimentsThatFoundStandard.get(expName)>=lowerThreshold)
          expThatAreFine.put(expName, expName);
      }
      InternalStandardStatistics stat = standStatistics.get(isName);
      float median = stat.getMedian();
      if (median>1){
        Hashtable<String,ResultAreaVO> isValues =  standResults.get(isName);
        Hashtable<String,Double> correctiveFactors = correctionFactors.get(isName);
        String currentExperiment = isValues.keySet().iterator().next();
        float closestToMedian = this.calculateRatioToMedianOneDir((float)correctAreaCorrespondingly(isValues.get(currentExperiment).getTotalArea(stat.getAmountIsotopesFoundInAll()), molGroupName, currentExperiment, stat.getAmountIsotopesFoundInAll(), correctiveFactors, respectDilution, standMethod,corrStandard),median);
        for (String experimentName : isValues.keySet()){
          if (expThatAreFine.containsKey(experimentName)){
            float area = (float)(correctAreaCorrespondingly(isValues.get(experimentName).getTotalArea(stat.getAmountIsotopesFoundInAll()), molGroupName, experimentName, stat.getAmountIsotopesFoundInAll(), correctiveFactors, respectDilution, standMethod,corrStandard));
            float ratio = this.calculateRatioToMedianOneDir(area,median);
            if (ratio<closestToMedian){
              closestToMedian = ratio;
              currentExperiment = experimentName;
            }
          }
        }
        expForStandard.put(isName, currentExperiment );
      }
    }
  }
  
  private void getApplicableESStandardsForExperiment(){
    inProperForBestComparableProbeES_ = new Hashtable<String,Hashtable<String,String>>();
    applicableStandardsES_ = new Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>>();
    this.getApplicableStandardsForExperiment(esResults_, allESNames_, esCorrectionFactors_, extStandStatistics_,
        inProperForBestComparableProbeES_, applicableStandardsES_,true,ResultCompVO.NO_STANDARD_CORRECTION);
    boolean isFound = isISAvailable();
    if (isFound){
      inProperForBestComparableProbeESISCorr_ = new Hashtable<String,Hashtable<String,String>>();
      applicableStandardsESISCorr_ = new Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>>();
      this.getApplicableStandardsForExperiment(esResults_, allESNames_, esCorrectionFactors_, extStandStatisticsISCorrected_,
          inProperForBestComparableProbeESISCorr_, applicableStandardsESISCorr_,true,ResultCompVO.STANDARD_CORRECTION_INTERNAL);
      inProperForBestComparableProbeESMedianCorr_ = new Hashtable<String,Hashtable<String,String>>();
      applicableStandardsESMedianCorr_ = new Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>>();
      this.getApplicableStandardsForExperiment(esResults_, allESNames_, esCorrectionFactors_, extStandStatisticsISMedianCorrected_,
          inProperForBestComparableProbeESMedianCorr_, applicableStandardsESMedianCorr_,true,ResultCompVO.STANDARD_CORRECTION_MEDIAN);
      applicableStandardsESSingleCorr_ = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>>>();
      inProperForBestComparableProbeESSingleCorr_ = new Hashtable<String,Hashtable<String,Hashtable<String,String>>>();
      for (String molGroupName : allResults_.keySet()){
        Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>> applicableStandardsCorrSingle = new Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>>();
        Hashtable<String,Hashtable<String,String>> inProperForBestComparableProbeSingleCorr = new Hashtable<String,Hashtable<String,String>>();
        for (String isName : allISNames_.get(molGroupName).keySet()){
          Hashtable<String,String> inProperForBestComparableProbe = new Hashtable<String,String>();
          Hashtable<String,Hashtable<String,Boolean>> applicableStandards = new Hashtable<String,Hashtable<String,Boolean>>();
          if (allESNames_.get(molGroupName).size()>0)
            this.getApplicableStandardsForExperimentGroup(molGroupName, esResults_.get(molGroupName), esCorrectionFactors_.get(molGroupName), 
              extStandStatisticsISSingleCorrected_.get(molGroupName).get(isName), inProperForBestComparableProbe,applicableStandards, true, 
              correctionTypeISLookup_.get(molGroupName).get(isName).intValue(), isName);
          applicableStandardsCorrSingle.put(isName, applicableStandards);
          inProperForBestComparableProbeSingleCorr.put(isName, inProperForBestComparableProbe);
          
        }
        applicableStandardsESSingleCorr_.put(molGroupName, applicableStandardsCorrSingle);
        inProperForBestComparableProbeESSingleCorr_.put(molGroupName, inProperForBestComparableProbeSingleCorr);
      }
    }
  }
  
  private void getApplicableISStandardsForExperiment(){
    inProperForBestComparableProbe_ = new Hashtable<String,Hashtable<String,String>>();
    applicableStandards_ = new Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>>();
    this.getApplicableStandardsForExperiment(isResults_, allISNames_, isCorrectionFactors_, intStandStatistics_,
        inProperForBestComparableProbe_, applicableStandards_,false,ResultCompVO.NO_STANDARD_CORRECTION);
  }
  
  private void getApplicableStandardsForExperiment(Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> standResults,
      Hashtable<String,Hashtable<String,String>> standNames, Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactors,
      Hashtable<String,Hashtable<String,InternalStandardStatistics>> standStatistics, Hashtable<String,Hashtable<String,String>> inproperForBestComparableProbe, 
      Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>> applicableStandards, boolean respectDilution, int standMethod){
    for (String moleculeGroup : standResults.keySet()){
      if (standNames.get(moleculeGroup).size()>0){
        Hashtable<String,Hashtable<String,Boolean>> groupApplic = new Hashtable<String,Hashtable<String,Boolean>>();
        Hashtable<String,String> inproper = new Hashtable<String,String>();
        this.getApplicableStandardsForExperimentGroup(moleculeGroup,standResults.get(moleculeGroup),correctionFactors.get(moleculeGroup),
            standStatistics.get(moleculeGroup), inproper, groupApplic, respectDilution, standMethod,null);
        inproperForBestComparableProbe.put(moleculeGroup, inproper);
        applicableStandards.put(moleculeGroup, groupApplic);
      }    
    }
  }
    
  private void getApplicableStandardsForExperimentGroup(String moleculeGroup, Hashtable<String,Hashtable<String,ResultAreaVO>> standResults,
      Hashtable<String,Hashtable<String,Double>> correctiveFactors, Hashtable<String,InternalStandardStatistics> standStatistics,
      Hashtable<String,String> inproperForBestComparableProbe,
      Hashtable<String,Hashtable<String,Boolean>> groupApplic, boolean respectDilution, int standMethod, String corrStandard){
    for (int i=1;i!=(maxIsotopesOfGroup_.get(moleculeGroup)+1);i++){
      for (File resultFile : resultFiles_){
        String fileName = StaticUtils.extractFileName(resultFile.getAbsolutePath());
        fileName = fileName.substring(0,fileName.lastIndexOf("."));
        //fileName = fileName.substring(charsToCutPrev_,fileName.length()-charsToCut_);
        fileName = fileName.substring(0,fileName.length()-charsToCut_);
        Hashtable<String,Boolean> isApplicable = new Hashtable<String,Boolean>();
        boolean foundOneStandard = false;
        for (String standardName: standResults.keySet()){
          Hashtable<String,ResultAreaVO> areasForStandard = standResults.get(standardName);
          if (areasForStandard.containsKey(fileName)){
            ResultAreaVO areaVO = areasForStandard.get(fileName);
            InternalStandardStatistics stat = standStatistics.get(standardName);
            int isosToTake = i;
            if (stat.getAmountIsotopesFoundInAll()<isosToTake)
              isosToTake = stat.getAmountIsotopesFoundInAll();
            float comparableArea = (float)correctAreaCorrespondingly(areaVO.getTotalArea(isosToTake), moleculeGroup, fileName, isosToTake,
                correctiveFactors.get(standardName), respectDilution, standMethod, corrStandard);
                  
            if (stat.getLowerOutlierThreshold()<=comparableArea&&comparableArea<=stat.getUpperOutlierThreshold()){
              isApplicable.put(standardName, true);
              foundOneStandard = true;
            }else{
              isApplicable.put(standardName, false);
            }  
          }else
            isApplicable.put(standardName, false);
        }
        if (!foundOneStandard){
          inproperForBestComparableProbe.put(fileName, fileName);
          for (String standardName: standResults.keySet()){
            Hashtable<String,ResultAreaVO> areasForStandard = standResults.get(standardName);
            if (areasForStandard.containsKey(fileName))
              isApplicable.put(standardName, true);
            else 
              isApplicable.put(standardName, false);
          }         
        }
        groupApplic.put(fileName, isApplicable);
        // here I deleted the calculation of a reference value, because it is not read anywhere
        // if there possibly occur problems -> check the log files
      }
    }
  }
  
  private void calculateESReferenceValueForExperiments(){
    referenceValuesES_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>();
    calculateReferenceValueForExperiments(esResults_,allESNames_,extStandStatistics_,correctionFactorsToBestES_,applicableStandardsES_,extstandsOrderedConcerningReliability_,
        bestExpForExtStandard_,referenceValuesES_,true,ResultCompVO.NO_STANDARD_CORRECTION);
    if (this.isISAvailable()){
      referenceValuesESISCorr_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>();
      calculateReferenceValueForExperiments(esResults_,allESNames_,extStandStatisticsISCorrected_,correctionFactorsToBestESISCorr_,applicableStandardsESISCorr_,extstandsISCorrOrderedConcerningReliability_,
          bestExpForExtStandardISCorr_,referenceValuesESISCorr_,true,ResultCompVO.STANDARD_CORRECTION_INTERNAL);
      referenceValuesESMedianCorr_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>();
      calculateReferenceValueForExperiments(esResults_,allESNames_,extStandStatisticsISMedianCorrected_,correctionFactorsToBestESMedianCorr_,applicableStandardsESMedianCorr_,extstandsMedianCorrOrderedConcerningReliability_,
          bestExpForExtStandardMedianCorr_,referenceValuesESMedianCorr_,true,ResultCompVO.STANDARD_CORRECTION_MEDIAN);     
      referenceValuesESSingleCorr_ = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>>();
      for (String molGroupName : allResults_.keySet()){
        Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> referenceValuesESGroup = new Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>();
        for (String isName : allISNames_.get(molGroupName).keySet()){
          Hashtable<String,Hashtable<Integer,Double>> refs = new Hashtable<String,Hashtable<Integer,Double>>();
          if (allESNames_.get(molGroupName).size()>0 && bestExpForExtStandardSingleCorr_.get(molGroupName)!=null && bestExpForExtStandardSingleCorr_.get(molGroupName).get(isName)!=null)
            this.calculateReferenceValueForExperimentsGroup(molGroupName, esResults_.get(molGroupName), extStandStatisticsISSingleCorrected_.get(molGroupName).get(isName), 
                correctionFactorsToBestESSingleCorr_.get(molGroupName).get(isName), applicableStandardsESSingleCorr_.get(molGroupName).get(isName), 
                extstandsSingleCorrOrderedConcerningReliability_.get(molGroupName).get(isName), bestExpForExtStandardSingleCorr_.get(molGroupName).get(isName), 
                refs, true, correctionTypeISLookup_.get(molGroupName).get(isName).intValue(),isName);
          referenceValuesESGroup.put(isName, refs);
        }
        referenceValuesESSingleCorr_.put(molGroupName, referenceValuesESGroup);
      }
    }
  }
  
  private void calculateISReferenceValueForExperiments(){
    referenceValues_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>();
    calculateReferenceValueForExperiments(isResults_,allISNames_,intStandStatistics_,correctionFactorsToBestIS_,applicableStandards_,standardsOrderedConcerningReliability_,
        bestExpForStandard_,referenceValues_,false,ResultCompVO.NO_STANDARD_CORRECTION);
  }
  
  private void calculateReferenceValueForExperiments(Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> standResults,
      Hashtable<String,Hashtable<String,String>> standNames,
      Hashtable<String,Hashtable<String,InternalStandardStatistics>> standStatistics,
      Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactorsGroup,
      Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>> applicableStandards,
      Hashtable<String,Vector<String>> standardsOrderedConcerningReliability, 
      Hashtable<String,Hashtable<String,String>> bestExpForStandard,
      Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> referenceValues,
      boolean respectDilution, int standMethod){
    
    for (String moleculeGroup : applicableStandards.keySet()){
      if (standNames.get(moleculeGroup).size()>0 && bestExpForStandard.get(moleculeGroup) !=null){
        Hashtable<String,Hashtable<Integer,Double>> refs = new Hashtable<String,Hashtable<Integer,Double>>();
        calculateReferenceValueForExperimentsGroup(moleculeGroup,standResults.get(moleculeGroup),standStatistics.get(moleculeGroup),
            correctionFactorsGroup.get(moleculeGroup), applicableStandards.get(moleculeGroup),
            standardsOrderedConcerningReliability.get(moleculeGroup),bestExpForStandard.get(moleculeGroup),
            refs, respectDilution, standMethod,null);
        referenceValues.put(moleculeGroup, refs);
      }
    }
  }  
      
  private void calculateReferenceValueForExperimentsGroup(String moleculeGroup, Hashtable<String,Hashtable<String,ResultAreaVO>> standResults,
      Hashtable<String,InternalStandardStatistics> standStatistics, Hashtable<String,Hashtable<String,Double>> correctionFactors,
      Hashtable<String,Hashtable<String,Boolean>> standsForGroup, Vector<String> standardsOrderedConcerningReliability,
      Hashtable<String,String> bestExpForStandard,Hashtable<String,Hashtable<Integer,Double>> refs,
      boolean respectDilution, int standMethod, String corrStandard){     
    if (standardsOrderedConcerningReliability.size()>0){
      String mostReliableStandard = standardsOrderedConcerningReliability.get(0);
      String expEverythingHasToBeComparedTo = bestExpForStandard.get(mostReliableStandard);
      InternalStandardStatistics stat = standStatistics.get(mostReliableStandard);
      for (int iso=1;iso!=(maxIsotopesOfGroup_.get(moleculeGroup)+1);iso++){
        int isoToUse = iso;
        if (iso>stat.getAmountIsotopesFoundInAll()){
          isoToUse = stat.getAmountIsotopesFoundInAll();
        }
        // the correction factor should be 1 here; it is just there for completeness
        double refValue =  (float)this.correctAreaCorrespondingly(standResults.get(mostReliableStandard).get(expEverythingHasToBeComparedTo).getTotalArea(isoToUse), moleculeGroup, expEverythingHasToBeComparedTo, isoToUse, correctionFactors.get(mostReliableStandard), respectDilution, standMethod,corrStandard); 
        if (iso>stat.getAmountIsotopesFoundInAll()){
//          String chemicalFormula = standResults.get(mostReliableStandard).get(expEverythingHasToBeComparedTo).getChemicalFormula();
          // the correction factor should be 1 here; it is just there for completeness
//          double zeroIso = (float)this.correctAreaCorrespondingly(standResults.get(mostReliableStandard).get(expEverythingHasToBeComparedTo).getTotalArea(1), moleculeGroup, expEverythingHasToBeComparedTo, 1, correctionFactors.get(mostReliableStandard), respectDilution, standMethod,corrStandard); 
//          refValue =  this.theoreticallyCorrectForAddIsotope(chemicalFormula,refValue,zeroIso,stat.getAmountIsotopesFoundInAll(),iso);
          refValue = correctAreaCorrespondingly(standResults.get(mostReliableStandard).get(expEverythingHasToBeComparedTo).getTheoreticalIsotopeValue(elementParser_, iso),moleculeGroup, expEverythingHasToBeComparedTo, 1,
              correctionFactors.get(mostReliableStandard), respectDilution, standMethod, corrStandard);
        }
        for (String expName : standsForGroup.keySet()){
          Hashtable<Integer,Double> isoValues = new Hashtable<Integer,Double>();
          if (refs.containsKey(expName))
            isoValues = refs.get(expName);
          String partner = null;
          for (int i=0;i!=standardsOrderedConcerningReliability.size();i++){
            if (partner==null){
              String standard = standardsOrderedConcerningReliability.get(i);
              String expToCompareTo = bestExpForStandard.get(standard);

              if (this.doesPartnerFit(expName,expToCompareTo,standsForGroup)){
                if (this.doesPartnerFit(expToCompareTo, expEverythingHasToBeComparedTo, standsForGroup)){
                  partner = expToCompareTo;
                  break;
                } 
              }
            }
          }
          double refToPut = Double.NaN;
          if (partner == null /*|| refValue<=0*/){
            System.out.println("!!!!!!!!!!!!!!!!!!!! The experiment "+expName+" cannot be compared to anything !!!!!!!!!!");
            refToPut = Double.NaN;
          }else if (expName.equalsIgnoreCase(expEverythingHasToBeComparedTo)){
            refToPut = refValue;
          } else if (partner.equalsIgnoreCase(expEverythingHasToBeComparedTo)){
            refToPut = refValue*this.getRelativeValue(standResults,expName, partner, moleculeGroup,iso,correctionFactors,standsForGroup,respectDilution, standMethod, corrStandard);
          } else{
            double relativePartnerToExpEverythingHasToBeComparedTo = this.getRelativeValue(standResults,partner, expEverythingHasToBeComparedTo, moleculeGroup,isoToUse,correctionFactors,standsForGroup,respectDilution, standMethod,corrStandard);
            refToPut = refValue*this.getRelativeValue(standResults,expName, partner, moleculeGroup,iso,correctionFactors,standsForGroup,respectDilution, standMethod,corrStandard)*relativePartnerToExpEverythingHasToBeComparedTo;
          }
          isoValues.put(iso-1, refToPut);
          refs.put(expName, isoValues);
        }
      }
    }
    
  }
  
  private double getRelativeValue(Hashtable<String,Hashtable<String,ResultAreaVO>> standResults,
      String exp,String partner,String moleculeGroup,int isotopes,
      Hashtable<String,Hashtable<String,Double>> correctionFactors,
      Hashtable<String,Hashtable<String,Boolean>> standsForGroup,boolean respectDilution, int standMethod, String corrStandard){
    Vector<Float> ratioOfIs = new Vector<Float>();
    for (String isName : standsForGroup.get(exp).keySet()){
      // asks if the standard is applicable for the experiment
      if (standsForGroup.get(exp).get(isName)){
     // asks if the standard is applicable for the partner
        if (standsForGroup.containsKey(partner)&&
            standsForGroup.get(partner).get(isName)){
          // exists the standard - I do not know why I have to ask this?
          if (standResults.get(isName).get(exp)!=null){
            float areaExp = (float)this.correctAreaCorrespondingly(standResults.get(isName).get(exp).getTotalArea(isotopes), moleculeGroup, exp, isotopes, correctionFactors.get(isName), respectDilution, standMethod, corrStandard);
            if (standResults.get(isName).get(partner)!=null){
              float areaPartner = (float)this.correctAreaCorrespondingly(standResults.get(isName).get(partner).getTotalArea(isotopes),moleculeGroup,partner,isotopes,correctionFactors.get(isName), respectDilution, standMethod, corrStandard);
              if (areaExp>0 && areaPartner>0)
                ratioOfIs.add(areaExp/areaPartner);
            }
          }
        }
      }
    }
    Float[] medians = ComparativeAnalysis.medianPlusUpperLowerValues(ratioOfIs);
    float[] lowerUpperThreshold = ComparativeAnalysis.getLowerUpperOutlierThresholds(medians,2f);
    double lowerThreshold = lowerUpperThreshold[0];
    double upperThreshold = lowerUpperThreshold[1];
    Vector<Float> ratioOfIs2 = new Vector<Float>();
    // this removes the outlier from the statistics
    for (Float value : ratioOfIs){
      if (lowerThreshold<=value&&value<=upperThreshold)
        ratioOfIs2.add(value);
    }

    if (ratioOfIs2.size()>0)
      return Calculator.median(ratioOfIs2);
    else
      return Double.NaN;
  }
  
  private boolean doesPartnerFit(String expName,String expEverythingHasToBeComparedTo,Hashtable<String,Hashtable<String,Boolean>> standsForGroup){
    if (expName.equalsIgnoreCase(expEverythingHasToBeComparedTo))
      return true;
    else{
      boolean foundCompareInExpEverythingHasToBeComparedTo = false;
      for (String isName :standsForGroup.get(expName).keySet()){
        if (standsForGroup.get(expName).get(isName)&&standsForGroup.containsKey(expEverythingHasToBeComparedTo)&&standsForGroup.get(expEverythingHasToBeComparedTo).get(isName)){
          foundCompareInExpEverythingHasToBeComparedTo = true;
        }
      }
      if (foundCompareInExpEverythingHasToBeComparedTo){
        return true;
      }else{
        System.out.println("Must look for a different parnter");
        return false;
      }
    }
  }
  
  private float calculateRatioToMedianOneDir(float value, float mean){
    float ratio = value/mean;
    if (ratio==0)
      return Float.MAX_VALUE;
    if (ratio<1)
      ratio = 1/ratio;
    return ratio;
  }
  
  protected void parseResultFile(File resultFile, String fileName) throws ExcelInputFileException, LipidCombinameEncodingException{
    Hashtable<String,Vector<LipidParameterSet>> results = new Hashtable<String,Vector<LipidParameterSet>>();
    Hashtable<String,Boolean> showMods = new Hashtable<String,Boolean>();
    QuantificationResult quantRes = LDAResultReader.readResultFile(resultFile.getAbsolutePath(), showMods);
    int nrChains;
    if (quantRes.getLcbHydroxyEncoding()!=null)
      this.lcbHydroxyEncoding_ = quantRes.getLcbHydroxyEncoding();
    if (quantRes.getFaHydroxyEncoding()!=null)
      this.faHydroxyEncoding_ = quantRes.getFaHydroxyEncoding();  
    results = quantRes.getIdentifications();

    Hashtable<String,Vector<Hashtable<String,ResultAreaVO>>> areaSheetVOs = new Hashtable<String,Vector<Hashtable<String,ResultAreaVO>>>();    
    ElementConfigParser elementParser = Settings.getElementParser();
    try{
      if (LipidomicsConstants.isotopicCorrection()) results = IsotopeCorrector.correctIsotopicPattern(elementParser_,results);
      for (String sheetName : results.keySet()){
        if (!chainsOfClass_.containsKey(sheetName))
          chainsOfClass_.put(sheetName, 0);
        Hashtable<String,Hashtable<String,ResultAreaVO>> areaVOs = new Hashtable<String,Hashtable<String,ResultAreaVO>>();
        Hashtable<String,String> modifications = new Hashtable<String,String>();
        Vector<String> orderedAnalytes = new Vector<String>();

        if (modifications_.containsKey(sheetName)) modifications = modifications_.get(sheetName);
        Vector<LipidParameterSet> params = results.get(sheetName);
        float highestZeroIsoArea = 0f;
        String recentModification = null;
        for (LipidParameterSet param: params){
          String rtDef = "";
          String analId = StaticUtils.generateLipidNameString(param.getName(), param.getDoubleBonds());
          boolean isInternalStandard = (isSelectionPrefix_!=null && param.getName().startsWith(isSelectionPrefix_));
          boolean isExternalStandard = (esSelectionPrefix_!=null && param.getName().startsWith(esSelectionPrefix_));;
          double retentionTime = -1d;
          if (param.getRt()!=null&&param.getRt().length()>0 ){
            retentionTime = new Double(param.getRt());
            if (expRtGroupingTime_>0 && !(isInternalStandard||isExternalStandard))rtDef = param.getRt();
          }
          
          //calculating the neutral mass
          double neutralMass = (double)param.Mz[0];
          double massDiffToNeutral = 0d;
          Hashtable<String,Integer> elements = StaticUtils.categorizeFormula(param.getModificationFormula());
          for (String element : elements.keySet()){
            massDiffToNeutral += elementParser.getElementDetails(element).getMonoMass()*elements.get(element);
          }
          
          neutralMass -= massDiffToNeutral;
          
        //TODO: this was when I read directly from the Excel file; might be useful
//        if (areaVO!=null){
//          if (!(retentionTime<0)){
//            areaVO.setRetentionTime(recentModification,retentionTime);
//            areaVO.setMoreThanOnePeak(recentModification,moreThanOnePeak);
//          }
//          Hashtable<String,ResultAreaVO> sameMoleculeDiffRet = new Hashtable<String,ResultAreaVO>();
//          if (areaVOs.containsKey(areaVO.getMoleculeNameWoRT())) sameMoleculeDiffRet = areaVOs.get(areaVO.getMoleculeNameWoRT());
//          sameMoleculeDiffRet.put(areaVO.getRt(),areaVO);
//          areaVOs.put(areaVO.getMoleculeNameWoRT(),sameMoleculeDiffRet);
//        }
          String formula = enterSpacesToFormula(param.getAnalyteFormula());
          ResultAreaVO areaVO = null;

          if (areaVOs.containsKey(analId) && param.getArea()>0){
            Hashtable<String,ResultAreaVO> sameMoleculeDiffRet = areaVOs.get(analId);
            if (sameMoleculeDiffRet.size()==1 && sameMoleculeDiffRet.keySet().iterator().next().equalsIgnoreCase(""))
              areaVO = sameMoleculeDiffRet.values().iterator().next();
            else if (rtDef!=null&&rtDef.length()>0){
              areaVO = hasAreaSameRt(rtDef,sameMoleculeDiffRet);
              if (areaVO==null) areaVO = new ResultAreaVO(param.getName(),param.getDoubleBonds(),rtDef,fileName,formula,param.getPercentalSplit(),neutralMass,
                  isInternalStandard,isExternalStandard);
              // Juergen: I am not sure if this "if" and the setting of the retention time is required; for various charge states it seems to be counterproductive
              else /*if (areaVO.hasModification(param.getModificationName()))*/{
                sameMoleculeDiffRet.remove(areaVO.getRt());
              }
            }  
          }else{
            if (param.getArea()>0f){
              areaVO = new ResultAreaVO(param.getName(),param.getDoubleBonds(),rtDef,fileName,formula,param.getPercentalSplit(),neutralMass,
                  isInternalStandard,isExternalStandard);
              if (orderedAnalytes.size()==0||!orderedAnalytes.get(orderedAnalytes.size()-1).equalsIgnoreCase(analId))
                orderedAnalytes.add(analId);
            }else{
              Hashtable<String,Hashtable<String,Boolean>> fileHash = new Hashtable<String,Hashtable<String,Boolean>>();
              if (this.isNullResult_.containsKey(fileName)) fileHash = isNullResult_.get(fileName);
              Hashtable<String,Boolean> sheetHash = new Hashtable<String,Boolean>();
              if (fileHash.containsKey(sheetName)) sheetHash = fileHash.get(sheetName);
              ResultAreaVO dummyVO = new ResultAreaVO(param.getName(),param.getDoubleBonds(),rtDef,fileName,formula,param.getPercentalSplit(),neutralMass,
                  isInternalStandard,isExternalStandard);
              sheetHash.put(dummyVO.getMoleculeName(), true);
              fileHash.put(sheetName, sheetHash);
              isNullResult_.put(fileName,fileHash);
            }
          }
          if (areaVO!=null){
            
            //some init variables
            highestZeroIsoArea = areaVO.getHighestZeroIsoArea(param.getModificationName());
            recentModification = param.getModificationName();
            Vector<CgProbe> zeroIsos = param.getIsotopicProbes().get(0);
            float sumMass = 0f;
            for (CgProbe probe : zeroIsos){
              sumMass += probe.Mz;
            }
            areaVO.addResultPart(recentModification, param.getModificationFormula(), param.Mz[0], sumMass/((double)zeroIsos.size()),
                param.getCharge(),param.getRt());
            Hashtable<Integer,Boolean> moreThanOnePeak = new Hashtable<Integer,Boolean>();
            modifications.put(recentModification, recentModification);
            
            for (Vector<CgProbe> probes : param.getIsotopicProbes()){
              int isotope = 0;
              if (probes.size()>0) isotope = probes.get(0).isotopeNumber;
              for (CgProbe probe : probes){
                Hashtable<Integer,Boolean> mtp = areaVO.addArea(recentModification,isotope,probe.Area);
                if (mtp!=null) moreThanOnePeak = mtp;
                if (isotope==0 && probe.Area>highestZeroIsoArea){
                  retentionTime = probe.Peak/60f;
                  highestZeroIsoArea = probe.Area;
                  areaVO.setRtOriginal(rtDef);
                }
                int iso = isotope;
                if (iso<0) iso = iso*-1;
                if (moreThanOnePeak.containsKey(iso)){
                  moreThanOnePeak.put(iso, true);
                }else
                  moreThanOnePeak.put(iso, false);
              }
              
            }
            areaVO.setRetentionTime(recentModification,retentionTime);
            areaVO.setMoreThanOnePeak(recentModification,moreThanOnePeak);
            //add the MSn information
            if (param instanceof LipidomicsMSnSet) {
              LipidomicsMSnSet msn = (LipidomicsMSnSet)param;
              areaVO.setMsnEvidence(true);
              //TODO: here I use only the ones where I find chain information!!!!
              if (msn.getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED) {
                Hashtable<String,Double> fullAreas = new Hashtable<String,Double>();
                for (String combi : msn.getChainCombinationRelativeAreas().keySet()) {
                  fullAreas.put(combi, msn.getChainCombinationRelativeAreas().get(combi)*((double)param.Area));
                }
                nrChains = areaVO.addChainInformation(recentModification,fullAreas);
                if (nrChains>this.chainsOfClass_.get(sheetName))
                  this.chainsOfClass_.put(sheetName, nrChains);
              }
               
            }
            Hashtable<String,ResultAreaVO> sameMoleculeDiffRet = new Hashtable<String,ResultAreaVO>();      
            if (areaVOs.containsKey(areaVO.getMoleculeNameWoRT())) sameMoleculeDiffRet = areaVOs.get(areaVO.getMoleculeNameWoRT());
            sameMoleculeDiffRet.put(areaVO.getRt(),areaVO);
            areaVOs.put(areaVO.getMoleculeNameWoRT(),sameMoleculeDiffRet);
          }
        }
        modifications_.put(sheetName, modifications);
        Vector<Hashtable<String,ResultAreaVO>> inCorrectOrder = new Vector<Hashtable<String,ResultAreaVO>>();
        for (String analyte : orderedAnalytes) inCorrectOrder.add(areaVOs.get(analyte));
        areaSheetVOs.put(sheetName, inCorrectOrder);

      }
    } catch (ChemicalFormulaException cfx){
      cfx.printStackTrace();
      throw new ExcelInputFileException(cfx);
    } catch (SpectrummillParserException spx){
      throw new ExcelInputFileException(spx);
    }
    unprocessedResults_.put(fileName, areaSheetVOs);
  }
  
  @SuppressWarnings("unchecked")
  protected void buildResultHashes(){
    if (expRtGroupingTime_>0){
      //first the retention times have to be normalized, for this purpose we have to know all the lipid classes and molecule names
      Hashtable<String,Hashtable<String,String>> groupAndMols = new Hashtable<String,Hashtable<String,String>>();
      Vector<String> fileNames = new Vector<String>();
      //first key: file name; second key: analyte group; third key: molecule name; third key retention time; fifth key: modification
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>> resultsInHash = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>>();
      for (String fileName : unprocessedResults_.keySet()){
        fileNames.add(fileName);
        Hashtable<String,Vector<Hashtable<String,ResultAreaVO>>> fileResult = unprocessedResults_.get(fileName);
        Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> fileResultInHash = new Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>(); 
        for (String groupName : fileResult.keySet()){
          Vector<Hashtable<String,ResultAreaVO>> mols = fileResult.get(groupName);
          Hashtable<String,String> molNames = new Hashtable<String,String>();
          Hashtable<String,Hashtable<String,ResultAreaVO>> molsInHash = new Hashtable<String,Hashtable<String,ResultAreaVO>>();
          if (groupAndMols.containsKey(groupName)) molNames = groupAndMols.get(groupName);
          for (Hashtable<String,ResultAreaVO> diffTime : mols){
            if (diffTime.size()>0){
              String name = diffTime.values().iterator().next().getMoleculeNameWoRT();
              molNames.put(name, name);
              molsInHash.put(name, diffTime);
            }
          }
          fileResultInHash.put(groupName, molsInHash);
          groupAndMols.put(groupName, molNames);
        }
        resultsInHash.put(fileName, fileResultInHash);
      }
      
      
      
      
      
      for (String groupName : groupAndMols.keySet()){
        for (String molName : groupAndMols.get(groupName).keySet()){
          //standards should not run through this process
          if ((isSelectionPrefix_!=null && molName.startsWith(isSelectionPrefix_))||(esSelectionPrefix_!=null && molName.startsWith(esSelectionPrefix_)))
            continue;
          int clusterId = 0;
          //first key: cluster id; second key: experiment name; third key: retention time of added VO; value; the assigned area VOs
          Hashtable<Integer,Hashtable<String,Hashtable<String,ResultAreaVO>>> valuesInClusters = new  Hashtable<Integer,Hashtable<String,Hashtable<String,ResultAreaVO>>>();
          Hashtable<String,Set<String>> usedRts = new Hashtable<String,Set<String>>();
          Hashtable<Integer,Double> rtClusters = new Hashtable<Integer,Double>();
          while (areThereHitsOutsideClusterRange(groupName, molName, fileNames, rtClusters, usedRts, resultsInHash)){
            String[] fileNameAndRt = findStrongestHitOutsideExistingClusters(groupName, molName, fileNames, rtClusters, usedRts, resultsInHash);
            String fileName = fileNameAndRt[0];
            String rt = fileNameAndRt[1];
            addClosestPeaksToCluster(clusterId, fileName, groupName, molName, rt, fileNames, usedRts, valuesInClusters, resultsInHash);
            addAreaWeightedMeanRt(rtClusters,clusterId,valuesInClusters.get(clusterId));
            clusterId++;
          }
          addRemainingPeaksToClosestClusters(rtClusters,valuesInClusters,usedRts,resultsInHash,groupName,molName,fileNames);
          rtClusters = new Hashtable<Integer,Double>();
          for (int cluster : valuesInClusters.keySet()){
            addAreaWeightedMeanRt(rtClusters,cluster,valuesInClusters.get(cluster));
          }
          while (clusterOverlap(rtClusters)){
            int[] overlapIds = detectClosestOverlap(rtClusters);
            uniteTwoClusters(overlapIds[0],overlapIds[1],valuesInClusters);
            rtClusters = new Hashtable<Integer,Double>();
            for (int cluster : valuesInClusters.keySet()){
              addAreaWeightedMeanRt(rtClusters,cluster,valuesInClusters.get(cluster));
            }
          }
          
          //unite the peaks in the unprocessedResults_
          for (int i=0;i!=valuesInClusters.size();i++){
            Hashtable<String,Hashtable<String,ResultAreaVO>> cluster = valuesInClusters.get(i);
            String newRt = Calculator.FormatNumberToString(rtClusters.get(i),2);
            for (String fileName : fileNames){
              if (!cluster.containsKey(fileName))
                continue;
              Hashtable<String,ResultAreaVO> clusterOfFile = cluster.get(fileName);
              //check for the strongest one
              ResultAreaVO highestArea = null;
              for (ResultAreaVO vo : clusterOfFile.values()){
                if (highestArea==null || vo.getTotalArea(Integer.MAX_VALUE)>highestArea.getTotalArea(Integer.MAX_VALUE))
                  highestArea = vo;
              }
              for (ResultAreaVO vo : clusterOfFile.values()){
                if (highestArea.getRt().equalsIgnoreCase(vo.getRt()))
                  continue;
                Vector<Hashtable<String,ResultAreaVO>> resInSequence = unprocessedResults_.get(vo.getExpName()).get(groupName);
                //now remove the vo that has to be united from the unprocessedResults_
                for (int j=0; j!=resInSequence.size(); j++){
                  ResultAreaVO oneVO = resInSequence.get(j).values().iterator().next();
                  if (!oneVO.getMoleculeNameWoRT().equalsIgnoreCase(vo.getMoleculeNameWoRT()))
                    continue;
                  resInSequence.get(j).remove(vo.getRt());
                  break;
                }
                highestArea.combineVOs(vo);
              }
              highestArea.setRt(newRt);
            }
          }          
        }
      }
      
    }
    // now build the conventional Hashes
    for (String fileName : unprocessedResults_.keySet()){
      Hashtable<String,Vector<Hashtable<String,ResultAreaVO>>> areaSheetVOs = unprocessedResults_.get(fileName);
      for (String sheetName : areaSheetVOs.keySet()){
        Vector<Hashtable<String,ResultAreaVO>> molecules = areaSheetVOs.get(sheetName);
        Hashtable<String,ResultAreaVO> molHash = new Hashtable<String,ResultAreaVO>();
        Hashtable<String,Vector<ResultAreaVO>> resultsMoleculeGroup = new Hashtable<String,Vector<ResultAreaVO>>();
        Hashtable<String,Hashtable<String,ResultAreaVO>> resultsMoleculeGroupHash = new Hashtable<String,Hashtable<String,ResultAreaVO>>();
        if (allResults_.containsKey(sheetName)){
          resultsMoleculeGroup = allResults_.get(sheetName);
          resultsMoleculeGroupHash = allResultsHash_.get(sheetName);
        }
        Vector<ResultAreaVO> resMols = new Vector<ResultAreaVO>();
        for (Hashtable<String,ResultAreaVO> mol : molecules){
          Vector<DoubleStringVO> forSort = new Vector<DoubleStringVO>();
          if (expRtGroupingTime_>0 && !mol.keySet().iterator().next().equalsIgnoreCase("")){
            for (String rt : mol.keySet()) forSort.add(new DoubleStringVO(rt,new Double(rt)));
            Collections.sort(forSort,new GeneralComparator("at.tugraz.genome.lda.vos.DoubleStringVO", "getValue", "java.lang.Double"));
            for (DoubleStringVO key : forSort){
              ResultAreaVO vo = mol.get(key.getKey());
              molHash.put(vo.getMoleculeName(), vo);
              resMols.add(vo);
            }
          }else{
            ResultAreaVO vo = mol.get(mol.keySet().iterator().next());
            molHash.put(vo.getMoleculeName(), vo);
            resMols.add(vo);            
          }
        } 
        resultsMoleculeGroup.put(fileName, resMols);
        resultsMoleculeGroupHash.put(fileName, molHash);
        allResults_.put(sheetName,resultsMoleculeGroup);
        allResultsHash_.put(sheetName, resultsMoleculeGroupHash);
      }
    }
    unprocessedResults_ = null;
    if (this.classCutoffs_!=null)
      filterOutSmallIntensities();
  }
  
  private void filterOutSmallIntensities(){
    for (String lipidClass : classCutoffs_.keySet()){
      double cutoff = classCutoffs_.get(lipidClass);
      if (cutoff>0d){
        Hashtable<String,Vector<ResultAreaVO>> oneClass = allResults_.get(lipidClass);
        Hashtable<String,Hashtable<String,ResultAreaVO>> oneClassHash = allResultsHash_.get(lipidClass);
        @SuppressWarnings("rawtypes")
        Vector<Hashtable> results = calculateCutoffAreas(cutoff,oneClass);
        @SuppressWarnings("unchecked")
        Hashtable<String,Double> cutoffAreas = (Hashtable<String,Double>)results.get(0);
        @SuppressWarnings("unchecked")
        Hashtable<String,String> moleculeNames = (Hashtable<String,String>)results.get(1);
        Hashtable<String,String> moleculesToRemove = findMoleculesBelowCutoff(oneClassHash,moleculeNames,cutoffAreas);
        //now remove the small areas
        for (String expName : oneClass.keySet()){
          Vector<ResultAreaVO> areasFiltered = new Vector<ResultAreaVO>();
          for (ResultAreaVO vo : oneClass.get(expName)){
            if (!moleculesToRemove.containsKey(vo.getMoleculeName()))
              areasFiltered.add(vo);
          }
          oneClass.put(expName, areasFiltered);
        }
        allResults_.put(lipidClass, oneClass);
        for (String expName : oneClassHash.keySet()){
          Hashtable<String,ResultAreaVO> areaHash = oneClassHash.get(expName);
          for (String molecule : moleculesToRemove.keySet()){
            if (areaHash.containsKey(molecule))
              areaHash.remove(molecule);
          }
          oneClassHash.put(expName, areaHash);
        }
        allResultsHash_.put(lipidClass, oneClassHash);
      }
    }
  }
  
  @SuppressWarnings("rawtypes")
  private Vector<Hashtable> calculateCutoffAreas(double cutoff,Hashtable<String,Vector<ResultAreaVO>> areas){
    Hashtable<String,Double> cutoffAreas = new Hashtable<String,Double>();
    Hashtable<String,String> moleculeNames = new Hashtable<String,String>();
    for (String expName : areas.keySet()){
      Vector<ResultAreaVO> areasOfExp = areas.get(expName);
      double highestIntensity = 0d;
      for (ResultAreaVO area : areasOfExp){
        // standards are excluded for the calculation of the highest area - and for the cutoff removal
        if ((isSelectionPrefix_==null || !area.getMoleculeName().startsWith(isSelectionPrefix_))&&(esSelectionPrefix_==null || !area.getMoleculeName().startsWith(esSelectionPrefix_))){
          moleculeNames.put(area.getMoleculeName(), area.getMoleculeName());
          double areaValue = area.getTotalArea(maxCutoffIsotope_+1);
          if (areaValue>highestIntensity) highestIntensity = areaValue;
        }
      }
      cutoffAreas.put(expName, highestIntensity*cutoff);
    }
    Vector<Hashtable> results = new Vector<Hashtable>();
    results.add(cutoffAreas);
    results.add(moleculeNames);
    return results;
  }
  
  private Hashtable<String,String> findMoleculesBelowCutoff(Hashtable<String,Hashtable<String,ResultAreaVO>> oneClassHash, Hashtable<String,String> moleculeNames, Hashtable<String,Double> cutoffAreas){
    Hashtable<String,String> belowCutoff = new Hashtable<String,String>();
    for (String molName : moleculeNames.keySet()){
      boolean oneAboveThreshold = false;
      for (String expName : oneClassHash.keySet()){
        if (oneClassHash.get(expName).containsKey(molName)){
          ResultAreaVO areaVO = oneClassHash.get(expName).get(molName);
          if (areaVO.getTotalArea(maxCutoffIsotope_+1)>cutoffAreas.get(expName)){
            oneAboveThreshold = true;
            break;
          }
        }
      }
      if (!oneAboveThreshold)
        belowCutoff.put(molName, molName);
    }
    return belowCutoff;
  }
  
  
  /**
   * checks whether an ResultAreaVO with the same retention time exists in the hash,
   * and if yes, this VO will be returned
   * @param rtDef the retention time string to look for 
   * @param sameMoleculeDiffRet the hash containing the retention time objects
   * @return the ResultAreaVO that has the same retention time, null otherwise
   */
  private ResultAreaVO hasAreaSameRt(String rtDef,Hashtable<String,ResultAreaVO> sameMoleculeDiffRet){
    ResultAreaVO found = null;
    if (rtDef==null||rtDef.length()<1)
      return found;
    for (String refRt : sameMoleculeDiffRet.keySet()){
    	  if (rtDef.equalsIgnoreCase(refRt))
    	    return sameMoleculeDiffRet.get(refRt);
    }
    return found;
  }
  
  public boolean isWithinRtGroupingBoundaries(double rt, double refTime){
    return StaticUtils.isWithinTolerance(expRtGroupingTime_, refTime, rt);
  }
  
  private String enterSpacesToFormula(String chemicalFormula){
    String finalFormula = "";
    char[] chars = chemicalFormula.toCharArray();
    boolean previousNumber = false;
    for (char oneChar: chars){
      if (Character.isLetter(oneChar)){
        if (previousNumber)
          finalFormula+=" ";
      }
      if (Character.isDigit(oneChar))
        previousNumber = true;
      else
        previousNumber = false;
      finalFormula+=String.valueOf(oneChar);
    }
    return finalFormula;
  }
  
  public static String removeChemicalFormula(String name){
    if (name != null && name.contains("_")){
      return name.substring(0,name.lastIndexOf("_"));
    } else if (name!=null && name.contains(":")){
      char[] nameChars = name.toCharArray();
      int idx = name.lastIndexOf(":")+1;
      String returnName = name.substring(0,idx);
      for (int i=idx;i!=nameChars.length;i++){
        if (Character.isDigit(nameChars[i]))
          returnName += String.valueOf(nameChars[i]);
        else
          break;
      }
      return returnName;
    }else
      return name;
  }
  
  private static float[] getLowerUpperOutlierThresholds(Float[] medians, float tolerance){
    float[] thresholds = new float[2];
    float multFactor = medians[0]/medians[1];
    thresholds[0] = medians[0]/(tolerance*multFactor);
    multFactor = medians[2]/medians[0];
    thresholds[1] = medians[0]*multFactor*tolerance;
    return thresholds;
  }
  
  private static Float[] medianPlusUpperLowerValues(Vector<Float> values)
  {
    Collections.sort(values);
    Float[] medians = new Float[3];
    Float median;
    Float lowerSectionMedian;
    Float upperSectionMedian;
    if (values.size() > 0) {
      if (values.size() % 2 == 0) {
        median = new Float(((values.get((values.size() / 2) - 1))
            .floatValue() + (values.get(values.size() / 2)).floatValue()) / 2);
      } else {
        median = new Float(values.get(((values.size() + 1) / 2) - 1));
      }
      int position = values.size() / 4;
      if (position==0){
        if (values.size()>1){
          lowerSectionMedian = new Float(((values.get((position))).floatValue()
              + (values.get(position+1)).floatValue()) / 2);
        }else{
          lowerSectionMedian = new Float(values.get(position));
        }
      }else{
        lowerSectionMedian = new Float(((values.get((position-1))).floatValue()
            + (values.get(position)).floatValue()) / 2);          
      }
      position = (3*values.size())/4;
      if (values.size()>1){
        upperSectionMedian = new Float(((values.get((position-1))).floatValue()
            + (values.get(position)).floatValue()) / 2);
      }else{
        upperSectionMedian = new Float(values.get(position));
      }        
        
    } else {
      median = new Float(0);
      lowerSectionMedian = new Float(0);
      upperSectionMedian = new Float(0);
      
    }
    medians[0] = median;
    medians[1] = lowerSectionMedian;
    medians[2] = upperSectionMedian;
    return medians;
  }
  
  @SuppressWarnings("unchecked")
  private void calculateRelativeComparisonValues() throws LipidCombinameEncodingException{
    allMoleculeNames_ = new Hashtable<String,Vector<String>>();
    Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> comparativeAreas = new Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>();
    Hashtable<String,Hashtable<String,Integer>> comparativeMaxIsotopes_ = new Hashtable<String,Hashtable<String,Integer>>();
    Hashtable<String,Hashtable<String,Hashtable<String,Vector<Double>>>> moleculeMassHash = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<Double>>>>();
    for (String groupName : this.allResults_.keySet()){
      Vector<String> moleculeNames = new Vector<String>();
      Hashtable<String,Hashtable<String,Vector<Double>>> moleculeMassHashGroup = new Hashtable<String,Hashtable<String,Vector<Double>>>();
      Hashtable<String,Vector<ResultAreaVO>> oneGroupResults = this.allResults_.get(groupName);
      Hashtable<String,Hashtable<String,ResultAreaVO>> comparativeAreasOneGroup = new Hashtable<String,Hashtable<String,ResultAreaVO>>();
      Hashtable<String,Integer> maxIsotopes = new Hashtable<String,Integer>();
      boolean foundPreviousForOrder = true;
      //for (String expName : oneGroupResults.keySet()){
      for (String expName : this.expNamesInSequence_){
        if (oneGroupResults.containsKey(expName)){
          Vector<ResultAreaVO> resVOs = oneGroupResults.get(expName);
          for (int i=0;i!=resVOs.size();i++){
            ResultAreaVO resVO = resVOs.get(i);            
            Hashtable<String,Vector<Double>> masses = new Hashtable<String,Vector<Double>>();
            // this if is to get an ordered list of results
            if (!moleculeMassHashGroup.containsKey(resVO.getMoleculeName())){
              if (moleculeNames.size()==0){
                moleculeNames.add(resVO.getMoleculeName());
              }else{
                if (i>0 ){
                  int position = findPositionOfItem(resVOs.get(i-1).getMoleculeName(),moleculeNames);
                  if (position==-1)
                    moleculeNames.add(resVO.getMoleculeName());
                  else
                    moleculeNames.add((position+1),resVO.getMoleculeName());
                }else{
                  if ((i+1)<resVOs.size()){
                    int position = findPositionOfItem(resVOs.get(i+1).getMoleculeName(),moleculeNames);
                    if (position==-1){
                      foundPreviousForOrder = false;
                      moleculeNames.add(resVO.getMoleculeName());
                    }else
                      moleculeNames.add(position,resVO.getMoleculeName());
                  }else  
                    moleculeNames.add(resVO.getMoleculeName());
                }
              }
            }else{
              masses = moleculeMassHashGroup.get(resVO.getMoleculeName());
            }
            masses.put(expName, resVO.getWeightedNeutralMass());           
            moleculeMassHashGroup.put(resVO.getMoleculeName(), masses);
            //here I am grouping the values to calculate the comparative values afterwards
            Hashtable<String,ResultAreaVO> resOfOneMolecule = new Hashtable<String,ResultAreaVO>();
            int isotopes = 1000;
            if (comparativeAreasOneGroup.containsKey(resVO.getMoleculeName())){
              resOfOneMolecule = comparativeAreasOneGroup.get(resVO.getMoleculeName());
              isotopes = maxIsotopes.get(resVO.getMoleculeName());
            }
            // TOTHINK: Here I am not quite sure if it is good to use several isotopes
            int resIsotopes = resVO.getMaxIsotope();
            if (resIsotopes>0 && resVO.getMaxIsotope()<isotopes)
              isotopes = resVO.getMaxIsotope();
            resOfOneMolecule.put(expName, resVO);
            comparativeAreasOneGroup.put(resVO.getMoleculeName(), resOfOneMolecule);
            maxIsotopes.put(resVO.getMoleculeName(), isotopes);
          }
        }
        moleculeMassHash.put(groupName, moleculeMassHashGroup);
      }
      if (this.correctAnalyteSequence_==null){
        if (!foundPreviousForOrder){
          List<String> molNameList = new ArrayList<String>(moleculeNames);
          Collections.sort(molNameList);
          moleculeNames = new Vector<String> (molNameList);
        }
        moleculeNames = sortLipidNames(moleculeNames);
      }
      allMoleculeNames_.put(groupName, moleculeNames);
      comparativeAreas.put(groupName, comparativeAreasOneGroup);
      comparativeMaxIsotopes_.put(groupName, maxIsotopes);
    }
    this.extractPotentialIsotopicLabels();
    // with this the order is brought into the one of the original input file
    if (correctAnalyteSequence_!=null){
      for (String groupName : correctAnalyteSequence_.keySet()){
        Vector<String> analytesInSequence = correctAnalyteSequence_.get(groupName);
        if (!allMoleculeNames_.containsKey(groupName)) continue;
        Vector<String> unorderedSequence = allMoleculeNames_.get(groupName);
        Vector<String> correctOrder = new Vector<String>();        
        for (String analPlain : analytesInSequence){
          Vector<String> analytesInclIsotopicLabels = new Vector<String>();
          analytesInclIsotopicLabels.add(analPlain);
          for (IsotopicLabelVO labelVO : isoLabels_) {
            for (String prefix : labelVO.getPrefixes().keySet())
              analytesInclIsotopicLabels.add(prefix+analPlain);
          }
          for (String anal : analytesInclIsotopicLabels) {
            Vector<DoubleStringVO> sameRt = new Vector<DoubleStringVO>();
            for (String isThereAnal : unorderedSequence){
              if (expRtGroupingTime_>0 && !isThereAnal.startsWith(isSelectionPrefix_) && !isThereAnal.startsWith(esSelectionPrefix_)){
                if (isThereAnal.startsWith(anal)&&isThereAnal.toCharArray()[anal.length()]=='_'){
                  try{
                    String rtString = isThereAnal.substring(anal.length()+1);
                    double rt = Double.parseDouble(rtString);
                    sameRt.add(new DoubleStringVO(rtString,rt));
                  } catch (NumberFormatException nfx){}
                }
              }else{
                if (anal.equalsIgnoreCase(isThereAnal)){
                  correctOrder.add(anal);
                  break;
                }
              }
            }
            if (expRtGroupingTime_>0){
              Collections.sort(sameRt,new GeneralComparator("at.tugraz.genome.lda.vos.DoubleStringVO", "getValue", "java.lang.Double"));
              for (DoubleStringVO rt : sameRt){
                correctOrder.add(anal+"_"+rt.getKey());
              }
            }
          }
        }
        //this is just for cross check if the quant file contains all analytes;
        //analytes that are not in the file are added at the end
        for (String isThereAnal : unorderedSequence){
          boolean found = false;
          for (String anal : correctOrder){
            if (anal.equalsIgnoreCase(isThereAnal)){
              found = true;
              break;
            }
          }
          if (!found){
            System.out.println("Warning: The molecule "+isThereAnal+" is not in the Quant file! Putting it to the end!!!");
            correctOrder.add(isThereAnal);
          }
        }
        allMoleculeNames_.put(groupName,correctOrder);
      }
    }
    
    comparativeRatios_ = new Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>>();
    for (String groupName : this.allResults_.keySet()){
      Hashtable<String,Hashtable<String,ResultAreaVO>> comparativeAreasOneGroup = comparativeAreas.get(groupName);
      Hashtable<String,Integer> isos = comparativeMaxIsotopes_.get(groupName);
      Vector<String> moleculeNames = allMoleculeNames_.get(groupName);
      Hashtable<String,Hashtable<String,ResultCompVO>> compForOneGroup = new Hashtable<String,Hashtable<String,ResultCompVO>>();
      String mostReliableStandard = null;
      String bestExp = null;
      ResultAreaVO areaOfBestStandard = null;
      
      String mostReliableExtStandardNoCorr = null;
      String bestExpESNoCorr = null;
      ResultAreaVO areaOfBestStandardESNoCorr = null;
      String mostReliableExtStandardInternalCorr = null;
      String bestExpESInternalCorr = null;
      ResultAreaVO areaOfBestStandardESInternalCorr = null;
      String mostReliableExtStandardMedianCorr = null;
      String bestExpESMedianCorr = null;
      ResultAreaVO areaOfBestStandardESMedianCorr = null;
      Hashtable<String,String> mostReliableExtStandardSingleCorr = new Hashtable<String,String>();
      Hashtable<String,String> bestExpESSingleCorr = new Hashtable<String,String>();
      Hashtable<String,ResultAreaVO> areaOfBestStandardESSingleCorr = new Hashtable<String,ResultAreaVO>();
      

      
      if (this.allISNames_.get(groupName).size()>0){
        mostReliableStandard = standardsOrderedConcerningReliability_.get(groupName).get(0);
        bestExp = bestExpForStandard_.get(groupName).get(mostReliableStandard);
        areaOfBestStandard = isResults_.get(groupName).get(mostReliableStandard).get(bestExp);
      }
      if (this.allESNames_.get(groupName).size()>0){
        mostReliableExtStandardNoCorr = extstandsOrderedConcerningReliability_.get(groupName).get(0);
        bestExpESNoCorr = bestExpForExtStandard_.get(groupName).get(mostReliableExtStandardNoCorr);
        areaOfBestStandardESNoCorr = esResults_.get(groupName).get(mostReliableExtStandardNoCorr).get(bestExpESNoCorr);
        if (this.isISAvailable() && bestExpForExtStandardISCorr_.get(groupName)!=null){
          mostReliableExtStandardInternalCorr = extstandsISCorrOrderedConcerningReliability_.get(groupName).get(0);
          bestExpESInternalCorr = bestExpForExtStandardISCorr_.get(groupName).get(mostReliableExtStandardInternalCorr);
          areaOfBestStandardESInternalCorr = esResults_.get(groupName).get(mostReliableExtStandardInternalCorr).get(bestExpESInternalCorr);
          mostReliableExtStandardMedianCorr = extstandsMedianCorrOrderedConcerningReliability_.get(groupName).get(0);
          bestExpESMedianCorr = bestExpForExtStandardMedianCorr_.get(groupName).get(mostReliableExtStandardMedianCorr);
          areaOfBestStandardESMedianCorr = esResults_.get(groupName).get(mostReliableExtStandardMedianCorr).get(bestExpESMedianCorr);
          
          for (String isName : allISNames_.get(groupName).keySet()){
            String mostReliable = extstandsSingleCorrOrderedConcerningReliability_.get(groupName).get(isName).get(0);
            if (bestExpForExtStandardSingleCorr_.get(groupName).get(isName)!=null){
              String bestExpSingle = bestExpForExtStandardSingleCorr_.get(groupName).get(isName).get(mostReliable);
              ResultAreaVO areaOfBestSingle = esResults_.get(groupName).get(mostReliable).get(bestExpSingle);
              mostReliableExtStandardSingleCorr.put(isName, mostReliable);
              bestExpESSingleCorr.put(isName, bestExpSingle);
              areaOfBestStandardESSingleCorr.put(isName,areaOfBestSingle);
            }
          }
        }
      }
      for (String molecule : moleculeNames){    
        Hashtable<String,ResultCompVO> relativeValues = new Hashtable<String,ResultCompVO>();
        int isoNr = 0;
        Double esVolumeInternalCorr = null;
        Double esConcentrationInternalCorr = null;
        Double esVolumeMedianCorr = null;
        Double esConcentrationMedianCorr = null;
        Hashtable<String,Vector<Double>> moleculeMasses = new Hashtable<String,Vector<Double>>();
        if (moleculeMassHash.get(groupName).containsKey(molecule)){
          moleculeMasses = moleculeMassHash.get(groupName).get(molecule);
        }
        Hashtable<Integer,Double> esVolumeSingleCorr = new Hashtable<Integer,Double>();
        Hashtable<Integer,Double> esConcentrationSingleCorr = new Hashtable<Integer,Double>();
        Double dilutionFactor = 1d;
        if (absSetting_!=null && this.allESNames_.get(groupName).size()>0){
          Hashtable<String,Hashtable<String,VolumeConcVO>> concVOs = absSetting_.getClassSettings().get(groupName).getEsStandards();
          if (isISAvailable()){
            esVolumeInternalCorr = concVOs.get(mostReliableExtStandardInternalCorr).get(bestExpESInternalCorr).getVolume();
            esConcentrationInternalCorr = concVOs.get(mostReliableExtStandardInternalCorr).get(bestExpESInternalCorr).getConcentration();
            esVolumeMedianCorr = concVOs.get(mostReliableExtStandardMedianCorr).get(bestExpESMedianCorr).getVolume();
            esConcentrationMedianCorr = concVOs.get(mostReliableExtStandardMedianCorr).get(bestExpESMedianCorr).getConcentration();
            for (String isName : allISNames_.get(groupName).keySet()){
              int isType = correctionTypeISLookup_.get(groupName).get(isName);
              esVolumeSingleCorr.put(isType, concVOs.get(mostReliableExtStandardSingleCorr.get(isName)).get(bestExpESSingleCorr.get(isName)).getVolume());
              esConcentrationSingleCorr.put(isType, concVOs.get(mostReliableExtStandardSingleCorr.get(isName)).get(bestExpESSingleCorr.get(isName)).getConcentration());
            }
          }
        }
        int resultType = ResultCompVO.ANALYTE_TYPE;
        if (this.allISNames_.get(groupName).size()>0 && this.allISNames_.get(groupName).containsKey(molecule))
          resultType = ResultCompVO.INTERNAL_STANDARD_TYPE;
        else if (this.allESNames_.get(groupName).size()>0 && this.allESNames_.get(groupName).containsKey(molecule))
          resultType = ResultCompVO.EXTERNAL_STANDARD_TYPE;
        if (comparativeAreasOneGroup.containsKey(molecule)){
          isoNr = isos.get(molecule);
          Hashtable<String,ResultAreaVO> resultsMolecule = comparativeAreasOneGroup.get(molecule);
          Hashtable<String,Vector<Double>> originalAreasHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Hashtable<String,Boolean>>> moreThanOnePeakHash = new Hashtable<String,Vector<Hashtable<String,Boolean>>>();
          
          Hashtable<String,Vector<Double>> correctionInternalISHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> correctionMedianISHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> bestISAreasHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> medianAreasHash = new Hashtable<String,Vector<Double>>();
          
          Hashtable<String,Vector<Double>> correctionInternalESNoISCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> bestESAreasNoISCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> correctionMedianESNoISCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> medianAreasESNoISCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> correctionInternalESISInternalCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> bestESAreasISInternalCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> correctionMedianESISInternalCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> medianAreasESISInternalCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> correctionInternalESISMedianCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> bestESAreasISMedianCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> correctionMedianESISMedianCorrHash = new Hashtable<String,Vector<Double>>();
          Hashtable<String,Vector<Double>> medianAreasESISMedianCorrHash = new Hashtable<String,Vector<Double>>();

          Hashtable<String,Hashtable<Integer,Vector<Double>>> correctionInternalESISSingleCorrHash = new Hashtable<String,Hashtable<Integer,Vector<Double>>>();
          Hashtable<String,Hashtable<Integer,Vector<Double>>> bestESAreasISSingleCorrHash = new Hashtable<String,Hashtable<Integer,Vector<Double>>>();
          Hashtable<String,Hashtable<Integer,Vector<Double>>> correctionMedianESISSingleCorrHash = new Hashtable<String,Hashtable<Integer,Vector<Double>>>();
          Hashtable<String,Hashtable<Integer,Vector<Double>>> medianAreasESISSingleCorrHash = new Hashtable<String,Hashtable<Integer,Vector<Double>>>();

          Hashtable<String,Hashtable<String,Double>> retentionTimes = new Hashtable<String,Hashtable<String,Double>>();
          Hashtable<String,Boolean> modsFound = new Hashtable<String,Boolean>();
          
          for (String expName : this.expNamesInSequence_){
            Vector<Double> originalAreas = new Vector<Double>();
            Vector<Hashtable<String,Boolean>> moreThanOnePeak = new Vector<Hashtable<String,Boolean>>();
            Vector<Double> correctionInternalIS = new Vector<Double>();
            Vector<Double> bestISAreas = new Vector<Double>();
            Vector<Double> correctionMedianIS = new Vector<Double>();
            Vector<Double> medianAreas = new Vector<Double>();
                        
            Vector<Double> correctionInternalESNoISCorr = new Vector<Double>();
            Vector<Double> bestESAreasNoISCorr = new Vector<Double>();
            Vector<Double> correctionMedianESNoISCorr = new Vector<Double>();
            Vector<Double> medianAreasESNoISCorr = new Vector<Double>();
            Vector<Double> correctionInternalESISInternalCorr = new Vector<Double>();
            Vector<Double> bestESAreasISInternalCorr = new Vector<Double>();
            Vector<Double> correctionMedianESISInternalCorr = new Vector<Double>();
            Vector<Double> medianAreasESISInternalCorr = new Vector<Double>();
            Vector<Double> correctionInternalESISMedianCorr = new Vector<Double>();
            Vector<Double> bestESAreasISMedianCorr = new Vector<Double>();
            Vector<Double> correctionMedianESISMedianCorr = new Vector<Double>();
            Vector<Double> medianAreasESISMedianCorr = new Vector<Double>();
            
            Hashtable<Integer,Vector<Double>> correctionInternalESISSingleCorr = new Hashtable<Integer,Vector<Double>>();
            Hashtable<Integer,Vector<Double>> bestESAreasISSingleCorr = new Hashtable<Integer,Vector<Double>>();
            Hashtable<Integer,Vector<Double>> correctionMedianESISSingleCorr = new Hashtable<Integer,Vector<Double>>();
            Hashtable<Integer,Vector<Double>> medianAreasESISSingleCorr = new Hashtable<Integer,Vector<Double>>();
            
            Hashtable<String,Double> retentionTime = new Hashtable<String,Double> ();
            
            for (int i=1;i!=(isoNr+1);i++){
              double correctionFactor = 1d;
              double medianCorrectionFactor = 1d;
              
              double correctionFactorESNoISCorr = 1d;
              double medianCorrectionESNoISCorr = 1d;
              double correctionFactorESISInternalCorr = 0d;
              double medianCorrectionESISInternalCorr = 0d;
              double correctionFactorESISMedianCorr = 0d;
              double medianCorrectionESISMedianCorr = 0d;
              
              if (this.allISNames_.get(groupName).size()>0){
                double bestArea = correctAreaCorrespondingly(areaOfBestStandard.getTotalArea(i),groupName,bestExp,i,correctionFactorsToBestIS_.get(groupName).get(mostReliableStandard),false,ResultCompVO.NO_STANDARD_CORRECTION,null);
                bestISAreas.add(bestArea);
                correctionFactor = bestArea/referenceValues_.get(groupName).get(expName).get(i-1);
                double bestMedianArea = medianOfRatios_.get(groupName).get(bestExp).get(i-1);
                medianAreas.add(bestMedianArea);
                medianCorrectionFactor = bestMedianArea/medianOfRatios_.get(groupName).get(expName).get(i-1);
              }
              correctionInternalIS.add(correctionFactor);
              correctionMedianIS.add(medianCorrectionFactor);

              if (this.allESNames_.get(groupName).size()>0){
                double bestArea = correctAreaCorrespondingly(areaOfBestStandardESNoCorr.getTotalArea(i),groupName,bestExpESNoCorr,i,correctionFactorsToBestES_.get(groupName).get(mostReliableExtStandardNoCorr),true,ResultCompVO.NO_STANDARD_CORRECTION,null);
                bestESAreasNoISCorr.add(bestArea);
                correctionFactorESNoISCorr = bestArea/referenceValuesES_.get(groupName).get(expName).get(i-1);
                double bestMedianArea = medianOfESRatios_.get(groupName).get(bestExpESNoCorr).get(i-1);
                medianAreasESNoISCorr.add(bestMedianArea);
                medianCorrectionESNoISCorr = bestMedianArea/medianOfESRatios_.get(groupName).get(expName).get(i-1);
                if (isISAvailable() && bestExpForExtStandardISCorr_.get(groupName)!=null){
                  bestArea = correctAreaCorrespondingly(areaOfBestStandardESInternalCorr.getTotalArea(i), groupName, bestExpESInternalCorr, i, correctionFactorsToBestESISCorr_.get(groupName).get(mostReliableExtStandardInternalCorr), true, ResultCompVO.STANDARD_CORRECTION_INTERNAL,null);
                  bestESAreasISInternalCorr.add(bestArea);
                  correctionFactorESISInternalCorr = bestArea/referenceValuesESISCorr_.get(groupName).get(expName).get(i-1);
                  bestMedianArea = medianOfESRatiosISCorr_.get(groupName).get(bestExpESInternalCorr).get(i-1);
                  medianAreasESISInternalCorr.add(bestMedianArea);
                  medianCorrectionESISInternalCorr = bestMedianArea/medianOfESRatiosISCorr_.get(groupName).get(expName).get(i-1);
                  
                  bestArea = correctAreaCorrespondingly(areaOfBestStandardESMedianCorr.getTotalArea(i),groupName,bestExpESMedianCorr,i,correctionFactorsToBestESMedianCorr_.get(groupName).get(mostReliableExtStandardMedianCorr),true,ResultCompVO.STANDARD_CORRECTION_MEDIAN,null);
                  bestESAreasISMedianCorr.add(bestArea);
                  correctionFactorESISMedianCorr = bestArea/referenceValuesESMedianCorr_.get(groupName).get(expName).get(i-1);
                  bestMedianArea = medianOfESRatiosMedianCorr_.get(groupName).get(bestExpESMedianCorr).get(i-1);
                  medianAreasESISMedianCorr.add(bestMedianArea);
                  medianCorrectionESISMedianCorr = bestMedianArea/medianOfESRatiosMedianCorr_.get(groupName).get(bestExpESMedianCorr).get(i-1);
                                    
                  for (String isName : allISNames_.get(groupName).keySet()){
                    int isType = correctionTypeISLookup_.get(groupName).get(isName);
                    Vector<Double> correctionInternalES = new Vector<Double>();
                    Vector<Double> areasESISSingleCorr = new Vector<Double>();
                    Vector<Double> correctionMedianES = new Vector<Double>();
                    Vector<Double> medianCorrectionESISSingleCorr = new Vector<Double>();
                    if (bestESAreasISSingleCorr.containsKey(isType)){
                      correctionInternalES = correctionInternalESISSingleCorr.get(isType);
                      areasESISSingleCorr = bestESAreasISSingleCorr.get(isType);
                      correctionMedianES = correctionMedianESISSingleCorr.get(isType);                      
                      medianCorrectionESISSingleCorr = medianAreasESISSingleCorr.get(isType);
                    }
                    if (areaOfBestStandardESSingleCorr.get(isName)!=null){
                      bestArea = correctAreaCorrespondingly(areaOfBestStandardESSingleCorr.get(isName).getTotalArea(i),groupName,bestExpESSingleCorr.get(isName),i,correctionFactorsToBestESSingleCorr_.get(groupName).get(isName).get(mostReliableExtStandardSingleCorr.get(isName)),true,isType,isName);
                      areasESISSingleCorr.add(bestArea);
                      double correctionFactorSingle = bestArea/referenceValuesESSingleCorr_.get(groupName).get(isName).get(expName).get(i-1);
                      correctionInternalES.add(correctionFactorSingle);
                      bestMedianArea = medianOfESRatiosSingleCorr_.get(groupName).get(isName).get(bestExpESSingleCorr.get(isName)).get(i-1);
                      medianCorrectionESISSingleCorr.add(bestMedianArea);
                      double correctionMedianSingle = bestMedianArea/medianOfESRatiosSingleCorr_.get(groupName).get(isName).get(expName).get(i-1);
                      correctionMedianES.add(correctionMedianSingle);
                    
                      correctionInternalESISSingleCorr.put(isType, correctionInternalES);
                      bestESAreasISSingleCorr.put(isType,areasESISSingleCorr);
                      correctionMedianESISSingleCorr.put(isType,correctionMedianES);                      
                      medianAreasESISSingleCorr.put(isType,medianCorrectionESISSingleCorr);
                    }
                  }
                }
              }            
              correctionInternalESNoISCorr.add(correctionFactorESNoISCorr);
              correctionMedianESNoISCorr.add(medianCorrectionESNoISCorr);
              correctionInternalESISInternalCorr.add(correctionFactorESISInternalCorr);
              correctionMedianESISInternalCorr.add(medianCorrectionESISInternalCorr);
              correctionInternalESISMedianCorr.add(correctionFactorESISMedianCorr);
              correctionMedianESISMedianCorr.add(medianCorrectionESISMedianCorr);
              

              if (resultsMolecule.containsKey(expName)){
                ResultAreaVO resVO = resultsMolecule.get(expName);
                double areaMeasured = resVO.getTotalArea(i);
                originalAreas.add(areaMeasured);
                moreThanOnePeak.add(resVO.getMoreThanOnePeak(i));
              }else{
                originalAreas.add(0d);
                moreThanOnePeak.add(new Hashtable<String,Boolean>());
              }
            }
            boolean allModsFound = true;
            if (resultsMolecule.containsKey(expName)){
              retentionTime = resultsMolecule.get(expName).getRetentionTimes();
              allModsFound = resultsMolecule.get(expName).containsAllModifications(this.modifications_.get(groupName));
            }
            retentionTimes.put(expName, retentionTime);
            modsFound.put(expName, allModsFound);
            
            originalAreasHash.put(expName, originalAreas);
            moreThanOnePeakHash.put(expName, moreThanOnePeak);
            correctionInternalISHash.put(expName, correctionInternalIS);
            bestISAreasHash.put(expName, bestISAreas);
            correctionMedianISHash.put(expName, correctionMedianIS);
            medianAreasHash.put(expName, medianAreas);
            
            correctionInternalESNoISCorrHash.put(expName, correctionInternalESNoISCorr);
            bestESAreasNoISCorrHash.put(expName, bestESAreasNoISCorr);
            correctionMedianESNoISCorrHash.put(expName, correctionMedianESNoISCorr);
            medianAreasESNoISCorrHash.put(expName, medianAreasESNoISCorr);
            correctionInternalESISInternalCorrHash.put(expName, correctionInternalESISInternalCorr);
            bestESAreasISInternalCorrHash.put(expName, bestESAreasISInternalCorr);
            correctionMedianESISInternalCorrHash.put(expName, correctionMedianESISInternalCorr);
            medianAreasESISInternalCorrHash.put(expName, medianAreasESISInternalCorr);
            correctionInternalESISMedianCorrHash.put(expName, correctionInternalESISMedianCorr);
            bestESAreasISMedianCorrHash.put(expName, bestESAreasISMedianCorr);
            correctionMedianESISMedianCorrHash.put(expName, correctionMedianESISMedianCorr);
            medianAreasESISMedianCorrHash.put(expName, medianAreasESISMedianCorr);
            
            correctionInternalESISSingleCorrHash.put(expName, correctionInternalESISSingleCorr);
            bestESAreasISSingleCorrHash.put(expName,bestESAreasISSingleCorr);
            correctionMedianESISSingleCorrHash.put(expName,correctionMedianESISSingleCorr);
            medianAreasESISSingleCorrHash.put(expName, medianAreasESISSingleCorr);
          }
                    
          for (String expName : this.expNamesInSequence_){
            Double endVolume = null;
            Double probeVolume = null;
            Double sampleWeight = null;
            Double proteinConcentration = null;
            Double neutralLipidConcentration = null;
            Vector<Double> moleculeMass = null;
            if (moleculeMasses.containsKey(expName))
              moleculeMass = moleculeMasses.get(expName);
            if (absSetting_!=null){
              dilutionFactor = absSetting_.getClassSettings().get(groupName).getDilutionFactors().get(expName);
              ProbeVolConcVO concVO = absSetting_.getVolumeSettings().get(expName);
              endVolume = concVO.getEndVolume();
              probeVolume = concVO.getProbeVolume();
              sampleWeight = concVO.getSampleWeight();
              proteinConcentration = concVO.getProteinConc();
              neutralLipidConcentration = concVO.getNeutralLipidConc();
            }
            Hashtable<Integer,Vector<Double>> isSingleCorrection = new Hashtable<Integer,Vector<Double>>();
            Hashtable<Integer,Vector<Double>> isSingleRefAreas = new Hashtable<Integer,Vector<Double>>();
            if (isSingleCorrectiveFactors_!=null && isSingleCorrectiveFactors_.get(groupName)!=null){
              isSingleCorrection = isSingleCorrectiveFactors_.get(groupName).get(expName);
              isSingleRefAreas = isSingleRefAreas_.get(groupName);
            }
            Hashtable<Integer,Vector<Double>> esSingleNoCorr = new Hashtable<Integer,Vector<Double>>();
            Hashtable<Integer,Vector<Double>> esSingleAreaNoCorr = new Hashtable<Integer,Vector<Double>>();
            if (esSingleCorrectiveFactorsNoCorr_!=null && esSingleCorrectiveFactorsNoCorr_.get(groupName)!=null){
              esSingleNoCorr = esSingleCorrectiveFactorsNoCorr_.get(groupName).get(expName);
              esSingleAreaNoCorr = esSingleRefAreasNoCorr_.get(groupName);
            }
            Hashtable<Integer,Vector<Double>> esSingleIntCorr = new Hashtable<Integer,Vector<Double>>();
            Hashtable<Integer,Vector<Double>> esSingleAreaIntCorr = new Hashtable<Integer,Vector<Double>>();
            Hashtable<Integer,Vector<Double>> esSingleMedCorr = new Hashtable<Integer,Vector<Double>>();
            Hashtable<Integer,Vector<Double>> esSingleAreaMedCorr = new Hashtable<Integer,Vector<Double>>();
            Hashtable<Integer,Hashtable<Integer,Vector<Double>>> esSingleSingleCorr = new Hashtable<Integer,Hashtable<Integer,Vector<Double>>>();
            Hashtable<Integer,Hashtable<Integer,Vector<Double>>> esSingleAreaSingleCorr = new Hashtable<Integer,Hashtable<Integer,Vector<Double>>>();
            if (esSingleCorrectiveFactorsIntCorr_!=null && esSingleCorrectiveFactorsIntCorr_.get(groupName)!=null){
              esSingleIntCorr = esSingleCorrectiveFactorsIntCorr_.get(groupName).get(expName);
              esSingleAreaIntCorr = esSingleRefAreasIntCorr_.get(groupName);
              esSingleMedCorr = esSingleCorrectiveFactorsMedCorr_.get(groupName).get(expName);
              esSingleAreaMedCorr = esSingleRefAreasMedCorr_.get(groupName);
              esSingleSingleCorr = esSingleCorrectiveFactorsSingleCorr_.get(groupName).get(expName);
              esSingleAreaSingleCorr = esSingleRefAreasSingleCorr_.get(groupName);
            }
            Vector<Double> originalAreas = originalAreasHash.get(expName);
            Vector<Hashtable<String,Boolean>> moreThanOnePeak =  moreThanOnePeakHash.get(expName);
            boolean existsInFile = false;
            if (allResultsHash_.containsKey(groupName)&&allResultsHash_.get(groupName).containsKey(expName)&&allResultsHash_.get(groupName).get(expName).containsKey(molecule)){
              existsInFile = true;
            }
            boolean isNullInFile = false;
            if (!existsInFile && isNullResult_.containsKey(expName) && isNullResult_.get(expName).containsKey(groupName) && isNullResult_.get(expName).get(groupName).containsKey(molecule)) isNullInFile = true;
            relativeValues.put(expName, new ResultCompVO(existsInFile,isNullInFile,resultType,moleculeMass,retentionTimes.get(expName),expNameToFile_.get(expName).getAbsolutePath(), isoNr,modsFound.get(expName), originalAreas, 
                correctionInternalISHash.get(expName), correctionMedianISHash.get(expName),  isSingleCorrection,
                bestISAreasHash.get(expName), medianAreasHash.get(expName), isSingleRefAreas, 
                correctionInternalESNoISCorrHash.get(expName), correctionMedianESNoISCorrHash.get(expName),esSingleNoCorr,
                correctionInternalESISInternalCorrHash.get(expName), correctionMedianESISInternalCorrHash.get(expName), esSingleIntCorr,
                correctionInternalESISMedianCorrHash.get(expName), correctionMedianESISMedianCorrHash.get(expName),esSingleMedCorr,
                correctionInternalESISSingleCorrHash.get(expName), correctionMedianESISSingleCorrHash.get(expName),esSingleSingleCorr,              
                bestESAreasNoISCorrHash.get(expName),medianAreasESNoISCorrHash.get(expName), esSingleAreaNoCorr,
                bestESAreasISInternalCorrHash.get(expName),medianAreasESISInternalCorrHash.get(expName), esSingleAreaIntCorr,
                bestESAreasISMedianCorrHash.get(expName),medianAreasESISMedianCorrHash.get(expName), esSingleAreaMedCorr,
                bestESAreasISSingleCorrHash.get(expName), medianAreasESISSingleCorrHash.get(expName), esSingleAreaSingleCorr,
                isAmountLookup_.get(groupName), dilutionFactor, esAmountLookup_.get(groupName), esVolumeInternalCorr,
                esConcentrationInternalCorr, esVolumeMedianCorr, esConcentrationMedianCorr, esVolumeSingleCorr, esConcentrationSingleCorr,
                endVolume, probeVolume,sampleWeight, proteinConcentration, neutralLipidConcentration,moreThanOnePeak,absSetting_!=null));

          }
        }else{
          for (String expName : this.expNamesInSequence_)
            relativeValues.put(expName, new ResultCompVO(false,false,resultType,new Vector<Double>(), new Hashtable<String,Double>(),expNameToFile_.get(expName).getAbsolutePath(), isoNr, false, new Vector<Double>(), 
                new Vector<Double>(), new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(), new Vector<Double>(),
                new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(), new Vector<Double>(), new Vector<Double>(), 
                new Hashtable<Integer,Vector<Double>>(), new Vector<Double>(), new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(),
                new Vector<Double>(), new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(), new Hashtable<Integer,Vector<Double>>(), 
                new Hashtable<Integer,Vector<Double>>(), new Hashtable<Integer,Hashtable<Integer,Vector<Double>>>(), new Vector<Double>(), 
                new Vector<Double>(),  new Hashtable<Integer,Vector<Double>>(), new Vector<Double>(), new Vector<Double>(),  
                new Hashtable<Integer,Vector<Double>>(), new Vector<Double>(), new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(),
                new Hashtable<Integer,Vector<Double>>(), new Hashtable<Integer,Vector<Double>>(), new Hashtable<Integer,Hashtable<Integer,Vector<Double>>>(),
                isAmountLookup_.get(groupName), dilutionFactor, esAmountLookup_.get(groupName), esVolumeInternalCorr,
              esConcentrationInternalCorr, esVolumeMedianCorr, esConcentrationMedianCorr,esVolumeSingleCorr,esVolumeSingleCorr,
              null, null, null,null, null, new Vector<Hashtable<String,Boolean>>(),false));          
        }
        compForOneGroup.put(molecule, relativeValues);
      }
      comparativeRatios_.put(groupName, compForOneGroup);
    }
    calculateSumAndHighestExp(expNamesInSequence_,comparativeRatios_);
    // this are the values for the grouped experiments
    if (groups_!=null && groups_.size()>0){
      comparativeRatiosGroups_ = new Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>>();
      for (String molGroup : comparativeRatios_.keySet()){
        Hashtable<String,Hashtable<String,ResultCompVO>> compForOneGroup = comparativeRatios_.get(molGroup);
        Hashtable<String,Hashtable<String,ResultCompVO>> groupedCompForOneGroup = new Hashtable<String,Hashtable<String,ResultCompVO>>();
        for (String molecule : compForOneGroup.keySet()){
          Hashtable<String,ResultCompVO> compForOneMolecule = compForOneGroup.get(molecule);
          Hashtable<String,ResultCompVO> groupedCompForOneMolecule = new Hashtable<String,ResultCompVO>();
          for (String groupName : expNamesOfGroup_.keySet()){
            Hashtable<String,ResultCompVO> participatingExps = new Hashtable<String,ResultCompVO>();
            for (String expName:expNamesOfGroup_.get(groupName)){
              participatingExps.put(expName,compForOneMolecule.get(expName));
            }
            ResultCompGroupVO groupVO = new ResultCompGroupVO(participatingExps);
            groupedCompForOneMolecule.put(groupName, groupVO);
          }
          
          groupedCompForOneGroup.put(molecule, groupedCompForOneMolecule);
        }
        comparativeRatiosGroups_.put(molGroup, groupedCompForOneGroup);
      }
      calculateRatiosGroup(groups_,comparativeRatiosGroups_);
    }
  }
  
  private void calculateSumAndHighestExp(Vector<String> expNames, Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> resultsHash){
    for (String expName : expNames){
      Hashtable<String,Vector<Double>> highestValuesGroup = new Hashtable<String,Vector<Double>>();
      Hashtable<String,Vector<Double>> sumValuesGroup = new Hashtable<String,Vector<Double>>();
      Hashtable<String,Vector<Double>> massValuesGroup = new Hashtable<String,Vector<Double>>();
      Vector<Double> highestTotalValue = new Vector<Double>();
      Vector<Double> totalSumValue = new Vector<Double>();
      int highestIso = 0;
      for (String groupName : resultsHash.keySet()){
        int maxIso = maxIsotopesOfGroup_.get(groupName);
        if (maxIso>highestIso)
          highestIso = maxIso;
        Hashtable<String,Hashtable<String,ResultCompVO>> results = resultsHash.get(groupName);
        Hashtable<Integer,Double> highestValues = new Hashtable<Integer,Double>();
        Hashtable<Integer,Double> sumValues = new Hashtable<Integer,Double>();
        Hashtable<Integer,Double> massValues = new Hashtable<Integer,Double>();
        for (int i=0;i!=maxIso;i++){
          double highestValue = 0;
          double sumValue = 0;
          double weightValue = 0;
          for (String molecule : results.keySet()){
            ResultCompVO resultVO = results.get(molecule).get(expName);          
            if (resultVO.getType() == ResultCompVO.ANALYTE_TYPE){
              double area = resultVO.getOriginalArea(resultVO.getAvailableIsotopeNr(i));
              if (area>highestValue)
                highestValue = area;
              sumValue+=area;
              weightValue+=(area*resultVO.getMass(resultVO.getAvailableIsotopeNr(i)));
            }
          }
          highestValues.put(i,highestValue);
          sumValues.put(i,sumValue);
          massValues.put(i,weightValue);
        }
      
        Vector<Double> highest = new Vector<Double>();
        Vector<Double> sum = new Vector<Double>();
        Vector<Double> mass = new Vector<Double>();
        for (int i=0;i!=maxIso;i++){
          highest.add(highestValues.get(i));
          sum.add(sumValues.get(i));
          mass.add(massValues.get(i));
        }
        highestValuesGroup.put(groupName, highest);
        sumValuesGroup.put(groupName, sum);
        massValuesGroup.put(groupName, mass);
      }
      for (int i=0;i!=highestIso;i++){
        double highest = 0;
        double sum = 0;
        for (String groupName : highestValuesGroup.keySet()){
          Vector<Double> valuesOfGroup = highestValuesGroup.get(groupName);
          double oneHighestValue = 0;
          double aSumValue = 0;
          if (i<valuesOfGroup.size()){
            oneHighestValue = valuesOfGroup.get(i);
            aSumValue = sumValuesGroup.get(groupName).get(i);
          }else if (valuesOfGroup.size()>0){
            oneHighestValue = valuesOfGroup.get(valuesOfGroup.size()-1);
            aSumValue = sumValuesGroup.get(groupName).get(valuesOfGroup.size()-1);
          }  
          if (oneHighestValue>highest)
            highest = oneHighestValue;
          sum += aSumValue;
        }
        highestTotalValue.add(highest);
        totalSumValue.add(sum);
      }
      for (String groupName : resultsHash.keySet()){
        Hashtable<String,Hashtable<String,ResultCompVO>> results = resultsHash.get(groupName);
        for (String molecule : results.keySet()){
          ResultCompVO resultVO = results.get(molecule).get(expName);
          resultVO.setHighestGroupIntensity(highestValuesGroup.get(groupName));
          resultVO.setTotalGroupIntensity(sumValuesGroup.get(groupName));
          resultVO.setTotalGroupMass(massValuesGroup.get(groupName));
          resultVO.setHighestFoundIntensity(highestTotalValue);
          resultVO.setTotalFoundIntensity(totalSumValue);
        }      
      }
    }    
  }
  
  private void calculateRatiosGroup(Vector<String> groupNames,Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> groups){   
    for (String groupName : groupNames){
      int highestIso = 0;
      Hashtable<String,Vector<Double>> sumMeans = new  Hashtable<String,Vector<Double>>();
      Hashtable<String,Vector<Double>> sumSds = new  Hashtable<String,Vector<Double>>();
      Vector<Double> sumTotalMeans = new  Vector<Double>();
      Vector<Double> sumTotalSds = new  Vector<Double>();
      
      Hashtable<Integer,Double> sumSumTotalMean = new Hashtable<Integer,Double>();
      Vector<Hashtable<Integer,Double>> sdsSumTotalMean = new Vector<Hashtable<Integer,Double>>();
      for (String className : groups.keySet()){
        int maxIso = maxIsotopesOfGroup_.get(className);
        if (maxIso>highestIso)
          highestIso = maxIso;
        Hashtable<String,Hashtable<String,ResultCompVO>> groupOfClass = groups.get(className);
        Hashtable<Integer,Double> sumSumMeanHash = new Hashtable<Integer,Double>();
        Vector<Hashtable<Integer,Double>> sdsSumMean = new Vector<Hashtable<Integer,Double>>();
        for (String molName : groupOfClass.keySet()){
          ResultCompGroupVO compVO = (ResultCompGroupVO)groupOfClass.get(molName).get(groupName);
          Hashtable<Integer,Double> sdsSumMeanMolecule = new Hashtable<Integer,Double>();
          Hashtable<Integer,Double> sdsSumTotalMeanMolecule = new Hashtable<Integer,Double>();
          if (compVO.getType() == ResultCompVO.ANALYTE_TYPE){
            for (int i=0;i!=maxIso;i++){
              double relativeSumMean = compVO.getMeanOfRatioToTotalIntensity(compVO.getAvailableIsotopeNr(i));
              double relativeSumSD = compVO.getSDOfRatioToTotalIntensity(compVO.getAvailableIsotopeNr(i));
              double relativeSumTotalMean = compVO.getMeanRatioToOverallGroupsIntensity(compVO.getAvailableIsotopeNr(i));
              double relativeSumTotalSD = compVO.getSDRatioToOverallGroupsIntensity(compVO.getAvailableIsotopeNr(i));
              
              double sumSumValue = 0;
              double sumSumTotalValue = 0;
              if (sumSumMeanHash.containsKey(i)){
                sumSumValue = sumSumMeanHash.get(i);
              }
              if (sumSumTotalMean.containsKey(i)){
                sumSumTotalValue = sumSumTotalMean.get(i);
              }
              
              sumSumValue += relativeSumMean;
              sumSumTotalValue += relativeSumTotalMean;
              sumSumMeanHash.put(i, sumSumValue);
              sumSumTotalMean.put(i, sumSumTotalValue);
              
              sdsSumMeanMolecule.put(i, relativeSumSD);
              sdsSumTotalMeanMolecule.put(i, relativeSumTotalSD);
            }
            sdsSumMean.add(sdsSumMeanMolecule);
            sdsSumTotalMean.add(sdsSumTotalMeanMolecule);
          }
        }
        Vector<Double> sumSumMean = new Vector<Double>();
        Vector<Double> sdSumMean = new Vector<Double>();
        for (int i=0;i!=maxIso;i++){
          sumSumMean.add(i,sumSumMeanHash.get(i));
          Vector<Double> stdevSumIso = new Vector<Double>();
          for (Hashtable<Integer,Double> stdevOfMolecule : sdsSumMean){
            stdevSumIso.add(stdevOfMolecule.get(i));
          }
          double stdevSumErrorPropagated = Calculator.calculateSumStdevErrorPropagated(stdevSumIso);
          sdSumMean.add(stdevSumErrorPropagated);
        }
        sumMeans.put(className, sumSumMean);
        sumSds.put(className, sdSumMean);
      }
      for (int i=0;i!=highestIso;i++){
        double aSumTotalMean = 0;
        if (i<sumSumTotalMean.size())
          aSumTotalMean = sumSumTotalMean.get(i);
        else if (sumSumTotalMean.size()>0)
          aSumTotalMean = sumSumTotalMean.get(sumSumTotalMean.size()-1);
        sumTotalMeans.add(aSumTotalMean);
        Vector<Double> stdevSumIso = new Vector<Double>();
        for (Hashtable<Integer,Double> stdevOfMolecule : sdsSumTotalMean){      
          double aStdevSumIso = 0;
          if (i<stdevOfMolecule.size())
            aSumTotalMean = stdevOfMolecule.get(i);
          else if (stdevOfMolecule.size()>0)
            aSumTotalMean = stdevOfMolecule.get(stdevOfMolecule.size()-1);
          stdevSumIso.add(aStdevSumIso);
        }        
        double stdevSumErrorPropagated = Calculator.calculateSumStdevErrorPropagated(stdevSumIso);
        sumTotalSds.add(stdevSumErrorPropagated);
      }
      
      for (String className : groups.keySet()){
        Hashtable<String,Hashtable<String,ResultCompVO>> groupOfClass = groups.get(className);
        for (String molName : groupOfClass.keySet()){
          ResultCompGroupVO compVO = (ResultCompGroupVO)groupOfClass.get(molName).get(groupName);
          compVO.setSumGroupMeans(sumMeans.get(className));
          compVO.setSumGroupSds(sumSds.get(className));
          compVO.setSumTotalMeans(sumTotalMeans);
          compVO.setSumTotalSds(sumTotalSds);
        }
      }  
    }
  }
  
  protected Vector<String> sortMainGroup(Vector<String> mainGroup){
    return mainGroup;
  }
  
  public Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> getResults(){
    return this.comparativeRatios_;
  }
  
  public Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> getGroupedResults(){
    return this.comparativeRatiosGroups_;
  }

  public Hashtable<String,Vector<String>> getAllMoleculeNames()
  {
    return allMoleculeNames_;
  }
  
  public Hashtable<String,Hashtable<String,String>> getAllISNames()
  {
    return allISNames_;
  }
  
  public Hashtable<String,Hashtable<String,String>> getAllESNames()
  {
    return allESNames_;
  }
  
  /**
   * 
   * @return the original class sequence if "Quant file" is entered, otherwise null
   */
  public LinkedHashMap<String,Integer> getClassSequence(){
    return this.classSequence_;
  }
  
  /**
   * 
   * @return the original quantification objects if "Quant file" is entered, otherwise null; first key is the class, second the analyte name, third, the modification
   */
  public Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> getQuantObjects(){
    return quantObjects_;
  }

  public Hashtable<String,Hashtable<String,Integer>> getCorrectionTypeISLookup()
  {
    return correctionTypeISLookup_;
  }

  public Hashtable<String,Hashtable<String,Integer>> getCorrectionTypeESLookup()
  {
    return correctionTypeESLookup_;
  }

  public int getMaxIsotopesOfGroup(String molGroupName)
  {
    return maxIsotopesOfGroup_.get(molGroupName);
  }
  
  private void calculteESRelativeCorrectionFactors(){
    esCorrectionFactors_ = new Hashtable<String,Hashtable<String,Hashtable<String,Double>>>();
    this.calculateRelativeCorrectionFactors(esResults_,allESNames_,esCorrectionFactors_,false);
  }
  
  private void calculteISRelativeCorrectionFactors(){
    isCorrectionFactors_ = new Hashtable<String,Hashtable<String,Hashtable<String,Double>>>();
    this.calculateRelativeCorrectionFactors(isResults_,allISNames_,isCorrectionFactors_,true);
  }
  
  private void calculateRelativeCorrectionFactors(Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> standResults,
      Hashtable<String,Hashtable<String,String>> standNames, Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactors,
      boolean isIS){
    for (String molGroup : standResults.keySet()){
      if (standNames.get(molGroup).size()>0){
      Hashtable<String,Hashtable<String,ResultAreaVO>> issOfGroup = standResults.get(molGroup);
      Hashtable<String,Hashtable<String,Double>> groupISCorrectionFactors = new Hashtable<String,Hashtable<String,Double>>();
      for (String isName : issOfGroup.keySet()){
        Hashtable<String,Double> relativeCorrectionFactors = new Hashtable<String,Double>();
        if (absSetting_ != null){
          LipidClassSettingVO settingVO = absSetting_.getClassSettings().get(molGroup);
          Hashtable<String,VolumeConcVO> standValues = null;
          if (isIS)
            standValues = settingVO.getIsStandards().get(isName);
          else
            standValues = settingVO.getEsStandards().get(isName);
          Vector<Double> amounts = new Vector<Double>();
          for (VolumeConcVO concVO : standValues.values()) amounts.add(concVO.getVolume()*concVO.getConcentration());
          double median = DoubleCalculator.median(amounts);
          for (String expName : this.expNamesInSequence_){
            VolumeConcVO isVO = standValues.get(expName);
            relativeCorrectionFactors.put(expName,median/(isVO.getVolume()*isVO.getConcentration()));
          }
        }else{
          for (String expName : this.expNamesInSequence_)
            relativeCorrectionFactors.put(expName, 1d);
        }
        groupISCorrectionFactors.put(isName, relativeCorrectionFactors);
      }
      correctionFactors.put(molGroup, groupISCorrectionFactors);
      }
    }    
  }
  
  private void calculateRelativeCorrectiveValuesComparedToBestES(){
    correctionFactorsToBestES_ = new Hashtable<String,Hashtable<String,Hashtable<String,Double>>>();
    calculateRelativeCorrectiveValuesComparedToBestStandard(esResults_,allESNames_,esCorrectionFactors_,extstandsOrderedConcerningReliability_,
        bestExpForExtStandard_,false,correctionFactorsToBestES_);
    boolean isFound = isISAvailable();
    if (isFound){
      correctionFactorsToBestESISCorr_ = new Hashtable<String,Hashtable<String,Hashtable<String,Double>>>();
      calculateRelativeCorrectiveValuesComparedToBestStandard(esResults_,allESNames_,esCorrectionFactors_,extstandsISCorrOrderedConcerningReliability_,
          bestExpForExtStandardISCorr_,false,correctionFactorsToBestESISCorr_);
      correctionFactorsToBestESMedianCorr_ = new Hashtable<String,Hashtable<String,Hashtable<String,Double>>>();
      calculateRelativeCorrectiveValuesComparedToBestStandard(esResults_,allESNames_,esCorrectionFactors_,extstandsMedianCorrOrderedConcerningReliability_,
          bestExpForExtStandardMedianCorr_,false,correctionFactorsToBestESMedianCorr_);
      
      correctionFactorsToBestESSingleCorr_ = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,Double>>>>();
      for (String molGroupName : allResults_.keySet()){
        Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactorsToBestGroup = new Hashtable<String,Hashtable<String,Hashtable<String,Double>>>();
        for (String isName : allISNames_.get(molGroupName).keySet()){
          Hashtable<String,Hashtable<String,Double>> correctionFactors = new Hashtable<String,Hashtable<String,Double>>();
          if (allESNames_.get(molGroupName).size()>0 && bestExpForExtStandardSingleCorr_.get(molGroupName)!=null)
            calculateRelativeCorrectiveValuesComparedToBestStandardGroup(molGroupName,esResults_.get(molGroupName),esCorrectionFactors_.get(molGroupName),
                extstandsSingleCorrOrderedConcerningReliability_.get(molGroupName).get(isName),
                bestExpForExtStandardSingleCorr_.get(molGroupName).get(isName),false,correctionFactors);
          correctionFactorsToBestGroup.put(isName, correctionFactors);
        }
        correctionFactorsToBestESSingleCorr_.put(molGroupName, correctionFactorsToBestGroup);
      }
    }
  }
  
  private void calculateRelativeCorrectiveValuesComparedToBestIS(){
    correctionFactorsToBestIS_ = new Hashtable<String,Hashtable<String,Hashtable<String,Double>>>();
    calculateRelativeCorrectiveValuesComparedToBestStandard(isResults_,allISNames_,isCorrectionFactors_,standardsOrderedConcerningReliability_,
        bestExpForStandard_,true, correctionFactorsToBestIS_);
  }
  
  private void calculateRelativeCorrectiveValuesComparedToBestStandard(Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> standResults,
      Hashtable<String,Hashtable<String,String>> standNames, Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactors, 
      Hashtable<String,Vector<String>> standardsOrderedConcerningReliability, Hashtable<String,Hashtable<String,String>> bestExpForStandard,  boolean isIS, 
      Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactorsToBestStandard){
    for (String groupName : standResults.keySet()){
      if (standNames.get(groupName).size()>0){
        Hashtable<String,Hashtable<String,Double>> groupCorrectionFactors = new Hashtable<String,Hashtable<String,Double>>();
        if (bestExpForStandard.get(groupName)!=null){
          calculateRelativeCorrectiveValuesComparedToBestStandardGroup(groupName,standResults.get(groupName),correctionFactors.get(groupName),
            standardsOrderedConcerningReliability.get(groupName),bestExpForStandard.get(groupName),isIS,groupCorrectionFactors);
          correctionFactorsToBestStandard.put(groupName, groupCorrectionFactors);
        }
      }
    }
  }
 
  private void calculateRelativeCorrectiveValuesComparedToBestStandardGroup(String groupName,Hashtable<String,Hashtable<String,ResultAreaVO>> standResults,
      Hashtable<String,Hashtable<String,Double>> correctionFactors,Vector<String> standardsOrderedConcerningReliability,
      Hashtable<String,String> bestExpForStandard,boolean isIS,Hashtable<String,Hashtable<String,Double>> groupCorrectionFactors){
    String mostReliableStandard = standardsOrderedConcerningReliability.get(0);
    if (absSetting_ != null){
      String bestExp = bestExpForStandard.get(mostReliableStandard);
      Hashtable<String,Hashtable<String,VolumeConcVO>> standValues = null;
      LipidClassSettingVO settingVO = absSetting_.getClassSettings().get(groupName);
      if (isIS)
        standValues = settingVO.getIsStandards();
      else
        standValues = settingVO.getEsStandards();
      double amountBestStand = standValues.get(mostReliableStandard).get(bestExp).getAmount();
      for (String isName : standResults.keySet()){
        Hashtable<String,Double> isCorrectionFactors = new Hashtable<String,Double>();
        for (String expName : this.expNamesInSequence_){
          isCorrectionFactors.put(expName,amountBestStand/standValues.get(isName).get(expName).getAmount());
        }
        groupCorrectionFactors.put(isName, isCorrectionFactors);
      }
    }else{
      groupCorrectionFactors.putAll(correctionFactors);
    }
    
  }
  
  
  
  private void calculateESMedians(){
    medianOfESRatios_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>();
    esSingleRefAreasNoCorr_ = new Hashtable<String,Hashtable<Integer,Vector<Double>>>();
    esSingleCorrectiveFactorsNoCorr_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<Double>>>>();
    correctionTypeESLookup_ = new Hashtable<String,Hashtable<String,Integer>>();
    calculateMedians(esResults_, allESNames_, correctionFactorsToBestES_, extStandStatistics_, applicableStandardsES_,medianOfESRatios_,
        esSingleRefAreasNoCorr_, esSingleCorrectiveFactorsNoCorr_,correctionTypeESLookup_,bestExpForExtStandard_, true, ResultCompVO.NO_STANDARD_CORRECTION);
    if (isISAvailable()){
      medianOfESRatiosISCorr_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>();
      esSingleRefAreasIntCorr_ = new Hashtable<String,Hashtable<Integer,Vector<Double>>>();
      esSingleCorrectiveFactorsIntCorr_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<Double>>>>();
      calculateMedians(esResults_, allESNames_, correctionFactorsToBestESISCorr_, extStandStatisticsISCorrected_, applicableStandardsESISCorr_,medianOfESRatiosISCorr_,
          esSingleRefAreasIntCorr_, esSingleCorrectiveFactorsIntCorr_, correctionTypeESLookup_, bestExpForExtStandardISCorr_,true, ResultCompVO.STANDARD_CORRECTION_INTERNAL);
      
      medianOfESRatiosMedianCorr_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>();
      esSingleRefAreasMedCorr_ = new Hashtable<String,Hashtable<Integer,Vector<Double>>>();
      esSingleCorrectiveFactorsMedCorr_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<Double>>>>();
      calculateMedians(esResults_, allESNames_, correctionFactorsToBestESISCorr_, extStandStatisticsISMedianCorrected_, applicableStandardsESMedianCorr_,medianOfESRatiosMedianCorr_,
          esSingleRefAreasMedCorr_, esSingleCorrectiveFactorsMedCorr_, correctionTypeESLookup_,
          bestExpForExtStandardMedianCorr_, true, ResultCompVO.STANDARD_CORRECTION_MEDIAN);
      medianOfESRatiosSingleCorr_ = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>>();
      esSingleRefAreasSingleCorr_ = new Hashtable<String,Hashtable<Integer,Hashtable<Integer,Vector<Double>>>>();
      esSingleCorrectiveFactorsSingleCorr_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Hashtable<Integer,Vector<Double>>>>>();

      for (String molGroupName : allResults_.keySet()){
        Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> medianOfESRatiosGroup = new Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>();
        Hashtable<Integer,Hashtable<Integer,Vector<Double>>> esSingleRefAreasGroup = new Hashtable<Integer,Hashtable<Integer,Vector<Double>>>();
        Hashtable<String,Hashtable<Integer,Hashtable<Integer,Vector<Double>>>> esSingleCorrectiveFactorsGroup = new Hashtable<String,Hashtable<Integer,Hashtable<Integer,Vector<Double>>>>();
        for (String isName : allISNames_.get(molGroupName).keySet()){
          int isType = correctionTypeISLookup_.get(molGroupName).get(isName).intValue();
          Hashtable<String,Hashtable<Integer,Double>> medianOfRatiosExps = new Hashtable<String,Hashtable<Integer,Double>>();
          Hashtable<Integer,Vector<Double>> groupRefAreas = new Hashtable<Integer,Vector<Double>>();
          Hashtable<String,Hashtable<Integer,Vector<Double>>> groupCorectiveFactors = new Hashtable<String,Hashtable<Integer,Vector<Double>>>();
          if (allESNames_.get(molGroupName).size()>0 && bestExpForExtStandardSingleCorr_.get(molGroupName)!=null && bestExpForExtStandardSingleCorr_.get(molGroupName).get(isName)!=null){

            this.calculateMediansGroup(molGroupName, esResults_.get(molGroupName), correctionFactorsToBestESSingleCorr_.get(molGroupName).get(isName), 
                extStandStatisticsISSingleCorrected_.get(molGroupName).get(isName), applicableStandardsESSingleCorr_.get(molGroupName).get(isName), 
                medianOfRatiosExps, groupRefAreas, groupCorectiveFactors, correctionTypeESLookup_.get(molGroupName), bestExpForExtStandardSingleCorr_.get(molGroupName).get(isName), 
                true, isType, isName);
          }  
          medianOfESRatiosGroup.put(isName, medianOfRatiosExps);
          esSingleRefAreasGroup.put(isType, groupRefAreas);
          for (String expName : groupCorectiveFactors.keySet()){
            Hashtable<Integer,Hashtable<Integer,Vector<Double>>> corrOfStandard = new Hashtable<Integer,Hashtable<Integer,Vector<Double>>>();
            if (esSingleCorrectiveFactorsGroup.containsKey(expName))
              corrOfStandard = esSingleCorrectiveFactorsGroup.get(expName);
            corrOfStandard.put(isType, groupCorectiveFactors.get(expName));           
            esSingleCorrectiveFactorsGroup.put(expName, corrOfStandard);
          }
        }
        medianOfESRatiosSingleCorr_.put(molGroupName, medianOfESRatiosGroup);
        esSingleRefAreasSingleCorr_.put(molGroupName, esSingleRefAreasGroup);
        esSingleCorrectiveFactorsSingleCorr_.put(molGroupName, esSingleCorrectiveFactorsGroup);
      }
    }    
  }
  
  private void calculateISMedians(){
    medianOfRatios_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>>();
    isSingleRefAreas_ = new Hashtable<String,Hashtable<Integer,Vector<Double>>>();
    isSingleCorrectiveFactors_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<Double>>>>();
    correctionTypeISLookup_ = new Hashtable<String,Hashtable<String,Integer>>();
    calculateMedians(isResults_, allISNames_, correctionFactorsToBestIS_, intStandStatistics_, applicableStandards_,medianOfRatios_,
        isSingleRefAreas_, isSingleCorrectiveFactors_,correctionTypeISLookup_,bestExpForStandard_, false, ResultCompVO.NO_STANDARD_CORRECTION);
  }
  
  private void calculateMedians(Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> standResults,
      Hashtable<String,Hashtable<String,String>> standNames, Hashtable<String,Hashtable<String,Hashtable<String,Double>>> correctionFactors,  
      Hashtable<String,Hashtable<String,InternalStandardStatistics>> standardsStatistics, Hashtable<String,Hashtable<String,Hashtable<String,Boolean>>> applicableStandards,
      Hashtable<String,Hashtable<String,Hashtable<Integer,Double>>> medianOfRatios, Hashtable<String,Hashtable<Integer,Vector<Double>>> singleRefAreas,
      Hashtable<String,Hashtable<String,Hashtable<Integer,Vector<Double>>>> singleCorrectiveFactors, 
      Hashtable<String,Hashtable<String,Integer>> correctionTypeLookup,Hashtable<String,Hashtable<String,String>> bestExpForStandard, 
      boolean respectDilution, int standMethod){  
    for (String moleculeGroup : standResults.keySet()){
      if (standNames.get(moleculeGroup).size()>0 && bestExpForStandard.get(moleculeGroup)!=null){
        Hashtable<String,Hashtable<Integer,Double>> medianOfRatiosExps = new Hashtable<String,Hashtable<Integer,Double>>();
        Hashtable<Integer,Vector<Double>> groupRefAreas = new Hashtable<Integer,Vector<Double>>();
        Hashtable<String,Hashtable<Integer,Vector<Double>>> groupCorectiveFactors = new Hashtable<String,Hashtable<Integer,Vector<Double>>>();
        Hashtable<String,Integer> correctionType = new Hashtable<String,Integer>();
        if (correctionTypeLookup.containsKey(moleculeGroup))
          correctionType = correctionTypeLookup.get(moleculeGroup);
        calculateMediansGroup(moleculeGroup,standResults.get(moleculeGroup),correctionFactors.get(moleculeGroup),standardsStatistics.get(moleculeGroup),
            applicableStandards.get(moleculeGroup),medianOfRatiosExps,groupRefAreas,groupCorectiveFactors,correctionType,
            bestExpForStandard.get(moleculeGroup),respectDilution,standMethod,null);      
        medianOfRatios.put(moleculeGroup, medianOfRatiosExps);
        singleRefAreas.put(moleculeGroup, groupRefAreas);
        singleCorrectiveFactors.put(moleculeGroup,groupCorectiveFactors);
        correctionTypeLookup.put(moleculeGroup, correctionType);

      }
    }
  }  
  
  private void calculateMediansGroup(String moleculeGroup, Hashtable<String,Hashtable<String,ResultAreaVO>> standResults, Hashtable<String,Hashtable<String,Double>> correctionFactors,
      Hashtable<String,InternalStandardStatistics> standardsStatistics, Hashtable<String,Hashtable<String,Boolean>> groupApplic,
      Hashtable<String,Hashtable<Integer,Double>> medianOfRatiosExps, Hashtable<Integer,Vector<Double>> groupRefAreas,
      Hashtable<String,Hashtable<Integer,Vector<Double>>> groupCorectiveFactors, Hashtable<String,Integer> correctionType,
      Hashtable<String,String> bestExpForStandard, boolean respectDilution, int standMethod, String corrStandard){
        for (int i=1;i!=(maxIsotopesOfGroup_.get(moleculeGroup)+1);i++){
          for (String expName : expNamesInSequence_){
            Hashtable<String,Boolean> isApplicable = groupApplic.get(expName);
            Vector<Float> standardValues = new  Vector<Float>();
            Hashtable<Integer,Double> medianRatioIsos = new Hashtable<Integer,Double>();
            if (medianOfRatiosExps.containsKey(expName)) medianRatioIsos = medianOfRatiosExps.get(expName);

            for (String standardName : isApplicable.keySet()){

              Hashtable<String,ResultAreaVO> oneIsSeveralExps = standResults.get(standardName);
              if (oneIsSeveralExps.containsKey(expName)&&isApplicable.containsKey(standardName) &&
                isApplicable.get(standardName)){
                int isosToTake = i;
                if (standardsStatistics.get(standardName).getAmountIsotopesFoundInAll()<isosToTake)
                  isosToTake = standardsStatistics.get(standardName).getAmountIsotopesFoundInAll();
                double refValue = this.correctAreaCorrespondingly(oneIsSeveralExps.get(expName).getTotalArea(isosToTake), moleculeGroup, expName, isosToTake, correctionFactors.get(standardName), respectDilution, standMethod, corrStandard);
                if (standardsStatistics.get(standardName).getAmountIsotopesFoundInAll()<i){
                  refValue = correctAreaCorrespondingly(oneIsSeveralExps.get(expName).getTheoreticalIsotopeValue(elementParser_, i),moleculeGroup, expName, 1,
                      correctionFactors.get(standardName), respectDilution, standMethod, corrStandard);

              }
              if (refValue>0)
                standardValues.add((float)refValue);
            }
          }
          Float median = Calculator.median(standardValues);
          if (median!=null&&!median.isNaN()&&!median.isInfinite()){
            medianRatioIsos.put(i-1, (double)median.floatValue());
            medianOfRatiosExps.put(expName, medianRatioIsos);
          }           
        }
      }
      int currentType = ResultCompVO.EXTERNAL_STANDARD_TYPE;
      boolean useExistingTypes = false;
      if (correctionType.size()>0)
        useExistingTypes = true;
      // this is for the calculation of the corrections by single standards
      for (String standName : standResults.keySet()){
        currentType++;
        if (useExistingTypes)
          currentType = correctionType.get(standName);
        String bestExp = bestExpForStandard.get(standName);
        Hashtable<String,ResultAreaVO> oneStandardSeveralExps = standResults.get(standName);
        Vector<Double> refAreas = new Vector<Double>();       
        for (int i=1;i!=(maxIsotopesOfGroup_.get(moleculeGroup)+1);i++){
          int isosToTake = i;
          if (standardsStatistics.get(standName).getAmountIsotopesFoundInAll()<isosToTake)
            isosToTake = standardsStatistics.get(standName).getAmountIsotopesFoundInAll();
          double refValue = this.correctAreaCorrespondingly(oneStandardSeveralExps.get(bestExp).getTotalArea(isosToTake), moleculeGroup, bestExp, isosToTake, correctionFactors.get(standName), respectDilution, standMethod, corrStandard);
          if (standardsStatistics.get(standName).getAmountIsotopesFoundInAll()<i){
            refValue = correctAreaCorrespondingly(oneStandardSeveralExps.get(bestExp).getTheoreticalIsotopeValue(elementParser_, i),moleculeGroup, bestExp, 1,
                correctionFactors.get(standName), respectDilution, standMethod, corrStandard);
          }
          // this correction is necessary since this standard should be now centric and not the best one any more
          refValue = refValue/correctionFactors.get(standName).get(bestExp);
          refAreas.add(refValue);
        }
        for (String expName : expNamesInSequence_){
          if (oneStandardSeveralExps.get(expName) !=null){
          Hashtable<Integer,Vector<Double>> standCorrect = new Hashtable<Integer,Vector<Double>>();
          if (groupCorectiveFactors.containsKey(expName))
            standCorrect = groupCorectiveFactors.get(expName);
          Vector<Double> correctiveFactors = new Vector<Double>();
          for (int i=1;i!=(maxIsotopesOfGroup_.get(moleculeGroup)+1);i++){
            int isosToTake = i;
            if (standardsStatistics.get(standName).getAmountIsotopesFoundInAll()<isosToTake)
              isosToTake = standardsStatistics.get(standName).getAmountIsotopesFoundInAll();
            double expValue = this.correctAreaCorrespondingly(oneStandardSeveralExps.get(expName).getTotalArea(isosToTake), moleculeGroup, expName, isosToTake, correctionFactors.get(standName), respectDilution, standMethod, corrStandard);
            if (standardsStatistics.get(standName).getAmountIsotopesFoundInAll()<i){
              expValue = correctAreaCorrespondingly(oneStandardSeveralExps.get(expName).getTheoreticalIsotopeValue(elementParser_, i),moleculeGroup, expName, 1,
                  correctionFactors.get(standName), respectDilution, standMethod, corrStandard);
            }
            // this correction is necessary since this standard should be now centric and not the best one any more
            expValue = expValue/correctionFactors.get(standName).get(bestExp);
            correctiveFactors.add(refAreas.get(i-1)/expValue);
          }
          standCorrect.put(currentType, correctiveFactors);
          groupCorectiveFactors.put(expName, standCorrect);
        }else{
          Vector<Double> correctiveFactors = new Vector<Double>();
          for (int i=1;i!=(maxIsotopesOfGroup_.get(moleculeGroup)+1);i++)correctiveFactors.add(Double.NaN);
          Hashtable<Integer,Vector<Double>> standCorrect = new Hashtable<Integer,Vector<Double>>();
          if (groupCorectiveFactors.containsKey(expName))
            standCorrect = groupCorectiveFactors.get(expName);
          standCorrect.put(currentType, correctiveFactors);
          groupCorectiveFactors.put(expName, standCorrect);
        }
        }
        groupRefAreas.put(currentType, refAreas);
        correctionType.put(standName, currentType);
      }
      
  }
  
  public Hashtable<String,Boolean> getISAvailability(){
    Hashtable<String,Boolean> availability = new Hashtable<String,Boolean>();
    for (String molGroupName : allResults_.keySet()){
      boolean available = false;
      if (allISNames_.containsKey(molGroupName) && allISNames_.get(molGroupName).size()>0) available = true;
      availability.put(molGroupName, available);
    }
    return availability;
  }
  
  public Hashtable<String,Boolean> getESAvailability(){
    Hashtable<String,Boolean> availability = new Hashtable<String,Boolean>();
    for (String molGroupName : allResults_.keySet()){
      boolean available = false;
      if (allESNames_.containsKey(molGroupName) && allESNames_.get(molGroupName).size()>0) available = true;
      availability.put(molGroupName, available);
    }
    return availability;
  }
  
  public Hashtable<String,Hashtable<String,Double>> getDilutionFactors(){
    Hashtable<String,Hashtable<String,Double>> dilutionFactors = new Hashtable<String,Hashtable<String,Double>>();
    for (String groupName : allResults_.keySet()){
      LipidClassSettingVO setVO = absSetting_.getClassSettings().get(groupName);
      dilutionFactors.put(groupName, setVO.getDilutionFactors());
    }
    
    return dilutionFactors;
  }  

  
  public boolean hasAbsoluteSettings(){
    if (absSetting_!=null)
      return true;
    else
      return false;
  }
  
  private void extractISAmountValues(){
    isAmountLookup_ = new Hashtable<String,Hashtable<Integer,VolumeConcVO>>();
    for (String groupName : this.allResults_.keySet()){
      Hashtable<Integer,VolumeConcVO> lookups = new Hashtable<Integer,VolumeConcVO>();
      if (absSetting_!=null && this.allISNames_.get(groupName).size()>0){
        extractAmountValues(lookups,absSetting_.getClassSettings().get(groupName).getIsStandards(),
            standardsOrderedConcerningReliability_.get(groupName), bestExpForStandard_.get(groupName),
            correctionTypeISLookup_.get(groupName));
      }
      isAmountLookup_.put(groupName, lookups);
    }
  }
  
  private void extractESAmountValues(){
    esAmountLookup_ = new Hashtable<String,Hashtable<Integer,VolumeConcVO>>();
    for (String groupName : this.allResults_.keySet()){
      Hashtable<Integer,VolumeConcVO> lookups = new Hashtable<Integer,VolumeConcVO>();
      if (absSetting_!=null && this.allESNames_.get(groupName).size()>0){
        extractAmountValues(lookups,absSetting_.getClassSettings().get(groupName).getEsStandards(),
            extstandsOrderedConcerningReliability_.get(groupName), bestExpForExtStandard_.get(groupName),
            correctionTypeESLookup_.get(groupName));
      }
      esAmountLookup_.put(groupName, lookups);
    }  
  }
  
  private void extractAmountValues(Hashtable<Integer,VolumeConcVO> lookups, Hashtable<String,Hashtable<String,VolumeConcVO>> concHash,
      Vector<String> standardsOrdered,Hashtable<String,String> bestExps,Hashtable<String,Integer> standToInt){
    String bestStand = standardsOrdered.get(0);
    String bestExp = bestExps.get(bestStand);
    lookups.put(ResultCompVO.STANDARD_CORRECTION_INTERNAL, concHash.get(bestStand).get(bestExp));
    lookups.put(ResultCompVO.STANDARD_CORRECTION_MEDIAN, concHash.get(bestStand).get(bestExp));
    for (String standard : standToInt.keySet()){
      lookups.put(standToInt.get(standard), concHash.get(standard).get(bestExps.get(standard)));
    }
  }

  public Hashtable<String,Vector<String>> getModifications()
  {
    Hashtable<String,Vector<String>> returnHash = new Hashtable<String,Vector<String>>();
    for (String groupName : modifications_.keySet()){
      Vector<String> modsAsVect = new Vector<String>();
      for (String modName : modifications_.get(groupName).keySet()) modsAsVect.add(modName);
      returnHash.put(groupName,modsAsVect);
    }
    return returnHash;
  }
  
  public double getModificationAreaOfAnalyte(String molGroup, String exp, String mol, String modification,int isotopes){
    double area = 0;
    if (this.allResultsHash_.get(molGroup).containsKey(exp)&&this.allResultsHash_.get(molGroup).get(exp).containsKey(mol))
      area = this.allResultsHash_.get(molGroup).get(exp).get(mol).getTotalAreaOfModification(modification, isotopes);
    return area;
  }

  public boolean neglectRtInformation (String analyteName){
    boolean neglectRt = false;
    if (getRtTolerance()==null||(isSelectionPrefix_!=null && analyteName.startsWith(isSelectionPrefix_))||(esSelectionPrefix_!=null && analyteName.startsWith(esSelectionPrefix_))) neglectRt = true;
    return neglectRt;
  }
  
  public Double getRtTolerance()
  {
    if (expRtGroupingTime_>0)
      return expRtGroupingTime_;
    else
      return null;
  }
  
  
  public void cleanup(){
    if (unprocessedResults_!=null)unprocessedResults_.clear();
    unprocessedResults_  = null;
    if (isNullResult_!=null)isNullResult_.clear();
    isNullResult_  = null;
    if (allResults_!=null)allResults_.clear();
    allResults_ = null;
    if (allResultsHash_!=null)allResultsHash_.clear();
    allResultsHash_ = null;
    if (isResults_!=null)isResults_.clear();
    isResults_ = null;
    if (isCorrectionFactors_!=null)isCorrectionFactors_.clear();
    isCorrectionFactors_ = null;
    if (esResults_!=null)esResults_.clear();
    esResults_ = null;
    if (esCorrectionFactors_!=null)esCorrectionFactors_.clear();
    esCorrectionFactors_ = null;
    if (correctionFactorsToBestIS_!=null)correctionFactorsToBestIS_.clear();
    correctionFactorsToBestIS_ = null;
    if (correctionFactorsToBestES_!=null)correctionFactorsToBestES_.clear();
    correctionFactorsToBestES_ = null;
    if (correctionFactorsToBestESISCorr_!=null)correctionFactorsToBestESISCorr_.clear();
    correctionFactorsToBestESISCorr_ = null;
    if (correctionFactorsToBestESMedianCorr_!=null)correctionFactorsToBestESMedianCorr_.clear();
    correctionFactorsToBestESMedianCorr_ = null;
    if (correctionFactorsToBestESSingleCorr_!=null)correctionFactorsToBestESSingleCorr_.clear();
    correctionFactorsToBestESSingleCorr_ = null;
    if (intStandStatistics_!=null)intStandStatistics_.clear();
    intStandStatistics_ = null;
    if (extStandStatistics_!=null)extStandStatistics_.clear();
    extStandStatistics_ = null;
    if (extStandStatisticsISCorrected_!=null)extStandStatisticsISCorrected_.clear();
    extStandStatisticsISCorrected_ = null;
    if (extStandStatisticsISMedianCorrected_!=null)extStandStatisticsISMedianCorrected_.clear();
    extStandStatisticsISMedianCorrected_ = null;
    if (extStandStatisticsISSingleCorrected_!=null)extStandStatisticsISSingleCorrected_.clear();
    extStandStatisticsISSingleCorrected_ = null;
    if (standardsOrderedConcerningReliability_!=null)standardsOrderedConcerningReliability_.clear();
    standardsOrderedConcerningReliability_ = null;
    if (extstandsOrderedConcerningReliability_!=null)extstandsOrderedConcerningReliability_.clear();
    extstandsOrderedConcerningReliability_ = null;
    if (extstandsISCorrOrderedConcerningReliability_!=null)extstandsISCorrOrderedConcerningReliability_.clear();
    extstandsISCorrOrderedConcerningReliability_ = null;
    if (extstandsMedianCorrOrderedConcerningReliability_!=null)extstandsMedianCorrOrderedConcerningReliability_.clear();
    extstandsMedianCorrOrderedConcerningReliability_ = null;
    if (extstandsSingleCorrOrderedConcerningReliability_!=null)extstandsSingleCorrOrderedConcerningReliability_.clear();
    extstandsSingleCorrOrderedConcerningReliability_ = null;
    if (modifications_!=null)modifications_.clear();
    modifications_ = null;
    if (applicableStandards_!=null)applicableStandards_.clear();
    applicableStandards_ = null;
    if (bestExpForStandard_!=null)bestExpForStandard_.clear();
    bestExpForStandard_ = null;
    if (inProperForBestComparableProbe_!=null)inProperForBestComparableProbe_.clear();
    inProperForBestComparableProbe_ = null;
    if (inProperForBestComparableProbeES_!=null)inProperForBestComparableProbeES_.clear();
    inProperForBestComparableProbeES_ = null;
    if (inProperForBestComparableProbeESISCorr_!=null)inProperForBestComparableProbeESISCorr_.clear();
    inProperForBestComparableProbeESISCorr_ = null;
    if (inProperForBestComparableProbeESMedianCorr_!=null)inProperForBestComparableProbeESMedianCorr_.clear();
    inProperForBestComparableProbeESMedianCorr_ = null;
    if (inProperForBestComparableProbeESSingleCorr_!=null)inProperForBestComparableProbeESSingleCorr_.clear();
    inProperForBestComparableProbeESSingleCorr_ = null;
    if (applicableStandardsES_!=null)applicableStandardsES_.clear();
    applicableStandardsES_ = null;
    if (applicableStandardsESISCorr_!=null)applicableStandardsESISCorr_.clear();
    applicableStandardsESISCorr_ = null;
    if (applicableStandardsESMedianCorr_!=null)applicableStandardsESMedianCorr_.clear();
    applicableStandardsESMedianCorr_ = null;
    if (applicableStandardsESSingleCorr_!=null)applicableStandardsESSingleCorr_.clear();
    applicableStandardsESSingleCorr_ = null;
    if (bestExpForExtStandard_!=null)bestExpForExtStandard_.clear();
    bestExpForExtStandard_ = null;
    if (bestExpForExtStandardISCorr_!=null)bestExpForExtStandardISCorr_.clear();
    bestExpForExtStandardISCorr_ = null;
    if (bestExpForExtStandardMedianCorr_!=null)bestExpForExtStandardMedianCorr_.clear();
    bestExpForExtStandardMedianCorr_ = null;
    if (bestExpForExtStandardSingleCorr_!=null)bestExpForExtStandardSingleCorr_.clear();
    bestExpForExtStandardSingleCorr_ = null;
    if (referenceValues_!=null)referenceValues_.clear();
    referenceValues_ = null;
    if (referenceValuesES_!=null)referenceValuesES_.clear();
    referenceValuesES_ = null;
    if (referenceValuesESISCorr_!=null)referenceValuesESISCorr_.clear();
    referenceValuesESISCorr_ = null;
    if (referenceValuesESMedianCorr_!=null)referenceValuesESMedianCorr_.clear();
    referenceValuesESMedianCorr_ = null;
    if (referenceValuesESSingleCorr_!=null)referenceValuesESSingleCorr_.clear();
    referenceValuesESSingleCorr_ = null;
    if (comparativeRatios_!=null)comparativeRatios_.clear();
    comparativeRatios_ = null;
    if (comparativeRatiosGroups_!=null)comparativeRatiosGroups_.clear();
    comparativeRatiosGroups_ = null;
    if (allMoleculeNames_!=null)allMoleculeNames_.clear();
    allMoleculeNames_ = null;   
    if (medianOfRatios_!=null)medianOfRatios_.clear();
    medianOfRatios_ = null;
    if (medianOfESRatios_!=null)medianOfESRatios_.clear();
    medianOfESRatios_ = null;
    if (medianOfESRatiosISCorr_!=null)medianOfESRatiosISCorr_.clear();
    medianOfESRatiosISCorr_ = null;
    if (medianOfESRatiosMedianCorr_!=null)medianOfESRatiosMedianCorr_.clear();
    medianOfESRatiosMedianCorr_ = null;
    if (medianOfESRatiosSingleCorr_!=null)medianOfESRatiosSingleCorr_.clear();
    medianOfESRatiosSingleCorr_ = null;
    if (expNamesOfGroup_!=null)expNamesOfGroup_.clear();
    expNamesOfGroup_ = null;
    if (maxIsotopesOfGroup_!=null)maxIsotopesOfGroup_.clear();
    maxIsotopesOfGroup_ = null;
    //maybe TODO: clean of elementParser
    elementParser_ = null;
    //maybe TODO: clean of absSetting
    absSetting_ = null;
    if (isSingleRefAreas_!=null)isSingleRefAreas_.clear();
    isSingleRefAreas_ = null;
    if (isSingleCorrectiveFactors_!=null)isSingleCorrectiveFactors_.clear();
    isSingleCorrectiveFactors_ = null;
    if (correctionTypeISLookup_!=null)correctionTypeISLookup_.clear();
    correctionTypeISLookup_ = null;
    if (esSingleRefAreasNoCorr_!=null)esSingleRefAreasNoCorr_.clear();
    esSingleRefAreasNoCorr_ = null;
    if (esSingleCorrectiveFactorsNoCorr_!=null)esSingleCorrectiveFactorsNoCorr_.clear();
    esSingleCorrectiveFactorsNoCorr_ = null;
    if (correctionTypeESLookup_!=null)correctionTypeESLookup_.clear();
    correctionTypeESLookup_ = null;
    if (esSingleRefAreasIntCorr_!=null)esSingleRefAreasIntCorr_.clear();
    esSingleRefAreasIntCorr_ = null;
    if (esSingleCorrectiveFactorsIntCorr_!=null)esSingleCorrectiveFactorsIntCorr_.clear();
    esSingleCorrectiveFactorsIntCorr_ = null;
    if (esSingleRefAreasMedCorr_!=null)esSingleRefAreasMedCorr_.clear();
    esSingleRefAreasMedCorr_ = null;
    if (esSingleCorrectiveFactorsMedCorr_!=null)esSingleCorrectiveFactorsMedCorr_.clear();
    esSingleCorrectiveFactorsMedCorr_ = null;
    if (esSingleRefAreasSingleCorr_!=null)esSingleRefAreasSingleCorr_.clear();
    esSingleRefAreasSingleCorr_ = null;
    if (esSingleCorrectiveFactorsSingleCorr_!=null)esSingleCorrectiveFactorsSingleCorr_.clear();
    esSingleCorrectiveFactorsSingleCorr_ = null;
    if (isAmountLookup_!=null)isAmountLookup_.clear();
    isAmountLookup_ = null;
    if (esAmountLookup_!=null)esAmountLookup_.clear();
    esAmountLookup_ = null;
    if (correctAnalyteSequence_!=null)correctAnalyteSequence_.clear();
    correctAnalyteSequence_ = null;
    if (quantObjects_!=null)quantObjects_.clear();
    quantObjects_ = null;
    this.isoLabels_ = null;
  }
  
  public Vector<String> getExpsOfGroup (String group){
    if (expNamesOfGroup_.containsKey(group)) return expNamesOfGroup_.get(group);
    else return new Vector<String>();
  }
  
  public ResultAreaVO getResultAreaVO (String molGroup, String molName, String expName){
    ResultAreaVO areaVO = null;
    if (this.allResultsHash_.containsKey(molGroup) && this.allResultsHash_.get(molGroup).containsKey(expName) && this.allResultsHash_.get(molGroup).get(expName).containsKey(molName))
      areaVO = allResultsHash_.get(molGroup).get(expName).get(molName);
    return areaVO;
  }
  
  
  /**
   * 
   * @return true when the RT grouping parameter is set
   */
  public boolean isRtGrouped(){
    if (this.expRtGroupingTime_>0)
      return true;
    else
      return false;
  }
  
  
  /**
   * checks whether there are still hits that have not been assigned to one of the existing clusters, and are outside the influence region of a cluster
   * @param groupName the analyte class
   * @param molName the analyte name
   * @param fileNames the names of the MS-runs
   * @param rtClusters hash containing the cluster id and weighted retention time center; key: cluster ID; value: retention time center
   * @param usedRts hits whose retention times were already added to clusters; key: file name; set: the retention time strings
   * @param resultsInHash the available results; first key: name of MS-run; second key: analyte group; third key: analyte name; fourth key: retention times of hits 
   * @return true when there are still outside the influence regions of existing clusters
   */
  private boolean areThereHitsOutsideClusterRange(String groupName, String molName, Vector<String> fileNames, Hashtable<Integer,Double> rtClusters,
      Hashtable<String,Set<String>> usedRts, Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>> resultsInHash){
    boolean anyFound = false;
    boolean outsideDetected = false;
    for (String fileName : fileNames){
      if (!resultsInHash.containsKey(fileName) || !resultsInHash.get(fileName).containsKey(groupName) ||
          !resultsInHash.get(fileName).get(groupName).containsKey(molName) ||
          resultsInHash.get(fileName).get(groupName).get(molName).size()==0)
        continue;
      anyFound = true;
      Hashtable<String,ResultAreaVO> results = resultsInHash.get(fileName).get(groupName).get(molName);
      for (String rtHitString : results.keySet()){
        if (outsideDetected) break;
        if (usedRts.containsKey(fileName) && usedRts.get(fileName).contains(rtHitString))
          continue;
        double rtHit = Double.parseDouble(rtHitString);
        boolean insideACluster = false;
        for (double rt : rtClusters.values()){
          if (isWithinRtGroupingBoundaries(rtHit,rt)){
            insideACluster = true;
            break;
          }
        }
        
        if (!insideACluster){
          outsideDetected = true;
        }
      }
      if (outsideDetected) break;
    }
    if (!anyFound)
      return false;
    else
      return outsideDetected;
  }
  
  
  /**
   * detects the strongest detection that is outside the influence regions of existing clusters
   * @param groupName the analyte class
   * @param molName the analyte name
   * @param fileNames the names of the MS-runs
   * @param rtClusters hash containing the cluster id and weighted retention time center; key: cluster ID; value: retention time center
   * @param usedRts hits whose retention times were already added to clusters; key: file name; set: the retention time strings
   * @param resultsInHash the available results; first key: name of MS-run; second key: analyte group; third key: analyte name; fourth key: retention times of hits
   * @return String[0]: fileName; String[1] retention time
   */
  private String[] findStrongestHitOutsideExistingClusters(String groupName, String molName, Vector<String> fileNames, Hashtable<Integer,Double> rtClusters,
      Hashtable<String,Set<String>> usedRts, Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>> resultsInHash){
    String[] fileNameAndRt = new String[2];
    double highestArea = 0d;
    for (String fileName : fileNames){
      if (!resultsInHash.containsKey(fileName) || !resultsInHash.get(fileName).containsKey(groupName) ||
          !resultsInHash.get(fileName).get(groupName).containsKey(molName) ||
          resultsInHash.get(fileName).get(groupName).get(molName).size()==0)
        continue;
      Hashtable<String,ResultAreaVO> results = resultsInHash.get(fileName).get(groupName).get(molName);
      for (String rtHitString : results.keySet()){
        if (usedRts.containsKey(fileName) && usedRts.get(fileName).contains(rtHitString))
          continue;
        double rtHit = Double.parseDouble(rtHitString);
        boolean inCluster = false;
        for (double rt : rtClusters.values()){
          if (isWithinRtGroupingBoundaries(rtHit,rt)){
            inCluster = true;
            break;
          }
        }
        if (inCluster)
          continue;      
        ResultAreaVO result = results.get(rtHitString);
        double area = result.getTotalArea(Integer.MAX_VALUE);
        if (area<highestArea)
          continue;
        highestArea = area;
        fileNameAndRt[0] = fileName;
        fileNameAndRt[1] = rtHitString;
      }
    }
    return fileNameAndRt;
  }
  

  /**
   * this method searches for hits which are inside the influence region of the newly created clusters, and adds them
   * only hits that are inside the region from the "strongest" hit to 1.2 times the median of the closest hits are added
   * @param clusterId an identifier for the cluster in integer format
   * @param strongestFile the name of the file for the strongest hit, i.e. the cluster founder
   * @param groupName analyte class
   * @param molName analyte name
   * @param strongestRt the retention time of the file for the strongest hit, i.e. the cluster founder
   * @param fileNames the names of the MS-runs
   * @param usedRts hits whose retention times were already added to clusters; key: file name; set: the retention time strings
   * @param clusters hash containing the clusters; first key: cluster id; second key: experiment name; third key: retention time of added VO; value; the assigned area VOs
   * @param resultsInHash the available results; first key: name of MS-run; second key: analyte group; third key: analyte name; fourth key: retention times of hits
   */
  private void addClosestPeaksToCluster(int clusterId, String strongestFile, String groupName, String molName,
      String strongestRt, Vector<String> fileNames, Hashtable<String,Set<String>> usedRts,
      Hashtable<Integer,Hashtable<String,Hashtable<String,ResultAreaVO>>> clusters,
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>> resultsInHash){
    ResultAreaVO strongestHit = resultsInHash.get(strongestFile).get(groupName).get(molName).get(strongestRt);
    double strongRt = Double.parseDouble(strongestRt);
    //first key: fileName; second key rt
    Hashtable<String,Hashtable<String,ResultAreaVO>> closest = new Hashtable<String,Hashtable<String,ResultAreaVO>>();
    //get hits within the tolerance and closest to strongest one
    for (String fileName : fileNames){
      if (fileName.equalsIgnoreCase(strongestFile))
        continue;
      if (!resultsInHash.containsKey(fileName) || !resultsInHash.get(fileName).containsKey(groupName) ||
          !resultsInHash.get(fileName).get(groupName).containsKey(molName) ||
          resultsInHash.get(fileName).get(groupName).get(molName).size()==0)
        continue;
      Hashtable<String,ResultAreaVO> results = resultsInHash.get(fileName).get(groupName).get(molName);
      String closestRt = null;
      for (String rtHitString : results.keySet()){
        if (usedRts.containsKey(fileName) && usedRts.get(fileName).contains(rtHitString))
          continue;
        double rtHit = Double.parseDouble(rtHitString);
        if (!isWithinRtGroupingBoundaries(rtHit,strongRt))
          continue;
        if (closestRt==null || (Math.abs(rtHit-strongRt)<Math.abs(Double.valueOf(closestRt)-strongRt))){
          closestRt = rtHitString;
        }
      }
      if (closestRt!=null){
        Hashtable<String,ResultAreaVO> closestHash = new Hashtable<String,ResultAreaVO>();
        closestHash.put(closestRt, results.get(closestRt));
        closest.put(fileName, closestHash);
      }
    }
    Hashtable<String,Hashtable<String,ResultAreaVO>> cluster = new Hashtable<String,Hashtable<String,ResultAreaVO>>();
    Set<String> rts = new HashSet<String>();
    Hashtable<String,ResultAreaVO> oneFile = new Hashtable<String,ResultAreaVO>();
    oneFile.put(strongestRt,strongestHit);
    rts.add(strongestRt);
    cluster.put(strongestFile, oneFile);
    usedRts.put(strongestFile, rts);
    //calculate the median of the closest peaks
    if (closest.size()>0){
      float[] differenceToStrongest = new float[closest.size()];
      int count = 0;
      for (String fileName : closest.keySet()){
        differenceToStrongest[count] = Math.abs(Float.parseFloat(closest.get(fileName).keySet().iterator().next())-(float)strongRt);
        count++;
      }
      float medianDifference = Calculator.median(differenceToStrongest);
      //add all hits to the cluster which are within the range of two times the median difference
      for (String fileName : closest.keySet()){
        String rt = closest.get(fileName).keySet().iterator().next();
        float diff = Math.abs(Float.parseFloat(rt)-(float)strongRt);
        if (diff>1.2d*medianDifference)
          continue;
        rts = new HashSet<String>();
        oneFile = new Hashtable<String,ResultAreaVO>();
        oneFile.put(rt,closest.get(fileName).get(rt));
        rts.add(rt);
        cluster.put(fileName, oneFile);
        usedRts.put(fileName, rts);
      }
    }
    clusters.put(clusterId, cluster);
  }

  
  /**
   * calculates an mean retention time that is weighted by the peak areas; this value is added to the rtClusters hash table
   * @param rtClusters hash containing the cluster id and weighted retention time center; key: cluster ID; value: area weighted retention time center
   * @param clusterId an identifier for the cluster in integer format
   * @param cluster hash containing the VOs of one cluster;first key key: experiment name; second key: retention time of added VO; value: the assigned area VOs
   */
  private void addAreaWeightedMeanRt(Hashtable<Integer,Double> rtClusters, int clusterId, Hashtable<String,Hashtable<String,ResultAreaVO>> cluster){
    double totalArea = 0d;
    double rtTimesArea = 0d;
    double rt;
    double area;
    for (Hashtable<String,ResultAreaVO> vos : cluster.values()){
      for (ResultAreaVO vo : vos.values()){
        rt = Double.parseDouble(vo.getRt());
        area = vo.getTotalArea(Integer.MAX_VALUE);
        totalArea += area;
        rtTimesArea += area*rt;
      }
    }
    rtClusters.put(clusterId, rtTimesArea/totalArea);
  }
  
  
  /**
   * when there no more peaks outside the influence regions of the clusters, the remaining ones are added to the clusters that are closest
   * @param rtClusters hash containing the cluster id and weighted retention time center; key: cluster ID; value: area weighted retention time center
   * @param valuesInClusters hash containing the clusters; first key: cluster id; second key: experiment name; third key: retention time of added VO; value; the assigned area VOs
   * @param usedRts hits whose retention times were already added to clusters; key: file name; set: the retention time strings
   * @param resultsInHash the available results; first key: name of MS-run; second key: analyte group; third key: analyte name; fourth key: retention times of hits
   * @param groupName analyte class
   * @param molName analyte name
   * @param fileNames the names of the MS-runs
   */
  private void addRemainingPeaksToClosestClusters(Hashtable<Integer,Double> rtClusters,
      Hashtable<Integer,Hashtable<String,Hashtable<String,ResultAreaVO>>> valuesInClusters, Hashtable<String,Set<String>> usedRts,
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>> resultsInHash,
      String groupName, String molName, Vector<String> fileNames){
    double rtCluster;
    double diff;
    for (String fileName : fileNames){
      if (!resultsInHash.containsKey(fileName) || !resultsInHash.get(fileName).containsKey(groupName) ||
          !resultsInHash.get(fileName).get(groupName).containsKey(molName) ||
          resultsInHash.get(fileName).get(groupName).get(molName).size()==0)
        continue;
      Hashtable<String,ResultAreaVO> results = resultsInHash.get(fileName).get(groupName).get(molName);
      for (String rtHitString : results.keySet()){
        if (usedRts.containsKey(fileName) && usedRts.get(fileName).contains(rtHitString))
          continue;
        double rtHit = Double.parseDouble(rtHitString);
        int bestClusterId = -1;
        double smallestDiff  = Double.MAX_VALUE;
        for (Integer clusterId : rtClusters.keySet()){
          rtCluster = rtClusters.get(clusterId);
          diff = Math.abs(rtCluster-rtHit);
          if (diff<smallestDiff){
            smallestDiff = diff;
            bestClusterId = clusterId;
          }
        }
        //add areaVO to hash
        ResultAreaVO areaVO = results.get(rtHitString);
        Hashtable<String,Hashtable<String,ResultAreaVO>> oneCluster = valuesInClusters.get(bestClusterId);
        Hashtable<String,ResultAreaVO> inCluster = new Hashtable<String,ResultAreaVO>();
        if (oneCluster.containsKey(fileName)) inCluster = oneCluster.get(fileName);
        inCluster.put(areaVO.getRt(),areaVO);
        oneCluster.put(fileName, inCluster);
        valuesInClusters.put(bestClusterId,oneCluster);
        //add the used RT
        Set<String> rts = new HashSet<String>();
        if (usedRts.containsKey(fileName))
          rts = usedRts.get(fileName);
        usedRts.put(fileName, rts);
      }
    }  
  }
  
  
  /**
   * checks whether cluster centers are within the other influence regions
   * @param rtClusters hash containing the cluster id and weighted retention time center; key: cluster ID; value: area weighted retention time center
   * @return true when cluster centers are overlapping with one another
   */
  private boolean clusterOverlap(Hashtable<Integer,Double> rtClusters){
    for (int i=0; i!=rtClusters.size(); i++){
      double firstRt = rtClusters.get(i);
      for (int j=(i+1); j<(rtClusters.size()); j++){
        double secondRt = rtClusters.get(j);
        if (isWithinRtGroupingBoundaries(firstRt, secondRt)){
          return true;
        }
      }
    }
    return false;
  }
  
  
  /**
   * returns two cluster ids which have the closest overlap
   * @param rtClusters hash containing the cluster id and weighted retention time center; key: cluster ID; value: area weighted retention time center
   * @return int[0] the cluster id of the first one of the two overlapping clusters; int[1] the cluster id of the seceond one of the two overlapping clusters 
   */
  private int[] detectClosestOverlap(Hashtable<Integer,Double> rtClusters){
    double lowestOverlap = this.expRtGroupingTime_;
    int[] clusterIds = new int[2];
    for (int i=0; i!=rtClusters.size(); i++){
      double firstRt = rtClusters.get(i);
      for (int j=(i+1); j<(rtClusters.size()); j++){
        double secondRt = rtClusters.get(j);
        if (!isWithinRtGroupingBoundaries(firstRt, secondRt))
          continue;
        double overlap = Math.abs(firstRt-secondRt);
        if (overlap>lowestOverlap)
          continue;
        clusterIds[0] = i;
        clusterIds[1] = j;
      }
    }
    return clusterIds;
  }

  
  /**
   * unites two overlapping clusters
   * @param id1 the cluster id of the first one of the two overlapping clusters
   * @param id2 the cluster id of the second one of the two overlapping clusters
   * @param valuesInClusters hash containing the clusters; first key: cluster id; second key: experiment name; third key: retention time of added VO; value; the assigned area VOs
   */
  private void uniteTwoClusters(int id1, int id2, Hashtable<Integer,Hashtable<String,Hashtable<String,ResultAreaVO>>> valuesInClusters){
    int lower = id1;
    int upper = id2;
    if (lower>upper){
      lower = id2;
      upper = id1;
    }
    Hashtable<String,Hashtable<String,ResultAreaVO>> strongerCluster = valuesInClusters.get(lower);
    Hashtable<String,Hashtable<String,ResultAreaVO>> weakerCluster = valuesInClusters.get(upper);
    for (String fileName : weakerCluster.keySet()){
      Hashtable<String,ResultAreaVO> toAdd = new Hashtable<String,ResultAreaVO>();
      if (strongerCluster.containsKey(fileName)) toAdd = strongerCluster.get(fileName);
      Hashtable<String,ResultAreaVO> toBeAdded = weakerCluster.get(fileName);
      for (String rt : toBeAdded.keySet()) toAdd.put(rt,toBeAdded.get(rt));
      strongerCluster.put(fileName, toAdd);
    }
    int count = upper;
    //reorganize the clusterIds
    while ((count+1)<valuesInClusters.size()){
      valuesInClusters.put(count, valuesInClusters.get(count+1));
      count++;
    }
    valuesInClusters.remove(valuesInClusters.size()-1);
  }
  
  protected void disableRtGrouping(){
    this.expRtGroupingTime_ = -1d;
  }
  
  
  @SuppressWarnings("unchecked")
  /**
   * extracts information about potential labels out of the provided data
   */
  private void extractPotentialIsotopicLabels() throws LipidCombinameEncodingException {
    char[] aChars;
    int stop;
    String lab;
    
    //first, check which prefixes are present, and which labels do they cause
    isoLabels_ = new ArrayList<IsotopicLabelVO>();
    Set<String> prefixes = new HashSet<String>(); 
    String pref;
    for(String lClass : this.allMoleculeNames_.keySet()) {
      for (String anal : this.allMoleculeNames_.get(lClass)) {
        aChars = anal.toCharArray();
        if (Character.isDigit(aChars[0]))
          continue;
        if (anal.startsWith(isSelectionPrefix_) && anal.length()>isSelectionPrefix_.length() && Character.isDigit(aChars[isSelectionPrefix_.length()]))
          continue;
        if (anal.startsWith(esSelectionPrefix_) && anal.length()>esSelectionPrefix_.length() && Character.isDigit(aChars[esSelectionPrefix_.length()]))
          continue;
        stop = 0;
        while (!Character.isDigit(aChars[stop]))
          stop++;
        pref = anal.substring(0,stop);
        if (!lcbHydroxyEncoding_.values().contains(pref) && !faHydroxyEncoding_.values().contains(pref))
          prefixes.add(pref);
      }
    }
    Set<String> singleLabels = new HashSet<String>();
    Hashtable<String,String> singleLabelLookup = new Hashtable<String,String>();
    Hashtable<String,Integer> labelToNrOfSingles = new Hashtable<String,Integer>();
    StaticUtils.extractIsoLabelInformation(prefixes, singleLabels, singleLabelLookup, labelToNrOfSingles);
    Hashtable<String,Set<String>> singlesToPrefixes = new Hashtable<String,Set<String>>();
    for (String prefix : singleLabelLookup.keySet()) {
      lab = singleLabelLookup.get(prefix);
      if (!singlesToPrefixes.containsKey(lab))
        singlesToPrefixes.put(lab, new HashSet<String>());
      singlesToPrefixes.get(lab).add(prefix);
    }
    
    String analOnly;
    Hashtable<String,Integer> diffsPerLabel;
    //label information is extracted - now look for partners
    for (String label : singleLabels) {
      Set<String> possPrefixes = singlesToPrefixes.get(label);
      // key: molecule name; value Hashtable containing the difference in the chemical formula: key: element sign; value: difference
      Hashtable<String,Hashtable<String,Integer>> partnerChemFormDiff = new Hashtable<String,Hashtable<String,Integer>>();
      for(String lClass : this.allMoleculeNames_.keySet()) {
        for (String anal : this.allMoleculeNames_.get(lClass)) {
          analOnly = new String(anal);
          if (expRtGroupingTime_>0 && !this.isISorES(lClass,anal))
            analOnly = analOnly.substring(0,analOnly.lastIndexOf("_"));
          aChars = analOnly.toCharArray();
          if (!Character.isDigit(aChars[0]))
            continue;
          for (String prefix : possPrefixes) {
            Hashtable<String,Integer> labelDiffs = checkForIsoPartnerMassDifference(lClass, analOnly, prefix+analOnly);
            if (labelDiffs == null)
              continue;
            boolean ok = true;
            for (String element : labelDiffs.keySet()) {
              if (!(labelDiffs.get(element)%labelToNrOfSingles.get(prefix)==0))
                ok = false;
            }
            if (ok) {
              diffsPerLabel = new Hashtable<String,Integer>();
              for (String element : labelDiffs.keySet()) diffsPerLabel.put(element, labelDiffs.get(element)/labelToNrOfSingles.get(prefix));
              partnerChemFormDiff.put(prefix+analOnly, diffsPerLabel);
            }
          }
        }
      }
      Hashtable<String,Integer> numberOfDifferences = new Hashtable<String,Integer>();
      Hashtable<String,Hashtable<String,Integer>> formulaLookup = new Hashtable<String,Hashtable<String,Integer>>();
      for (String mol : partnerChemFormDiff.keySet()) {
        String formula = StaticUtils.getFormulaInHillNotation(partnerChemFormDiff.get(mol),false);
        if (!numberOfDifferences.containsKey(formula)) {
          numberOfDifferences.put(formula, 0);
          formulaLookup.put(formula, partnerChemFormDiff.get(mol));
        }
        numberOfDifferences.put(formula, numberOfDifferences.get(formula)+1);
      }

      String mostOftenChange = null;
      int highestNumber = 0;
      int totalNumber = 0;
      for (String formula : numberOfDifferences.keySet()) {
        if (numberOfDifferences.get(formula)>highestNumber) {
          mostOftenChange = formula;
          highestNumber = numberOfDifferences.get(formula);
        }
        totalNumber += numberOfDifferences.get(formula);
      }
      //if the same formula difference has been detected in 80% of the cases, this label is accepted
      if (mostOftenChange!=null && numberOfDifferences.get(mostOftenChange)>=((totalNumber*4)/5)) {
        List<String> prefixesOnly = new ArrayList<String>(singlesToPrefixes.get(label));
        Collections.sort(prefixesOnly);
        LinkedHashMap<String,Integer> prefixesAndNr = new LinkedHashMap<String,Integer>();
        for (String prefix : prefixesOnly)
          prefixesAndNr.put(prefix, labelToNrOfSingles.get(prefix));
        isoLabels_.add(new IsotopicLabelVO(label,formulaLookup.get(mostOftenChange),prefixesAndNr));
      }
    }
    Collections.sort(isoLabels_,new GeneralComparator("at.tugraz.genome.lda.vos.IsotopicLabelVO", "getLabelId", "java.lang.String"));
    
    // now extract the best matching retention times
    if (isoLabels_.size()==0 || !(expRtGroupingTime_>0))
      return;
    // the RT differences for each label
    Hashtable<String,Vector<Float>> labelRtDiffs = new Hashtable<String,Vector<Float>>();
    for (IsotopicLabelVO label : isoLabels_)
      labelRtDiffs.put(label.getLabelId(), new Vector<Float>());
    Vector<Object> isoAnal;
    Hashtable<String,Integer> elements = new Hashtable<String,Integer>();
    Hashtable<String,Integer> chainFrequencies;
    int mostOftenDetection;
    int timesFound;
    String name;
    String mostOftenFound;
//    int nrIsos = 0;
    DoubleStringVO isoPartner;
    double diff;
    for (String lClass : allMoleculeNames_.keySet()) {
      Hashtable<String,Vector<ResultAreaVO>> resultsOfClass = this.allResults_.get(lClass);
      Hashtable<String,String> usedAnalytes = new Hashtable<String,String>();
      for (String analWithRt : this.allMoleculeNames_.get(lClass)) {
        String anal = analWithRt.substring(0,analWithRt.lastIndexOf("_"));        
        if (usedAnalytes.containsKey(anal) || anal.startsWith(isSelectionPrefix_) || anal.startsWith(esSelectionPrefix_))
          continue;
        usedAnalytes.put(anal,anal);
        // first, we have to select the analyte only if it is not an isotopic label
        boolean mightBeIsoLabel = false;
        for (IsotopicLabelVO isoLabel : isoLabels_) {
          if (anal.startsWith(isoLabel.getLabelId())) {
            mightBeIsoLabel = true;
            break;
          }
        }
        if (mightBeIsoLabel)
          continue;
        //now, get the RTs of all the unlabeled hits of this analyte in the sample
        List<DoubleStringVO> foundRTs = new ArrayList<DoubleStringVO>();
        for (String otherAnal : this.allMoleculeNames_.get(lClass)) {
          if (otherAnal.startsWith(anal) && otherAnal.lastIndexOf("_")==anal.length()) {
            foundRTs.add(new DoubleStringVO(otherAnal,new Double(otherAnal.substring(otherAnal.lastIndexOf("_")+1))));
          }
        }
        Collections.sort(foundRTs,new GeneralComparator("at.tugraz.genome.lda.vos.DoubleStringVO", "getValue", "java.lang.Double"));
        
        //now get the elemental composition of the unlabeled sample, and read the strongest matching chain information if any
        Hashtable<String,Integer> unlabeledElements = new Hashtable<String,Integer>();
        Hashtable<String,String> mostOftenChainCombination = new Hashtable<String,String>();
        Hashtable<String,Hashtable<String,String>> mostOftenChainCombiInEachFile = new Hashtable<String,Hashtable<String,String>>();
        for (DoubleStringVO vo : foundRTs)
          mostOftenChainCombiInEachFile.put(vo.getKey(), new Hashtable<String,String>());
        for (String file : resultsOfClass.keySet()) {
          Vector<ResultAreaVO> resultsOfFile = resultsOfClass.get(file);
          for (ResultAreaVO areaVO : resultsOfFile) {
            if (!areaVO.getMoleculeNameWoRT().equalsIgnoreCase(anal))
              continue;
            if (unlabeledElements.size()==0)
              unlabeledElements = areaVO.getChemicalFormulaElements();
            if (areaVO.getStrongestChainIdentification()!=null)
              mostOftenChainCombiInEachFile.get(areaVO.getMoleculeName()).put(file, areaVO.getStrongestChainIdentification());
          }
        }
        // calculate mostOftenChainIdentification
        for (String analWRT : mostOftenChainCombiInEachFile.keySet()) {
          chainFrequencies = new Hashtable<String,Integer>();
          mostOftenDetection = 0;
          Hashtable<String,String> eachFile = mostOftenChainCombiInEachFile.get(analWRT);
          mostOftenFound = "";
          for (String chainCombi : eachFile.values()) {
            timesFound = 0;
            name = chainCombi;
            if (chainFrequencies.containsKey(name))
              timesFound = chainFrequencies.get(name);
            else {
              for (String otherCombi : chainFrequencies.keySet()) {
                if (StaticUtils.isAPermutedVersion(name, otherCombi, LipidomicsConstants.CHAIN_COMBI_SEPARATOR)) {
                  name = otherCombi;
                  timesFound = chainFrequencies.get(name);
                }
              }
            }
            timesFound++;
            chainFrequencies.put(name, timesFound);
            if (timesFound>mostOftenDetection) {
              mostOftenDetection = timesFound;
              mostOftenFound = name;
            }
          }
          if (mostOftenDetection>0)
            mostOftenChainCombination.put(analWRT, mostOftenFound);
        }
        
        //now, search for the labeled species and check which one comes closest to the unlabeled ones
        for (IsotopicLabelVO label : isoLabels_) {
          //now, we use the original data to get an more accurate comparison
          Hashtable<String,Integer> labeledSpecies = new Hashtable<String,Integer>();
          Hashtable<String,Hashtable<String,Integer>> labeledSpeciesElements = new Hashtable<String,Hashtable<String,Integer>>();
          for (String prefix : label.getPrefixes().keySet()) {
            labeledSpecies.put(prefix+anal, label.getPrefixes().get(prefix));
            elements = new Hashtable<String,Integer>(unlabeledElements);
            addLabelElements(elements,label.getLabelElements(),label.getPrefixes().get(prefix));
            labeledSpeciesElements.put(prefix+anal, elements);
          }
          for (String file : resultsOfClass.keySet()) {
            Vector<ResultAreaVO> resultsOfFile = resultsOfClass.get(file);
            for (ResultAreaVO areaVO : resultsOfFile) {
              if (!labeledSpecies.containsKey(areaVO.getMoleculeNameWoRT()))
                continue;
              if (!StaticUtils.isChemicalFormulaTheSame(areaVO.getChemicalFormulaElements(),labeledSpeciesElements.get(areaVO.getMoleculeNameWoRT())))
                continue;
              double isoRt = Double.parseDouble(areaVO.getRt());
              Vector<DoubleStringVO> closestInformation = findClosestPartners(isoRt,foundRTs);
              String[] chainIdentifications = getChainIdentifications(closestInformation,file,mostOftenChainCombiInEachFile,mostOftenChainCombination);
//              if (!lClass.equalsIgnoreCase("PI"))
//                continue;
//              System.out.println(areaVO.getMoleculeName());
              isoAnal = checkChainInformationOrSelectClosestPartner(label, labeledSpecies.get(areaVO.getMoleculeNameWoRT()), closestInformation, areaVO.getStrongestChainIdentification(), chainIdentifications);
              if (isoAnal.size()==0)
                continue;
              isoPartner = (DoubleStringVO)isoAnal.get(0);
              diff = isoPartner.getValue();
              if (closestInformation.get(1).getKey()!=null && closestInformation.get(1).getKey().equalsIgnoreCase(isoPartner.getKey()))
                diff = diff*-1d;
              diff = diff/((double)labeledSpecies.get(areaVO.getMoleculeNameWoRT()));            
              if ((Boolean)isoAnal.get(1))
//                System.out.println("diff: "+diff);
                labelRtDiffs.get(label.getLabelId()).add((float)diff);
              
              
//              The next step is to look for the closest identifications in the one or the other direction. If there is information about molecular species -> the one is taken that matches, if both match the closer one shall be selected
//              Then the RT difference is written into the list for each isotopic label (if the molecular species matches, it is written twice, otherwise only once)
              
              
              //if (lClass.equalsIgnoreCase("PI")) {
                //System.out.println(areaVO.getMoleculeName()+": "+(closestPartners[0]!=null ? closestPartners[0] : "")+" "+(closestPartners[1]!=null ? closestPartners[1] : ""));
                //System.out.println(areaVO.getMoleculeName()+" "+(areaVO.getStrongestChainIdentification()!=null?areaVO.getStrongestChainIdentification():"")+": "+(chainIdentifications[0]!=null ? chainIdentifications[0] : "")+" "+(chainIdentifications[1]!=null ? chainIdentifications[1] : ""));
 //               nrIsos++;
              //}  
            }
          }
        }
        
//      if (lClass.equalsIgnoreCase("PI"))
//      System.out.println(anal+": "+otherAnal);

        
//        if (lClass.equalsIgnoreCase("PI"))
//          System.out.println("anal: "+anal);
      }
    }
    //calculate the median of the retention time and add it to the label
    for (IsotopicLabelVO label : isoLabels_) {
      if (labelRtDiffs.get(label.getLabelId()).size()>0) {
        List<Float> values = new ArrayList<Float>(labelRtDiffs.get(label.getLabelId()));
        Collections.sort(values);
        label.setRtShift(Calculator.median(labelRtDiffs.get(label.getLabelId())));
        //System.out.println(label.getLabelId()+": "+label.getRtShift()+"     "+labelRtDiffs.get(label.getLabelId()).size());
      }
    }
    //some omega-positions are not that frequent in nature - thus, the prediction should rather be based on the more frequent ones
    //a prediction is done on the number of labels, where linearity is assumed (a single isotopic label causes the same shift)
    //first, group together the labels that are caused by the same isotope
    Hashtable<String,Vector<IsotopicLabelVO>> labelsOfSameIsotope = new Hashtable<String,Vector<IsotopicLabelVO>>();
    List<String> isotopes;
    String isotopeId;
    //first, group together the labels that are caused by the same isotope
    for (IsotopicLabelVO label : isoLabels_) {
      isotopes = new ArrayList<String>();
      for (String element : label.getLabelElements().keySet()) {
        if (label.getLabelElements().get(element)>0)
          isotopes.add(element);
      }
      Collections.sort(isotopes);
      isotopeId = "";
      for (String iso : isotopes) 
        isotopeId += iso;
      if (!labelsOfSameIsotope.containsKey(isotopeId))
        labelsOfSameIsotope.put(isotopeId, new Vector<IsotopicLabelVO>());
      labelsOfSameIsotope.get(isotopeId).add(label);
    }
    //now calculate the retention time values weighted by the frequency of their occurrence
    for (Vector<IsotopicLabelVO> labels : labelsOfSameIsotope.values()) {
      isotopes = new ArrayList<String>();
      for (String element : labels.get(0).getLabelElements().keySet()) {
        if (labels.get(0).getLabelElements().get(element)>0)
          isotopes.add(element);
      }
      float totalShift = 0f;
      int totalDivisor = 0;
      int nrElements;
      for (IsotopicLabelVO label : labels) {
        if (label.getRtShift()==null)
          continue;
        nrElements = 0;
        for (String element : isotopes)
          nrElements += label.getLabelElements().get(element);
        totalShift += ((float)labelRtDiffs.get(label.getLabelId()).size())*(label.getRtShift()/((float)nrElements));
        totalDivisor += labelRtDiffs.get(label.getLabelId()).size();
      }
      float shiftEachLabelElement = totalShift/((float)totalDivisor);
//      System.out.println("shiftEachLabelElement: "+shiftEachLabelElement);
      for (IsotopicLabelVO label : labels) {
        nrElements = 0;
        for (String element : isotopes)
          nrElements += label.getLabelElements().get(element);
        label.setRtShift(shiftEachLabelElement*((float)nrElements));
      }
    }
    // the results
//    for (IsotopicLabelVO label : isoLabels_) {
//      System.out.println(label.getLabelId()+": "+label.getRtShift());
//    }
    
//    System.out.println("nrIsos: "+nrIsos);
  }
  
  
  /**
   * searches for elemental composition of the difference between labeled partners
   * @param lClass the lipid class
   * @param unlabeled the name of the unlabeled version
   * @param labeled the name of the labeled version
   * @return the difference in elements for the labeling
   */
  private Hashtable<String,Integer> checkForIsoPartnerMassDifference(String lClass, String unlabeled, String labeled){
    Hashtable<String,Integer> result = null;
    String analOnly;
    int diff;
    for (String anal : this.allMoleculeNames_.get(lClass)) {
      if (!anal.startsWith(labeled))
        continue;
      analOnly = new String(anal);
      if (expRtGroupingTime_>0)
        analOnly = analOnly.substring(0,analOnly.lastIndexOf("_"));
      if (!analOnly.equalsIgnoreCase(labeled))
        continue;
      Hashtable<String,Integer> formulaUnlabeled = getFormulaFromResultAreaVOHash(lClass, unlabeled);
      Hashtable<String,Integer> formulaLabeled = getFormulaFromResultAreaVOHash(lClass, labeled);
      if (formulaUnlabeled==null || formulaLabeled==null)
        continue;
      result = new Hashtable<String,Integer>();
      for (String element : formulaLabeled.keySet()) {
        diff = formulaLabeled.get(element)-(formulaUnlabeled.containsKey(element) ? formulaUnlabeled.get(element) : 0);
        //TODO: I am not sure whether I should display the subtracted chemical elements
////        if (diff<1)
////          continue;
        result.put(element, diff);
      }
    }
    return result;
  }
  
  
  /**
   * searches in the allResultsHash_ for a specific analyte and returns the chemical formula as hash table
   * @param lClass the class name of the analyte
   * @param anal the name of the analyte
   * @return the chemical formula as hash table
   */
  private Hashtable<String,Integer> getFormulaFromResultAreaVOHash(String lClass, String anal){
    Hashtable<String,Integer> formula = null;
    for (Hashtable<String,ResultAreaVO> hitsPerFile : allResultsHash_.get(lClass).values()) {
      for (String otherAnal : hitsPerFile.keySet()) {
        String other = (expRtGroupingTime_>0 ? otherAnal.substring(0,otherAnal.lastIndexOf("_")) : otherAnal);
        if (!anal.contentEquals(other))
          continue;
        formula = hitsPerFile.get(otherAnal).getChemicalFormulaElements();
        break;        
      }
      if (formula!=null)
        break;
    }
    return formula;
  }

  /**
   * 
   * @return list of potential isotopic labels
   */
  public List<IsotopicLabelVO> getIsoLabels()
  {
    return isoLabels_;
  }
  
  /**
   * adds the number of elements to a chemical formula in form of a hashtable
   * @param elements the chemical formula
   * @param labelElements the label elements that have to be added
   * @param nrOfLabels how often has this label been applied to this species
   */
  private void addLabelElements(Hashtable<String,Integer> elements, Hashtable<String,Integer> labelElements, int nrOfLabels) {
    for (String element : labelElements.keySet()) {
      int nrOfElements = 0;
      if (elements.containsKey(element))
        nrOfElements = elements.get(element);
      nrOfElements += nrOfLabels*labelElements.get(element);
      elements.put(element, nrOfElements);
    }
  }
  
  /**
   * searches for the closest unlabeled partners out of list if identifications - the partner eluting sooner and later are returned
   * @param isoRt the retention time of the isotopically labeled species
   * @param foundRTs the retention times of the detected unlabeled species
   * @return a vector containing two species: the first one is the unlabeled partner that elutes sooner, the second one elutes later; if no partner is detected, the DoubleStringVO.getKey() returns null 
   */
  private Vector<DoubleStringVO> findClosestPartners(double isoRt, List<DoubleStringVO> foundRTs) {
    Vector<DoubleStringVO> closestInformation = new Vector<DoubleStringVO>();
    String[] partners = new String[2];
    partners[0] = null;
    partners[1] = null;
    Double[] timeDiffs = new Double[2];
    timeDiffs[0] = Double.MAX_VALUE;
    timeDiffs[1] = Double.MAX_VALUE;
    for (DoubleStringVO hit : foundRTs) {
      if (hit.getValue()<=isoRt) {
        if ((isoRt-hit.getValue())<timeDiffs[0]) {
          timeDiffs[0] = isoRt-hit.getValue();
          partners[0] = hit.getKey();
        }
      }else{
        if ((hit.getValue()-isoRt)<timeDiffs[1]) {
          timeDiffs[1] = (hit.getValue()-isoRt);
          partners[1] = hit.getKey();
        }
      }
    }
    closestInformation.add(new DoubleStringVO(partners[0],timeDiffs[0]));
    closestInformation.add(new DoubleStringVO(partners[1],timeDiffs[1]));
    return closestInformation;
  }
  
  /**
   * detects the strongest chain identification of closely eluting unlabeled species
   * @param closestInfo a vector containing two species: the first one is the unlabeled partner that elutes sooner, the second one elutes later; if no partner is detected, the DoubleStringVO.getKey() returns null; the DoubleStringVO.getValue() contains the retention time difference
   * @param file the file where the detection had been found
   * @param mostOftenChainCombiInEachFile a hash table containing information about the strongest chain identification in each search; first key: the molecule name inclusive retention time; second key: the file name; value: the chain identification
   * @param mostOftenChainCombination a hash table containing the strongest, most often detected chain identification detected in all of the provided MS-runs; key: the molecule name inclusive retention time; value: the chain identification
   * @return the strongest chain identification of closely eluting unlabeled species; String[0] is the sooner eluting partner; String[1] is the later eluting partner; the corresponding values are null if there was no partner identified
   */
  private String[] getChainIdentifications(Vector<DoubleStringVO> closestInfo, String file, Hashtable<String,Hashtable<String,String>> mostOftenChainCombiInEachFile,
      Hashtable<String,String> mostOftenChainCombination) {
    String[] partnerChains = new String[2];
    partnerChains[0]=null;
    partnerChains[1]=null;
    
    if (closestInfo.get(0).getKey()!=null && mostOftenChainCombiInEachFile.get(closestInfo.get(0).getKey()).containsKey(file))
      partnerChains[0] = mostOftenChainCombiInEachFile.get(closestInfo.get(0).getKey()).get(file);
    else if (closestInfo.get(0).getKey()!=null && mostOftenChainCombination.containsKey(closestInfo.get(0).getKey()))
      partnerChains[0] = mostOftenChainCombination.get(closestInfo.get(0).getKey());
    
    if (closestInfo.get(1).getKey()!=null && mostOftenChainCombiInEachFile.get(closestInfo.get(1).getKey()).containsKey(file))
      partnerChains[1] = mostOftenChainCombiInEachFile.get(closestInfo.get(1).getKey()).get(file);
    else if (closestInfo.get(1).getKey()!=null && mostOftenChainCombination.containsKey(closestInfo.get(1).getKey()))
      partnerChains[1] = mostOftenChainCombination.get(closestInfo.get(1).getKey());
    
    return partnerChains;
  }
  
  
  @SuppressWarnings("unused")
  /**
   * detects the partner that has the same chain identification as the isotopically labeled species; or if there are no chain identifications, the closest partner
   * @param label value object containing information about the isotopic label
   * @param nrOfLabels the number of labels applied to to this isotopically labeled species
   * @param closestInfo a vector containing two species: the first one is the unlabeled partner that elutes sooner, the second one elutes later; if no partner is detected, the DoubleStringVO.getKey() returns null; the DoubleStringVO.getValue() contains the retention time difference
   * @param labeledChain the chain identification of the isotopically labeled species
   * @param chainIdentifications the strongest chain identification of closely eluting unlabeled species; String[0] is the sooner eluting partner; String[1] is the later eluting partner; the corresponding values are null if there was no partner identified
   * @return vector of two values; the first one is a DoubleStringVO with the closest partner; the second is a boolean that is true when the closest partner has the same chain identification
   * @throws LipidCombinameEncodingException thrown if the encoded lipid name does not follow the LDA syntax
   */
  private Vector<Object> checkChainInformationOrSelectClosestPartner(IsotopicLabelVO label, int nrOfLabels, Vector<DoubleStringVO> closestInfo, String labeledChain, String[] chainIdentifications) throws LipidCombinameEncodingException {
    Vector<Object> partnerValues = new Vector<Object>();
    DoubleStringVO partner = null;
    boolean isChainIdent = false;
    if (labeledChain!=null && (chainIdentifications[0]!=null || chainIdentifications[1]!=null)) {
//      System.out.println(labeledChain+";"+closestInfo.get(0).getKey()+" "+chainIdentifications[0]+";"+closestInfo.get(1).getKey()+" "+chainIdentifications[1]);
      isChainIdent = true;
      boolean containsFirstChain = containsLabeledChain(label, nrOfLabels, labeledChain, chainIdentifications[0]);
      boolean containsSecondChain = containsLabeledChain(label, nrOfLabels, labeledChain, chainIdentifications[1]);
      if (containsFirstChain && containsSecondChain) {
        if (closestInfo.get(1).getValue()<closestInfo.get(0).getValue())
          partner = closestInfo.get(1);
        else
          partner = closestInfo.get(0);
      }else if (containsFirstChain)
        partner = closestInfo.get(0);
      else if (containsSecondChain)
        partner = closestInfo.get(1);
    } else if (labeledChain==null && (chainIdentifications[0]==null || chainIdentifications[1]==null)) {
      if (closestInfo.get(1).getValue()<closestInfo.get(0).getValue())
        partner = closestInfo.get(1);
      else
        partner = closestInfo.get(0);
    }
    if (partner!=null) {
      partnerValues.add(partner);
      partnerValues.add(isChainIdent);
    }
    return partnerValues;
  }
  
  /**
   * checks whether the unlabeled chain identification is the same as the one carrying the isotopic label
   * @param label value object containing information about the isotopic label
   * @param nrOfLabels the number of labels applied to to this isotopically labeled species
   * @param labeledChain the identification of the labeled chain in its LDA encoded stale
   * @param chainIdentification the identification of the unlabeled chain in its LDA encoded style
   * @return true when the unlabeled chain and the labeled version are principally the same (except for the additional label)
   * @throws LipidCombinameEncodingException thrown if the encoded lipid name does not follow the LDA syntax
   */
  private boolean containsLabeledChain(IsotopicLabelVO label, int nrOfLabels, String labeledChain, String chainIdentification) throws LipidCombinameEncodingException {
    if (chainIdentification==null)
      return false;
    Vector<String> permutedLabels = createPermutedLabeledChains(label,nrOfLabels,chainIdentification);
    for (String permutedLabel : permutedLabels) {
      if (StaticUtils.isAPermutedVersion(permutedLabel, labeledChain, LipidomicsConstants.CHAIN_COMBI_SEPARATOR))
        return true;
    }
    return false;
  }
  
  /**
   * creates a list of possible chain combinations out of the unlabeled species carrying the label(s) at different chain positions 
   * @param label value object containing information about the isotopic label
   * @param nrOfLabels the number of labels applied to to this isotopically labeled species
   * @param chainIdentification the identification of the unlabeled chain in its LDA encoded style
   * @return a list of possible chain combinations species carrying the label(s) at different chain positions 
   * @throws LipidCombinameEncodingException thrown if the encoded lipid name does not follow the LDA syntax
   */
  private Vector<String> createPermutedLabeledChains(IsotopicLabelVO label, int nrOfLabels, String chainIdentification) throws LipidCombinameEncodingException {
    Vector<FattyAcidVO> fas = StaticUtils.decodeLipidNamesFromChainCombi(chainIdentification);
    Vector<Vector<Boolean>> labelPositions = getAllPossibleCombinations(fas.size());
    Vector<Vector<Boolean>> correctLabelNr = new Vector<Vector<Boolean>>();
    for (Vector<Boolean> combi : labelPositions) {
      int labels = 0;
      for (Boolean labelThere : combi) {
        if (labelThere) labels++;
      }
      if (labels==nrOfLabels)
        correctLabelNr.add(combi);
    }
    Vector<String> chainCombis = new Vector<String>();
    //this was added that only species with at least two chains are used for the retention time shift detection
    if (fas.size()<2)
      return chainCombis;
    for (Vector<Boolean> correct : correctLabelNr) {
      Vector<FattyAcidVO> fasNew = new Vector<FattyAcidVO>();
      for (int i=0;i!=correct.size();i++) {
        if (correct.get(i)) {
          FattyAcidVO old = fas.get(i);
          //it does not matter if I take the wrong old masses
          fasNew.add(new FattyAcidVO(old.getChainType(), label.getLabelId(), old.getcAtoms(), old.getDoubleBonds(), old.getOhNumber(), old.getMass(), old.getFormula()));
        }else {
          fasNew.add(fas.get(i));
        }
      }
      chainCombis.add(StaticUtils.encodeLipidCombi(fasNew));
    }
    return chainCombis;
  }
  
  Vector<Vector<Boolean>> getAllPossibleCombinations(int size){
    Vector<Vector<Boolean>> combis = new Vector<Vector<Boolean>>();
    if (size>1) {
      Vector<Vector<Boolean>> otherCombis = getAllPossibleCombinations(size-1);
      for (Vector<Boolean> combi : otherCombis) {
        Vector<Boolean> newCombi = new Vector<Boolean>(combi);
        newCombi.add(0,true);
        combis.add(newCombi);
      }
      for (Vector<Boolean> combi : otherCombis) {
        Vector<Boolean> newCombi = new Vector<Boolean>(combi);
        newCombi.add(0,false);
        combis.add(newCombi);
      }
    }else {
      Vector<Boolean> newCombi = new Vector<Boolean>();
      newCombi.add(true);
      combis.add(newCombi);
      newCombi = new Vector<Boolean>();
      newCombi.add(false);
      combis.add(newCombi);

    }
    return combis;
  }

  /**
   * 
   * @return the character encoding of the number of hydroxylation sites for the FA
   */
  public HydroxyEncoding getFaHydroxyEncoding()
  {
    return faHydroxyEncoding_;
  }

  /**
   * 
   * @return the character encoding of the number of hydroxylation sites for the LCB
   */
  public HydroxyEncoding getLcbHydroxyEncoding()
  {
    return lcbHydroxyEncoding_;
  }

  /**
   * 
   * @return a hash table containing the number of chains for each lipid class
   */
  public Hashtable<String,Integer> getNrOfChainsOfClass()
  {
    return chainsOfClass_;
  }
  
  /**
   * checks whether this species is an internal or external standard
   * @param className name of the analyte class
   * @param analyteName name of the analyte
   * @return true if the analyte is an internal or external standard
   */
  public boolean isISorES(String className, String analyteName) {
    if (this.allISNames_.containsKey(className) && this.allISNames_.get(className).containsKey(analyteName))
      return true;
    if (this.allESNames_.containsKey(className) && this.allESNames_.get(className).containsKey(analyteName))
      return true;
    return false;
  }
  
  
}
