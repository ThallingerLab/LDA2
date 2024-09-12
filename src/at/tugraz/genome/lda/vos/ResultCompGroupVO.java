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

package at.tugraz.genome.lda.vos;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.analysis.exception.CalculationNotPossibleException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ResultCompGroupVO extends ResultCompVO
{
  private Hashtable<String,ResultCompVO> oneGroup_;
  
//  private Vector<Double> highestGroupMeans_;
//  private Vector<Double> highestGroupSds_;
  private Vector<Double> sumGroupMeans_;
  @SuppressWarnings("unused")
  private Vector<Double> sumGroupSds_;
//  private Vector<Double> highestTotalMeans_;
//  private Vector<Double> highestTotalSds_;
  private Vector<Double> sumTotalMeans_;
  @SuppressWarnings("unused")
  private Vector<Double> sumTotalSds_;
  
  @SuppressWarnings("unused")
  private Vector<Double> sumPercentualMeans_;
  @SuppressWarnings("unused")
  private Vector<Double> sumPercentualSds_;
//  private Vector<Double> relativeMedianAreaSD_;
  
  public ResultCompGroupVO(Hashtable<String,ResultCompVO> oneGroup){
    oneGroup_ = oneGroup;

    super.usedIsotpes_ = 0;
    super.dilutionFactor_ = 1d;
    for (ResultCompVO compVO : oneGroup_.values()){
      type_ = compVO.type_;
      if (usedIsotpes_ == 0)
        usedIsotpes_ = compVO.getUsedIsotpes();
      if (compVO.getUsedIsotpes()<usedIsotpes_)
        usedIsotpes_ = compVO.getUsedIsotpes();
      if (compVO.dilutionFactor_!=null)
        dilutionFactor_ = compVO.dilutionFactor_;
    }
//    highestGroupMeans_ = null;
//    highestGroupSds_ = null;
    sumGroupMeans_ = null;
    sumGroupSds_ = null;
//    highestTotalMeans_ = null;
//    highestTotalSds_ = null;
    sumTotalMeans_ = null;
    sumTotalSds_ = null;
    sumPercentualMeans_ = null;
    sumPercentualSds_ = null;
  }
  
  public double getOriginalArea(int maxIsotope)
  {
    return getMeanValue(this.extractOriginalAreas(maxIsotope));
  }
  public double getOriginalAreaSD(int maxIsotope)
  {
    return Calculator.stddeviation(this.extractOriginalAreas(maxIsotope));
  }
  public double getOriginalAreaSE(int maxIsotope)
  {
    return this.getSE(this.extractOriginalAreas(maxIsotope));
  }
  private double getSE(Vector<Double> values){
    return Calculator.stddeviation(values)/Math.sqrt(values.size());
  } 
  private Vector<Double> extractOriginalAreas(int maxIsotope){
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getOriginalArea(maxIsotope);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getStandardizedArea(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor)
  {
    
    return getMeanValue(this.extractStandardizedAreas(maxIsotope, standMethod, esStandMethod, dilutionFactor));
  }
  public double getStandardizedAreaSD(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor)
  {
    return Calculator.stddeviation(this.extractStandardizedAreas(maxIsotope, standMethod, esStandMethod, dilutionFactor));
  }
  public double getStandardizedAreaSE(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor)
  {
    return this.getSE(this.extractStandardizedAreas(maxIsotope, standMethod, esStandMethod, dilutionFactor));
  }
  private Vector<Double> extractStandardizedAreas(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor){
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getStandardizedArea(maxIsotope, standMethod,esStandMethod, dilutionFactor);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getRelativeValue(int maxIsotope, ResultDisplaySettingsVO settingVO)
  {
    if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT)||settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT)){
      return super.getRelativeValue(maxIsotope, settingVO);
    }else{
      return getMeanValue(this.extractRelativeValues(maxIsotope, settingVO));
    }
  }
  public double getRelativeValueSD(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException
  {
      return getAreaSD(maxIsotope, settingVO)/relativeMedianAreas_.get(SUM_COMPOSITION).get(maxIsotope);
  }
  public double getRelativeValueSE(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException
  {
    Vector<Double> values = extractRelativeValues(maxIsotope, settingVO);
    double stdev = getRelativeValueSD(maxIsotope, settingVO);
    return stdev/Math.sqrt(values.size());
  }
  private Vector<Double> extractRelativeValues(int maxIsotope, ResultDisplaySettingsVO settingVO){
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getRelativeValue(maxIsotope,settingVO);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
//  public double getRelativeToMedianOfISRatios(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor)
//  {
//    return getMeanValue(this.extractRelativeToMedianOfISRatioss(maxIsotope, standMethod, esStandMethod, dilutionFactor));
//  } 
//  public double getRelativeToMedianOfISRatiosSD(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor)
//  {
//    return Calculator.stddeviation(this.extractRelativeToMedianOfISRatioss(maxIsotope, standMethod, esStandMethod, dilutionFactor));
//  } 
//  public double getRelativeToMedianOfISRatiosSE(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor)
//  {
//    return this.getSE(this.extractRelativeToMedianOfISRatioss(maxIsotope, standMethod, esStandMethod, dilutionFactor));
//  }
//  private Vector<Double> extractRelativeToMedianOfISRatioss(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor){
//    Vector<Double> values = new Vector<Double>();
//    ResultDisplaySettingsVO settingVO = new ResultDisplaySettingsVO(null,standMethod, esStandMethod, dilutionFactor,false);
//    for (ResultCompVO groupVO : oneGroup_.values()){
//      double value = groupVO.getRelativeToMedian(maxIsotope, settingVO);
//      if (value>0)
//        values.add(value);
//    }
//    return values;
//  }
  
//  public double getRelativeToISMedian(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor)
//  {
//    return getMeanValue(this.extractRelativeToISMedians(maxIsotope, standMethod, esStandMethod, dilutionFactor));
//  } 
//  public double getRelativeToISMedianSD(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor)
//  {
//    return Calculator.stddeviation(this.extractRelativeToISMedians(maxIsotope, standMethod, esStandMethod, dilutionFactor));
//  }
//  public double getRelativeToISMedianSE(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor)
//  {
//    return this.getSE(this.extractRelativeToISMedians(maxIsotope, standMethod, esStandMethod, dilutionFactor));
//  }
//  private Vector<Double> extractRelativeToISMedians(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor){
//    Vector<Double> values = new Vector<Double>();
//    for (ResultCompVO groupVO : oneGroup_.values()){
//      ResultDisplaySettingsVO settingVO = new ResultDisplaySettingsVO(null,standMethod, esStandMethod, dilutionFactor,false);
//      double value = groupVO.getRelativeToMedian(maxIsotope, settingVO);
//      if (value>0)
//        values.add(value);
//    }
//    return values;
//  }
  
  public double getAmountInEndVolume(int maxIsotope, int standMethod) throws CalculationNotPossibleException 
  {
    return getMeanValue(this.extractAmountInEndVolumes(maxIsotope, standMethod));
  } 
  public double getAmountInEndVolumeSD(int maxIsotope, int standMethod) throws CalculationNotPossibleException 
  {
    return Calculator.stddeviation(this.extractAmountInEndVolumes(maxIsotope, standMethod));
  }
  public double getAmountInEndVolumeSE(int maxIsotope, int standMethod) throws CalculationNotPossibleException 
  {
    return this.getSE(this.extractAmountInEndVolumes(maxIsotope, standMethod));
  }
  private Vector<Double> extractAmountInEndVolumes(int maxIsotope, int standMethod) throws CalculationNotPossibleException {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getAmountInEndVolume(maxIsotope, standMethod);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getConcentrationInEndVolume(int maxIsotope, int standMethod) throws CalculationNotPossibleException 
  {
    return getMeanValue(this.extractConcentrationInEndVolumes(maxIsotope, standMethod));
  } 
  public double getConcentrationInEndVolumeSD(int maxIsotope, int standMethod) throws CalculationNotPossibleException 
  {
    return Calculator.stddeviation(this.extractConcentrationInEndVolumes(maxIsotope, standMethod));
  }
  public double getConcentrationInEndVolumeSE(int maxIsotope, int standMethod) throws CalculationNotPossibleException 
  {
    return this.getSE(this.extractConcentrationInEndVolumes(maxIsotope, standMethod));
  }
  private Vector<Double> extractConcentrationInEndVolumes(int maxIsotope, int standMethod) throws CalculationNotPossibleException {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getConcentrationInEndVolume(maxIsotope, standMethod);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getWeightInEndVolume(int maxIsotope, int standMethod) throws CalculationNotPossibleException 
  {
    return getMeanValue(this.extractWeightInEndVolumes(maxIsotope, standMethod));
  } 
  public double getWeightInEndVolumeSD(int maxIsotope, int standMethod) throws CalculationNotPossibleException 
  {
    return Calculator.stddeviation(this.extractWeightInEndVolumes(maxIsotope, standMethod));
  }
  public double getWeightInEndVolumeSE(int maxIsotope, int standMethod) throws CalculationNotPossibleException 
  {
    return this.getSE(this.extractWeightInEndVolumes(maxIsotope, standMethod));
  }
  private Vector<Double> extractWeightInEndVolumes(int maxIsotope, int standMethod) throws CalculationNotPossibleException {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getWeightInEndVolume(maxIsotope, standMethod);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getAmountInProbeVolume(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return getMeanValue(this.extractAmountInProbeVolumes(maxIsotope, standMethod, esStandMethod));
  } 
  public double getAmountInProbeVolumeSD(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return Calculator.stddeviation(this.extractAmountInProbeVolumes(maxIsotope, standMethod, esStandMethod));
  }
  public double getAmountInProbeVolumeSE(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return this.getSE(this.extractAmountInProbeVolumes(maxIsotope, standMethod, esStandMethod));
  }
  private Vector<Double> extractAmountInProbeVolumes(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getAmountInProbeVolume(maxIsotope, standMethod, esStandMethod);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getConcentrationInProbeVolume(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return getMeanValue(this.extractConcentrationInProbeVolumes(maxIsotope, standMethod, esStandMethod));
  } 
  public double getConcentrationInProbeVolumeSD(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return Calculator.stddeviation(this.extractConcentrationInProbeVolumes(maxIsotope, standMethod, esStandMethod));
  }
  public double getConcentrationInProbeVolumeSE(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return this.getSE(this.extractConcentrationInProbeVolumes(maxIsotope, standMethod, esStandMethod));
  }
  private Vector<Double> extractConcentrationInProbeVolumes(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getConcentrationInProbeVolume(maxIsotope, standMethod, esStandMethod);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getWeightInProbeVolume(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return getMeanValue(this.extractWeightInProbeVolumes(maxIsotope, standMethod, esStandMethod));
  } 
  public double getWeightInProbeVolumeSD(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return Calculator.stddeviation(this.extractWeightInProbeVolumes(maxIsotope, standMethod, esStandMethod));
  }
  public double getWeightInProbeVolumeSE(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return this.getSE(this.extractWeightInProbeVolumes(maxIsotope, standMethod, esStandMethod));
  }
  private Vector<Double> extractWeightInProbeVolumes(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getWeightInProbeVolume(maxIsotope, standMethod, esStandMethod);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getAnalyteInRelationToProtein(int maxIsotope, int standMethod, int esStandMethod, boolean useAU) throws CalculationNotPossibleException 
  {
    return getMeanValue(this.extractAnalyteInRelationToProteins(maxIsotope, standMethod, esStandMethod, useAU));
  } 
  public double getAnalyteInRelationToProteinSD(int maxIsotope, int standMethod, int esStandMethod, boolean useAU) throws CalculationNotPossibleException 
  {
    return Calculator.stddeviation(this.extractAnalyteInRelationToProteins(maxIsotope, standMethod, esStandMethod, useAU));
  }
  public double getAnalyteInRelationToProteinSE(int maxIsotope, int standMethod, int esStandMethod, boolean useAU) throws CalculationNotPossibleException 
  {
    return this.getSE(this.extractAnalyteInRelationToProteins(maxIsotope, standMethod, esStandMethod, useAU));
  }
  private Vector<Double> extractAnalyteInRelationToProteins(int maxIsotope, int standMethod, int esStandMethod, boolean useAU) throws CalculationNotPossibleException {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getAnalyteInRelationToProtein(maxIsotope, standMethod, esStandMethod, useAU);
      if (value>0)
        values.add(value);
    }
    return values;
  }

  public double getAnalyteInRelationToSampleWeight(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return getMeanValue(this.extractAnalyteInRelationToSampleWeights(maxIsotope, standMethod, esStandMethod));
  } 
  public double getAnalyteInRelationToSampleWeightSD(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return Calculator.stddeviation(this.extractAnalyteInRelationToSampleWeights(maxIsotope, standMethod, esStandMethod));
  }
  public double getAnalyteInRelationToSampleWeightSE(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return this.getSE(this.extractAnalyteInRelationToSampleWeights(maxIsotope, standMethod, esStandMethod));
  }
  private Vector<Double> extractAnalyteInRelationToSampleWeights(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getAnalyteInRelationToSampleWeight(maxIsotope, standMethod, esStandMethod);
      if (value>0)
        values.add(value);
    }
    return values;
  }

  
  public double getAnalyteInRelationToNeutralLipidContent(int maxIsotope, int standMethod, int esStandMethod, boolean useAU) throws CalculationNotPossibleException 
  {
    return getMeanValue(this.extractAnalyteInRelationToNeutralLipidContents(maxIsotope, standMethod, esStandMethod,useAU));
  } 
  public double getAnalyteInRelationToNeutralLipidContentSD(int maxIsotope, int standMethod, int esStandMethod, boolean useAU) throws CalculationNotPossibleException 
  {
    return Calculator.stddeviation(this.extractAnalyteInRelationToNeutralLipidContents(maxIsotope, standMethod, esStandMethod,useAU));
  }
  public double getAnalyteInRelationToNeutralLipidContentSE(int maxIsotope, int standMethod, int esStandMethod, boolean useAU) throws CalculationNotPossibleException 
  {
    return this.getSE(this.extractAnalyteInRelationToNeutralLipidContents(maxIsotope, standMethod, esStandMethod,useAU));
  }
  private Vector<Double> extractAnalyteInRelationToNeutralLipidContents(int maxIsotope, int standMethod, int esStandMethod, boolean useAU) throws CalculationNotPossibleException {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getAnalyteInRelationToNeutralLipidContent(maxIsotope, standMethod, esStandMethod,useAU);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getAnalyteInRelationToMeasuredNeutralLipidContent(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return getMeanValue(this.extractAnalyteInRelationToMeasuredNeutralLipidContents(maxIsotope, standMethod, esStandMethod));
  } 
  public double getAnalyteInRelationToMeasuredNeutralLipidContentSD(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return Calculator.stddeviation(this.extractAnalyteInRelationToMeasuredNeutralLipidContents(maxIsotope, standMethod, esStandMethod));
  }
  public double getAnalyteInRelationToMeasuredNeutralLipidContentSE(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException 
  {
    return this.getSE(this.extractAnalyteInRelationToMeasuredNeutralLipidContents(maxIsotope, standMethod, esStandMethod));
  }
  private Vector<Double> extractAnalyteInRelationToMeasuredNeutralLipidContents(int maxIsotope, int standMethod, int esStandMethod) throws CalculationNotPossibleException {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      // I am leaving there the mass calculation as it is, because I do
      double value = groupVO.getAnalyteInRelationToMeasuredNeutralLipidContent(maxIsotope, standMethod, esStandMethod);
//      double value = -1;
//      if (totalGroupMass_!=null && totalGroupMass_.size()>0){
//        double amount = groupVO.getAmountInProbeVolume(maxIsotope, standMethod, esStandMethod);
//        double totalAmountGramm = groupVO.getAmountInProbeVolume(totalGroupMass_.get(maxIsotope), maxIsotope, standMethod, esStandMethod);
//        value = amount/totalAmountGramm;
//      }
      if (value>0)
        values.add(value);
    }
    return values;
  }

  public double getRatioToHighestPeak(int maxIsotope) 
  {
    return getMeanValue(this.extractRatioToHighestPeak(maxIsotope));
  } 
  public double getRatioToHighestPeakSD(int maxIsotope) 
  {
    return Calculator.stddeviation(this.extractRatioToHighestPeak(maxIsotope));
  }
  public double getRatioToHighestPeakSE(int maxIsotope) 
  {
    return this.getSE(this.extractRatioToHighestPeak(maxIsotope));
  }
  private Vector<Double> extractRatioToHighestPeak(int maxIsotope){
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getRatioToHighestPeak(maxIsotope);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  
  public double getRatioToTotalIntensity(int maxIsotope){
    double analyteInt = getMeanValue(extractRatioToTotalIntensity(maxIsotope));
    
    // this is for the normalisation on 100%
    return analyteInt/this.sumGroupMeans_.get(maxIsotope);
  }
  
  public double getRatioToTotalIntensitySD(int maxIsotope){
    Vector<Double> values = extractRatioToTotalIntensity(maxIsotope);
//    double analyteMean = getMeanValue(values);
    double analyteStdev = Calculator.stddeviation(values);
//    return Calculator.calculateRatioStdevErrorPropagated(analyteMean,Math.pow(analyteStdev,2), sumGroupMeans_.get(maxIsotope), Math.pow(sumGroupSds_.get(maxIsotope),2), 0);
    return analyteStdev;
  }
  
  public double getRatioToTotalIntensitySE(int maxIsotope) 
  {
    Vector<Double> values = extractRatioToTotalIntensity(maxIsotope);
    double stdev = getRatioToTotalIntensitySD(maxIsotope);
    return stdev/Math.sqrt(values.size());
  }
  
  public double getRatioToOverallGroupsIntensity(int maxIsotope){
    double analyteInt = getMeanValue(extractRatioToOverallGroupsIntensity(maxIsotope));
    // this is for the normalisation on 100%
    return analyteInt/sumTotalMeans_.get(maxIsotope);
  }
  
  public double getRatioToOverallGroupsIntensitySD(int maxIsotope){
    Vector<Double> values = extractRatioToOverallGroupsIntensity(maxIsotope);
//    double analyteMean = getMeanValue(values);
    double analyteStdev = Calculator.stddeviation(values);
    return analyteStdev;
//    return Calculator.calculateRatioStdevErrorPropagated(analyteMean, Math.pow(analyteStdev,2), sumTotalMeans_.get(maxIsotope), Math.pow(sumTotalSds_.get(maxIsotope),2), 0);
  }
  
  public double getRatioToOverallGroupsIntensitySE(int maxIsotope) 
  {
    Vector<Double> values = extractRatioToOverallGroupsIntensity(maxIsotope);
    double stdev = getRatioToOverallGroupsIntensitySD(maxIsotope);
    return stdev/Math.sqrt(values.size());
  }
  
  public double getMeanOfRatioToTotalIntensity(int maxIsotope) 
  {
    return getMeanValue(this.extractRatioToTotalIntensity(maxIsotope));
  } 
  public double getSDOfRatioToTotalIntensity(int maxIsotope) 
  {
    return Calculator.stddeviation(this.extractRatioToTotalIntensity(maxIsotope));
  }
  private Vector<Double> extractRatioToTotalIntensity(int maxIsotope){
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getRatioToTotalIntensity(maxIsotope);
      if (type_ == CLASS_TYPE)
        value = groupVO.getOriginalArea(maxIsotope);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getRelativeToMeasuredNeutralLipidContent(int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException
  {
    return getMeanValue(this.extractRelativeToMeasuredNeutralLipidContent(maxIsotope, standMethod, esMethod));
  } 
  public double getRelativeToMeasuredNeutralLipidContentSD(int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException
  {
    return Calculator.stddeviation(this.extractRelativeToMeasuredNeutralLipidContent(maxIsotope, standMethod, esMethod));
  }
  public double getRelativeToMeasuredNeutralLipidContentSE(int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException
  {
    return this.getSE(this.extractRelativeToMeasuredNeutralLipidContent(maxIsotope, standMethod, esMethod));
  }
  private Vector<Double> extractRelativeToMeasuredNeutralLipidContent(int maxIsotope, int standMethod, int esMethod)throws CalculationNotPossibleException{
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getRelativeToMeasuredNeutralLipidContent(maxIsotope, standMethod, esMethod);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
//  public double getArea(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException{
//    return getMeanValue(this.extractArea(maxIsotope, settingVO));
//  }
//  
//  public double getAreaSD(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException{
//    return Calculator.stddeviation(this.extractArea(maxIsotope, settingVO));
//  }
//  
//  public double getAreaSE(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException{
//    return this.getSE(this.extractArea(maxIsotope, settingVO));
//  }

//private Vector<Double> extractArea(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException {
//Vector<Double> values = new Vector<Double>();
//for (ResultCompVO groupVO : oneGroup_.values()){
//  double value = groupVO.getArea(maxIsotope, settingVO);
//  if (value>0)
//    values.add(value);
//}
//return values;
//}

  public double getAreaSD(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException{
    double area = 0;
    if (settingVO.isPercent()){
      area = getRatioToPercentualValueSD(maxIsotope,settingVO);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_VALUE)){
      area = getStandardizedAreaSD(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod(), settingVO.considerDilution());
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_BASE_PEAK)){
      area = getRatioToHighestPeakSD(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT)){
      area = getRatioToTotalIntensitySD(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_HIGHEST_TOTAL_PEAK)){
      area = getRatioToHighestFoundPeakSD(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT)){
      area = this.getRatioToOverallGroupsIntensitySD(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase("amount end-volume")){
      area = getAmountInEndVolumeSD(maxIsotope, settingVO.getISStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("conc. end-volume")){
      area = getConcentrationInEndVolumeSD(maxIsotope, settingVO.getISStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("weight end-volume")){
      area = getWeightInEndVolumeSD(maxIsotope, settingVO.getISStandMethod());      
    } else if (settingVO.getType().equalsIgnoreCase("amount sample-volume")){
      area = getAmountInProbeVolumeSD(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("conc. sample-volume")){
      area = getConcentrationInProbeVolumeSD(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("weight sample-volume")){
      area = getWeightInProbeVolumeSD(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod()); 
    } else if (settingVO.getType().equalsIgnoreCase("relative to sample weight")){
      area = getAnalyteInRelationToSampleWeightSD(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod()); 
    } else if (settingVO.getType().equalsIgnoreCase("relation to protein content")){
      area = getAnalyteInRelationToProteinSD(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod(),settingVO.isAu());
    } else if (settingVO.getType().equalsIgnoreCase("relation to neutral lipid content")){
      area = getAnalyteInRelationToNeutralLipidContentSD(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod(),settingVO.isAu());
    } else if (settingVO.getType().equalsIgnoreCase("relation to measured neutral lipid")){
      area = this.getAnalyteInRelationToMeasuredNeutralLipidContentSD(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    }else if (settingVO.getType().equalsIgnoreCase("percentual value")){
      area = getRatioToPercentualValueSD(maxIsotope,settingVO);
    }
    if (area>0 && !Double.isInfinite(area) && !Double.isNaN(area)){
      if (settingVO.getType().equalsIgnoreCase("conc. end-volume")||settingVO.getType().equalsIgnoreCase("conc. sample-volume")||
          settingVO.getType().equalsIgnoreCase("relative to sample weight")||settingVO.getType().equalsIgnoreCase("relation to protein content")||
          settingVO.getType().equalsIgnoreCase("relation to neutral lipid content")||settingVO.getType().equalsIgnoreCase("relation to measured neutral lipid"))
        area = StaticUtils.getValueDividedByUnit(area,settingVO.getDivisorMagnitude());
      return area;
    }
    return area;
  }
  
  public double getAreaSE(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException{
    double area = 0;
    if (settingVO.isPercent()){
      area = getRatioToPercentualValueSE(maxIsotope,settingVO);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_VALUE)){
      area = getStandardizedAreaSE(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod(), settingVO.considerDilution());
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_BASE_PEAK)){
      area = getRatioToHighestPeakSE(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT)){
      area = getRatioToTotalIntensitySE(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_HIGHEST_TOTAL_PEAK)){
      area = getRatioToHighestFoundPeakSE(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT)){
      area = this.getRatioToOverallGroupsIntensitySE(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase("amount end-volume")){
      area = getAmountInEndVolumeSE(maxIsotope, settingVO.getISStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("conc. end-volume")){
      area = getConcentrationInEndVolumeSE(maxIsotope, settingVO.getISStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("amount sample-volume")){
      area = getAmountInProbeVolumeSE(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("conc. sample-volume")){
      area = getConcentrationInProbeVolumeSE(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("relative to sample weight")){
      area = getAnalyteInRelationToSampleWeightSE(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());  
    } else if (settingVO.getType().equalsIgnoreCase("relation to protein content")){
      area = getAnalyteInRelationToProteinSE(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod(),settingVO.isAu());
    } else if (settingVO.getType().equalsIgnoreCase("relation to neutral lipid content")){
      area = getAnalyteInRelationToNeutralLipidContentSE(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod(),settingVO.isAu());
    } else if (settingVO.getType().equalsIgnoreCase("relation to measured neutral lipid")){
      area = this.getAnalyteInRelationToMeasuredNeutralLipidContentSE(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    }else if (settingVO.getType().equalsIgnoreCase("percentual value")){
      area = getRatioToPercentualValueSE(maxIsotope,settingVO);
    }
    if (area>0 && !Double.isInfinite(area) && !Double.isNaN(area)){
      if (settingVO.getType().equalsIgnoreCase("conc. end-volume")||settingVO.getType().equalsIgnoreCase("conc. sample-volume")||
          settingVO.getType().equalsIgnoreCase("relative to sample weight")||settingVO.getType().equalsIgnoreCase("relation to protein content")||
          settingVO.getType().equalsIgnoreCase("relation to neutral lipid content")||settingVO.getType().equalsIgnoreCase("relation to measured neutral lipid"))
        area = StaticUtils.getValueDividedByUnit(area,settingVO.getDivisorMagnitude());
      return area;
    }  
    return area;
  }
  
  
  public double getRelativeToMedian(int maxIsotope, ResultDisplaySettingsVO settingVO) {
    return getMeanValue(extractRelativeToMedian(maxIsotope, settingVO));
  }
  
  public double getRelativeToMedianSD(int maxIsotope, ResultDisplaySettingsVO settingVO) {
    return Calculator.stddeviation(extractRelativeToMedian(maxIsotope, settingVO));
  }
  
  public double getRelativeToMedianSE(int maxIsotope, ResultDisplaySettingsVO settingVO) {
    return getSE(extractRelativeToMedian(maxIsotope, settingVO));
  }
  
  private Vector<Double> extractRelativeToMedian(int maxIsotope, ResultDisplaySettingsVO settingVO) {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getRelativeToMedian(maxIsotope, settingVO);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public void addRelativeMedianArea(ArrayList<Double> relativeMedianArea)
  {
  	super.addRelativeMedianArea(relativeMedianArea);
  	for (ResultCompVO groupVO : oneGroup_.values()){
      groupVO.addRelativeMedianArea(relativeMedianArea);
    }
  }
  
//  public void setRelativeMedianAreaSD(Vector<Double> relativeMedianAreaSD)
//  {
//    this.relativeMedianAreaSD_ = relativeMedianAreaSD;    
//  }
  
  public double getRetentionTime(String modName){
    Vector<Double> values = extractRetentionTimes(modName);
    if (values.size()>0)
      return getMeanValue(extractRetentionTimes(modName));
    else
      return -1;
  }
  
  public double getRetentionTimeSD(String modName){
    Vector<Double> values = extractRetentionTimes(modName);
    if (values.size()>0)
      return Calculator.stddeviation(extractRetentionTimes(modName));
    else
      return -1;
  }
  
  private Vector<Double> extractRetentionTimes(String modName) {
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getRetentionTime(modName);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  private double getMeanValue(Vector<Double> areas){
    double mean = Calculator.mean(areas);
    if (!Double.isNaN(mean) && !Double.isInfinite(mean)){
      return mean;
    }else
      return 0d;
  }
  
  public boolean getMoreThanOnePeak(int maxIsotope)
  {
    return false;
  }
  
  public double getRatioToHighestFoundPeak(int maxIsotope) 
  {
    return getMeanValue(this.extractRatioToHighestFoundPeak(maxIsotope));
  } 
  public double getRatioToHighestFoundPeakSD(int maxIsotope) 
  {
    return Calculator.stddeviation(this.extractRatioToHighestFoundPeak(maxIsotope));
  }
  public double getRatioToHighestFoundPeakSE(int maxIsotope) 
  {
    return this.getSE(this.extractRatioToHighestFoundPeak(maxIsotope));
  }
  private Vector<Double> extractRatioToHighestFoundPeak(int maxIsotope){
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getRatioToHighestFoundPeak(maxIsotope);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public double getMeanRatioToOverallGroupsIntensity(int maxIsotope) 
  {
    return getMeanValue(this.extractRatioToOverallGroupsIntensity(maxIsotope));
  } 
  public double getSDRatioToOverallGroupsIntensity(int maxIsotope) 
  {
    return Calculator.stddeviation(this.extractRatioToOverallGroupsIntensity(maxIsotope));
  }
//  public double getRatioToOverallGroupsIntensitySE(int maxIsotope) 
//  {
//    return this.getSE(this.extractRatioToOverallGroupsIntensity(maxIsotope));
//  }
  private Vector<Double> extractRatioToOverallGroupsIntensity(int maxIsotope){
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getRatioToOverallGroupsIntensity(maxIsotope);
//      double value = -1;
//      if (totalFoundIntensity_!=null && totalFoundIntensity_.get(maxIsotope)>0)
//        value = groupVO.getOriginalArea(getAvailableIsotopeNr(maxIsotope))/totalFoundIntensity_.get(maxIsotope);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  @Override
  public double getMass(int maxIsotope)
  {
  	Collection<ResultCompVO> vos = oneGroup_.values();
  	int count = 0;
  	Double sum = 0d;
  	for (ResultCompVO vo : vos)
  	{
  		Double value = vo.getMass(maxIsotope);
  		if (value > 0)
  		{
  			sum += value;
  			count++;
  		}
  	}
  	return sum/count;
  }
  
  @Override
  public boolean isMSnVerifiedOrStandard()
  {
  	for (ResultCompVO compVO : oneGroup_.values())
  	{
  		if (compVO.isMSnVerifiedOrStandard())
  			return true;
  	}
  	return false;
  }
  
  public Hashtable<String,ResultCompVO> getGroupingPartners(){
    return oneGroup_;
  }

//  public void setHighestGroupMeans(Vector<Double> highestGroupMeans)
//  {
//    this.highestGroupMeans_ = highestGroupMeans;
//  }
//
//  public void setHighestGroupSds(Vector<Double> highestGroupSds)
//  {
//    this.highestGroupSds_ = highestGroupSds;
//  }

  public void setSumGroupMeans(Vector<Double> sumGroupMeans)
  {
    this.sumGroupMeans_ = sumGroupMeans;
  }

  public void setSumGroupSds(Vector<Double> sumGroupSds)
  {
    this.sumGroupSds_ = sumGroupSds;
  }

//  public void setHighestTotalMeans(Vector<Double> highestTotalMeans)
//  {
//    this.highestTotalMeans_ = highestTotalMeans;
//  }
//
//  public void setHighestTotalSds(Vector<Double> highestTotalSds)
//  {
//    this.highestTotalSds_ = highestTotalSds;
//  }

  public void setSumTotalMeans(Vector<Double> sumTotalMeans)
  {
    this.sumTotalMeans_ = sumTotalMeans;
  }

  public void setSumTotalSds(Vector<Double> sumTotalSds)
  {
    this.sumTotalSds_ = sumTotalSds;
  }
  
  public double getRatioToPercentualValue(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException{
    double analyteInt = getMeanValue(extractRatioToPercentualValue(maxIsotope,settingVO));
    
    // this is for the normalisation on 100%
//    return analyteInt/this.sumPercentualMeans_.get(maxIsotope);
    return analyteInt;
  }
  
  public double getRatioToPercentualValueSD(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException{
    Vector<Double> values = extractRatioToPercentualValue(maxIsotope,settingVO);
//    double analyteMean = getMeanValue(values);
    double analyteStdev = Calculator.stddeviation(values);
//    return analyteStdev/this.sumPercentualMeans_.get(maxIsotope);
//    return Calculator.calculateRatioStdevErrorPropagated(analyteMean,Math.pow(analyteStdev,2), sumPercentualMeans_.get(maxIsotope), Math.pow(sumPercentualSds_.get(maxIsotope),2), 0);
    return analyteStdev;
  }
  
  public double getRatioToPercentualValueSE(int maxIsotope, ResultDisplaySettingsVO settingVO)  throws CalculationNotPossibleException
  {
    Vector<Double> values = extractRatioToPercentualValue(maxIsotope,settingVO);
    double stdev = getRatioToPercentualValueSD(maxIsotope,settingVO);
    return stdev/Math.sqrt(values.size());
  }
  
  private Vector<Double> extractRatioToPercentualValue(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException{
    Vector<Double> values = new Vector<Double>();
    for (ResultCompVO groupVO : oneGroup_.values()){
      double value = groupVO.getRatioToPercentualValue(maxIsotope,settingVO);
//      if (type_ == CLASS_TYPE)
//        value = groupVO.getOriginalArea(maxIsotope);
      if (value>0)
        values.add(value);
    }
    return values;
  }
  
  public void setSumPercentualMeans(Vector<Double> sumGroupMeans)
  {
    this.sumPercentualMeans_ = sumGroupMeans;
  }

  public void setSumPercentualSds(Vector<Double> sumPercentualSds)
  {
    this.sumPercentualSds_ = sumPercentualSds;
  }
  
  public boolean hasAllModsFound(){
    return true;
  }
  
}
