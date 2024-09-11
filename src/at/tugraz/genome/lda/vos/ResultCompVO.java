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
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.analysis.exception.CalculationNotPossibleException;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ResultCompVO
{
  public final static int NO_STANDARD_CORRECTION = 0;
  public final static int STANDARD_CORRECTION_INTERNAL = 1;
  public final static int STANDARD_CORRECTION_MEDIAN = 2;
  
  public final static int ANALYTE_TYPE = 0;
  public final static int INTERNAL_STANDARD_TYPE = 1;
  public final static int EXTERNAL_STANDARD_TYPE = 2;
  public final static int CLASS_TYPE = 3;
  
  private ResultAreaVO resultMolecule_;
  
  protected int type_;
  protected Vector<Double> mass_;
  protected Hashtable<String,Double> retentionTime_;
  /** true if all this are is constructed of all found modifications*/
  protected boolean allModsFound_;
  
  private Vector<Double> originalArea_;
  private Vector<Double> correctionFactorInternalIS_;
  private Vector<Double> correctionFactorMedianIS_;
  private Hashtable<Integer,Vector<Double>> correctionFactorSingleIS_;
  protected Hashtable<String,ArrayList<Double>> relativeMedianAreas_ = new Hashtable<String,ArrayList<Double>>();
  public final static String SUM_COMPOSITION = "Sum composition";
  private Vector<Double> areaISInternalComparison_;
  private Vector<Double> areaISMedianComparison_;
  private Hashtable<Integer,Vector<Double>> areaISSingleComparison_;
  
  protected Double dilutionFactor_;
  
  protected Vector<Double> correctionFactorESNoISCorrInternal_;
  protected Vector<Double> correctionFactorESNoISCorrMedian_;
  private Hashtable<Integer,Vector<Double>> correctionFactorESNoISCorrSingle_;
  protected Vector<Double> correctionFactorESISInternalCorrInternal_;
  protected Vector<Double> correctionFactorESISInternalCorrMedian_;
  private Hashtable<Integer,Vector<Double>> correctionFactorESISInternalCorrSingle_;
  protected Vector<Double> correctionFactorESISMedianCorrInternal_;
  protected Vector<Double> correctionFactorESISMedianCorrMedian_;
  private Hashtable<Integer,Vector<Double>> correctionFactorESISMedianCorrSingle_ ;
  private Hashtable<Integer,Vector<Double>> correctionFactorESISSingleCorrInternal_;
  private Hashtable<Integer,Vector<Double>> correctionFactorInternalESISSingleCorrMedian_;
  private Hashtable<Integer,Hashtable<Integer,Vector<Double>>> correctionFactorESISSingleCorrSingle_;

  // the area of the external standards are already multiplied with the dilution factor
  protected Vector<Double> areaESNoISCorrInternalComparison_;
  protected Vector<Double> areaESNoISCorrMedianComparison_;  
  private Hashtable<Integer,Vector<Double>> areaESNoISCorrSingleComparison_;
  protected Vector<Double> areaESISInternalCorrInternalComparison_;
  protected Vector<Double> areaESISInternalCorrMedianComparison_;
  private Hashtable<Integer,Vector<Double>> areaESISInternalCorrSingleComparison_;
  protected Vector<Double> areaESISMedianCorrInternalComparison_;
  protected Vector<Double> areaESISMedianCorrMedianComparison_;
  private Hashtable<Integer,Vector<Double>> areaESISMedianCorrSingleComparison_;
  private Hashtable<Integer,Vector<Double>> areaESISSingleCorrInternalComparison_;
  private Hashtable<Integer,Vector<Double>> areaESISSingleCorrMedianComparison_;
  private Hashtable<Integer,Hashtable<Integer,Vector<Double>>> areaESISSingleCorrSingleComparison_;
  
  protected Double relativeValue_;
  
  protected int usedIsotpes_;
  protected String absoluteFilePath_;
  
  protected Hashtable<Integer,VolumeConcVO> isAmounts_;
  protected Hashtable<Integer,VolumeConcVO> esAmounts_;
  protected Double esStandardVolumeInternalCorr_;
  protected Double esStandardConcentrationInternalCorr_;
  protected Double esStandardVolumeMedianCorr_;
  protected Double esStandardConcentrationMedianCorr_;
  protected Hashtable<Integer,Double> esStandardVolumeSingleCorr_;
  protected Hashtable<Integer,Double> esStandardConcentrationSingleCorr_;
  protected Double endVolume_;
  protected Double probeVolume_;
  protected Double sampleWeight_;
  protected Double proteinConcentration_;
  protected Double neutralLipidConcentration_;
  protected Vector<Double> measuredNeutralLipidConcentration_;
  
  protected Vector<Double> highestGroupIntensity_;
  protected Vector<Double> totalGroupIntensity_;  
  protected Vector<Double> totalGroupMass_;
  protected Vector<Double> highestFoundIntensity_;
  protected Vector<Double> totalFoundIntensity_;
  protected Vector<Double> sumValueForPercentage_;  
  
  protected Vector<Hashtable<String,Boolean>> moreThanOnePeak_ = new Vector<Hashtable<String,Boolean>>();
  
  protected boolean hasAbs_;
  protected boolean existsInFile_;
  protected boolean isNullInFile_;
  
  
  // this should be just for CLASS_TYPE only
  public ResultCompVO (int type, int usedIsotopes, Vector<Double> originalArea){
    this.setConstructorValues(true,false,type, usedIsotopes, originalArea);
  }
  
  public ResultCompVO(ResultAreaVO resultMolecule, boolean existsInFile, boolean isNullInFile, int type,Vector<Double> mass,Hashtable<String,Double> retentionTime, String absoluteFilePath, int usedIsotopes, boolean allModsFound, Vector<Double> originalArea, Vector<Double> correctionFactorInternalIS, Vector<Double> correctionFactorMedianIS,
      Hashtable<Integer,Vector<Double>> correctionFactorSingleIS, Vector<Double> areaISInternalComparision,
      Vector<Double> areaISMedianComparison, Hashtable<Integer,Vector<Double>> areaISSingleComparison,
      Vector<Double> correctionFactorESNoISCorrInternal,Vector<Double> correctionFactorESNoISCorrMedian, Hashtable<Integer,Vector<Double>> correctionFactorESNoISCorrSingle,
      Vector<Double> correctionFactorESISInternalCorrInternal,Vector<Double> correctionFactorESISInternalCorrMedian, Hashtable<Integer,Vector<Double>> correctionFactorESISInternalCorrSingle,
      Vector<Double> correctionFactorESISMedianCorrInternal,Vector<Double> correctionFactorESISMedianCorrMedian, Hashtable<Integer,Vector<Double>>  correctionFactorESISMedianCorrSingle,
      Hashtable<Integer,Vector<Double>> correctionFactorESISSingleCorrInternal, Hashtable<Integer,Vector<Double>> correctionFactorInternalESISSingleCorrMedian,
      Hashtable<Integer,Hashtable<Integer,Vector<Double>>> correctionFactorESISSingleCorrSingle,
      Vector<Double> areaESNoISCorrInternalComparison, Vector<Double> areaESNoISCorrMedianComparison,  Hashtable<Integer,Vector<Double>> areaESNoISCorrSingleComparison,
      Vector<Double> areaESISInternalCorrInternalComparison, Vector<Double> areaESISInternalCorrMedianComparison,  Hashtable<Integer,Vector<Double>> areaESISInternalCorrSingleComparison,
      Vector<Double> areaESISMedianCorrInternalComparison, Vector<Double> areaESISMedianCorrMedianComparison, Hashtable<Integer,Vector<Double>> areaESISMedianCorrSingleComparison,
      Hashtable<Integer,Vector<Double>> areaESISSingleCorrInternalComparison, Hashtable<Integer,Vector<Double>> areaESISSingleCorrMedianComparison,
      Hashtable<Integer,Hashtable<Integer,Vector<Double>>> areaESISSingleCorrSingleComparison,
      Hashtable<Integer,VolumeConcVO> isAmount, Double dilutionFactor, 
      Hashtable<Integer,VolumeConcVO> esAmount, Double esStandardVolumeInternalCorr,
      Double esStandardConcentrationInternalCorr, Double esStandardVolumeMedianCorr, Double esStandardConcentrationMedianCorr,
      Hashtable<Integer,Double> esStandardVolumeSingleCorr, Hashtable<Integer,Double> esStandardConcentrationSingleCorr,
      Double endVolume, Double probeVolume, Double sampleWeight, Double proteinConcentration, Double neutralLipidConcentration, 
      Vector<Hashtable<String,Boolean>> moreThanOnePeak,
      boolean hasAbs){    
    super();
    this.setConstructorValues(resultMolecule,existsInFile,isNullInFile,type,mass,retentionTime,absoluteFilePath,usedIsotopes, allModsFound, originalArea, correctionFactorInternalIS, correctionFactorMedianIS, correctionFactorSingleIS,
        areaISInternalComparision, areaISMedianComparison, areaISSingleComparison,
        correctionFactorESNoISCorrInternal, correctionFactorESNoISCorrMedian, correctionFactorESNoISCorrSingle, correctionFactorESISInternalCorrInternal,
        correctionFactorESISInternalCorrMedian, correctionFactorESISInternalCorrSingle, correctionFactorESISMedianCorrInternal, correctionFactorESISMedianCorrMedian,
        correctionFactorESISMedianCorrSingle, correctionFactorESISSingleCorrInternal, correctionFactorInternalESISSingleCorrMedian,
        correctionFactorESISSingleCorrSingle,
        areaESNoISCorrInternalComparison, areaESNoISCorrMedianComparison, areaESNoISCorrSingleComparison, areaESISInternalCorrInternalComparison,
        areaESISInternalCorrMedianComparison, areaESISInternalCorrSingleComparison, areaESISMedianCorrInternalComparison, areaESISMedianCorrMedianComparison, areaESISMedianCorrSingleComparison,
        areaESISSingleCorrInternalComparison, areaESISSingleCorrMedianComparison, areaESISSingleCorrSingleComparison,
        isAmount, dilutionFactor, esAmount,esStandardVolumeInternalCorr, esStandardConcentrationInternalCorr,
        esStandardVolumeMedianCorr, esStandardConcentrationMedianCorr, esStandardVolumeSingleCorr, esStandardConcentrationSingleCorr,
        endVolume, probeVolume, sampleWeight, proteinConcentration, neutralLipidConcentration, moreThanOnePeak,hasAbs);    
  }
  
  protected void setConstructorValues(boolean existsInFile, boolean isNullInFile, int type, int usedIsotopes, Vector<Double> originalArea){
    existsInFile_ = existsInFile;
    this.isNullInFile_ = isNullInFile;
    this.type_ = type;
    this.usedIsotpes_ = usedIsotopes;
    this.originalArea_ = originalArea;
  }
  
  protected void setConstructorValues(ResultAreaVO resultMolecule, boolean existsInFile, boolean isNullInFile, int type, Vector<Double> mass, Hashtable<String,Double> retentionTime,String absoluteFilePath, int usedIsotopes, boolean allModsFound, Vector<Double> originalArea, Vector<Double> correctionFactorInternalIS, Vector<Double> correctionFactorMedianIS,
      Hashtable<Integer,Vector<Double>> correctionFactorSingleIS,
      Vector<Double> areaISInternalComparision, Vector<Double> areaISMedianComparison, Hashtable<Integer,Vector<Double>> areaISSingleComparison,
      Vector<Double> correctionFactorESNoISCorrInternal,Vector<Double> correctionFactorESNoISCorrMedian, Hashtable<Integer,Vector<Double>> correctionFactorESNoISCorrSingle,
      Vector<Double> correctionFactorESISInternalCorrInternal,Vector<Double> correctionFactorESISInternalCorrMedian, Hashtable<Integer,Vector<Double>> correctionFactorESISInternalCorrSingle,
      Vector<Double> correctionFactorESISMedianCorrInternal,Vector<Double> correctionFactorESISMedianCorrMedian, Hashtable<Integer,Vector<Double>>  correctionFactorESISMedianCorrSingle,
      Hashtable<Integer,Vector<Double>> correctionFactorESISSingleCorrInternal, Hashtable<Integer,Vector<Double>> correctionFactorInternalESISSingleCorrMedian,
      Hashtable<Integer,Hashtable<Integer,Vector<Double>>> correctionFactorESISSingleCorrSingle,
      Vector<Double> areaESNoISCorrInternalComparison, Vector<Double> areaESNoISCorrMedianComparison, Hashtable<Integer,Vector<Double>> areaESNoISCorrSingleComparison,
      Vector<Double> areaESISInternalCorrInternalComparison, Vector<Double> areaESISInternalCorrMedianComparison, Hashtable<Integer,Vector<Double>> areaESISInternalCorrSingleComparison,
      Vector<Double> areaESISMedianCorrInternalComparison, Vector<Double> areaESISMedianCorrMedianComparison, Hashtable<Integer,Vector<Double>> areaESISMedianCorrSingleComparison,
      Hashtable<Integer,Vector<Double>> areaESISSingleCorrInternalComparison, Hashtable<Integer,Vector<Double>> areaESISSingleCorrMedianComparison,
      Hashtable<Integer,Hashtable<Integer,Vector<Double>>> areaESISSingleCorrSingleComparison,
      Hashtable<Integer,VolumeConcVO> isAmount, Double dilutionFactor, 
      Hashtable<Integer,VolumeConcVO> esAmount, Double esStandardVolumeInternalCorr, 
      Double esStandardConcentrationInternalCorr, Double esStandardVolumeMedianCorr, Double esStandardConcentrationMedianCorr,
      Hashtable<Integer,Double> esStandardVolumeSingleCorr, Hashtable<Integer,Double> esStandardConcentrationSingleCorr,
      Double endVolume, Double probeVolume, Double sampleWeight, Double proteinConcentration, Double neutralLipidConcentration, Vector<Hashtable<String,Boolean>> moreThanOnePeak, boolean hasAbs){
    this.setConstructorValues(existsInFile, isNullInFile,type, usedIsotopes, originalArea);
    this.resultMolecule_ = resultMolecule;
    this.mass_ = mass;
    if (mass_==null){
      mass_ = new Vector<Double>();
      for (int i=0;i!=usedIsotpes_;i++)mass_.add(0d);
    }
    this.retentionTime_ = retentionTime;
    this.allModsFound_ = allModsFound;
    this.absoluteFilePath_ = absoluteFilePath;
    this.correctionFactorInternalIS_ = correctionFactorInternalIS;
    this.correctionFactorMedianIS_ = correctionFactorMedianIS;
    this.correctionFactorSingleIS_ = correctionFactorSingleIS;
    this.areaISInternalComparison_ = areaISInternalComparision;
    this.areaISMedianComparison_ = areaISMedianComparison;
    this.areaISSingleComparison_ = areaISSingleComparison;
    
    highestGroupIntensity_ = null;
    totalGroupIntensity_ = null;
    totalGroupMass_ = null;
    correctionFactorESNoISCorrInternal_ = correctionFactorESNoISCorrInternal;
    correctionFactorESNoISCorrMedian_ = correctionFactorESNoISCorrMedian;
    correctionFactorESNoISCorrSingle_ = correctionFactorESNoISCorrSingle;
    correctionFactorESISInternalCorrInternal_ = correctionFactorESISInternalCorrInternal;
    correctionFactorESISInternalCorrMedian_ = correctionFactorESISInternalCorrMedian;
    correctionFactorESISInternalCorrSingle_ = correctionFactorESISInternalCorrSingle;
    correctionFactorESISMedianCorrInternal_ = correctionFactorESISMedianCorrInternal;
    correctionFactorESISMedianCorrMedian_ = correctionFactorESISMedianCorrMedian;
    correctionFactorESISMedianCorrSingle_ = correctionFactorESISMedianCorrSingle;
    correctionFactorESISSingleCorrInternal_ = correctionFactorESISSingleCorrInternal;
    correctionFactorInternalESISSingleCorrMedian_ = correctionFactorInternalESISSingleCorrMedian;
    correctionFactorESISSingleCorrSingle_ = correctionFactorESISSingleCorrSingle;
    
    areaESNoISCorrInternalComparison_ = areaESNoISCorrInternalComparison;
    areaESNoISCorrMedianComparison_ = areaESNoISCorrMedianComparison;
    areaESNoISCorrSingleComparison_ = areaESNoISCorrSingleComparison;
    areaESISInternalCorrInternalComparison_ = areaESISInternalCorrInternalComparison;
    areaESISInternalCorrMedianComparison_ = areaESISInternalCorrMedianComparison;
    areaESISInternalCorrSingleComparison_ = areaESISInternalCorrSingleComparison;
    areaESISMedianCorrInternalComparison_ = areaESISMedianCorrInternalComparison;
    areaESISMedianCorrMedianComparison_ = areaESISMedianCorrMedianComparison;
    areaESISMedianCorrSingleComparison_ = areaESISMedianCorrSingleComparison;
    areaESISSingleCorrInternalComparison_ = areaESISSingleCorrInternalComparison;
    areaESISSingleCorrMedianComparison_ = areaESISSingleCorrMedianComparison;
    areaESISSingleCorrSingleComparison_ = areaESISSingleCorrSingleComparison;    
    
    
    this.dilutionFactor_ = dilutionFactor;
//    this.correctionFactorInternalES_ = correctionFactorInternalES;
//    this.correctionFactorMedianES_ = correctionFactorMedianES;
//    this.areaESInternalComparison_ = areaESInternalComparison;
//    this.areaESMedianComparison_ = areaESMedianComparison;
    
//    this.isStandardVolume_ = isStandardVolume;
//    this.isStandardConcentration_ = isStandardConcentration;
//    this.esStandardVolumeNoCorr_ = esStandardVolumeNoCorr;
//    this.esStandardConcentrationNoCorr_ = esStandardConcentrationNoCorr;
    this.isAmounts_ = isAmount;
    this.esAmounts_ = esAmount;
    this.esStandardVolumeInternalCorr_ = esStandardVolumeInternalCorr;
    this.esStandardConcentrationInternalCorr_ = esStandardConcentrationInternalCorr;
    this.esStandardVolumeMedianCorr_ = esStandardVolumeMedianCorr;
    this.esStandardConcentrationMedianCorr_ = esStandardConcentrationMedianCorr;
    this.esStandardVolumeSingleCorr_ = esStandardVolumeSingleCorr;
    this.esStandardConcentrationSingleCorr_ = esStandardConcentrationSingleCorr;
    this.endVolume_ = endVolume;
    this.probeVolume_ = probeVolume;
    this.sampleWeight_ = sampleWeight;
    this.proteinConcentration_ = proteinConcentration;
    this.neutralLipidConcentration_ = neutralLipidConcentration;
    this.moreThanOnePeak_ = moreThanOnePeak;
    this.hasAbs_ = hasAbs;
  }
  
  public ResultCompVO()
  {
    this(null, false,false,ResultCompVO.ANALYTE_TYPE, null, new Hashtable<String,Double>(),"", 0, false, new Vector<Double>(), 
        new Vector<Double>(), new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(), new Vector<Double>(), 
        new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(), new Vector<Double>(),
        new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(), new Vector<Double>(), new Vector<Double>(), 
        new Hashtable<Integer,Vector<Double>>(),new Vector<Double>(), new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(),
        new Hashtable<Integer,Vector<Double>>(), new Hashtable<Integer,Vector<Double>>(), new Hashtable<Integer,Hashtable<Integer,Vector<Double>>>(),
        new Vector<Double>(), new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(), new Vector<Double>(), new Vector<Double>(), 
        new Hashtable<Integer,Vector<Double>>(), new Vector<Double>(), new Vector<Double>(), new Hashtable<Integer,Vector<Double>>(),
        new Hashtable<Integer,Vector<Double>>(), new Hashtable<Integer,Vector<Double>>(),new Hashtable<Integer,Hashtable<Integer,Vector<Double>>>(),
        new Hashtable<Integer,VolumeConcVO>(), null, new Hashtable<Integer,VolumeConcVO>(), null, null, null, null, null, null, null, null, null, null,null,
        new Vector<Hashtable<String,Boolean>>(),false);
  }
  
  public ResultAreaVO getResultMolecule()
  {
  	return this.resultMolecule_;
  }
  
  public boolean isMSnVerifiedOrStandard()
  {
  	return this.resultMolecule_ != null && (this.resultMolecule_.isMSnVerified() || this.resultMolecule_.isAStandard());
  }
  
  public int getType()
  {
    return type_;
  }

  public double getOriginalArea(int maxIsotope)
  {
    if (maxIsotope>=usedIsotpes_ || maxIsotope<0)
      return 0;
    else
      return originalArea_.get(maxIsotope);
  }

  public double getStandardizedArea(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor){
    return getStandardizedArea(null,maxIsotope, standMethod, esStandMethod, dilutionFactor);
  }
  
   private double getStandardizedArea(Double area, int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor){
    if (area==null && maxIsotope>=usedIsotpes_ || maxIsotope<0)
      return 0;
    else{
      Double value = area;
      if (value==null)
        value = getOriginalArea(maxIsotope);
      if (dilutionFactor && dilutionFactor_!=null)
        value = value * dilutionFactor_;
      if (standMethod == NO_STANDARD_CORRECTION){
        if (esStandMethod == STANDARD_CORRECTION_INTERNAL){
          value = value * this.correctionFactorESNoISCorrInternal_.get(maxIsotope);
        } else if (esStandMethod == STANDARD_CORRECTION_MEDIAN){
          value = value * this.correctionFactorESNoISCorrMedian_.get(maxIsotope);
        } else if (esStandMethod > 0){
          value = value * this.correctionFactorESNoISCorrSingle_.get(esStandMethod).get(maxIsotope);
        }
      } else if (standMethod == STANDARD_CORRECTION_INTERNAL){
        value = value * correctionFactorInternalIS_.get(maxIsotope);
        if (esStandMethod == STANDARD_CORRECTION_INTERNAL){
          value = value * this.correctionFactorESISInternalCorrInternal_.get(maxIsotope);
        } else if (esStandMethod == STANDARD_CORRECTION_MEDIAN){
          value = value * this.correctionFactorESISInternalCorrMedian_.get(maxIsotope);
        } else if (esStandMethod > 0){
          value = value * this.correctionFactorESISInternalCorrSingle_.get(esStandMethod).get(maxIsotope);
        }
      } else if (standMethod == STANDARD_CORRECTION_MEDIAN){
        value = value * correctionFactorMedianIS_.get(maxIsotope);
        if (esStandMethod == STANDARD_CORRECTION_INTERNAL){
          value = value * this.correctionFactorESISMedianCorrInternal_.get(maxIsotope);
        } else if (esStandMethod == STANDARD_CORRECTION_MEDIAN){
          value = value * this.correctionFactorESISMedianCorrMedian_.get(maxIsotope);
        } else if (esStandMethod > 0){
          value = value * this.correctionFactorESISMedianCorrSingle_.get(esStandMethod).get(maxIsotope);
        }
      } else {
        if (correctionFactorSingleIS_ !=null && correctionFactorSingleIS_.containsKey(standMethod)){
          value = value * correctionFactorSingleIS_.get(standMethod).get(maxIsotope);
          if (esStandMethod == STANDARD_CORRECTION_INTERNAL){
            value = value * this.correctionFactorESISSingleCorrInternal_.get(standMethod).get(maxIsotope);
          } else if (esStandMethod == STANDARD_CORRECTION_MEDIAN){
            value = value * this.correctionFactorInternalESISSingleCorrMedian_.get(standMethod).get(maxIsotope);
          } else if (esStandMethod > 0){
            if (correctionFactorESISSingleCorrSingle_.containsKey(standMethod)&& correctionFactorESISSingleCorrSingle_.get(standMethod).containsKey(esStandMethod)){
              value = value * this.correctionFactorESISSingleCorrSingle_.get(standMethod).get(esStandMethod).get(maxIsotope);
            }else
              return 0;
          }
        }else
          return 0;
      }
      if (value>0 && !value.isInfinite() && !value.isNaN())
        return value;
      else
        return 0;
    }
  }
   
  
  public double getRelativeValue(int maxIsotope, ResultDisplaySettingsVO settingVO)
  {
  	return getRelativeValue(maxIsotope, settingVO, SUM_COMPOSITION, 1.0);
  }
  
  /**
   * Attention: the median areas have to be calculated before and set via addRelativeMedianArea(ArrayList<Double> relativeMedianArea),
   * before the getRelativeValue can be called!!!
   * @param maxIsotope
   * @param settingVO
   * @param molecularSpecies
   * @param molecularSpeciesContribution
   * @return
   */
  public double getRelativeValue(int maxIsotope, ResultDisplaySettingsVO settingVO, String molecularSpecies, double molecularSpeciesContribution)
  {
  	if (maxIsotope>=usedIsotpes_ || maxIsotope<0)
      return -1;
    else{
      
      Double standardizedArea = null;
      try {
        standardizedArea = getArea(maxIsotope, settingVO) * molecularSpeciesContribution;
      }
      catch (CalculationNotPossibleException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      Double relativeMedianArea = relativeMedianAreas_.get(molecularSpecies).get(maxIsotope);
//      if (dilutionFactor && dilutionFactor_ != null)
//        relativeMedianArea = relativeMedianArea/dilutionFactor_;
      return standardizedArea/relativeMedianArea;
    }
  }
  
  /**
   * Adds a relative median area. 
   * @param molecularSpecies				the molecular species this relativeMedianArea belongs to
   * @param relativeMedianArea			the relative median area per isotope		
   */
  public void addRelativeMedianArea(String molecularSpecies, ArrayList<Double> relativeMedianArea)
  {
  	this.relativeMedianAreas_.put(molecularSpecies, relativeMedianArea);
  }
  
  public void addRelativeMedianArea(ArrayList<Double> relativeMedianArea)
  {
  	addRelativeMedianArea(SUM_COMPOSITION, relativeMedianArea);
  }

  public int getUsedIsotpes()
  {
    return usedIsotpes_;
  }

  public String getAbsoluteFilePath()
  {
    return absoluteFilePath_;
  }
  
  public double getRelativeToMedian(int maxIsotope, ResultDisplaySettingsVO settingVO){
    return this.getRelativeToMedian(null, maxIsotope, settingVO);
  }
  
  private double getRelativeToMedian(Double area, int maxIsotope, ResultDisplaySettingsVO settingVO){
    if (area==null && (maxIsotope>=usedIsotpes_ || maxIsotope<0))
      return 0;
    else{
      Double standardizedArea = this.getStandardizedArea(area, maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod(), settingVO.considerDilution());
      if (standardizedArea>0){
        double standardArea = getStandardArea(settingVO, maxIsotope);
        return standardizedArea/standardArea;
      }else
        return 0;
    }
  }
  
  private double getCorrectionFactor(ResultDisplaySettingsVO settingVO, int maxIsotope){   
    ResultDisplaySettingsVO esAreaSet = new ResultDisplaySettingsVO(settingVO.getType(),NO_STANDARD_CORRECTION,settingVO.getESStandMethod(),false,true);
    ResultDisplaySettingsVO isAreaSet = new ResultDisplaySettingsVO(settingVO.getType(),settingVO.getISStandMethod(),NO_STANDARD_CORRECTION,false,true);
    double correctionFactor = getStandardArea(esAreaSet,maxIsotope)/getStandardArea(isAreaSet,maxIsotope);
    return correctionFactor;
  }
  
  private double getStandardArea(ResultDisplaySettingsVO settingVO, int maxIsotope){
    Double standardArea = 1d;
    if (settingVO.getISStandMethod() != NO_STANDARD_CORRECTION && settingVO.getESStandMethod() != NO_STANDARD_CORRECTION){
      if (settingVO.getISStandMethod() == STANDARD_CORRECTION_INTERNAL){
        if (settingVO.getESStandMethod() == STANDARD_CORRECTION_INTERNAL)
          standardArea = areaESISInternalCorrInternalComparison_.get(maxIsotope);
        else if (settingVO.getESStandMethod() == STANDARD_CORRECTION_MEDIAN)
          standardArea = areaESISInternalCorrMedianComparison_.get(maxIsotope);
        else
          standardArea = areaESISInternalCorrSingleComparison_.get(settingVO.getESStandMethod()).get(maxIsotope);
      } else if (settingVO.getISStandMethod() == STANDARD_CORRECTION_MEDIAN){
        if (settingVO.getESStandMethod() == STANDARD_CORRECTION_INTERNAL)
          standardArea = areaESISMedianCorrInternalComparison_.get(maxIsotope);
        else if (settingVO.getESStandMethod() == STANDARD_CORRECTION_MEDIAN)
          standardArea = areaESISMedianCorrMedianComparison_.get(maxIsotope);
        else
          standardArea = areaESISMedianCorrSingleComparison_.get(settingVO.getESStandMethod()).get(maxIsotope);
      } else {
        if (settingVO.getESStandMethod() == STANDARD_CORRECTION_INTERNAL)
          standardArea = areaESISSingleCorrInternalComparison_.get(settingVO.getISStandMethod()).get(maxIsotope);
        else if (settingVO.getESStandMethod() == STANDARD_CORRECTION_MEDIAN)
          standardArea = areaESISSingleCorrMedianComparison_.get(settingVO.getISStandMethod()).get(maxIsotope);
        else
          standardArea = areaESISSingleCorrSingleComparison_.get(settingVO.getISStandMethod()).get(settingVO.getESStandMethod()).get(maxIsotope);
      }
    } else if (settingVO.getISStandMethod() != NO_STANDARD_CORRECTION){
      if (settingVO.getISStandMethod() == STANDARD_CORRECTION_INTERNAL)
        standardArea = areaISInternalComparison_.get(maxIsotope);
      else if (settingVO.getISStandMethod() == STANDARD_CORRECTION_MEDIAN)
        standardArea = areaISMedianComparison_.get(maxIsotope);
      else
        standardArea = areaISSingleComparison_.get(settingVO.getISStandMethod()).get(maxIsotope);
    } else if (settingVO.getESStandMethod() != NO_STANDARD_CORRECTION){
      if (settingVO.getESStandMethod() == STANDARD_CORRECTION_INTERNAL)
        standardArea = areaESNoISCorrInternalComparison_.get(maxIsotope);
      else if (settingVO.getESStandMethod() == STANDARD_CORRECTION_MEDIAN)
        standardArea = areaESNoISCorrMedianComparison_.get(maxIsotope);
      else
        standardArea = areaESNoISCorrSingleComparison_.get(settingVO.getESStandMethod()).get(maxIsotope);
    }
    // the ES is already corrected with the dilution factor, thus we have to correct the dilution again
    // if we do not select it
    if (settingVO.getESStandMethod() != NO_STANDARD_CORRECTION && !settingVO.considerDilution())
      standardArea = standardArea/this.dilutionFactor_;
    return standardArea;
  }

  
//  public double getRelativeToMedianOfISRatios(int maxIsotope, int standMethod, int esMethod, boolean dilutionFactor)
//  {
//    if (maxIsotope>=usedIsotpes_)
//      return 0;
//    else{
//      Double standardizedArea = this.getStandardizedArea(maxIsotope, standMethod, esMethod, dilutionFactor);
//      if (standardizedArea>0){
//        Double standardArea = areaISMedianComparison_.get(maxIsotope);
//        return standardizedArea/standardArea;
//      }else
//        return 0;
//    }  
//  }
//
//  public double getRelativeToISMedian(int maxIsotope, int standMethod, int esMethod, boolean dilutionFactor)
//  {
//    if (maxIsotope>=usedIsotpes_)
//      return 0;
//    else{
//      Double standardizedArea = this.getStandardizedArea(maxIsotope, standMethod, esMethod, dilutionFactor);
//      if (standardizedArea>0){
//        Double standardArea = areaISInternalComparison_.get(maxIsotope);
//        return standardizedArea/standardArea;
//      }else
//        return 0;
//    }
//  }
  public double getAmountInEndVolume(int maxIsotope, int standMethod) throws CalculationNotPossibleException {
    return this.getAmountInEndVolume(null,maxIsotope, standMethod);
  }
  
  private double getAmountInEndVolume(Double area, int maxIsotope, int standMethod) throws CalculationNotPossibleException {
    if (area==null && (maxIsotope>=usedIsotpes_ || maxIsotope<0))
      return 0;
    else{
      if (standMethod == NO_STANDARD_CORRECTION)
        throw new CalculationNotPossibleException("For the calculation of the amount of the probe in the end volume an internal standard has to be used");
      double amountStandard = this.getAmountISStandard(standMethod);
      double amount = 0;
      ResultDisplaySettingsVO settingVO = new ResultDisplaySettingsVO(null,standMethod, NO_STANDARD_CORRECTION,false,false);
      amount = amountStandard*this.getRelativeToMedian(area,maxIsotope, settingVO);
//      if (standMethod == NO_STANDARD_CORRECTION){
//        throw new CalculationNotPossibleException("The method \"NO_STANDARD_CORRECTION\" is not possible for the calculation of the amount in the end volume");
//      } else if (standMethod == STANDARD_CORRECTION_INTERNAL){
//        amount = amountStandard*this.getRelativeToISMedian(maxIsotope, standMethod, NO_STANDARD_CORRECTION,false);
//        
//      } else if (standMethod == STANDARD_CORRECTION_MEDIAN){
//        amount = amountStandard*this.getRelativeToMedianOfISRatios(maxIsotope, standMethod, NO_STANDARD_CORRECTION, false);
//      }else{
//        amount = amountStandard*this.getRelativeToMedianOfISRatios(maxIsotope, standMethod, NO_STANDARD_CORRECTION, false);
//      }
      return amount;
    }
  }
  
  public double getConcentrationInEndVolume(int maxIsotope, int standMethod) throws CalculationNotPossibleException {
    if (maxIsotope>=usedIsotpes_ || maxIsotope<0)
      return 0;
    else{
      double amountInEndVolume = this.getAmountInEndVolume(maxIsotope, standMethod);
      if (endVolume_ == null)
        throw new CalculationNotPossibleException("For the calculation of the concentration of the probe in the end volume, the end volume has to be known");
      return amountInEndVolume/endVolume_;
    }
  }
  
  public double getWeightInEndVolume(int maxIsotope, int standMethod) throws CalculationNotPossibleException {
    if (maxIsotope>=usedIsotpes_ || maxIsotope<0)
      return 0;
    else{
      double amountInEndVolume = this.getAmountInEndVolume(maxIsotope, standMethod);
      if (endVolume_ == null)
        throw new CalculationNotPossibleException("For the calculation of the concentration of the probe in the end volume, the end volume has to be known");
      return amountInEndVolume*getMass(maxIsotope);
    }
  }
  
  public double getAmountInProbeVolume(int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException {
    return this.getAmountInProbeVolume(null,maxIsotope, standMethod, esMethod);
  }

  public double getAmountInProbeVolume(Double area, int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException {
    if (area==null && (maxIsotope>=usedIsotpes_ || maxIsotope<0))
      return 0;
    else{
      if (standMethod == NO_STANDARD_CORRECTION && esMethod == NO_STANDARD_CORRECTION)
        throw new CalculationNotPossibleException("For the calculation of absolute values the external or the internal standard must be respected; best both of them");
      double amountInProbeVolume = 0;
      if (standMethod != NO_STANDARD_CORRECTION && esMethod != NO_STANDARD_CORRECTION){
//        double amountInEndVolume = this.getAmountInEndVolume(area, maxIsotope, standMethod);
//        double idealAmountESEndVolume = 1d;
//        if (esMethod == STANDARD_CORRECTION_INTERNAL){
//          idealAmountESEndVolume = esStandardVolumeInternalCorr_*esStandardConcentrationInternalCorr_;
//        }else if (esMethod == STANDARD_CORRECTION_MEDIAN){
//          idealAmountESEndVolume = esStandardVolumeMedianCorr_*esStandardConcentrationMedianCorr_;
//        } else {
//          idealAmountESEndVolume = esStandardVolumeSingleCorr_.get(esMethod)*esStandardConcentrationSingleCorr_.get(esMethod);
//        }
//        idealAmountESEndVolume=idealAmountESEndVolume/dilutionFactor_;
//        double idealAmountISInEndVolume = isStandardVolume_*isStandardConcentration_;
//        double idealRatio =  idealAmountESEndVolume/idealAmountISInEndVolume;
//        double isAreaInEndVolume = 1d;
//        double esAreaInEndVolume = 1d;
//        if ( standMethod == STANDARD_CORRECTION_INTERNAL){
//          isAreaInEndVolume = areaISInternalComparison_.get(maxIsotope);
//          if (esMethod == STANDARD_CORRECTION_INTERNAL)
//            esAreaInEndVolume = areaESISInternalCorrInternalComparison_.get(maxIsotope);
//          else if (esMethod == STANDARD_CORRECTION_MEDIAN)
//            esAreaInEndVolume = areaESISInternalCorrMedianComparison_.get(maxIsotope);
//          else
//            esAreaInEndVolume = areaESISInternalCorrSingleComparison_.get(esMethod).get(maxIsotope);
//        } else if (standMethod == STANDARD_CORRECTION_MEDIAN){
//          isAreaInEndVolume = areaISMedianComparison_.get(maxIsotope);
//          if (esMethod == STANDARD_CORRECTION_INTERNAL)
//            esAreaInEndVolume = areaESISMedianCorrInternalComparison_.get(maxIsotope);
//          else if (esMethod == STANDARD_CORRECTION_MEDIAN)
//            esAreaInEndVolume = areaESISMedianCorrMedianComparison_.get(maxIsotope);
//          else
//            esAreaInEndVolume = areaESISMedianCorrSingleComparison_.get(esMethod).get(maxIsotope);
//        } else {
//          isAreaInEndVolume = areaISSingleComparison_.get(standMethod).get(maxIsotope);
//          if (esMethod == STANDARD_CORRECTION_INTERNAL)
//            esAreaInEndVolume = areaESISSingleCorrInternalComparison_.get(standMethod).get(maxIsotope);
//          else if (esMethod == STANDARD_CORRECTION_MEDIAN)
//            esAreaInEndVolume = areaESISSingleCorrMedianComparison_.get(standMethod).get(maxIsotope);
//          else
//            esAreaInEndVolume = areaESISSingleCorrSingleComparison_.get(standMethod).get(esMethod).get(maxIsotope);
//        }
//        esAreaInEndVolume = esAreaInEndVolume/dilutionFactor_;
//        double realRatio = esAreaInEndVolume/isAreaInEndVolume;
//        double correctiveFactor = idealRatio/realRatio;
//        amountInProbeVolume = amountInEndVolume*dilutionFactor_*correctiveFactor;
        double amountStandard = this.getAmountISStandard(standMethod);
//////        double amountStandard = this.getAmountESStandard(standMethod);
        ResultDisplaySettingsVO settingVO = new ResultDisplaySettingsVO(null,standMethod, esMethod,true,false);
        //this correction is necessary, because we continue our calculation with the amount of the internal standard,
        //while RelativeToMedian is calculated in reference to the external standard; there the correction of  
        //a different amount in internal standard volume is neglected
        if (standMethod>STANDARD_CORRECTION_MEDIAN && esMethod!=NO_STANDARD_CORRECTION)
          amountStandard = amountStandard*(isAmounts_.get(STANDARD_CORRECTION_INTERNAL).getAmount()/isAmounts_.get(settingVO.getISStandMethod()).getAmount());
        double correctiveFactor = getCorrectionFactor(settingVO, maxIsotope);
        amountInProbeVolume = amountStandard*correctiveFactor*this.getRelativeToMedian(area,maxIsotope, settingVO)*dilutionFactor_;
      }else if (standMethod == NO_STANDARD_CORRECTION && esMethod != NO_STANDARD_CORRECTION){
        Double relationAreaProbeToES = 1d;
        Double value = this.getStandardizedArea(area,maxIsotope, NO_STANDARD_CORRECTION, esMethod, true);
        if (esMethod == STANDARD_CORRECTION_INTERNAL){
          relationAreaProbeToES = value/areaESNoISCorrInternalComparison_.get(maxIsotope);
        } else if (esMethod == STANDARD_CORRECTION_MEDIAN){
          relationAreaProbeToES = value/areaESNoISCorrMedianComparison_.get(maxIsotope);
        } else {
          relationAreaProbeToES = value/areaESNoISCorrSingleComparison_.get(esMethod).get(maxIsotope);
        }
        double amountESProbeVolume = this.getAmountESStandard(esMethod);//esStandardVolumeNoCorr_*esStandardConcentrationNoCorr_;
        amountInProbeVolume = amountESProbeVolume*relationAreaProbeToES;
      }else if (standMethod != NO_STANDARD_CORRECTION && esMethod == NO_STANDARD_CORRECTION){
        amountInProbeVolume = this.getAmountInEndVolume(area,maxIsotope, standMethod)*dilutionFactor_;
      }
      return amountInProbeVolume;
//      return amountInEndVolume*this.getCorrectionFactorForProbeExtractionLosses(maxIsotope, standMethod)*this.dilutionFactor_;
    }
  }
  
  public double getConcentrationInProbeVolume(int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException {
    if (this.probeVolume_ ==null)
      throw new CalculationNotPossibleException("For the calculation of the concentration of the probe in the probe volume, the probe volume has to be known");
    double amountInProbeVolume = this.getAmountInProbeVolume(maxIsotope, standMethod, esMethod);
    return amountInProbeVolume/this.probeVolume_;
  }
  
  public double getWeightInProbeVolume(int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException {
    double amountInProbeVolume = this.getAmountInProbeVolume(maxIsotope, standMethod, esMethod);
    return amountInProbeVolume*this.getMass(maxIsotope);
  }
  
  public double getAnalyteInRelationToProtein(int maxIsotope, int standMethod, int esMethod, boolean useAU) throws CalculationNotPossibleException {
    double amountInProbeVolume = 0d;
    if ((standMethod != NO_STANDARD_CORRECTION || esMethod != NO_STANDARD_CORRECTION) && (!useAU))
      amountInProbeVolume = getAmountInProbeVolume(maxIsotope, standMethod, esMethod);
    else
      amountInProbeVolume = getStandardizedArea(maxIsotope, standMethod, esMethod, true);
    return amountInProbeVolume/getProteinAmount();
  }
  
  public double getAnalyteInRelationToNeutralLipidContent(int maxIsotope, int standMethod, int esMethod, boolean useAU) throws CalculationNotPossibleException {
    double amountInProbeVolume = 0d;
    if ((standMethod != NO_STANDARD_CORRECTION || esMethod != NO_STANDARD_CORRECTION) && (!useAU))
      amountInProbeVolume = getAmountInProbeVolume(maxIsotope, standMethod, esMethod);
    else
      amountInProbeVolume = getStandardizedArea(maxIsotope, standMethod, esMethod, true);
    return amountInProbeVolume/getNeutralLipidAmount();
  }
  
  public double getAnalyteInRelationToSampleWeight(int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException {
    double amountInProbeVolume = getAmountInProbeVolume(maxIsotope, standMethod, esMethod);
    return amountInProbeVolume/getRelativeSampleWeightValue();
  }
  
  private double getRelativeSampleWeightValue(){
    return sampleWeight_;
  }
  
  private double getProteinAmount() throws CalculationNotPossibleException {
    if (this.probeVolume_ ==null)
      throw new CalculationNotPossibleException("For the calculation of the relation to the protein, the probe volume has to be known");
    if (this.proteinConcentration_ ==null)
      throw new CalculationNotPossibleException("For the calculation of the relation to the protein, the protein concentration has to be known");
    return probeVolume_ * proteinConcentration_;
  }
  
  private double getNeutralLipidAmount() throws CalculationNotPossibleException {
    if (this.probeVolume_ ==null)
      throw new CalculationNotPossibleException("For the calculation of the relation to the protein, the probe volume has to be known");
    if (this.neutralLipidConcentration_ ==null)
      throw new CalculationNotPossibleException("For the calculation of the relation to the neutral lipid contents, the neutral lipid contents has to be known");
    return probeVolume_ * neutralLipidConcentration_;
  }
  
//  private double getCorrectionFactorForProbeExtractionLosses(int maxIsotope, int standMethod) throws CalculationNotPossibleException  {
//    double correctionFactor = 1d;
//    double shouldBeRatioEStoIS = this.getAmountESStandardInEndVolume()/this.getAmountISStandard(); 
//    if (maxIsotope<usedIsotpes_){
//      if (standMethod == NO_STANDARD_CORRECTION){
//        throw new CalculationNotPossibleException("The method \"NO_STANDARD_CORRECTION\" is not possible for the calculation of correction factor for probe losses");
//      } else if (standMethod == STANDARD_CORRECTION_INTERNAL){
//        double measuredRatioEStoIS = areaESInternalComparison_.get(maxIsotope)/areaISInternalComparison_.get(maxIsotope);
//        correctionFactor = shouldBeRatioEStoIS/measuredRatioEStoIS;
//      } else if (standMethod == STANDARD_CORRECTION_MEDIAN){
//        double measuredRatioEStoIS = areaESMedianComparison_.get(maxIsotope)/areaESMedianComparison_.get(maxIsotope);
//        correctionFactor = shouldBeRatioEStoIS/measuredRatioEStoIS;
//      }else{
//        throw new CalculationNotPossibleException("The method "+standMethod+" is not applicable for the calculation of the correction factor for probe losses");
//      }
//    }
//    return correctionFactor;
//  }
  
  private double getAmountISStandard(int method) throws CalculationNotPossibleException {
    if (hasAbs_){
      VolumeConcVO concVO = isAmounts_.get(method);
      if (concVO.getVolume() == null)
        throw new CalculationNotPossibleException("For the calculation of the amount of the probe in the volume of the internal standard has to be entered");
      if (concVO.getConcentration() == null)
        throw new CalculationNotPossibleException("For the calculation of the amount of the probe in the concentration of the internal standard has to be entered");
      return concVO.getAmount();
    }
    else return 1;
  }
  
  private double getAmountESStandard(int method) throws CalculationNotPossibleException {
    if (hasAbs_){
      VolumeConcVO concVO = esAmounts_.get(method);
      if (concVO.getVolume() == null)
        throw new CalculationNotPossibleException("For the calculation of the amount of the probe in the volume of the internal standard has to be entered");
      if (concVO.getConcentration() == null)
        throw new CalculationNotPossibleException("For the calculation of the amount of the probe in the concentration of the internal standard has to be entered");
      return concVO.getAmount();
    }else
      return 1;
  }
  
//  private double getAmountESStandardInEndVolume() throws CalculationNotPossibleException {
//    if (esStandardVolume_ == null)
//      throw new CalculationNotPossibleException("For the calculation of the external standard amount the volume of the standard has to be known");
//    if (esStandardConcentration_ == null)
//      throw new CalculationNotPossibleException("For the calculation of the external standard amount the concentration of the standard has to be known");
//    return (esStandardVolume_*esStandardConcentration_)/this.dilutionFactor_;
//  }
  
  //TODO: result comp is composed of one result area, which is composed of multiple params - we need to compute rel. contribution of each mol. species to a result area!
  public double getArea(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException{
    double area = 0;
    if (settingVO.isPercent()){
      area = getRatioToPercentualValue(maxIsotope,settingVO);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_VALUE)){
      area = getStandardizedArea(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod(), settingVO.considerDilution());
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_BASE_PEAK)){
      area = getRatioToHighestPeak(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT)){
      area = getRatioToTotalIntensity(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_HIGHEST_TOTAL_PEAK)){
      area = getRatioToHighestFoundPeak(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT)){
      area = this.getRatioToOverallGroupsIntensity(maxIsotope);
    } else if (settingVO.getType().equalsIgnoreCase("amount end-volume")){
      area = getAmountInEndVolume(maxIsotope, settingVO.getISStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("conc. end-volume")){
      area = getConcentrationInEndVolume(maxIsotope, settingVO.getISStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("weight end-volume")){
      area = getWeightInEndVolume(maxIsotope, settingVO.getISStandMethod());      
    } else if (settingVO.getType().equalsIgnoreCase("amount sample-volume")){
      area = getAmountInProbeVolume(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("conc. sample-volume")){
      area = getConcentrationInProbeVolume(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("weight sample-volume")){
      area = getWeightInProbeVolume(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("relative to sample weight")){
      area = getAnalyteInRelationToSampleWeight(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("relation to protein content")){
      area = getAnalyteInRelationToProtein(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod(),settingVO.isAu());
    } else if (settingVO.getType().equalsIgnoreCase("relation to neutral lipid content")){
      area = getAnalyteInRelationToNeutralLipidContent(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod(),settingVO.isAu());
    } else if (settingVO.getType().equalsIgnoreCase("relation to measured neutral lipid")){
      area = this.getAnalyteInRelationToMeasuredNeutralLipidContent(maxIsotope, settingVO.getISStandMethod(), settingVO.getESStandMethod());
    } else if (settingVO.getType().equalsIgnoreCase("percentual value")){
      area = getRatioToPercentualValue(maxIsotope,settingVO);
    }  
    if (area>0 && !Double.isInfinite(area) && !Double.isNaN(area)){
      if (settingVO.getType().equalsIgnoreCase("conc. end-volume")||settingVO.getType().equalsIgnoreCase("conc. sample-volume")||
          settingVO.getType().equalsIgnoreCase("relative to sample weight")||settingVO.getType().equalsIgnoreCase("relation to protein content")||
          settingVO.getType().equalsIgnoreCase("relation to neutral lipid content")||settingVO.getType().equalsIgnoreCase("relation to measured neutral lipid"))
        area = StaticUtils.getValueDividedByUnit(area,settingVO.getDivisorMagnitude());
      return area;
    }  
    else return 0;
  }

  public double getRatioToHighestPeak(int maxIsotope){
    if (highestGroupIntensity_!=null && maxIsotope>=0 && highestGroupIntensity_.get(maxIsotope)>0){
      return originalArea_.get(maxIsotope)/highestGroupIntensity_.get(maxIsotope);
    }else
      return 0;
  }
  
  public double getRatioToTotalIntensity(int maxIsotope){
    if (totalGroupIntensity_!=null && maxIsotope>=0 && totalGroupIntensity_.get(maxIsotope)>0){
      return originalArea_.get(getAvailableIsotopeNr(maxIsotope))/totalGroupIntensity_.get(maxIsotope);
    } else
      return 0;
  }
  
  public double getAnalyteInRelationToMeasuredNeutralLipidContent(int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException{
    if (totalGroupMass_!=null && maxIsotope>=0 && totalGroupMass_.size()>0){
      double amount = getAmountInProbeVolume(maxIsotope, standMethod, esMethod);
      double totalAmountGramm = this.getAmountInProbeVolume(totalGroupMass_.get(maxIsotope), maxIsotope, standMethod, esMethod);
      return amount/totalAmountGramm;
    }else
      return 0;
  }
  
  public double getRelativeToMeasuredNeutralLipidContent(int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException{
    return this.getAnalyteInRelationToMeasuredNeutralLipidContent(maxIsotope, standMethod, esMethod)*mass_.get(maxIsotope);
  }
  
//  public double get
  
  public int getAvailableIsotopeNr(int desiredIsotopeNr){
    int isotopesToTake = desiredIsotopeNr;
    if (isotopesToTake>(usedIsotpes_-1))
      isotopesToTake = (usedIsotpes_-1);
    return isotopesToTake;
  }

  public double getMass(int maxIsotope)
  {
    if (maxIsotope<0||isEmptyObject())
      return 0d;
    else
      return mass_.get(maxIsotope);
  }
  
  public double getRetentionTime(String modName)
  {
    double value = -1;
    if (retentionTime_.containsKey(modName))
      value = retentionTime_.get(modName);
    return value;
  }

  public void setHighestGroupIntensity(Vector<Double> highestGroupIntensity)
  {
    this.highestGroupIntensity_ = highestGroupIntensity;
  }


  public void setTotalGroupIntensity(Vector<Double> totalGroupIntensity)
  {
    this.totalGroupIntensity_ = totalGroupIntensity;
  }


  public void setTotalGroupMass(Vector<Double> totalGroupMass)
  {
    this.totalGroupMass_ = totalGroupMass;
  }


  public double getTotalGroupMass(int maxIsotope, int standMethod, int esMethod) throws CalculationNotPossibleException
  {
    return this.getAmountInProbeVolume(totalGroupMass_.get(maxIsotope), maxIsotope, standMethod, esMethod);
  }



  public boolean getMoreThanOnePeak(int maxIsotope)
  {
    if (maxIsotope>=usedIsotpes_ || maxIsotope<0){
      return false;
    }else{
      boolean moreThanOnePeak = false;
      for (boolean more : moreThanOnePeak_.get(maxIsotope).values()){
        if (more) moreThanOnePeak = more;
      }
      return moreThanOnePeak;
    }  
  }
  
  public boolean getMoreThanOnePeak(int maxIsotope, String modName)
  {
    if (maxIsotope>=usedIsotpes_ || maxIsotope<0){
      return false;
    }else{
      Hashtable<String,Boolean> hash = moreThanOnePeak_.get(maxIsotope);
      if (hash.containsKey(modName))
        return hash.get(modName);
      else 
        return false;
    }  
  }


  public void setHighestFoundIntensity(Vector<Double> highestFoundIntensity)
  {
    this.highestFoundIntensity_ = highestFoundIntensity;
  }



  public void setTotalFoundIntensity(Vector<Double> totalFoundIntensity)
  {
    this.totalFoundIntensity_ = totalFoundIntensity;
  }

  public double getRatioToHighestFoundPeak(int maxIsotope){
    if (highestFoundIntensity_!=null && maxIsotope>=0 && highestFoundIntensity_.get(maxIsotope)>0){
      return originalArea_.get(maxIsotope)/highestFoundIntensity_.get(maxIsotope);
    }else
      return 0;
  }
  
  public double getRatioToOverallGroupsIntensity(int maxIsotope){
    if (totalFoundIntensity_!=null && maxIsotope>=0 && totalFoundIntensity_.get(maxIsotope)>0){
      return originalArea_.get(getAvailableIsotopeNr(maxIsotope))/totalFoundIntensity_.get(maxIsotope);
    } else
      return 0;
  }
  
  public double getRatioToPercentualValue(int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException{
    if (sumValueForPercentage_!=null && maxIsotope>=0 && sumValueForPercentage_.get(maxIsotope)>0){
      double area = originalArea_.get(getAvailableIsotopeNr(maxIsotope));
      if (type_!=CLASS_TYPE){
        ResultDisplaySettingsVO newSetting = new ResultDisplaySettingsVO(settingVO);
        newSetting.setPercent(false);
        area = this.getArea(maxIsotope,newSetting);
      }
      
      return area/sumValueForPercentage_.get(maxIsotope);
    } else
      return 0;
  }
  
  public void setSumValueForPercentage(Vector<Double> sumValueForPercentage)
  {
    this.sumValueForPercentage_ = sumValueForPercentage;
  }
  
//  public double getTotalGroupIntensityStandardized(int maxIsotope, int standMethod, int esStandMethod, boolean dilutionFactor, boolean relativeToStandard){
//    if (totalGroupIntensity_!=null && totalGroupIntensity_.get(maxIsotope)>0){
//      if (relativeToStandard){
//        ResultDisplaySettingsVO settingVO = new ResultDisplaySettingsVO(null,standMethod, esStandMethod,dilutionFactor,true);
//        int isotopes = maxIsotope;
//        if (standMethod !=  NO_STANDARD_CORRECTION && correctionFactorInternalIS_!=null && correctionFactorInternalIS_.size()>0 && maxIsotope>(correctionFactorInternalIS_.size()-1))
//          isotopes = correctionFactorInternalIS_.size()-1;
//        else if (esStandMethod !=  NO_STANDARD_CORRECTION && correctionFactorESNoISCorrInternal_!=null && correctionFactorESNoISCorrInternal_.size()>0 && maxIsotope> correctionFactorESNoISCorrInternal_.size()-1)
//          isotopes = correctionFactorESNoISCorrInternal_.size()-1;
//        return this.getRelativeToMedian(totalGroupIntensity_.get(maxIsotope),isotopes, settingVO);
//      }else
//        return this.getStandardizedArea(totalGroupIntensity_.get(maxIsotope),maxIsotope, standMethod, esStandMethod, dilutionFactor);
//    }else
//      return -1;
//  }
  
  public boolean isEmptyObject(){
    return (retentionTime_.size()==0);
  }
  
  public boolean containsMod(String modName){
    if (isEmptyObject() && existsInFile()) return true;
    return retentionTime_.containsKey(modName);
  }
  
  public boolean existsInFile(){
    if (isEmptyObject())
      return (this.existsInFile_||this.isNullInFile_);
    else
      return true;
  }
  
  /**
   * 
   * @return true if all this are is constructed of all found modifications
   */
  public boolean hasAllModsFound(){
    return allModsFound_;
  }
  
}
