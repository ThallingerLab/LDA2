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

package at.tugraz.genome.lda.vos;

import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class QuantVO
{
  private String analyteClass_;
  private String analyteName_;
  private int dbs_;
  private String analyteFormula_;
  private double analyteMass_;
  private int charge_;
  private String modName_;
  private String modFormula_;
  private float retTime_;
  private float isobaricRetTime_;
  private float usedMinusTime_;
  private float isobaricMinusTime_;
  private float usedPlusTime_;
  private float isobaricPlusTime_;
  private Vector<Double> mustMatchProbabs_;
  private Vector<Double> probabs_;
  private int negStartValue_;
  /** other isobaric quantitation objects */
  private Vector<QuantVO> isobaricSpecies_;
  /** will this QuantVO be quantified by another isobar?*/
  private boolean quantifiedByOtherIsobar_;
  
  public QuantVO(String analyteClass, String analyteName, int dbs,
      String analyteFormula, double analyteMass, int charge,
      String modName, String modFormula, float retTime,
      float usedMinusTime, float usedPlusTime,
      Vector<Double> mustMatchProbabs, Vector<Double> probabs,
      int negativeStartValue)
  {
    super();
    this.analyteClass_ = analyteClass;
    this.analyteName_ = analyteName;
    this.dbs_ = dbs;
    this.analyteFormula_ = analyteFormula;
    this.analyteMass_ = analyteMass;
    this.charge_ = charge;
    this.modName_ = modName;
    this.modFormula_ = modFormula;
    this.retTime_ = retTime;
    this.isobaricRetTime_ = this.retTime_;
    this.usedMinusTime_ = usedMinusTime;
    this.isobaricMinusTime_ = this.usedMinusTime_;
    this.usedPlusTime_ = usedPlusTime;
    this.isobaricPlusTime_ = this.usedPlusTime_;
    this.mustMatchProbabs_ = mustMatchProbabs;
    this.probabs_ = probabs;
    this.negStartValue_ = negativeStartValue;
    this.isobaricSpecies_ = new Vector<QuantVO>();
    this.quantifiedByOtherIsobar_ = false;
  }
  
  public String getAnalyteClass()
  {
    return analyteClass_;
  }
  public String getAnalyteName()
  {
    return analyteName_;
  }
  public int getDbs()
  {
    return dbs_;
  }
  public String getAnalyteFormula()
  {
    return analyteFormula_;
  }
  public double getAnalyteMass()
  {
    return analyteMass_;
  } 
//  public double getAnalyteStartMass(){
//    return (analyteMass_+massReduction_);
//  }
  
  public int getCharge()
  {
    return charge_;
  }
  public String getModName()
  {
    return modName_;
  }
  public String getModFormula()
  {
    return modFormula_;
  }
  public float getRetTime()
  {
    return retTime_;
  }
  /**
   * setter method for the retention time
   * @param rt the retention time
   */
  public void setRetTime(float rt)
  {
    retTime_ = rt;
  }
  public float getUsedMinusTime()
  {
    return usedMinusTime_;
  }
  public float getUsedPlusTime()
  {
    return usedPlusTime_;
  }
  public Vector<Double> getMustMatchProbabs()
  {
    if (LipidomicsConstants.removeIfDistriDoesNotFit())
      return mustMatchProbabs_;
    else
      return new Vector<Double>();
  }
  public Vector<Double> getProbabs()
  {
    return probabs_;
  }
  public int getNegativeStartValue(){
    return this.negStartValue_;
  }
  
  /**
   * 
   * @return original string for the analyte name
   */
  public String getIdString(){
    return StaticUtils.generateLipidNameString(analyteName_,dbs_);
  }
  
  /**
   * adds a QuantVO of another isobaric species
   * @param isobar the other isobar
   */
  public void addIsobaricSpecies(QuantVO isobar){
    boolean alreadyThere = false;
    for (QuantVO quant : isobaricSpecies_){
      if (isobar.equals(quant)){
        alreadyThere = true;
        break;
      }
    }
    if (!alreadyThere){
      isobaricSpecies_.add(isobar);
      if (this.retTime_>0){
        float startTime = retTime_-this.usedMinusTime_;
        float stopTime = retTime_+this.usedPlusTime_;
        for (QuantVO vo : this.isobaricSpecies_){
          if (vo.retTime_>0){
            if ((vo.retTime_-vo.usedMinusTime_)<startTime) startTime = vo.retTime_-vo.usedMinusTime_;
            if ((vo.retTime_+vo.usedPlusTime_)>stopTime) stopTime = vo.retTime_+vo.usedPlusTime_;
          }else{
            startTime = -1f;
            stopTime = -1f;
            break;
          }
        }
        if (startTime>0 && stopTime>0){
          this.isobaricRetTime_ = -1f;
          this.isobaricPlusTime_ = -1f;
          this.isobaricMinusTime_ = -1f;
        } else {
          this.isobaricRetTime_ = (startTime+stopTime)/2f;
          this.isobaricMinusTime_ = this.isobaricMinusTime_-startTime;
          this.isobaricPlusTime_ = this.isobaricRetTime_-stopTime;
        }
      }
    }
  }

  /**
   * 
   * @return true if another isobar will quantify this object on the MS1 level
   */
  public boolean isQuantifiedByOtherIsobar()
  {
    return quantifiedByOtherIsobar_;
  }

  /**
   * 
   * @param quantifiedByOtherIsobar true, if this species shall be quantified by anohter isobar
   */
  public void setQuantifiedByOtherIsobar(boolean quantifiedByOtherIsobar)
  {
    this.quantifiedByOtherIsobar_ = quantifiedByOtherIsobar;
  }
  
  /**
   * 
   * @return true if there are other isobaric species present
   */
  public boolean hasOtherIsobars(){
    return this.isobaricSpecies_.size()>0;
  }

  /**
   * 
   * @return the other QuantVOs for other isobaric species
   */
  public Vector<QuantVO> getOtherIsobaricSpecies(){
    return this.isobaricSpecies_;
  }
  
  /**
   * removes the information about other isobars
   */
  public void removeOtherIsobaricSpecies(){
    this.isobaricSpecies_.clear();
  }

  public float getIsobaricRetTime_()
  {
    return isobaricRetTime_;
  }

  public float getIsobaricMinusTime_()
  {
    return isobaricMinusTime_;
  }

  public float getIsobaricPlusTime_()
  {
    return isobaricPlusTime_;
  }
  
  
  
}
