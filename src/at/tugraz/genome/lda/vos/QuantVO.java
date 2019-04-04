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
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class QuantVO
{
  protected String analyteClass_;
  protected String prefixOrName_;
  protected int carbons_;
  protected int dbs_;
  protected int oh_;
  protected String analyteFormula_;
  protected double analyteMass_;
  protected int charge_;
  protected String modName_;
  protected String modFormula_;
  protected float retTime_;
  protected float isobaricRetTime_;
  protected float usedMinusTime_;
  protected float isobaricMinusTime_;
  protected float usedPlusTime_;
  protected float isobaricPlusTime_;
  protected Vector<Double> mustMatchProbabs_;
  protected Vector<Double> probabs_;
  protected int negStartValue_;
  /** other isobaric quantitation objects */
  protected Vector<QuantVO> isobaricSpecies_;
  /** will this QuantVO be quantified by another isobar?*/
  protected boolean quantifiedByOtherIsobar_;
  
  /**
   * constructor for an object holding necessary information for 
   * @param analyteClass the name of the lipid class
   * @param analyteName the parsed name of the analyte
   * @param dbs the number of double bonds
   * @param ohNumber the number of hydroxylation sites
   * @param analyteFormula the neutral chemical formula
   * @param analyteMass the mass for this search
   * @param charge the charge
   * @param modName the name of the modification
   * @param modFormula the chemical formula for the modification
   * @param retTime the presumable retention time to search for the analyte
   * @param usedMinusTime the time tolerance (from retTime) in negative direction
   * @param usedPlusTime the time tolerance (from retTime) in positive direction
   * @param mustMatchProbabs the relative intensities of the probabilities of isotopes that must match
   * @param probabs the relative intensities of the probabilities of isotopes to be searched for
   * @param negativeStartValue true when a negative isotopic distribution is assumed
   * @throws HydroxylationEncodingException when there is no encoding for the provided ohNumber 
   */
  public QuantVO(String analyteClass, String analyteName, int dbs,
      int ohNumber, String analyteFormula, double analyteMass, int charge,
      String modName, String modFormula, float retTime,
      float usedMinusTime, float usedPlusTime,
      Vector<Double> mustMatchProbabs, Vector<Double> probabs,
      int negativeStartValue) throws HydroxylationEncodingException
  {
    super();
    this.analyteClass_ = analyteClass;
    //the next line must always come before splitInCarbonNumberAndPrefix to make TargetlistEntry work
    this.dbs_ = dbs;
    Object[] prefixAndC = splitInCarbonNumberAndPrefix(this.analyteClass_, analyteName);
    this.prefixOrName_ = (String)prefixAndC[0];
    this.carbons_ = (Integer)prefixAndC[1];
    this.oh_ = ohNumber;
    //this is a check whether an encoding for this OH number exists
    if (oh_>0)
      Settings.getLcbHydroxyEncoding().getEncodedPrefix((short)oh_);
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
    try {
      return prefixOrName_+(oh_>0 ? Settings.getLcbHydroxyEncoding().getEncodedPrefix((short)oh_) : "")
          +(carbons_>=0 ? String.valueOf(carbons_) : "");
    }catch (HydroxylationEncodingException e) {// this was caught before - this error cannot happen
    }
    return null;
  }
  public int getDbs()
  {
    return dbs_;
  }
  
  public int getOhNumber()
  {
    return oh_;
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
   * @return original string for the analyte name
   */
  public String getIdString(){
    return StaticUtils.generateLipidNameString(getAnalyteName(),dbs_);
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
  
  /**
   * this method splits an arbitrary input name to a string part, and the number of carbons that come at the end
   * @param analyteClass the class of the analyte
   * @param analyteName the name of the analyte (without ":" and followed by the double bonds)
   * @return [0] String containing prefix or analyte name; [1] Integer containing the number of carbons;
   */
  protected Object[] splitInCarbonNumberAndPrefix(String analyteClass, String analyteName) {
    Object[] prefixAndC = new Object[2];
    int cs = LipidomicsConstants.EXCEL_NO_OH_INFO;
    char[] chars = analyteName.toCharArray();
    int splitChar = chars.length;
    while (splitChar>0 && Character.isDigit(chars[splitChar-1])) 
      splitChar--;
    String prefix = analyteName.substring(0,splitChar);
    if (splitChar<chars.length)
      cs = Integer.parseInt(analyteName.substring(splitChar));
    prefixAndC[0] = prefix;
    prefixAndC[1] = cs;
    return prefixAndC;
  }
  
}
