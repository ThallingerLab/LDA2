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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ResultAreaVO
{
  private String name_;
  private Integer dbs_;
  private String expName_;
  private String rtOriginal_;
  private String rt_;
  private Hashtable<String,Integer> chemicalFormula_;
  private Hashtable <String,Hashtable<String,Integer>> modFormulas_; 
  private Double neutralMass_;
  private Hashtable<String,Double> mass_;
  private Hashtable<String,Double> expMass_;
  private Hashtable<String,Integer> charge_;
  private Hashtable<String,Double> retentionTime_;
  private Hashtable<String,Vector<Double>> areas_;
  private Hashtable<String,Vector<Boolean>> moreThanOnePeak_;
  /** if there is a split according to MSn intensities - a percentage of the usable peak intensity is stored*/
  private float percentalSplit_;
  
  public ResultAreaVO(String name, Integer dbs, String rt, String expName, String chemicalFormula, float percentalSplit/*,String modName, String modFormula*/,double neutralMass) throws ChemicalFormulaException
  {
    this(name,dbs,expName,chemicalFormula,percentalSplit,neutralMass);
    this.rtOriginal_ = rt;
    this.rt_ = rt;
  }
    
  public ResultAreaVO(String name, Integer dbs, String expName, String chemicalFormula, float percentalSplit/*,String modName, String modFormula*/,double neutralMass) throws ChemicalFormulaException
  {
    super();
    name_ = name;
    dbs_ = dbs;
    expName_ = expName;
    chemicalFormula_ = StaticUtils.categorizeFormula(chemicalFormula);
    modFormulas_ = new Hashtable <String,Hashtable<String,Integer>>();
    neutralMass_ = neutralMass;
    mass_ = new Hashtable<String,Double>();
    expMass_ = new Hashtable<String,Double>();
    charge_ = new Hashtable<String,Integer>();
    retentionTime_ = new Hashtable<String,Double>();
    areas_ = new Hashtable<String,Vector<Double>>();
    moreThanOnePeak_ = new Hashtable<String,Vector<Boolean>>();
    rtOriginal_ = null;
    rt_ = null;
    percentalSplit_ = 1f;
    if(percentalSplit>=1) percentalSplit_ = percentalSplit/100f;
  }
  
  public void addResultPart(String modName, String modFormula,double mass,double expMass,int charge) throws ChemicalFormulaException{
    String modToAdd = getNotNullName(modName);
    // this is truly a new modification, otherwise this hit is just a peak with a different retention time; no new modification
    if (!modFormulas_.containsKey(modToAdd)){
      String modFormToAdd = getNotNullName(modFormula);
      modFormulas_.put(modToAdd, StaticUtils.categorizeFormula(modFormToAdd));
      mass_.put(modToAdd, mass);
      expMass_.put(modToAdd, expMass);
      charge_.put(modToAdd, charge);
      Vector<Double> isoAreas = new Vector<Double>();
      areas_.put(modToAdd, isoAreas);
      Vector<Boolean> isoMoreThanOne = new Vector<Boolean>();
      moreThanOnePeak_.put(modToAdd, isoMoreThanOne);
    }
  }
  
  private String getNotNullName(String name){
    if (name!=null&&name.length()>0) return name;
    else return "";
  }
  
  public String getMoleculeName()
  {
    return StaticUtils.generateLipidNameString(name_, dbs_, rt_);
  } 
  
  public String getMoleculeNameWoRT()
  {
    return StaticUtils.generateLipidNameString(name_, dbs_);
  }
  
  public String getExpName()
  {
    return expName_;
  }  
//  public String getChemicalFormula()
//  {
//    return chemicalFormula;
//  }
  public double getTheoreticalIsotopeValue(ElementConfigParser elementParser_, int isosDesired){
    double result = 0;
    for (String modString : areas_.keySet()){
      String chemicalFormula = getChemicalFormula(modString);
//      double refValue = currentValue;
      try {
        Vector<Double> distris = StaticUtils.calculateChemicalFormulaIntensityDistribution(elementParser_, chemicalFormula, isosDesired, false);;
        Vector<Double> areas = areas_.get(modString);
        Vector<Double> isoValues = new Vector<Double>();
        for (int i=0; i!=isosDesired;i++){
          if (i<areas.size())
            isoValues.add(areas.get(i));
          else{
            isoValues.add(isoValues.get(0)*distris.get(i));
          }  
          //refValue += zeroIsoValue*distris.get(i);
        }
        result+=getTotalArea(isoValues,isosDesired);
      }
      catch (SpectrummillParserException e) {
        e.printStackTrace();
      }
    }    
    return result;
  }
  
  private String getChemicalFormula(String modString){
    String formula = "";
    Hashtable<String,Integer> modFormula = new Hashtable<String,Integer>();
    if (modFormulas_.containsKey(modString)) modFormula = modFormulas_.get(modString);
    Vector<String> elementsOrder = new Vector<String>(); 
    for (String element: chemicalFormula_.keySet()){
      if (element.equalsIgnoreCase("C")) elementsOrder.add(0,element);
      else elementsOrder.add(element);
    }
    for (String element: elementsOrder){
      if (formula.length()>0) formula+=" ";
      int amount = chemicalFormula_.get(element);
      if (modFormula.containsKey(element)) amount+=modFormula.get(element);
      formula+=element+String.valueOf(amount);
    }
    for (String element : modFormula.keySet()){
      if (!chemicalFormula_.containsKey(element)){
        if (formula.length()>0) formula+=" ";
        formula+=element+String.valueOf(modFormula.get(element));
      }
    }
    return formula;
  }
  
  public String getChemicalFormulaBase(){
    return getChemicalFormula("");
  }
  
  public Vector<Double> getWeightedNeutralMass()
  {
    int maxIsotpe = getMaxIsotope();
    Vector<Double> weightedMasses = new Vector<Double>();
    for (int i=0;i!=maxIsotpe;i++){
      double totalArea = 0;
      double totalMassArea = 0;
//      for (String modName : areas_.keySet()){
//        double area = getTotalArea(areas_.get(modName),i+1);
//        double mass = mass_.get(modName);
//        totalArea += area;
//        totalMassArea += area*mass;
//      }
      for (String modName : areas_.keySet()){
        Vector<Double> areas = areas_.get(modName);
        for (int j=0; j!=(i+1)&&j!=areas.size();j++){
          double area = areas.get(j);
          double mass = /****mass_.get(modName)*/neutralMass_+LipidomicsConstants.getNeutronMass()*j;
          totalArea += area;
          totalMassArea += area*mass;
        }
      }
      weightedMasses.add(totalMassArea/totalArea);
    }
    return weightedMasses;
  }

//  public Vector<Double> getAreas()
//  {
//    return areas;
//  }
  
  public Hashtable<Integer,Boolean> addArea(String modificationName, int iso, double area){
    String modName = getNotNullName(modificationName);
    Vector<Double> areas = areas_.get(modName);
    Vector<Boolean> moreThanOnePeak = moreThanOnePeak_.get(modName);
    boolean isAdditionalPeak = false;
    int isotope = iso;
    if (isotope<0) isotope = isotope*-1;
    if (areas.size()<=isotope){
      areas.add(area*percentalSplit_);
      moreThanOnePeak.add(false);
    }else{
      double totalArea = areas.get(isotope)+area*percentalSplit_;
      areas.remove(isotope);
      areas.add(isotope, totalArea);
      moreThanOnePeak.remove(isotope);
      moreThanOnePeak.add(isotope,true);
      isAdditionalPeak = true;
    }
    areas_.put(modName, areas);
    moreThanOnePeak_.put(modName, moreThanOnePeak);
    if (isAdditionalPeak){
      Hashtable<Integer,Boolean> mtp = new Hashtable<Integer,Boolean>();
      for (int i=0;i!=moreThanOnePeak.size();i++) mtp.put(i, moreThanOnePeak.get(i));
      return mtp;
    }
    return null;
  }
  
  public double getTotalArea(int amountOfIsotopes){
    double totalArea = 0;
    for (Vector<Double> areas : areas_.values())
      totalArea+=getTotalArea(areas,amountOfIsotopes);
    return totalArea;
  }
  
  public double getTotalAreaOfModification(String mod, int amountOfIsotopes){
    return this.getTotalArea(areas_.get(mod),amountOfIsotopes);
  }
  
  private double getTotalArea(Vector<Double> areas, int amountOfIsotopes){
    double totalArea = 0;
    if (areas!=null){
      for (int i=0; i!=areas.size()&&i<amountOfIsotopes; i++){
        totalArea += areas.get(i);
      }
    }
    return totalArea;
  }

  public double getRetentionTime(String modificationName)
  {
    return retentionTime_.get(getNotNullName(modificationName));
  }

  public void setRetentionTime(String modificationName, double retentionTime)
  {
    retentionTime_.put(getNotNullName(modificationName), retentionTime);
  }
  
  public Hashtable<String,Double>  getRetentionTimes()
  {
    return retentionTime_;
  }


  public Hashtable<String,Boolean> getMoreThanOnePeak(int amountOfIsotopes)
  {
//    boolean moreThanOnePeak = false;
//    for (Vector<Boolean> moreThanOnePeaks : moreThanOnePeak_.values())
//      for (int i=0; i!=moreThanOnePeaks.size()&&i<amountOfIsotopes; i++){
//        if (moreThanOnePeaks.get(i))
//          moreThanOnePeak = true;
//    }
//    return moreThanOnePeak;
    Hashtable<String,Boolean> moreThanOnePeakHash = new Hashtable<String,Boolean>();
    for (String modName : moreThanOnePeak_.keySet()){
      boolean moreThanOnePeak = false;
      Vector<Boolean> moreThanOnePeaks = moreThanOnePeak_.get(modName);
      for (int i=0; i!=moreThanOnePeaks.size()&&i<amountOfIsotopes; i++){
        if (moreThanOnePeaks.get(i))
          moreThanOnePeak = true;
      }
      moreThanOnePeakHash.put(modName, moreThanOnePeak);
    }
    return moreThanOnePeakHash;
  }
  
  public void setMoreThanOnePeak(String modificationName, Hashtable<Integer,Boolean> moreThanOnePeak){
    List<Integer> isotopes = new ArrayList<Integer>(moreThanOnePeak.keySet());
    Vector<Boolean> moreThanOnePeaks = new Vector<Boolean>();
    Collections.sort(isotopes);
    for (int i=0;i!=isotopes.size();i++){
      moreThanOnePeaks.add(moreThanOnePeak.get(isotopes.get(i)));
    }
    moreThanOnePeak_.put(modificationName, moreThanOnePeaks);
  }
  
  public int getMaxIsotope(){
    int maxIsotope = 0;
    for (Vector<Double> areas : areas_.values()){
      int size = areas.size();
      if (maxIsotope<size)maxIsotope=size;
    }
    return maxIsotope;
  }

  public void setRt(String rt)
  {
    this.rt_ = rt;
  }

  public String getRt()
  {
    return rt_;
  }

  public String getRtOriginal()
  {
    return rtOriginal_;
  }
  
  public void setRtOriginal(String rt){
    this.rtOriginal_ = rt;
    this.rt_ = rt;
  }
  
  public float getHighestZeroIsoArea(String modification){
    if (areas_.get(modification)!=null && areas_.get(modification).get(0)!=null)
        return areas_.get(modification).get(0).floatValue();
    else
      return 0f;
  }
  
  public boolean hasModification(String mod){
    return modFormulas_.containsKey(mod);
  }
  
  public String getModificationFormula(String mod){
    String modFormula = "";
    if (hasModification(mod)){
      Hashtable<String,Integer> formula = modFormulas_.get(mod);

      for (String element : formula.keySet()){
        String amount = String.valueOf(formula.get(element));
        if (amount.startsWith("-")){
          amount = amount.substring(1);
          modFormula+="-";
        }
        modFormula+=element+amount;
      }
    }
    return modFormula;
  }
  
  public void combineVOs(ResultAreaVO other){
    double highestArea = 0;
    for (String mod : areas_.keySet()){
      if (areas_.get(mod).size()>0 && areas_.get(mod).get(0)>highestArea){
        highestArea = areas_.get(mod).get(0);
      }
    }

    for (String mod : other.mass_.keySet()){
      float highestZeroArea = getHighestZeroIsoArea(mod);
      if (!this.mass_.containsKey(mod)){
        try {
          this.addResultPart(mod, other.getChemicalFormula(mod), other.mass_.get(mod), other.expMass_.get(mod), other.charge_.get(mod));
        }
        catch (ChemicalFormulaException e) {
          e.printStackTrace();
        }
      }
      Vector<Double> areas = other.areas_.get(mod);
      Hashtable<Integer,Boolean> moreThanOnePeak = new Hashtable<Integer,Boolean>();
      for (int i=0; i!=areas.size(); i++){
        Hashtable<Integer,Boolean> mtp = addArea(mod,i,areas.get(i));
        if (mtp!=null) moreThanOnePeak = mtp;
        if (moreThanOnePeak.containsKey(i))mtp.put(i, true);
        else moreThanOnePeak.put(i, false);
        if (i==0 && areas.get(i)>highestZeroArea){
          setRetentionTime(mod,other.getRetentionTime(mod));
        }
        if (i==0 && areas.get(i)>highestArea){
          setRtOriginal(other.rtOriginal_);
        }
      }
      setMoreThanOnePeak(mod,moreThanOnePeak);
    }
  }

  public Integer getCharge(String modification){
    if (hasModification(modification)) return charge_.get(modification);
    else return null;
  }
  
  public Double getTheoreticalMass(String modification){
    if (hasModification(modification)) return mass_.get(modification);
    else return null;
  }
  
  public Double getExperimentalMass(String modification){
    if (hasModification(modification)) return expMass_.get(modification);
    else return null;
  }
  
  /**
   * checks if all of the modification names in the hash table are really present in this VO
   * @param mods the modification names stored in a hash
   * @return true if all modifications are present
   */
  public boolean containsAllModifications(Hashtable<String,String> mods){
    boolean all = true;
    for (String mod : mods.keySet()){
      if (!this.areas_.containsKey(mod)){
        all = false;
        break;
      }
    }
    return all;
  }
  
}
