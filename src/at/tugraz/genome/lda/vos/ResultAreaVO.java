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
import java.util.concurrent.ConcurrentHashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.lda.quantification.LipidParameterSet;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ResultAreaVO
{
	/** 
	 * If thread safety is a concern, make sure to replace the HashSet with a thread safe Set. 
	 */
	private Set<LipidParameterSet> lipidParameterSets_ = ConcurrentHashMap.newKeySet();
  private String name_;
  private Integer dbs_;
  private String expName_;
  /** the original retention time values to have a lookup to the objects in the result files*/
  //TODO: this is probably not needed anymore since we have the params...
//  private Hashtable<String,Hashtable<String,String>> allOriginalRts_ = new Hashtable<String,Hashtable<String,String>>();
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
  /** true when this value is an internal standard*/
  private boolean internalStandard_;
  /** true when this value is an external standard*/
  private boolean externalStandard_;
  
  private String oxState_;
  
  private MolecularSpeciesAreaVO molecularSpeciesContributions_;
  
  
  /**
   * public constructor for creating a ResultAreaVO
   * @param lipidParameterSet 					the lipidParameterSet informing about the name, the amount of double bonds and the oxidation state of the analyte as well as the percental split if present
   * @param rt 													an identifier String for the retention time
   * @param expName 										the abbreviated name of the experiment 
   * @param chemicalFormula 						the chemical formula of the analyte
   * @param neutralMass 								the neutral mass of the analyte
   * @param internalStandard 						true when this value is an internal standard
   * @param externalStandard 						true when this value is an external standard
   * @throws ChemicalFormulaException 	thrown when an element is missing in the elementconfig.xml
   */
  public ResultAreaVO(LipidParameterSet lipidParameterSet, String rt, String expName, String chemicalFormula, double neutralMass,
      boolean internalStandard, boolean externalStandard) throws ChemicalFormulaException
  {
  	lipidParameterSets_.add(lipidParameterSet);
    name_ = lipidParameterSet.getName();
    dbs_ = lipidParameterSet.getDoubleBonds();
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
    percentalSplit_ = 1f;
    if(lipidParameterSet.getPercentalSplit()>=1) percentalSplit_ = lipidParameterSet.getPercentalSplit()/100f;
    internalStandard_ = internalStandard;
    externalStandard_ = externalStandard;
    rt_ = rt;
    oxState_ = lipidParameterSet.getOxState();
  }
  
  /**
   * adds a result belonging to the same analyte species (e.g. another adduct/modification)
   * @param lipidParameterSets 					the lipidParameterSets this result belongs to
   * @param modName 										the modification/adduct name
   * @param modFormula 									the chemical formula of the modification/adduct
   * @param mass 												the theoretical m/z value
   * @param expMass 										the measured m/z value
   * @param charge 											the presumed charge
   * @param rt 													the retention time
   * @throws ChemicalFormulaException
   */
  public void addResultPart(Set<LipidParameterSet> lipidParameterSets, String modName, String modFormula,double mass,double expMass,int charge, String rt) throws ChemicalFormulaException{
  	lipidParameterSets_.addAll(lipidParameterSets);
  	String modToAdd = getNotNullName(modName);
    // this is truly a new modification, otherwise this hit is just a peak with a different retention time; no new modification
    if (!modFormulas_.containsKey(modToAdd)){
      modFormulas_.put(modToAdd, StaticUtils.categorizeFormula(getNotNullName(modFormula)));
      mass_.put(modToAdd, mass);
      expMass_.put(modToAdd, expMass);
      charge_.put(modToAdd, charge);
      areas_.put(modToAdd, new Vector<Double>());
      moreThanOnePeak_.put(modToAdd, new Vector<Boolean>());
    }
  }
  
  
  private String getNotNullName(String name){
    if (name!=null&&name.length()>0) return name;
    else return "";
  }
  
  public String getMoleculeName()
  {
    return StaticUtils.generateLipidNameString(name_, dbs_, rt_, oxState_);
  } 
  
  public String getMoleculeNameWoRT()
  {
    return StaticUtils.generateLipidNameString(name_, dbs_,-1, oxState_);
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
  
  /**
   * @return the chemical formula as a hash table
   */
  public Hashtable<String,Integer> getChemicalFormulaElements(){
    return this.chemicalFormula_;
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
          double mass = neutralMass_+LipidomicsConstants.getNeutronMass()*j;
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
  
  /**
   * Adds an area. Make sure the lipidParameterSet this area is taken from is already part of this object, 
   * either via the constructor, the method 'addResultPart', or combineVOs
   * @param lipidParameterSets
   * @param modificationName
   * @param iso
   * @param area
   * @return
   */
  public Hashtable<Integer,Boolean> addArea(Set<LipidParameterSet> lipidParameterSets, String modificationName, int iso, double area){
  	for (LipidParameterSet param : lipidParameterSets)
  	{
  		
  	}
  	lipidParameterSets_.addAll(lipidParameterSets);
  	
  	//TODO: this is the only part, where Areas are actually added, thus, ensure we also get the molecular species and their rel. contributions here. 
  	//TODO: what do we do about different modifications?? C=C??? like... I guess I need modification, human readable, pos insensitive and area all in one object?
  	
    String modName = getNotNullName(modificationName);
    Vector<Double> areas = areas_.get(modName);
    Double totalAreaBefore = getTotalArea(areas.size());
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
    getMolecularSpeciesContributions();
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
  	LipidParameterSet param = lipidParameterSets_.iterator().next();
  	if (param.getNameStringWithoutRt().equals("32:2") 
  			&& param.getChemicalFormula().contains("C41")
  			)
  	{
  		System.out.println(param.getChemicalFormula());
  		System.out.println("hi");
  	}
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
          this.addResultPart(other.getLipidParameterSets(), mod, other.getChemicalFormula(mod), other.mass_.get(mod), other.expMass_.get(mod), other.charge_.get(mod), null);
        }
        catch (ChemicalFormulaException e) {
          e.printStackTrace();
        }
      }
      Vector<Double> areas = other.areas_.get(mod);
      Hashtable<Integer,Boolean> moreThanOnePeak = new Hashtable<Integer,Boolean>();
      for (int i=0; i!=areas.size(); i++){
        Hashtable<Integer,Boolean> mtp = addArea(other.getLipidParameterSets(),mod,i,areas.get(i));
        if (mtp!=null) moreThanOnePeak = mtp;
        if (moreThanOnePeak.containsKey(i))mtp.put(i, true);
        else moreThanOnePeak.put(i, false);
        if (i==0 && areas.get(i)>highestZeroArea){
          setRetentionTime(mod,other.getRetentionTime(mod));
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

  
  /**
   * checks whether this retention time belongs to this ResultAreaVO (one rectangle in the heat map)
   * @param rt the retention time to check
   * @param mod the modification belonging to this retention time
   * @return true when this retention time belongs to this ResultAreaVO
   */
  public boolean belongsRtToThisAreaVO(String rt, String mod)
  {
  	for (LipidParameterSet param : lipidParameterSets_)
  	{
  		if (param.getModificationName().equals(mod) && param.getRt().equals(rt))
  		{
  			return true;
  		}
  	}
  	return false;
  }
  
  
  /**
   * 
   * @return true when this belongs to an internal or an external standard
   */
  public boolean isAStandard(){
    return (internalStandard_ || this.externalStandard_);
  }
  
  protected Set<LipidParameterSet> getLipidParameterSets()
  {
  	return lipidParameterSets_;
  }
  
  private MolecularSpeciesAreaVO getMolecularSpeciesContributions()
  {
  	for (LipidParameterSet param : lipidParameterSets_)
  	{
  		if (param instanceof LipidomicsMSnSet)
  		{
  			LipidomicsMSnSet paramMSn = (LipidomicsMSnSet) param;
  			if (paramMSn.getStatus() > LipidomicsMSnSet.HEAD_GROUP_DETECTED)
  			{
  				Set<String> names = paramMSn.getHumanReadableNameSet();
  				if (names.size() > 1)
  				{
  					Vector<String> chains = paramMSn.getValidChainCombinations();
  					System.out.println("hi");
  					paramMSn.getHumanReadableNameSet();
  				}
  				String mod = paramMSn.getModificationName();
  			}
  		}
  	}
  	return new MolecularSpeciesAreaVO(expName_, expName_, expName_, neutralMass_);
  }
  
  /**
   * Computes a table of position insensitive molecular species names and their relative contribution 
   * to the overall isotope areas of the params that contributed to this object.
   * @return
   */
  public Hashtable<String,Vector<Double>> computeMolecularSpeciesContributions()
  {
  	Hashtable<String,Vector<Double>> contributions = new Hashtable<String,Vector<Double>>();
  	Set<String> names = ConcurrentHashMap.newKeySet();
  	for (LipidParameterSet param : lipidParameterSets_)
  	{
  		if (param instanceof LipidomicsMSnSet)
  		{
  			LipidomicsMSnSet paramMSn = (LipidomicsMSnSet) param;
  			if (paramMSn.getStatus() > LipidomicsMSnSet.HEAD_GROUP_DETECTED)
  			{
  				//in the algorithm, the areas of one modification are only contributed to by one LipidParameterSet 
  				Vector<Double> isoAreas = areas_.get(param.getModificationName());
  				//TODO: get valid chain combinations.
  				Hashtable<String,Vector<Double>> areas = areas_;
  				Vector<Vector<CgProbe>> probes = param.getIsotopicProbes();
  				System.out.println("yo!");
  			}
  		}
  	}
  	return contributions;
  }
  
  private class MolecularSpeciesAreaVO
  {
  	private String humanReadable_;
  	private String positionAndCCInsensitive_;
  	private String modification_;
  	private Double areaContribution_;
  	
  	private MolecularSpeciesAreaVO(String humanReadable, String positionAndCCInsensitive, String modification, Double areaContribution)
  	{
  		this.humanReadable_ = humanReadable;
  		this.positionAndCCInsensitive_ = positionAndCCInsensitive;
  		this.modification_ = modification;
  		this.areaContribution_ = areaContribution;
  	}
  	
  	
  }
  
}
