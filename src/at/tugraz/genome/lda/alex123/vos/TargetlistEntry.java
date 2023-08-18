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

package at.tugraz.genome.lda.alex123.vos;

import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.exception.AlexTargetlistParserException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.QuantVO;

/**
 * VO containing the Alex123 target list entries, extending the typical LDA object for quantitation
 * @author Juergen Hartler
 *
 */
public class TargetlistEntry extends QuantVO
{
  /** Alex123 class name*/
  private String originalClassName_;
  /** name of the detector*/
  private String detector_ = null;
  /** the polarity of this ion*/
  private String polarity_ = null;
  /** the msLevel of this observance*/
  private int msLevel_ = -1;
  /** the name of the fragment*/
  private String fragment_ = null;
  /** the verified lipid structure*/
  private String structure_ = null;
  /** the name of this species*/
  private String species_ = null;
  /** the original name of the lipid molecular species*/
  private String originalMolecularSpecies_ = null;
  /** the name of the lipid molecular species in LDA encoding*/
  private String molecularSpecies_ = null;
  /** the m/z of the MS2 precursor*/
  private String ms2Precursor_ = null;
  /** the activation mode for MS2*/
  private String ms2Activation_ = null;
  /** the MS3 precursor*/
  private String ms3Precursor_ = null;
  /** the activation mode for MS3*/
  private String ms3Activation_ = null;
  /** the Alex123 ID*/
  private String id_ = null;
  /** the lipid category*/
  private String category_ = null;
  /** is this a conflicting ion*/
  private String conflicts_ = null;
  /** the number of carbon atoms*/
  private int carbonNumber_ = -1;
  /** the number of double bonds*/
  private int dbNumber_ = -1;
  /** the sum composition of this lipid species*/
  private String sumComposition_ = null;
  /** the sum formula as stored in the target list*/
  private String originalSumFormula_ = null;
  /** the number of carbon atoms of the fragment*/
  private int fragmentCarbonNumber_ = -1;
  /** the number of double bonds of the fragment*/
  private int fragmentDbNumber_ = -1;
  /** the number of OH groups of the fragment*/
  private int fragmentOhNumber_ = -1;
  /** the sum composition of the fragment*/
  private String fragmentSumComposition_ = null;
  /** the chemical formula of the fragment*/
  private String fragmentFormula_ = null;
  /** has this class alex123 fragment information*/
  private boolean alex123FragmentsForClass_ = false;
  /** the MSn fragments - only possibly for MS1 species - first key: msLevel; second key: fragment name; third key molecular species name*/
  private Hashtable<Integer,Hashtable<String,Hashtable<String,TargetlistEntry>>> msnFragments_;
  
  
  /**
   * constructor for the value object 
   * @param detector name of the detector
   * @param polarity the polarity of this ion
   * @param msLevel the msLevel of this observance
   * @param mz the m/z of this ion
   * @param fragment the name of the fragment
   * @param structure the verified lipid structure
   * @param species the name of this species
   * @param molecularSpecies the name of the lipid molecular species
   * @param lipidClass the name of the lipid class
   * @param ms2Precursor  the m/z of the MS2 precursor 
   * @param ms2Activation the activation mode for MS2
   * @param ms3Precursor the MS3 precursor
   * @param ms3Activation the activation mode for MS3
   * @param adduct the name of the modification
   * @param adductFormula the chemical formula of the modification
   * @param id the Alex123 ID
   * @param category the lipid category
   * @param conflicts is this a conflicting ion
   * @param charge the charge of the ion
   * @param carbonNumber the number of carbon atoms
   * @param dbNumber the number of double bonds
   * @param ohNumber the number of OH groups
   * @param sumComposition the sum composition of this lipid species
   * @param formula the chemical formula of this species
   * @param fragmentCarbonNumber the number of carbon atoms of the fragment
   * @param fragmentDbNumber the number of double bonds of the fragment
   * @param fragmentOhNumber the number of OH groups of the fragment
   * @param fragmentSumComposition the sum composition of the fragment
   * @param fragmentFormula the chemical formula of the fragment
   * @throws AlexTargetlistParserException thrown whenever there is something wrong with the entries
   * @throws HydroxylationEncodingException when there is no encoding for the provided ohNumber
   */
  public TargetlistEntry(String detector, String polarity, int msLevel,
      double mz, String fragment, String structure, String species,
      String molecularSpecies, String lipidClass, String ms2Precursor,
      String ms2Activation, String ms3Precursor, String ms3Activation,
      String adduct, String adductFormula, String id, String category, String conflicts,
      int charge, int carbonNumber, int dbNumber, int ohNumber,
      String sumComposition, String formula, String originalSumFormula, int fragmentCarbonNumber,
      int fragmentDbNumber, int fragmentOhNumber,
      String fragmentSumComposition, String fragmentFormula) throws AlexTargetlistParserException, HydroxylationEncodingException
  {
    super(lipidClass, species, -1, ohNumber, formula, mz, charge, adduct, adductFormula, -1f, 0f, 0f,
        new Vector<Double>(), new Vector<Double>(),0,"");
	this.originalClassName_ = lipidClass;
    this.detector_ = detector;
    this.polarity_ = polarity;
    this.msLevel_ = msLevel;
    this.fragment_ = fragment;
    this.structure_ = structure;
//    if (/*lipidClass.endsWith(" O-") || lipidClass.endsWith(" P-") || lipidClass.equalsIgnoreCase("Cer")*/originalClassName_.startsWith("IS ")) 
//      System.out.println(originalClassName_+" ; "+species+" ; "+molecularSpecies);
    this.species_ = species;
    this.originalMolecularSpecies_ = molecularSpecies;
    try {
      Object[] ldaEncoding = StaticUtils.translateAlexNomenclatureToLDA(originalClassName_,species_,originalMolecularSpecies_);
      this.analyteClass_ = (String)ldaEncoding[0];
      FattyAcidVO speciesVO = (FattyAcidVO)ldaEncoding[1];
      this.prefixOrName_ = speciesVO.getPrefix();
      this.carbons_ = speciesVO.getcAtoms();
      this.dbs_ = speciesVO.getDoubleBonds();
      this.molecularSpecies_ = (String)ldaEncoding[2];
//      if (/*lipidClass.endsWith(" O-") || lipidClass.endsWith(" P-") || lipidClass.equalsIgnoreCase("Cer")*/originalClassName_.startsWith("IS ")) 
//        System.out.println("After: "+(String)ldaEncoding[0]+" ; "+getIdString()+" ; "+ldaEncoding[2]);
    } catch (LipidCombinameEncodingException e) {
      throw new AlexTargetlistParserException(e);
    }


    this.ms2Precursor_ = ms2Precursor;
    this.ms2Activation_ = ms2Activation;
    this.ms3Precursor_ = ms3Precursor;
    this.ms3Activation_ = ms3Activation;
    this.id_ = id;
    this.category_ = category;
    this.conflicts_ = conflicts;
    this.carbonNumber_ = carbonNumber;
    this.dbNumber_ = dbNumber;
    this.sumComposition_ = sumComposition;
    this.originalSumFormula_ = originalSumFormula;
    this.fragmentCarbonNumber_ = fragmentCarbonNumber;
    this.fragmentDbNumber_ = fragmentDbNumber;
    this.fragmentOhNumber_ = fragmentOhNumber;
    this.fragmentSumComposition_ = fragmentSumComposition;
    this.fragmentFormula_ = fragmentFormula;
    this.msnFragments_ = new Hashtable<Integer,Hashtable<String,Hashtable<String,TargetlistEntry>>>();
    this.alex123FragmentsForClass_ = false;
  }

  /**
   * @return Alex123 class name
   */
  public String getOriginalClassName()
  {
    return originalClassName_;
  }
  
  /**
   * @return name of the detector
   */
  public String getDetector()
  {
    return detector_;
  }

  /**
   * @return the polarity of this ion
   */
  public String getPolarity()
  {
    return polarity_;
  }

  /**
   * @return the msLevel of this observance
   */
  public int getMsLevel()
  {
    return msLevel_;
  }
  
  
  /**
   * sets the msLevel
   * @param msLevel the MS level
   */
  public void setMsLevel_(int msLevel)
  {
    this.msLevel_ = msLevel;
  }

  /**
   * @return the name of the fragment
   */
  public String getFragment()
  {
    return fragment_;
  }

  /**
   * @return the verified lipid structure
   */
  public String getStructure()
  {
    return structure_;
  }

  /**
   * @return the name of this species
   */
  public String getSpecies()
  {
    return species_;
  }

  /**
   * @return the name of the lipid molecular species
   */
  public String getMolecularSpecies()
  {
    return molecularSpecies_;
  }
  
  
  /**
   * @return the name of the lipid molecular species
   */
  public String getOriginalMolecularSpecies()
  {
    return originalMolecularSpecies_;
  }
  
  
  /**
   * sets the name of the molecular species
   * @param molecularSpecies name of the molecular species
   * @throws AlexTargetlistParserException thrown whenever there is something wrong with the entries
   */
  public void setOriginalMolecularSpecies(String originalMolecularSpecies) throws AlexTargetlistParserException
  {
    this.originalMolecularSpecies_ = originalMolecularSpecies;
    Object[] ldaEncoding;
    try {
      ldaEncoding = StaticUtils.translateAlexNomenclatureToLDA(originalClassName_,species_,originalMolecularSpecies_);
      this.analyteClass_ = (String)ldaEncoding[0];
      FattyAcidVO speciesVO = (FattyAcidVO)ldaEncoding[1];
      this.prefixOrName_ = speciesVO.getPrefix();
      this.carbons_ = speciesVO.getcAtoms();
      this.dbs_ = speciesVO.getDoubleBonds();
      this.molecularSpecies_ = (String)ldaEncoding[2];
    }
    catch (LipidCombinameEncodingException e) {
      throw new AlexTargetlistParserException(e);
    }
  }
  
  
  /**
   * @return the m/z of the MS2 precursor 
   */
  public String getMs2Precursor()
  {
    return ms2Precursor_;
  }
  
  
  /**
   * sets the value of the MS2 precursor
   * @param ms2Precursor the m/z of the MS2 precursor
   */
  public void setMs2Precursor(String ms2Precursor)
  {
    this.ms2Precursor_ = ms2Precursor;
  }

  /**
   * @return the activation mode for MS2
   */
  public String getMs2Activation()
  {
    return ms2Activation_;
  }
    
  /**
   * sets the MS2 activation type
   * @param ms2Activation MS2 activation type
   */
  public void setMs2Activation(String ms2Activation)
  {
    this.ms2Activation_ = ms2Activation;
  }

  /**
   * @return the MS3 precursor
   */
  public String getMs3Precursor()
  {
    return ms3Precursor_;
  }
 
  /**
   * sets the value of the MS3 precursor
   * @param ms3Precursor the m/z of the MS3 precursor
   */
  public void setMs3Precursor(String ms3Precursor)
  {
    this.ms3Precursor_ = ms3Precursor;
  }


  /**
   * @return the activation mode for MS3
   */
  public String getMs3Activation()
  {
    return ms3Activation_;
  }
  
  /**
   * sets the MS3 activation type
   * @param ms3Activation MS3 activation type
   */
  public void setMs3Activation(String ms3Activation)
  {
    this.ms3Activation_ = ms3Activation;
  }

  /**
   * @return the Alex123 ID
   */
  public String getId()
  {
    return id_;
  }

  /**
   * @return the lipid category
   */
  public String getCategory()
  {
    return category_;
  }

  /**
   * @return is this a conflicting ion
   */
  public String getConflicts()
  {
    return conflicts_;
  }

  /**
   * @return the number of carbon atoms
   */
  public int getCarbonNumber()
  {
    return carbonNumber_;
  }

  /**
   * @return the number of double bonds
   */
  public int getDbNumber()
  {
    return dbNumber_;
  }

  /**
   * @return the sum composition of this lipid species
   */
  public String getSumComposition()
  {
    return sumComposition_;
  }
  
  
  /**
   * @return the sum formula as stored in the target list
   */
  public String getOriginalSumFormula()
  {
    return originalSumFormula_;
  }

  /**
   * @return the number of carbon atoms of the fragment
   */
  public int getFragmentCarbonNumber()
  {
    return fragmentCarbonNumber_;
  }

  /**
   * @return the number of double bonds of the fragment
   */
  public int getFragmentDbNumber()
  {
    return fragmentDbNumber_;
  }

  /**
   * @return the number of OH groups of the fragment
   */
  public int getFragmentOhNumber()
  {
    return fragmentOhNumber_;
  }

  /**
   * @return the sum composition of the fragment
   */
  public String getFragmentSumComposition()
  {
    return fragmentSumComposition_;
  }

  /**
   * @return the chemical formula of the fragment
   */
  public String getFragmentFormula()
  {
    return fragmentFormula_;
  }
  
  /**
   * adds putative MSn fragments to an MS1 target
   * @param fragment the MSn fragment (target)
   */
  public void addFragment(TargetlistEntry fragment){
    Hashtable<String,Hashtable<String,TargetlistEntry>> fragments = new Hashtable<String,Hashtable<String,TargetlistEntry>>();
    if (msnFragments_.containsKey(fragment.getMsLevel())) fragments = msnFragments_.get(fragment.getMsLevel());
    Hashtable<String,TargetlistEntry> sameFragmentDiffMolSpecies = new Hashtable<String,TargetlistEntry>();
    if (fragments.containsKey(fragment.getFragment())) sameFragmentDiffMolSpecies = fragments.get(fragment.getFragment());
    if (sameFragmentDiffMolSpecies.containsKey(fragment.getMolecularSpecies())){
      System.out.println("The fragment "+fragment.getFragment()+" with the adduct \""+fragment.getModName()+"\" is defined more than once for the species \""+species_+"\" and molecular species \""+fragment.getMolecularSpecies()+"\"!");
      return;
    }
    sameFragmentDiffMolSpecies.put(fragment.getMolecularSpecies(), fragment);
    fragments.put(fragment.getFragment(), sameFragmentDiffMolSpecies);
    msnFragments_.put(fragment.getMsLevel(), fragments);
  }
  
  /**
   * sets the retention time constraints for quantitation
   * @param retTime the putative retention time of the analyte
   * @param minusTime the negative retention time tolerance
   * @param plusTime the positive retention time tolerance
   */
  public void setTimeConstraints(float retTime, float minusTime, float plusTime){
    super.isobaricRetTime_ = retTime;
    super.isobaricMinusTime_ = minusTime;
    super.isobaricPlusTime_ = plusTime;
  }
  
  /**
   * sets the required isotopic distribution values for this target
   * @param mustMatchProbabs the probabilities of the isotopes that must match
   * @param probabs probabilities of all isotopes to quantify
   * @param negativeStartValue if the distribution goes in the negative direction - how negative is the lowest isotope
   */
  public void setDistributionValues(Vector<Double> mustMatchProbabs, Vector<Double> probabs, int negativeStartValue){
    mustMatchProbabs_ = mustMatchProbabs;
    probabs_ = probabs;
    negStartValue_ = negativeStartValue;
  }
  
  /**
   * 
   * @return the MSn fragments - only possibly for MS1 species - first key: msLevel; second key: fragment name; third key molecular species name
   */
  public Hashtable<Integer,Hashtable<String,Hashtable<String,TargetlistEntry>>> getMsnFragments(){
    return this.msnFragments_;
  }

  /**
   * 
   * @return has this class alex123 fragment information
   */
  public boolean hasAlex123FragmentsForClass()
  {
    return alex123FragmentsForClass_;
  }

  /**
   * 
   * @param alex123FragmentsForClass has this class alex123 fragment information
   */
  public void setAlex123FragmentsForClass(boolean alex123FragmentsForClass)
  {
    this.alex123FragmentsForClass_ = alex123FragmentsForClass;
  }
  
  
  //do nothing; this should be handled at a different location for the TargetlistEntry
  protected Object[] splitInCarbonNumberAndPrefix(String analyteClass, String analyteName) {
    Object[] prefixAndC = new Object[2];
    prefixAndC[0] = "";
    prefixAndC[1] = 0;
    return prefixAndC;
  }
  
  public Hashtable<String,String> getMolSpeciesLookup(){
    Hashtable<String,String> molSpecies = new Hashtable<String,String>();
    for (Hashtable<String,Hashtable<String,TargetlistEntry>> targetsOfLevel : this.msnFragments_.values()) {
      for (Hashtable<String,TargetlistEntry> targetsOfFragment : targetsOfLevel.values()) {
        for (String molName : targetsOfFragment.keySet()) {
          if (molSpecies .containsKey(molName))
            continue;
          molSpecies.put(molName, targetsOfFragment.get(molName).getOriginalMolecularSpecies());
        }
      }
    }
    return molSpecies;
  }
  
}
