/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.tugraz.genome.lda.export.vos;

import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Vector;

import at.tugraz.genome.lda.vos.ExportOptionsVO;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * class containing cumulative information across the various aspects
 * @author Juergen Hartler
 *
 */

public class SummaryVO
{
  /** a unique identifier that must be different between the detected FeatureVOs*/
  private Integer id_;
  /** a unique identifier for this species*/
  private String speciesId_;
  /** ids referring to FeatureVOs that belong to this summary*/
  private Vector<Integer> featureRefs_;
  /** structural information of this SummaryVO; if none present, this value is null*/
  private String molecularId_;
  /** the chemical formula of the neutral molecule*/
  private String chemFormula_;
  /** the theoretical mass of this object*/
  private Double neutralMass_;
  /** sorted (by abundance) vector of modifications/adducts*/
  private Vector<String> mods_;
  /** a reliability score specific to mzTab*/
  private int mzTabReliability_;
  /** the totally detected values/areas for the experiments; key: experiment name; value: area*/
  private Hashtable<String,Double> areas_;
  /** the mean area for each selected group (heat map); key: group name; value: mean area*/
  private LinkedHashMap<String,Double> groupMeans_;
  /** the coefficient of variation of each selected group (heat map); key: group name; value: coefficient of variation*/
  private LinkedHashMap<String,Double> groupCoeffVar_;
  /** the retention times of the highest peaks of every modification; first key modification; second key: experiment*/
  private Hashtable<String,Hashtable<String,Double>> rtsOfMods_;
  /** the reliability of the evidence of every modification; first key modification; second key: experiment*/
  private Hashtable<String,Hashtable<String,Short>> evidenceReliabilityOfMods_;
  /** the mean retention time for each selected group (heat map); first key: modification; key: group name; value: mean area*/
  private Hashtable<String,Hashtable<String,Double>> groupRts_;
  /** the standard deviations for the retention times of each selected group (heat map); first key: modification; key: group name; value: mean area*/  
  private Hashtable<String,Hashtable<String,Double>> groupRtStdevs_;
  
  /** there is no ms2 evidence found for the hit*/
  public final static short EVIDENCE_MS1_ONLY = 0;
  /** ms2 found and no overlap*/
  public final static short EVIDENCE_MS2_UNAMBIGUOUS = 1;
  /** ms2 is found, but there is an overlap with an isobar from another species, but the two hits can be separated*/ 
  public final static short EVIDENCE_MS2_SPLIT = 2;
  /** ms2 is found, but there is an overlap with an isobar from another species, and the two hits cannot be separated*/
  public final static short EVIDENCE_MS2_NO_SPLIT_POSSIBLE = 3;
  
  /**
   * constructor for creating cumulative class containing information across the various aspects
   * @param id a unique identifier that must be different between the detected SummaryVOs
   * @param speciesId a unique identifier for this species
   * @param molecularId structural information of this SummaryVO; if none present, this value is null
   * @param featureRefs ids referring to FeatureVOs that belong to this summary
   * @param chemFormula the chemical formula of the neutral molecule
   * @param neutralMass the theoretical mass of this object
   * @param rt the retention time of the strongest detection peak across all experiments
   * @param mods sorted (by abundance) vector of modifications/adducts
   * @param mzTabReliability a reliability score specific to mzTab
   * @param areas the mean area for each selected group (heat map); key: group name; value: mean area
   * @param rtsOfMods retention times of the highest peaks of every modification; first key modification; second key: experiment
   * @param evidenceReliabilityOfMods what kind of evidence supports this area value; first key modification; second key: experiment
   * @param expsOfGroup key: group name; value: experiments belonging to this group
   */
  public SummaryVO(Integer id, String speciesId, String molecularId, Vector<Integer> featureRefs, String chemFormula,
      Double neutralMass, Float rt, Vector<String> mods, int mzTabReliability, Hashtable<String,Double> areas,
      Hashtable<String,Hashtable<String,Double>> rtsOfMods, Hashtable<String,Hashtable<String,Short>> evidenceReliabilityOfMods,
      LinkedHashMap<String,Vector<String>> expsOfGroup)
  {
    this.id_ = id;
    this.speciesId_ = speciesId;
    this.featureRefs_ = featureRefs;
    this.molecularId_ = molecularId;
    this.chemFormula_ = chemFormula;
    this.neutralMass_ = neutralMass;
    this.mods_ = mods;
    this.mzTabReliability_ = mzTabReliability;
    this.areas_ = areas;
    this.evidenceReliabilityOfMods_ = evidenceReliabilityOfMods;
    this.rtsOfMods_ = rtsOfMods;
    Vector<Double> areasOfGroup;
    Vector<Double> rtsOfGroup;
    if (expsOfGroup.size()>0){
      groupMeans_ = new LinkedHashMap<String,Double>();
      groupCoeffVar_ = new LinkedHashMap<String,Double>();
      groupRts_ = new Hashtable<String,Hashtable<String,Double>>();
      groupRtStdevs_ = new Hashtable<String,Hashtable<String,Double>>();
      for (String groupName : expsOfGroup.keySet()){
        areasOfGroup = new Vector<Double>();
        for (String exp : expsOfGroup.get(groupName)){
          if (areas.containsKey(exp) && areas.get(exp)>0)
            areasOfGroup.add(areas.get(exp));
        }
        if (areasOfGroup.size()==0){
          groupMeans_.put(groupName, 0d);
          groupCoeffVar_.put(groupName, Double.NaN);
          continue;
        }
        double[] doubleArray = new double[areasOfGroup.size()];
        for (int i=0; i!=areasOfGroup.size(); i++) doubleArray[i] = areasOfGroup.get(i);  
        double mean = Calculator.mean(doubleArray);
        groupMeans_.put(groupName, mean);
        groupCoeffVar_.put(groupName, Calculator.stddeviation(areasOfGroup)/mean);
        
        //for calculating the Rt values;
        for (String mod : rtsOfMods.keySet()){
          Hashtable<String,Double> rts = rtsOfMods.get(mod);
          rtsOfGroup = new Vector<Double>();
          for (String exp : expsOfGroup.get(groupName)){
            if (rts.containsKey(exp))
              rtsOfGroup.add(rts.get(exp));
          }
          if (!groupRts_.containsKey(mod)){
            groupRts_.put(mod, new Hashtable<String,Double>());
            groupRtStdevs_.put(mod, new Hashtable<String,Double>());
          }
          doubleArray = new double[rtsOfGroup.size()];
          for (int i=0; i!=rtsOfGroup.size(); i++) doubleArray[i] = rtsOfGroup.get(i);  
          groupRts_.get(mod).put(groupName, Calculator.mean(doubleArray));
          groupRtStdevs_.get(mod).put(groupName, Calculator.stddeviation(rtsOfGroup));
        }
      }
    }
  }
  
  /**
   * constructor for creating cumulative class containing information across the various aspects
   * @param speciesId a unique identifier for this species
   * @param molecularId structural information of this SummaryVO; if none present, this value is null
   * @param featureRefs ids referring to FeatureVOs that belong to this summary
   * @param chemFormula the chemical formula of the neutral molecule
   * @param neutralMass the theoretical mass of this object
   * @param rt the retention time of the strongest detection peak across all experiments
   * @param mods sorted (by abundance) vector of modifications/adducts
   * @param mzTabReliability a reliability score specific to mzTab
   * @param areas the mean area for each selected group (heat map); key: group name; value: mean area
   * @param rtsOfMods retention times of the highest peaks of every modification; first key modification; second key: experiment
   * @param evidenceReliabilityOfMods what kind of evidence supports this area value; first key modification; second key: experiment
   * @param expsOfGroup key: group name; value: experiments belonging to this group
   */
  public SummaryVO(String speciesId, String molecularId, Vector<Integer> featureRefs, String chemFormula,
      Double neutralMass, Float rt, Vector<String> mods, int mzTabReliability, Hashtable<String,Double> areas,
      Hashtable<String,Hashtable<String,Double>> rtsOfMods, Hashtable<String,Hashtable<String,Short>> evidenceReliabilityOfMods,
      LinkedHashMap<String,Vector<String>> expsOfGroup)
  {
    this(null,speciesId,molecularId,featureRefs,chemFormula,neutralMass,rt,mods,mzTabReliability,areas,rtsOfMods,evidenceReliabilityOfMods,
        expsOfGroup);
  }

  
  /**
   * 
   * @return the unique identifier for this summary
   */
  public Integer getId()
  {
    return id_;
  }

  
  /**
   * sets the unique identifier (in cases where sorting is required)
   * @param id the he unique identifier for this summary
   */
  public void setId(Integer id)
  {
    this.id_ = id;
  }


  /**
   * 
   * @return a unique identifier for this species
   */
  public String getSpeciesId()
  {
    return speciesId_;
  }

  
  /**
   * 
   * @return ids referring to FeatureVOs that belong to this summary
   */
  public Vector<Integer> getFeatureRefs()
  {
    return featureRefs_;
  }

  
  /**
   * 
   * @return structural information of this SummaryVO; if none present, this value is null
   */
  public String getMolecularId()
  {
    return molecularId_;
  }
  
  
  /**
   * 
   * @return the chemical formula of the neutral molecule
   */
  public String getChemFormula()
  {
    return chemFormula_;
  }

  
  /**
   * 
   * @return the theoretical mass of this object
   */
  public Double getNeutralMass()
  {
    return neutralMass_;
  }

  
  /**
   * 
   * @return sorted (by abundance) vector of modifications/adducts
   */
  public Vector<String> getModifications()
  {
    return mods_;
  }
  
  
  /**
   * 
   * @return a reliability score specific to mzTab
   */
  public int getMzTabReliability()
  {
    return mzTabReliability_;
  }

  
  /**
   * returns the intensity value of a certain experiment; if there is no detection, 0 is returned
   * @param exp the experiment name
   * @return the intensity value of a certain experiment; if there is no detection, 0 is returned
   */
  public Double getArea(String exp){
    if (areas_.containsKey(exp))
      return areas_.get(exp);
    else
      return 0d;
  }
  
  /**
   * returns the reliability of this identified area - for details see "EVIDENCE_" specifications in this object
   * @param exp
   * @return
   */
  public short getEvidenceReliabilty(String exp){
    short reliability = EVIDENCE_MS1_ONLY;
    for (Hashtable<String,Short> relies : evidenceReliabilityOfMods_.values()){
      if (!relies.containsKey(exp))
        continue;
      // returns the corresponding reliability, where MS2 evidence, and more uncertain MS2 evidence is preceding
      if (relies.get(exp)>reliability)
        reliability = relies.get(exp);
    }
    return reliability;
  }
  
  /**
   * returns the mean intensity value of a certain experiment group; if there is no detection, 0 is returned
   * @param groupName the name of the group
   * @return the mean intensity value of a certain experiment group; if there is no detection, 0 is returned
   */
  public Double getMeanArea(String groupName){
    Double area = 0d;
    if (groupMeans_!=null && groupMeans_.containsKey(groupName))
      area = groupMeans_.get(groupName);
    return area;
  }
  
  
  /**
   * returns the coefficient of variation of a certain experiment group; if there is no detection, 0 is returned
   * @param groupName the name of the group
   * @return the coefficient of variation of a certain experiment group; if there is no detection, 0 is returned
   */
  public Double getCoeffVar(String groupName){
    Double stdev = Double.NaN;
    if (groupCoeffVar_.containsKey(groupName))
      stdev = groupCoeffVar_.get(groupName);
    return stdev;
  }
  
  /**
   * calculates a deviation; the expOptions specifies which deviation value has to be calculated 
   * @param expOptions value object specifying the type of deviation value
   * @param expsOfGroup the names of the experiments belonging to this sample group
   * @return the calculated deviation value
   */
  public double calculateDeviationValue(ExportOptionsVO expOptions, Vector<String> expsOfGroup){
    double sdValue = Double.NaN;
    if (expOptions!=null&&expOptions.getExportType()!=ExportOptionsVO.EXPORT_NO_DEVIATION){
      Vector<Double> areasOfGroup = new Vector<Double>();
      double area;
      for (String exp : expsOfGroup){
        area = getArea(exp);
        if (area>0)
          areasOfGroup.add(area);
      }
      sdValue = Calculator.stddeviation(areasOfGroup);
      if (expOptions.getExportType() == ExportOptionsVO.EXPORT_SD_DEVIATION || expOptions.getExportType() == ExportOptionsVO.EXPORT_SD_DEV_AND_ERROR)
        sdValue = sdValue*Double.parseDouble(expOptions.getSdValue());
      if (expOptions.getExportType() == ExportOptionsVO.EXPORT_SD_ERROR || expOptions.getExportType() == ExportOptionsVO.EXPORT_SD_DEV_AND_ERROR)
        sdValue = sdValue/Math.sqrt(areasOfGroup.size());
    }
    return sdValue;
  }
  
  /**
   * returns the retention time of a certain modification of an experiment
   * @param mod the adduct name
   * @param exp the experiment name
   * @return the retention time of a certain modification of an experiment
   */
  public Double getRetentionTime(String mod, String exp){
    if (rtsOfMods_.containsKey(mod) && rtsOfMods_.get(mod).containsKey(exp))
      return rtsOfMods_.get(mod).get(exp);
    else
      return null;
  }
  
  /**
   * returns the mean retention time of a certain modification of a sample group
   * @param mod the adduct name
   * @param group the sample group
   * @return mean retention time of a certain modification of a sample group
   */
  public Double getMeanRetentionTime(String mod, String group){
    if (groupRts_.containsKey(mod) && groupRts_.get(mod).containsKey(group))
      return groupRts_.get(mod).get(group);
    else
      return null;
  }

  /**
   * returns the standard deviation of the retention time of a certain modification of a sample group
   * @param mod the adduct name
   * @param group the sample group
   * @return standard deviation of a certain modification of a sample group
   */
  public Double getStdevRetentionTime(String mod, String group){
    if (groupRtStdevs_.containsKey(mod) && groupRtStdevs_.get(mod).containsKey(group))
      return groupRtStdevs_.get(mod).get(group);
    else
      return null;
  }
  
  
}
