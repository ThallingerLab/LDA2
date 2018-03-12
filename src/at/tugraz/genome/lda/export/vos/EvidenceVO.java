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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;

/**
 * Class for exporting lipid evidence information. 
 * 
 * @author Juergen Hartler
 *
 */
public class EvidenceVO extends EvidenceBase
{
  
  /** a unique identifier that must be different between the EvidenceVOs*/
  private int id_;
  /** the assumed chemical formula of this identification (including the modification/adduct)*/
  private String chemFormula_;
  /** the actually detected m/z value*/
  private double expMz_;
  /** the theoretical m/z value of the assumed identification*/
  private double theorMz_;
  /** the (recommended) LDA original identification*/
  private String ldaId_;
  /** the evidence of a sample containing the highest structural information; key: experiment name; value; (molecular) species name*/
  private Hashtable<String,String> highestLDAStructuralEv_;
  /** the MS-level this evidence is based on*/
  private int msLevel_;
  /** the found scan numbers; key: experiment name; value: scan numbers*/
  private Hashtable<String,List<Integer>> scanNrs_;

  /**
   * constructor for EvidenceVO 
   * @param base an EvidenceBase object holding basic information that is common to the individual evidence value objects
   * @param id a unique identifier that must be different between the EvidenceVOs
   * @param ldaId the (recommended) LDA original identification
   * @param highestLDAStructuralEv the evidence of a sample containing the highest structural information; key: experiment name; value; (molecular) species name
   * @param chemFormula the assumed chemical formula of this identification (including the modification/adduct)
   * @param expMz the actually detected m/z value
   * @param theorMz the theoretical m/z value of the assumed identification
   * @param msLevel the MS-level this evidence is based on
   * @param scanNrs found scan numbers; key: experiment name; value: scan numbers
   */
  public EvidenceVO(EvidenceBase base, int id, String ldaId, Hashtable<String,String> highestLDAStructuralEv,
      String chemFormula, double expMz, double theorMz, int msLevel, Hashtable<String,Set<Integer>> scanNrs){
    super(base.getEvidenceGroupingId(),base.getModification(),base.getCharge());
    this.id_ = id;
    this.ldaId_ = ldaId;
    this.highestLDAStructuralEv_ = highestLDAStructuralEv;
    this.chemFormula_ = chemFormula;
    this.expMz_ = expMz;
    this.theorMz_ = theorMz;
    this.msLevel_ = msLevel;
    scanNrs_ = new Hashtable<String,List<Integer>>();
    for (String exp : scanNrs.keySet()){
      List<Integer> list = new ArrayList<Integer>(scanNrs.get(exp));
      Collections.sort(list);
      scanNrs_.put(exp, list);
    }
  }

  
  /**
   * 
   * @return a unique identifier for this EvidenceVO
   */
  public int getId()
  {
    return id_;
  }
  
  
  /**
   * 
   * @return the (recommended) LDA original identification
   */
  public String getLdaId()
  {
    return ldaId_;
  }

  
  /**
   * 
   * @return the assumed chemical formula of this identification (including the modification/adduct)
   */
  public String getChemFormula()
  {
    return chemFormula_;
  }

  
  /**
   * 
   * @return the actually detected m/z value
   */
  public double getExpMz()
  {
    return expMz_;
  }

  
  /**
   * 
   * @return the theoretical m/z value of the assumed identification
   */
  public double getTheorMz()
  {
    return theorMz_;
  }
  
  
  /**
   * 
   * @return the MS-level this evidence is based on
   */
  public int getMsLevel()
  {
    return msLevel_;
  }
  
  
  /**
   * returns the found MSn scan numbers of an experiment
   * @param exp the experiment name
   * @return the found MSn scan numbers of an experiment
   */
  public List<Integer> getScanNrs(String exp){
    if (scanNrs_.containsKey(exp))
      return scanNrs_.get(exp);
    else
      return null;
  }
  
  /**
   * the highest structural information detected for this experiment
   * @param exp experiment name
   * @return the highest structural information detected for this experiment
   */
  public String getBestIdentification(String exp){
    if (highestLDAStructuralEv_.containsKey(exp))
      return highestLDAStructuralEv_.get(exp);
    else
      return null;
  }
  
}
