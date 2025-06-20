/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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
  /** the MS-run identifier*/
  private String expName_;
  /** a unique identifier that must be different between the EvidenceVOs*/
  private int id_;
  /** an identifier for evidence originating from the same spectra*/
  private int evidenceGroupingId_;
  /** the assumed chemical formula of this identification (including the modification/adduct)*/
  private String chemFormula_;
  /** the actually detected m/z value*/
  private double expMz_;
  /** the theoretical m/z value of the assumed identification*/
  private double theorMz_;
  /** the identifier for the lipid species*/
  private String speciesId_;
  /** the (recommended) LDA original identification*/
  private String ldaStructure_;
  /** the MS-level this evidence is based on*/
  private int msLevel_;
  /** the found scan numbers*/
  private List<Integer> scanNrs_;

  /**
   * constructor for EvidenceVO 
   * @param expName the experiment identifier
   * @param base an EvidenceBase object holding basic information that is common to the individual evidence value objects
   * @param id a unique identifier that must be different between the EvidenceVOs
   * @param evidenceGroupingId an  identifier for evidence originating from the same spectra
   * @param speciesId the identifier for the lipid species
   * @param ldaStructure the (recommended) LDA original identification
   * @param chemFormula the assumed chemical formula of this identification (including the modification/adduct)
   * @param expMz the actually detected m/z value
   * @param theorMz the theoretical m/z value of the assumed identification
   * @param msLevel the MS-level this evidence is based on
   * @param scanNrs found scan numbers
   */
  public EvidenceVO(String expName, EvidenceBase base, int id, int evidenceGroupingId, String speciesId, String ldaStructure, String chemFormula, double expMz,
      double theorMz, int msLevel, Set<Integer> scanNrs){
    super(base.getModification(),base.getCharge());
    this.expName_ = expName;
    this.id_ = id;
    this.evidenceGroupingId_ = evidenceGroupingId;
    this.speciesId_ = speciesId;
    this.ldaStructure_ = ldaStructure;
    this.chemFormula_ = chemFormula;
    this.expMz_ = expMz;
    this.theorMz_ = theorMz;
    this.msLevel_ = msLevel;
    this.scanNrs_ = null;
    if (scanNrs!=null){
      scanNrs_ = new ArrayList<Integer>(scanNrs);
      Collections.sort(scanNrs_);
    }
  }

  
  /**
   * 
   * @return the experiment identifier
   */
  public String getExpName()
  {
    return expName_;
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
  * @return  an  identifier for evidence originating from the same spectra
  */
  public int getEvidenceGroupingId(){
    return evidenceGroupingId_;
  }

  
  /**
   * 
   * @return the identifier for the lipid species
   */
  public String getSpeciesId()
  {
    return speciesId_;
  }
  
  /**
   * 
   * @return the (recommended) LDA original identification
   */
  public String getLdaStructure()
  {
    return ldaStructure_;
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
  public List<Integer> getScanNrs(){
    if (scanNrs_!=null)
      return scanNrs_;
    else
      return null;
  }
  
}
