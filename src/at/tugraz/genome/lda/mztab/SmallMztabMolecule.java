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

package at.tugraz.genome.lda.mztab;

import java.util.Vector;

import de.isas.mztab1_1.model.SmallMoleculeEvidence;
import de.isas.mztab1_1.model.SmallMoleculeFeature;
import de.isas.mztab1_1.model.SmallMoleculeSummary;

/**
 * Class containing exportable information in mzTab-specific format
 * 
 * @author Juergen Hartler
 *
 */
public class SmallMztabMolecule
{
  /** the next unique summary id to use*/
  private int currentSummaryId_;
  /** exportable mzTab-specific SML objects*/
  private Vector<SmallMoleculeSummary> summary_;
  /** the next unique feature id to use*/
  private int currentFeatureId_;
  /** exportable mzTab-specific SMF objects*/
  private Vector<SmallMoleculeFeature> features_;
  /** the next unique identifier for an EvidenceVO*/
  private int currentEvidenceId_;
  /** the next unique identifier for an evidence group, i.e. an identifier for evidence originating from the same spectra*/
  private int currentEvGroupingId_;
  /** the detected evidence information*/
  private Vector<SmallMoleculeEvidence> evidence_;
  
  
  /**
   * Constructor for generating object containing exportable information in mzTab-specific format
   * @param currentSummaryId the next unique summary id to use
   * @param summary exportable mzTab-specific SML objects
   * @param currentFeatureId the next unique feature id to use
   * @param features exportable mzTab-specific SMF objects
   * @param currentEvidenceId the next unique identifier for an EvidenceVO
   * @param currentEvGroupingId the next unique identifier for an evidence group, i.e. an identifier for evidence originating from the same spectra
   * @param evidence exportable mzTab-specific SME objects
   */
  public SmallMztabMolecule(int currentSummaryId, Vector<SmallMoleculeSummary> summary, int currentFeatureId,
      Vector<SmallMoleculeFeature> features, int currentEvidenceId, int currentEvGroupingId,
      Vector<SmallMoleculeEvidence> evidence){
      this.currentSummaryId_ = currentSummaryId;
     summary_ = summary;
     this.currentFeatureId_ = currentFeatureId;
     this.features_ = features;
     this.currentEvidenceId_ = currentEvidenceId;
     this.currentEvGroupingId_ = currentEvGroupingId;
     this.evidence_ = evidence;
  }

  
  /**
   * 
   * @return the next unique summary id to use
   */
  public int getCurrentSummaryId()
  {
    return currentSummaryId_;
  }


  /**
   * 
   * @return exportable mzTab-specific SML objects
   */
  public Vector<SmallMoleculeSummary> getSummary()
  {
    return summary_;
  }

  
  /**
   * 
   * @return the next unique feature id to use
   */
  public int getCurrentFeatureId()
  {
    return currentFeatureId_;
  }

  
  /**
   * 
   * @return exportable mzTab-specific SMF objects
   */
  public Vector<SmallMoleculeFeature> getFeatures()
  {
    return features_;
  }

  
  /**
   * 
   * @return the next unique identifier for an EvidenceVO
   */
  public int getCurrentEvidenceId()
  {
    return currentEvidenceId_;
  }

  
  /**
   * 
   * @return the next unique identifier for an evidence group, i.e. an identifier for evidence originating from the same spectra
   */
  public int getCurrentEvGroupingId()
  {
    return currentEvGroupingId_;
  }

  
  /**
   * 
   * @return exportable mzTab-specific SME objects
   */
  public Vector<SmallMoleculeEvidence> getEvidence()
  {
    return evidence_;
  }  
  
}
