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

import java.util.Vector;

/**
 * Class containing all the different information regarding the composite export of several experiments
 * 
 * @author Juergen Hartler
 *
 */

public class SpeciesExportVO
{
  /** the next unique feature id to use*/ 
  private int currentSummaryId_;
  /** sorted vector containing cumulative information across the various aspects*/
  private Vector<SummaryVO> summaries_;
  /** the next unique feature id to use*/ 
  private int currentFeatureId_;
  /** the detected feature information*/
  private Vector<FeatureVO> features_;
  /** the next unique identifier for an EvidenceVO*/
  private int currentEvidenceId_;
  /** the next unique identifier for an evidence group, i.e. an identifier for evidence originating from the same spectra*/
  private int currentEvGroupingId_;
  /** the detected evidence information*/
  private Vector<EvidenceVO> evidence_;
  
  /**
   * constructor for creating the composite SpeciesExportVO
   * @param currentSummaryId the next unique summary id to use
   * @param summaries sorted vector containing cumulative information across the various aspects
   * @param currentFeatureId the next unique feature id to use
   * @param features the detected feature information
   * @param currentEvidenceId the next unique identifier for an EvidenceVO
   * @param currentEvGroupingId the next unique identifier for an evidence group, i.e. an identifier for evidence originating from the same spectra
   * @param evidence the detected evidence information
   */
  public SpeciesExportVO(int currentSummaryId, Vector<SummaryVO> summaries, int currentFeatureId, Vector<FeatureVO> features,
      int currentEvidenceId, int currentEvGroupingId, Vector<EvidenceVO> evidence)
  {
    super();
    this.currentSummaryId_ = currentSummaryId;
    this.summaries_ = summaries;
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
   * @return sorted vector containing cumulative information across the various aspects
   */
  public Vector<SummaryVO> getSummaries()
  {
    return summaries_;
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
   * @return the detected feature information
   */
  public Vector<FeatureVO> getFeatures()
  {
    return features_;
  }
  
  
  /**
   * 
   * @return the next unique feature id to use
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
   * @return the detected evidence information
   */
  public Vector<EvidenceVO> getEvidence()
  {
    return evidence_;
  }
  
}
