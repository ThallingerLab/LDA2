/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2021 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

import at.tugraz.genome.lda.vos.IsotopicLabelVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class AnalyteOmegaInfoVO
{
  
  /** the name of the unlabeled analyte*/
  private String analyteName_;
  /** the name of the label detection including retention time*/
  private String labelDetection_;
  /** the information about the identified labels*/
  private Vector<IsotopicLabelVO> labels_;
  /** has this molecule only a single double bond*/
  private boolean singelDoubleBond_;
  
  
  /**
   * standard constructor for this value object
   * @param analyteName name of the unlabeled analyte
   * @param labelDetection name of the label detection including retention time
   * @param labels information about the identified labels
   * @param singleDoubleBond does this analyte have only one double bond
   */
  public AnalyteOmegaInfoVO(String analyteName, String labelDetection, Vector<IsotopicLabelVO> labels, boolean singleDoubleBond) {
    this.analyteName_ = analyteName;
    this.labelDetection_ = labelDetection;
    this.labels_ = labels;
    this.singelDoubleBond_ = singleDoubleBond;
  }

  
  /**
   * 
   * @return unlabeled analyte
   */
  public String getAnalyteName() {
    return analyteName_;
  }

  
  /**
   * 
   * @return label detection including retention time
   */
  public String getLabelDetection() {
    return labelDetection_;
  }

  
  /**
   * 
   * @return information about the identified labels
   */
  public Vector<IsotopicLabelVO> getLabels() {
    return labels_;
  }
  
  
  /**
   * 
   * @return the labels that were added to this molecule
   */
  public String getAppliedLabelId(){
    String labelPrefix = "";
    for (IsotopicLabelVO label : labels_)
      labelPrefix += label.getLabelId();
    return labelPrefix;
  }


  /**
   * 
   * @return has this analyte only one double bond
   */
  public boolean isSingelDoubleBond()
  {
    return singelDoubleBond_;
  }


  /**
   * calculates the expected RT of the unlabeled species
   * @param rtOfLabeledSpecies the rt of the labeled detection
   * @return the expected RT of the unlabeled species
   */
  public double calculateRtIncludingTheShift(double rtOfLabeledSpecies) {
    for (IsotopicLabelVO labelVO : this.labels_) {
      //the following line is for using absolute retention time shifts; however, results seem to be better using relative retention time
      //rtOfLabeledSpecies -= (double)labelVO.getRtShift();
      rtOfLabeledSpecies = rtOfLabeledSpecies/(1+(((double)labelVO.getRrtShift())/100d));
    }
    return rtOfLabeledSpecies;
  }
    
}
