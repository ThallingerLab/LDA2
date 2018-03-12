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

import java.util.Vector;


/**
 * Class containing export information regarding one detected MS1-feature 
 * 
 * @author Juergen Hartler
 *
 */
public class FeatureVO
{
  /** a unique identifier that must be different between the detected FeatureVOs*/
  private int id_;
  /** the presumed adduct/modification of this hit*/
  private String adduct_;
  /** the actually detected m/z value*/
  private double expMz_;
  /** the presumed charge of this feature*/
  private int charge_;
  /** the retention time of the peak apex*/
  private double rtApex_;
  /** the retention time of the peak start*/
  private double rtStart_;
  /** the retention time of the peak end*/
  private double rtEnd_;
  /** the areas detected for this feature (sorted in the sequence of the experiments - 0 is entered for no detection)*/
  private Vector<Double> areas_;
  /** ids refering to EvidenceVOs that belong to this feature*/
  private Vector<Integer> evidenceRefs_;
  
  /**
   * construct for creating a feature object which contains information regarding one detected MS1-feature
   * @param id a unique identifier that must be different between the detected FeatureVOs
   * @param adduct the presumed adduct/modification of this hit
   * @param expMz the actually detected m/z value
   * @param charge the presumed charge of this feature
   * @param rtApex the retention time of the peak apex
   * @param rtStart the retention time of the peak start
   * @param rtEnd the retention time of the peak end
   * @param areas the areas detected for this feature (sorted in the sequence of the experiments - 0 is entered for no detection)
   */
  public FeatureVO(int id, String adduct, double expMz, int charge, double rtApex, double rtStart,
      double rtEnd, Vector<Double> areas)
  {
    this.id_ = id;
    this.adduct_ = adduct;
    this.expMz_ = expMz;
    this.charge_ = charge;
    this.rtApex_ = rtApex;
    this.rtStart_ = rtStart;
    this.rtEnd_ = rtEnd;
    this.areas_ = areas;
    this.evidenceRefs_ = new Vector<Integer>();
  }
  
  
  /**
   * 
   * @return the unique identifier for this feature
   */
  public int getId()
  {
    return id_;
  }
  
  
  /**
   * 
   * @return the presumed adduct/modification of this hit
   */
  public String getAdduct()
  {
    return adduct_;
  }


  /**
   * 
   * @return actually detected m/z value
   */
  public double getExpMz()
  {
    return expMz_;
  }


  /**
   * 
   * @return the presumed charge of this feature
   */
  public int getCharge()
  {
    return charge_;
  }


  /**
   * 
   * @return the retention time of the peak apex
   */
  public double getRtApex()
  {
    return rtApex_;
  }


  /**
   * 
   * @return the retention time of the peak start
   */
  public double getRtStart()
  {
    return rtStart_;
  }


  /**
   * 
   * @return the retention time of the peak end
   */
  public double getRtEnd()
  {
    return rtEnd_;
  }

  
  /**
   * 
   * @return the areas detected for this feature (sorted in the sequence of the experiments - 0 is entered for no detection)
   */
  public Vector<Double> getAreas()
  {
    return areas_;
  }
  
  
  /**
   * adds a reference to an EvidenceVO
   * @param evidenceRef the unique identifier for the EvidenceVO
   */
  public void addEvidenceRef(int evidenceRef){
    this.evidenceRefs_.add(evidenceRef);
  }


  /**
   * 
   * @return ids refering to EvidenceVOs belonging to this feature 
   */
  public Vector<Integer> getEvidenceRefs()
  {
    return evidenceRefs_;
  }
  
  
}
