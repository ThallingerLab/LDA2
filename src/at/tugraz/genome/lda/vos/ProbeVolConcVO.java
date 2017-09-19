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

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ProbeVolConcVO
{
  
  private Double probeVolume_;
  private Double endVolume_;
  private Double sampleWeight_;
  private Double proteinConc_;
  private Double neutralLipidConc_;
  
  public ProbeVolConcVO(Double probeVolume_, Double endVolume_,
      Double sampleWeight, Double proteinConc_, Double neutralLipidConc_)
  {
    super();
    this.probeVolume_ = probeVolume_;
    this.endVolume_ = endVolume_;
    this.sampleWeight_ = sampleWeight;
    this.proteinConc_ = proteinConc_;
    this.neutralLipidConc_ = neutralLipidConc_;
  }
  
  public Double getProbeVolume()
  {
    return probeVolume_;
  }
  public Double getEndVolume()
  {
    return endVolume_;
  }
  public Double getSampleWeight()
  {
    return sampleWeight_;
  }
  public Double getProteinConc()
  {
    return proteinConc_;
  }
  public Double getNeutralLipidConc()
  {
    return neutralLipidConc_;
  }
  
  
  
}
