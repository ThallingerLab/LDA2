/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp 
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
public class ExpVolumeSettingVO
{
  
  private String expName;
  private Double probeVolume;
  private Double sampleWeight;
  private Double endVolume;
  private Double proteinConc;
  private Double neutralLipidConc;
  
  public ExpVolumeSettingVO(String expName, Double probeVolume, Double sampleWeight, 
      Double endVolume, Double proteinConc, Double neutralLipidConc)
  {
    this.expName = expName;
    this.probeVolume = probeVolume;
    this.sampleWeight = sampleWeight;
    this.endVolume = endVolume;
    this.proteinConc = proteinConc;
    this.neutralLipidConc = neutralLipidConc;
  } 
  
  public String getExpName()
  {
    return expName;
  }

  public Double getProbeVolume()
  {
    return probeVolume;
  }
  
  public Double getSampleWeight()
  {
    return sampleWeight;
  }

  public void setSampleWeight(Double sampleWeight)
  {
    this.sampleWeight = sampleWeight;
  }
  
  public Double getEndVolume()
  {
    return endVolume;
  }

  public Double getProteinConc()
  {
    return proteinConc;
  }

  public Double getNeutralLipidConc()
  {
    return neutralLipidConc;
  }
  
  
}
