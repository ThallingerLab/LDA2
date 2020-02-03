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
package at.tugraz.genome.vos;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LdaDialStandardsEvidence
{
  private String ms1Name_;
  private boolean faDetectable_;
  private boolean ms1DetectedLDA_;
  private boolean faDetectedLDAWoFPs_;
  private boolean faDetectedLDA_;
  private boolean ms1DetectedLB_;
  private boolean faDetectedLBWoFPs_;
  private boolean faDetectedLB_;

  public LdaDialStandardsEvidence(String ms1Name){
    this.ms1Name_ = ms1Name;
    this.faDetectable_ = false;
    ms1DetectedLDA_ = false;
    faDetectedLDAWoFPs_ = false;
    faDetectedLDA_ = false;
    ms1DetectedLB_ = false;
    faDetectedLBWoFPs_ = false;
    faDetectedLB_ = false;
  }

 
  public boolean isFaDetectable()
  {
    return faDetectable_;
  }

  public void setFaDetectable(boolean faDetectable)
  {
    this.faDetectable_ = faDetectable;
  }
  
  public boolean isMs1DetectedLDA()
  {
    return ms1DetectedLDA_;
  }

  public void setMs1DetectedLDA(boolean ms1DetectedLDA)
  {
    this.ms1DetectedLDA_ = ms1DetectedLDA;
  }

  public boolean isFaDetectedLDAWoFPs()
  {
    return faDetectedLDAWoFPs_;
  }

  public void setFaDetectedLDAWoFPs(boolean faDetectedLDAWoFPs)
  {
    this.faDetectedLDAWoFPs_ = faDetectedLDAWoFPs;
  }

  public boolean isFaDetectedLDA()
  {
    return faDetectedLDA_;
  }

  public void setFaDetectedLDA(boolean faDetectedLDA)
  {
    this.faDetectedLDA_ = faDetectedLDA;
  }

  public boolean isMs1DetectedLB()
  {
    return ms1DetectedLB_;
  }

  public void setMs1DetectedLB(boolean ms1DetectedLB)
  {
    this.ms1DetectedLB_ = ms1DetectedLB;
  }

  public boolean isFaDetectedLBWoFPs()
  {
    return faDetectedLBWoFPs_;
  }

  public void setFaDetectedLBWoFPs(boolean faDetectedLBWoFPs)
  {
    this.faDetectedLBWoFPs_ = faDetectedLBWoFPs;
  }

  public boolean isFaDetectedLB()
  {
    return faDetectedLB_;
  }

  public void setFaDetectedLB(boolean faDetectedLB)
  {
    this.faDetectedLB_ = faDetectedLB;
  }

  public String getMs1Name()
  {
    return ms1Name_;
  }

  
}
