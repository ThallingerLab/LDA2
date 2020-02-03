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

import java.util.Hashtable;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ReferenceInfoVO
{
  
  private String ms2Name_;
  private boolean useInEvaluation_;
  private boolean positionAvailable_;
  private double[] correctRts_;
  private Hashtable<String,Double> adducts_ = new Hashtable<String,Double>();
  
  

  public ReferenceInfoVO(String ms2Name, double[] correctRts, boolean useInEvaluation, boolean positionAvailable) {
    this(ms2Name, useInEvaluation, positionAvailable);
    this.correctRts_ = correctRts;
  }
  
  public ReferenceInfoVO(String ms2Name, double correctRt, boolean useInEvaluation, boolean positionAvailable) {
    this(ms2Name, useInEvaluation, positionAvailable);
    this.correctRts_ = new double[1];
    this.correctRts_[0] = correctRt;
  }

  
  private ReferenceInfoVO(String ms2Name, boolean useInEvaluation, boolean positionAvailable)
  {
    this.ms2Name_ = ms2Name;
    this.useInEvaluation_ = useInEvaluation;
    this.positionAvailable_ = positionAvailable;
    this.adducts_ = new Hashtable<String,Double>();
  }
  
  public String getMS2Name(){
    return ms2Name_;
  }
  public boolean useInEvaluation()
  {
    return useInEvaluation_;
  }
  public boolean isPositionAvailable()
  {
    return positionAvailable_;
  }
  public double[] getCorrectRts()
  {
    return correctRts_;
  }
  public double[] setCorrectRts(double[] rts)
  {
    return this.correctRts_=rts;
  }
  public void addAdduct(String adduct, double mass) {
    this.adducts_.put(adduct, mass);
  }

  public Hashtable<String,Double> getAdducts()
  {
    return adducts_;
  }
  
}
