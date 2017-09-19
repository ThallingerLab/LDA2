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
public class AddAnalyteVO
{
  private String name_;
  private String formula_;
  private String modName_;
  private String modFormula_;
  private String mzTolerance_;
  private String exactMass_;
  private String charge_;
  private String rt_;
  
  
  public AddAnalyteVO(String name, String formula, String modName, String modFormula, String mzTolerance,
      String exactMass, String charge, String rt)
  {
    super();
    this.name_ = name;
    this.formula_ = formula;
    this.modName_ = modName;
    this.modFormula_ = modFormula;
    this.mzTolerance_ = mzTolerance;
    this.exactMass_ = exactMass;
    this.charge_ = charge;
    this.rt_ = rt;
  }


  public String getName()
  {
    if (name_.lastIndexOf(":")!=-1)
      return name_.substring(0,name_.lastIndexOf(":"));
    else
      return name_;
  }

  public Integer getDoubleBonds(){
    if (name_.lastIndexOf(":")!=-1){
      return new Integer(name_.substring(name_.lastIndexOf(":")+1));
    } else
      return -1;
  }

  public String getFormula()
  {
    return formula_;
  }


  public String getMzTolerance()
  {
    return mzTolerance_;
  }


  public String getExactMass()
  {
    return exactMass_;
  }


  public String getModName()
  {
    return modName_;
  }


  public String getModFormula()
  {
    return modFormula_;
  }


  public String getCharge()
  {
    return charge_;
  }


  public String getRt()
  {
    return rt_;
  }
  
  
}
