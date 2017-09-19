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
public class ResultDisplaySettingsVO
{
  private String type_;
  private int isStand_;
  private int esStand_;
  private boolean dilution_;
  private boolean au_;
  private boolean percent_;
  private String divisorMagnitude_;
 
  public ResultDisplaySettingsVO(String type, int isStand, int esStand,
      boolean dilution, boolean au)
  {
    this(type,isStand,esStand,dilution,au,"");
  }
  
  public ResultDisplaySettingsVO(String type, int isStand, int esStand,
      boolean dilution, boolean au, String divisorMagnitude)
  {
    super();
    this.type_ = type;
    this.isStand_ = isStand;
    this.esStand_ = esStand;
    this.dilution_ = dilution;
    this.au_ = au;
    this.percent_ = false;
    this.divisorMagnitude_ = divisorMagnitude;
  }
  
  public String getType()
  {
    return type_;
  }
  
  public void setType(String type_)
  {
    this.type_ = type_;
  }

  public int getISStandMethod()
  {
    return isStand_;
  }
  public int getESStandMethod()
  {
    return esStand_;
  }
  public boolean considerDilution()
  {
    return dilution_;
  }
  public boolean isAu()
  {
    return au_;
  }

  public boolean isPercent()
  {
    return percent_;
  }

  public void setPercent(boolean percent_)
  {
    this.percent_ = percent_;
  }

  public String getDivisorMagnitude()
  {
    return divisorMagnitude_;
  }

  public void setDivisorMagnitude(String divisorMagnitude)
  {
    this.divisorMagnitude_ = divisorMagnitude;
  }

  public ResultDisplaySettingsVO(ResultDisplaySettingsVO from){
    this(from.type_,from.isStand_,from.esStand_,from.dilution_, from.au_, from.divisorMagnitude_);
  }
}
