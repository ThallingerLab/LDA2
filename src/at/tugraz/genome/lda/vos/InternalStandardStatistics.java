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
public class InternalStandardStatistics
{
  String standardName;
  int howOftenFound;
  int amountIsotopesFoundInAll;
  float mean;
  float median;
  float stdev;
  float coefficentOfVariation;
  float lowerOutlierThreshold;
  float upperOutlierThreshold;
  
  public InternalStandardStatistics(String standardName,int howOftenFound,
      int amountIsotopesFoundInAll, float mean, float median, float stdev,
      float coefficentOfVariation, float lowerOutlierThreshold, float upperOutlierThreshold)
  {
    super();
    this.standardName = standardName;
    this.howOftenFound = howOftenFound;
    this.amountIsotopesFoundInAll = amountIsotopesFoundInAll;
    this.mean = mean;
    this.median = median;
    this.stdev = stdev;
    this.coefficentOfVariation = coefficentOfVariation;
    this.lowerOutlierThreshold = lowerOutlierThreshold;
    this.upperOutlierThreshold = upperOutlierThreshold;
  }
  
  
  
  public String getStandardName()
  {
    return standardName;
  }

  public int getHowOftenFound()
  {
    return howOftenFound;
  }
  public int getAmountIsotopesFoundInAll()
  {
    return amountIsotopesFoundInAll;
  }
  
  public float getMean()
  {
    return mean;
  }
  public float getMedian()
  {
    return median;
  }
  public float getStdev()
  {
    return stdev;
  }
  public float getCoefficentOfVariation()
  {
    return coefficentOfVariation;
  }
  public float getLowerOutlierThreshold()
  {
    return lowerOutlierThreshold;
  }
  public float getUpperOutlierThreshold()
  {
    return upperOutlierThreshold;
  }
  
}
