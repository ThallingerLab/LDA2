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

package at.tugraz.genome;

import java.math.BigDecimal;

import at.tugraz.genome.maspectras.utils.Calculator;
import JSci.maths.statistics.TDistribution;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class IndependentTwoSamplesTTest
{
  double[] results_;
  
  public static void main(String[] args)
  {
    if (args.length!=5){
      if (args.length!=0){
        System.out.println("ERROR: you must provide 5 arguements! You provided just: "+args.length);
        System.out.println();
      }
      printUsage();
      return;
    }
    int sampleSize = 0;
    double mean1 = 0;
    double mean2 = 0;
    double stdev1 = 0;
    double stdev2 = 0;
    try{ sampleSize= Integer.parseInt(args[0]);}catch (NumberFormatException nfx){
      System.out.println("sampleSize is not integer format");
      System.out.println();
      return;
    }
    try{ mean1=Double.parseDouble(args[1].replace(",", "."));}catch (NumberFormatException nfx){
      System.out.println("mean1 is not double format");
      System.out.println();
      return;
    }
    try{ mean2=Double.parseDouble(args[2].replace(",", "."));}catch (NumberFormatException nfx){
      System.out.println("mean2 is not double format");
      System.out.println();
      return;
    }
    try{ stdev1=Double.parseDouble(args[3].replace(",", "."));}catch (NumberFormatException nfx){
      System.out.println("stdev1 is not double format");
      System.out.println();
      return;
    }
    try{ stdev2=Double.parseDouble(args[4].replace(",", "."));}catch (NumberFormatException nfx){
      System.out.println("stdev2 is not double format");
      System.out.println();
      return;
    }
    new IndependentTwoSamplesTTest(sampleSize,mean1,mean2,stdev1,stdev2);
  }

  public IndependentTwoSamplesTTest(int sampleSize, double mean1, double mean2, double stdev1, double stdev2){
    this(sampleSize, sampleSize,mean1,mean2,stdev1,stdev2);
  }
  
  public IndependentTwoSamplesTTest(int sampleSize1, int sampleSize2, double mean1, double mean2, double stdev1, double stdev2){
    results_ = calculatePValueForIndependentTwoSamplesTTest(sampleSize1, sampleSize2,mean1,mean2,stdev1,stdev2);
  }
  
  private double[] calculatePValueForIndependentTwoSamplesTTest(int s1Size, int s2Size, double mean1, double mean2, double stdev1, double stdev2){
    //this is for equal variance and sample size
//    double mixedStdev = Math.sqrt(0.5d*(Math.pow(stdev1, 2)+Math.pow(stdev2, 2)));
    double[] values = new double[3];
    double mixedStdev = Math.sqrt(Math.pow(stdev1, 2)/s1Size+Math.pow(stdev2, 2)/s2Size);
    double tValue = (mean1-mean2)/(mixedStdev);
    if (tValue<0) tValue = tValue*-1d;
//    int df = 2*sampleSize-2;
//    double dfDouble = ((sampleSize-1)*(Math.pow(stdev1, 4)+2*Math.pow(stdev1, 2)*Math.pow(stdev2, 2)+Math.pow(stdev2, 4)))/(Math.pow(stdev1, 4)+Math.pow(stdev2, 4));
    double dfDouble = Math.pow(Math.pow(stdev1, 2)/s1Size+Math.pow(stdev2, 2)/s2Size, 2)/(Math.pow(Math.pow(stdev1, 2)/s1Size, 2)/(s1Size-1)+Math.pow(Math.pow(stdev2, 2)/s2Size, 2)/(s2Size-1));
    int df = (int)Calculator.roundDBL(dfDouble, 0, BigDecimal.ROUND_HALF_UP);
    System.out.println("dfDouble: "+dfDouble+" "+df);
    TDistribution dist = new TDistribution(df);
    float pValue = (float)((1-dist.cumulative(tValue))*2);
    System.out.println("pValue: "+pValue);
    values[0] = pValue;
    values[1] = dfDouble;
    values[2] = df;
    return values;
  }

  
  private static void printUsage()
  {
    System.out.println("\nTwo Independent Sample T-test usage for equal sample sizes:");
    System.out.println();
    System.out.println("twoSampleTTest $sampleSize $mean1 $mean2 $stdev1 $stdev2");
  }

  public double[] getResults(){
    return this.results_;
  }
}
