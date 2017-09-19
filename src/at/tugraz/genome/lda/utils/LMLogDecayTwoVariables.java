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

package at.tugraz.genome.lda.utils;

import at.tugraz.genome.lda.exception.LMException;
import at.tugraz.genome.util.FloatMatrix;

/**
 * Implements the LevenbergMarquardtOptimizer for a logarithmic function of variable one
 * and an Euler decay function for variable two with slowed decrease depending on variable one
 * this function worked best for the lipid chromatography problem.
 * Equation:
 * f(x,y) = A*ln(B*x) + C*e^(-D*y + E*x) + F
 * @author Juergen Hartler
 *
 */
public class LMLogDecayTwoVariables extends LevenbergMarquardtOptimizer
{

  /** the values of the input parameters (a m x 2 matrix)*/
  private float[][] values_;
  /** the vector of equation results (a m x 1 matrix) */
  private float[][] observations_;
  
  /**
   * constructor setting the input values and the measured observations
   * @param values the values of the input parameters (a m x n matrix)
   * @param observations the measured results for the input values
   * @param parameters proposed starting parameters
   */
  public LMLogDecayTwoVariables(float[][] xValues, float[] observations, FloatMatrix parameters){
    values_ = xValues;
    observations_ = new float[observations.length][1];
    for (int i=0; i!=observations.length;i++){
      observations_[i][0] = observations[i];
    }
    if (parameters!=null) resultParams_ = parameters;
  }

  
  protected FloatMatrix calculateEquationResults(float[][] values,
      FloatMatrix paramsVector)
  {
    float[][] results = new float[values.length][1];
    for (int i=0;i!=values.length;i++){
      float equationResult = paramsVector.A[0][0]*(float)Math.log(values[i][0]*paramsVector.A[1][0]);
      equationResult += (float)(((double)paramsVector.A[2][0])*Math.exp((double)(paramsVector.A[4][0]*values[i][0]-paramsVector.A[3][0]*values[i][1])));
      equationResult += paramsVector.A[5][0];
      results[i][0] = equationResult;
    }
    return new FloatMatrix(results);
  }

  protected FloatMatrix calculateJacobianMatrix(float[][] values,
      float[][] parameters)
  {
    float[][] jacobian = new float[values.length][6];
    for (int i=0; i!=values.length; i++){      
      jacobian[i][0] = (float)Math.log(values[i][0]*parameters[1][0]);
      jacobian[i][1] = parameters[0][0]/parameters[1][0];
      jacobian[i][2] = (float)Math.exp((double)(parameters[4][0]*values[i][0]*values[i][1]-parameters[3][0]*values[i][1]));
      jacobian[i][3] = -1f*parameters[2][0]*values[i][1]*(float)Math.exp((double)(parameters[4][0]*values[i][0]-parameters[3][0]*values[i][1]));
      jacobian[i][4] = parameters[2][0]*values[i][0]*(float)Math.exp((double)(parameters[4][0]*values[i][0]-parameters[3][0]*values[i][1]));
      jacobian[i][5] = 1f;
    }
    return new FloatMatrix(jacobian);
  }

  protected float[][] getValues()
  {
    return values_;
  }

  protected float[][] getObservations()
  {
    return observations_;
  }

  public void fit() throws LMException{
    float[][] parameters = null;
    if (this.values_.length<6) throw new LMException("There are too less observations for a model fit!");
    if (resultParams_ == null){
      parameters = new float[6][1];
      parameters[0][0] = 0.5f;
      parameters[1][0] = 100f;
      parameters[2][0] = 10f;
      parameters[3][0] = 0.3f;
      parameters[4][0] = 0.001f;
      parameters[5][0] = 1f;
    } else parameters = resultParams_.A;
    fit(parameters);
  }

  public float calculateFitValue(float[] input) throws LMException
  {
    if (input.length!=2) throw new LMException("The input must consist of two variables");
    float[][] values = new float[1][2];
    values[0][0] = input[0];
    values[0][1] = input[1];
    float result = calculateEquationResults(values, resultParams_).A[0][0];
    return result;
  }

  protected float getLambdaStartValue(){
    return 5f;
  }
  
  protected float getLambdaMultiplierStartValue(){
    return 1.5f;
  }
  
}
