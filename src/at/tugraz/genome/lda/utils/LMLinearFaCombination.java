/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2021 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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
 * Implements the LevenbergMarquardtOptimizer for a variable set of fatty acid combinations
 * @author Juergen Hartler
 *
 */
public class LMLinearFaCombination extends LevenbergMarquardtOptimizer
{
  
  /** the values of the input parameters (a m x n matrix)*/
  private float[][] values_;
  /** the vector of equation results (a m x 1 matrix) */
  private float[][] observations_;
  /** the sum of all intensities*/
  private float totalInt_;

  
  public LMLinearFaCombination(float[][] values, float[] observations){
    values_ = values;
    observations_ = new float[observations.length][1];
    totalInt_ = 0f;
    for (int i=0; i!=observations.length;i++){
      observations_[i][0] = observations[i];
      totalInt_ += observations_[i][0];
    }
    totalInt_ = totalInt_/3f;
  }
  
  
  protected FloatMatrix calculateEquationResults(float[][] values,
      FloatMatrix paramsVector)
  {
    float[][] results = new float[values.length][1];
    for (int i=0;i!=values.length;i++){
      float equationResult = 0f;
      for (int j=0; j!=values[i].length;j++){
        equationResult += paramsVector.A[j][0]*values[i][j];
      }
      results[i][0] = equationResult;
    }
    return new FloatMatrix(results);
  }

  
  protected FloatMatrix calculateJacobianMatrix(float[][] values,
      float[][] parameters)
  {
    //in the linear case, the jacobian matrix is the same as the values matrix*/
    float[][] jacobian = new float[values.length][values[0].length];
    for (int i=0; i!=values.length; i++){
      for (int j=0; j!=values[i].length; j++){
        jacobian[i][j] = values[i][j];
      }
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

  
  public void fit() throws LMException
  {
    float[][] parameters = null;
    if (resultParams_ == null){
      int parametersNumber = values_[0].length;
      parameters = new float[parametersNumber][1];
      float relativeShare = totalInt_/((float)parametersNumber);
      for (int i=0; i!=parametersNumber; i++)
        parameters[i][0] = relativeShare;
    } else parameters = resultParams_.A;
    fit(parameters);
  }

  
  public float calculateFitValue(float[] input) throws LMException
  {
    float[][] values = new float[1][input.length];
    for (int i=0; i!=input.length; i++){
      values[0][i] = input[i];
    }
    float result = calculateEquationResults(values, resultParams_).A[0][0];
    return result;
  }

}
