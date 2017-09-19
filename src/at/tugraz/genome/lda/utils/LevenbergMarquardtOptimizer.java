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

import java.util.Vector;

import at.tugraz.genome.lda.exception.LMException;
import at.tugraz.genome.util.FloatMatrix;

/**
 * Abstract class for fitting an arbitrary curve to measured values by Levenberg-Marquardt algorithm
 * @author Juergen Hartler
 *
 */
public abstract class LevenbergMarquardtOptimizer
{

  /** the parameters of the equation returned by Levenberg-Marquardt algorithm*/
  protected FloatMatrix resultParams_;
  /** the chi-squared value for the equation*/
  protected double resultChiSqr_;
  /** the lambda that was used for the final step of Levenberg-Marquardt algorithm*/
  protected float resultLambda_;
  /** if the deviation from the fitted curve should not allow for a higher deviation*/
  protected Float maxDev_;
  
  /**
   * executes the fitting operation of an arbitrary curve (defined by the implementation class)
   * @param initParameters parameters for the parameters vector - it is a n x 1 matrix
   * @throws LMException exception if model adaption is not possible
   */
  protected void fit(float[][] initParameters) throws LMException{
    
    maxDev_ = null;
    FloatMatrix paramsVector = new FloatMatrix(initParameters);
    FloatMatrix resultVector = new FloatMatrix(getObservations());

//    int df = getDegreesOfFreedom(paramsVector,resultVector);
//    ChiSqrDistribution dist = new ChiSqrDistribution(df);
    //double pValue = 1e-12d;
    double previousChiSqr = Double.POSITIVE_INFINITY;
    
    float lambda = getLambdaStartValue();
    float vLambda = getLambdaMultiplierStartValue();
    // algorithm iterates until a pre-defined number of iterations or if chi squared does not change anymore 
    for (int i=0; i!=getMaximumOfIterations(); i++){
      @SuppressWarnings("rawtypes")
      Vector paramsAndLambda = levenbergMarquadtIteration(paramsVector,getValues(),resultVector,lambda,vLambda);
      paramsVector = (FloatMatrix)paramsAndLambda.get(0);
      lambda = (Float)paramsAndLambda.get(1);
      double chiSquared = calculateChiSquared(paramsVector,getValues(),resultVector);
      // this is for stopping the iterations if the result does not change anymore
      if (/*chiSquared<chiSqr ||*/ previousChiSqr == chiSquared && !Double.isInfinite(chiSquared)){
        previousChiSqr = chiSquared;
        //System.out.println("Break at: "+i+" ; "+dist.cumulative(chiSquared));
        break;
      }
      previousChiSqr = chiSquared;
    }
    //pValue = 5e-7;
    //double chiSqr = dist.inverse(pValue);
    //if (previousChiSqr>chiSqr) throw new LMException("The fit did not converge - no fitting possible!");
    
    resultParams_ = paramsVector;
    resultChiSqr_ = previousChiSqr;
    resultLambda_ = lambda;
  }
  
  /**
   * calculates the degrees of freedom for chi squared calculation
   * @param params the current parameters vector
   * @param resultVector the vector of equation results (a m x 1 matrix)
   * @return degrees of freedom for chi squared calculation
   */
  private int getDegreesOfFreedom(FloatMatrix params,FloatMatrix resultVector){
    int df = resultVector.m-params.m;
    return df;
  }
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  /**
   * performs one Levenberg-Marquardt iteration 
   * @param params the current parameters vector (a m x 1 matrix)
   * @param values the values of the input parameters (a m x n matrix)
   * @param resultVector the measured results of the input parameters (a m x 1 matrix)
   * @param lambdaBefore the current lambda which was used for the last LM cycle - lambda is the damping parameter of the equation
   * @param vLambda multiplication factor to change lambda
   * @return a Vector containing the corrected parameters of this cycle and the used lambda for this correction
   * @throws LMException exception if model adaption is not possible
   */
  private Vector levenbergMarquadtIteration(FloatMatrix params, float[][] values,FloatMatrix resultVector, float lambdaBefore, float vLambda) throws LMException{    
    Vector paramsAndLambda = new Vector();
    // for each LM iteration - various versions of the Jacobian matrix have to be calculated
    FloatMatrix jacobian = calculateJacobianMatrix(values,params.A);    
    FloatMatrix jacobianTransposed = jacobian.transpose();
    FloatMatrix jacobianProduct = jacobianTransposed.times(jacobian);
    FloatMatrix jacobianDiag = getDiagonalMatrix(jacobianProduct);
    // the lambda changes in each LM cycle - this method finds the best one to use
    float lambda = detectBestLambda(params, values,resultVector, lambdaBefore, vLambda,jacobianTransposed,jacobianProduct,
        jacobianDiag);
    // with the lambda, we can correct the parameters (the outcome of each LM iteration)
    FloatMatrix correctedParams = solveOneLMCycle(params,values,resultVector,lambda,jacobianTransposed,jacobianProduct,jacobianDiag);
    paramsAndLambda.add(correctedParams);
    paramsAndLambda.add(lambda);
    return paramsAndLambda;

  }
  
  /**
   * detects the lambda that should be used for this cycle according to Marquardts recommendation:
   * the LM equation is solved with lambda and lambda/vLambda;
   * the sum of squares is calculated for the current model, the one corrected with lambda and the one with lambda/vLambda
   * if lambda/vLambda leads to an improvement lambda/vLambda is the next proposed lambda,
   * else if lambda leads to a reduced sum of squares, lambda is taken,
   * if none of these improve the sum of squares, lambda is increased by successive multiplication of vLambda until the sum of squares are reduced
   * @param params the current parameters vector (a m x 1 matrix)
   * @param values the values of the input parameters (a m x n matrix)
   * @param resultVector the measured results of the input parameters (a m x 1 matrix)
   * @param lambdaBefore the current lambda which was used for the last LM cycle - lambda is the damping parameter of the equation
   * @param vLambda multiplication factor to change lambda
   * @param jacobianTransposed the transposed Jacobian matrix
   * @param jacobianProduct the product of transposed Jacobian matrix and the Jacobian matrix
   * @param jacobianDiag diagonal matrix of the Jacobian product
   * @return the proposed lambda value to be used
   * @throws LMException exception if model adaption is not possible
   */
  private float detectBestLambda(FloatMatrix params, float[][] values,FloatMatrix resultVector, float lambdaBefore, float vLambda,
      FloatMatrix jacobianTransposed, FloatMatrix jacobianProduct, FloatMatrix jacobianDiag) throws LMException{
    FloatMatrix currentResidues = calculateResidues(resultVector,values,params);
    float sumOfSquaresCurrent = calculateSumOfSquares(currentResidues);
    FloatMatrix paramsLambda = solveOneLMCycle(params,values,resultVector,lambdaBefore,jacobianTransposed,jacobianProduct,jacobianDiag);
    float sumOfSquaresLambda = Float.NaN;
    if (paramsLambda!=null)
      sumOfSquaresLambda = calculateSumOfSquares(calculateResidues(resultVector,values,paramsLambda));
    FloatMatrix paramsLambdaSmaller = solveOneLMCycle(params,values,resultVector,lambdaBefore/vLambda,jacobianTransposed,jacobianProduct,jacobianDiag);
    float sumOfSquaresLambdaSmaller = Float.NaN;
    if (paramsLambdaSmaller!=null)
      sumOfSquaresLambdaSmaller = calculateSumOfSquares(calculateResidues(resultVector,values,paramsLambdaSmaller));
    if ((!Float.isNaN(sumOfSquaresLambda) && sumOfSquaresLambda<sumOfSquaresCurrent) || (!Float.isNaN(sumOfSquaresLambdaSmaller) && sumOfSquaresLambdaSmaller<sumOfSquaresCurrent)){
      if (Float.isNaN(sumOfSquaresLambdaSmaller)) return lambdaBefore;
      else if (Float.isNaN(sumOfSquaresLambda)) return (lambdaBefore/vLambda);
      else if (sumOfSquaresLambdaSmaller<sumOfSquaresCurrent) return (lambdaBefore/vLambda);
      else return lambdaBefore;
    }else{
      float lambda = lambdaBefore;
      while (Float.isNaN(sumOfSquaresLambda) || sumOfSquaresLambda>sumOfSquaresCurrent){
        lambda = lambda*vLambda;
        if (Float.isInfinite(lambda)) throw new LMException("The curve cannot be fitted - singular matrix");
        paramsLambda = solveOneLMCycle(params,values,resultVector,lambda,jacobianTransposed,jacobianProduct,jacobianDiag);
        sumOfSquaresLambda = Float.NaN;
        if (paramsLambda!=null)
          sumOfSquaresLambda = calculateSumOfSquares(calculateResidues(resultVector,values,paramsLambda));
      }
      return lambda;
    }
  }

  /**
   * solves the Levenberg Marquardt equation to calculate delta:
   * (JtJ + lambda*diag(JtJ))*delta = Jt[y - f(beta))]
   * J ...           jacobian matrix
   * Jt ...          transposed jacobian matrix
   * lambda ...      the damping parameter
   * diag () ...     diagonal matrix
   * delta ...       delta of parameters (the values we want to find out)
   * y ...           the results vector (measured values)
   * f (beta) ...    the results calculated by solving the equation
   * @param params the current parameters vector (a m x 1 matrix)
   * @param values the values of the input parameters (a m x n matrix)
   * @param resultVector the measured results of the input parameters (a m x 1 matrix)
   * @param lambda the damping factor of the equation
   * @param jacobianTransposed the transposed Jacobian matrix
   * @param jacobianProduct the product of transposed Jacobian matrix and the Jacobian matrix
   * @param jacobianDiag diagonal matrix of the Jacobian product
   * @return the corrected parameters vector of this LM cycle
   */
  private FloatMatrix solveOneLMCycle(FloatMatrix params, float[][] values,FloatMatrix resultVector, float lambda, FloatMatrix jacobianTransposed,
      FloatMatrix jacobianProduct, FloatMatrix jacobianDiag) {
    // calculates the residus (deviation from measured values - y-f(beta))) 
    FloatMatrix residues = calculateResidues(resultVector,values,params);
    // calculates the right side of the equation: Jt[y - f(beta))
    FloatMatrix rightSide = jacobianTransposed.times(residues);
    // calculates the lambda matrix: lambda*diag(JtJ)
    FloatMatrix lambdaMatrix = new FloatMatrix(jacobianDiag.A.clone());
    lambdaMatrix = lambdaMatrix.times(lambda);
    // the left side of the equation: (JtJ + lambda*diag(JtJ))
    FloatMatrix leftSide = jacobianProduct.plus(lambdaMatrix);
    FloatMatrix delta = null;
    // the parameter corrective values we want to find out: delta = (JtJ + lambda*diag(JtJ))^-1 * Jt[y - f(beta))]
    delta = leftSide.inverse().times(rightSide);
    // correct the parameters by delta
    FloatMatrix correctedParams = params.plus(delta);
    return correctedParams;
  }

  private FloatMatrix calculateResidues(FloatMatrix resultVector, float[][] values, FloatMatrix paramsVector){
    float[][] residues = new float[values.length][1];
    FloatMatrix results = calculateEquationResults(values, paramsVector);
    for (int i=0;i!=values.length;i++){
      residues[i][0] = resultVector.A[i][0]-results.A[i][0];
    }
    return new FloatMatrix(residues);
  }

  /**
   * calculates the diagonal matrix of a given matrix
   * @param in the given matrix
   * @return the diagonal matrix
   */
  private FloatMatrix getDiagonalMatrix(FloatMatrix in){
    FloatMatrix mx = in.copy();
    for (int i=0;i!=mx.A.length;i++){
      for (int j=0; j!=mx.A[i].length;j++){
        if (i!=j)mx.A[i][j] = 0f;
      }
    }
    return mx;
  }
  
  /**
   * calculates the sum of squares of the calculated residues (a m x 1 matrix)
   * @param res the residues (a m x 1 matrix)
   * @return the sum of squares of this matrix
   */
  private float calculateSumOfSquares(FloatMatrix res){
    float sumOfSquares = 0f;
    for (int i=0; i!=res.m; i++){
      sumOfSquares += Math.pow(res.A[i][0],2d);
    }
    return sumOfSquares;
  }
  
  /**
   * calculates the chi squared value for the equation
   * @param params the current parameters vector (a m x 1 matrix)
   * @param values the values of the input parameters (a m x n matrix)
   * @param resultVector the measured results of the input parameters (a m x 1 matrix)
   * @return chi squared value
   */
  private float calculateChiSquared(FloatMatrix params, float[][] values,FloatMatrix resultVector){
    FloatMatrix residues = calculateResidues(resultVector,values,params);
    FloatMatrix results = calculateEquationResults(values,params);
    float chiSquared = 0f;
    for (int i=0;i!=values.length;i++){
      chiSquared +=   Math.pow(residues.A[i][0],2d)/Math.abs(results.A[i][0]);
    }
    int df = getDegreesOfFreedom(params,resultVector);
    chiSquared = chiSquared/df;
    return chiSquared;
  }
  
  /**
   * 
   * @return the mean deviation from the calculated curve - if a max Deviation value was set, this value is returned instead
   * @throws LMException exception if model adaption is not possible
   */
  public float getMeanDeviation() throws LMException{
    if (maxDev_!=null) return maxDev_;
    else{
      FloatMatrix res = calculateResidues(new FloatMatrix(getObservations()), getValues(), resultParams_);
      float sumOfSquares = calculateSumOfSquares(res);
      return (float)Math.sqrt(sumOfSquares/((float)res.m));
    }
  }
  
  /**
   * sets a max deviation value instead of the mean deviation
   * @param maxDev value for maximum deviation
   */
  public void setMaxDeviation(float maxDev){
    maxDev_ = maxDev;
  }
  
  /**
   * 
   * @return the lambda value to start the LM algorithm
   */
  protected float getLambdaStartValue(){
    return 10f;
  }
  
  /**
   * 
   * @return the multiplication factor vor lambda
   */
  protected float getLambdaMultiplierStartValue(){
    return 1.5f;
  }

  /**
   * 
   * @return the maximum amount of LM iterations
   */
  protected int getMaximumOfIterations(){
    return 10000;
  }
  
  /**
   * calculates an results vector based on the fitted equation
   * @param values the values of the input parameters (a m x n matrix)
   * @param paramsVector the measured results of the input parameters (a m x 1 matrix)
   * @return equation results
   */
  protected abstract FloatMatrix calculateEquationResults(float[][] values, FloatMatrix paramsVector);
  
  /**
   * calculates the Jacobian matrix based on the used function and the input values
   * @param values the values of the input parameters (a m x n matrix)
   * @param parameters the current parameters vector (a m x 1 matrix)
   * @return Jacobian matrix
   */
  protected abstract FloatMatrix calculateJacobianMatrix(float[][] values, float[][] parameters);
  
  /**
   * 
   * @return matrix of input values
   */
  protected abstract float[][] getValues();
  
  /**
   * 
   * @return the measured results for the input values
   */
  protected abstract float[][] getObservations();
  
  /**
   * initiates the fitting of the function
   * @throws LMException exception if model adaption is not possible
   */
  public abstract void fit() throws LMException;
  
  /**
   * calculates the equation result of a fitted model for (a) certain input value/s
   * @param input the input values for which the model outcome should be calculated
   * @return the equation result of a fitted model for (a) certain input value/s
   * @throws LMException exception if model adaption is not possible
   */
  public abstract float calculateFitValue(float[] input) throws LMException;

  /**
   * 
   * @return the fitted parameters of the equation
   */
  public FloatMatrix getResultParams()
  {
    return resultParams_;
  }

  /**
   * 
   * @return the chi squared value where the fit did stop
   */
  public double getResultChiSqr()
  {
    return resultChiSqr_;
  }

  /**
   * 
   * @return the lambda used for the final fitting step
   */
  public float getResultLambda()
  {
    return resultLambda_;
  }
  
    
}
