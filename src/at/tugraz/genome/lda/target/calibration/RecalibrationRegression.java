/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2023 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.target.calibration;

import java.util.ArrayList;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;

import javafx.util.Pair;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class RecalibrationRegression
{
	private ArrayList<AnchorPointVO> differences_;
	private double grouping_;
	private ArrayList<Pair<Double,Double>> clustered_;
	private PolynomialSplineFunction function_;
	private String dataType_ = CalibrationGraphPanel.PLOT_ALL;
	private String lipidClass_ = CalibrationGraphPanel.PLOT_ALL;
	
	/**
	 * @param differences
	 * @param grouping
	 * @param dataType
	 * @param lipidClass
	 */
	protected RecalibrationRegression(ArrayList<AnchorPointVO> differences, double grouping, String dataType, String lipidClass)
	{
		this(differences, grouping);
		this.dataType_ = dataType;
		this.lipidClass_ = lipidClass;
	}
	
	/**
	 * @param differences
	 * @param grouping
	 */
	private RecalibrationRegression(ArrayList<AnchorPointVO> differences, double grouping)
	{
		this.differences_ = differences;
		this.grouping_ = grouping;
		initRegression();
	}
	
	/**
	 * Splits the data into windows on the x-axis, then does a fit with apache splineinterpolator.
	 * If there are too few data points the function is set to null.
	 * 
	 */
	protected void initRegression()
	{
		initRegression(clusterKeys(differences_));
	}
	
	/**
	 * Splits the data into windows on the x-axis, then does a fit with apache splineinterpolator.
	 * If there are too few data points the function is set to null.
	 * @param clustered		the clustered data points
	 * 
	 */
	protected void initRegression(ArrayList<Pair<Double,Double>> clustered)
	{
		this.clustered_ = clustered;
		SplineInterpolator interpolator = new SplineInterpolator();
		double[] xValues = new double[this.clustered_.size()];
		double[] yValues = new double[this.clustered_.size()];
		for (int i=0; i<this.clustered_.size(); i++)
		{
			xValues[i] = this.clustered_.get(i).getKey();
			yValues[i] = this.clustered_.get(i).getValue();
		}
		if (this.clustered_.size()>2) //less than two data points cannot be interpolated
		{
			this.function_ = interpolator.interpolate(xValues, yValues);
		}
		else
		{
			this.function_ = null;
		}
	}
	
  private ArrayList<Pair<Double,Double>> clusterKeys(ArrayList<AnchorPointVO> differences)
  {
  	ArrayList<Pair<Double,Double>> clustered = new ArrayList<Pair<Double,Double>>();
  	ArrayList<Double> xValues = new ArrayList<Double>();
  	ArrayList<Double> intervals = new ArrayList<Double>();
  	for (int i=0; i<=(getMax(differences)-getMin(differences))/grouping_+1;i++)
  	{
  		intervals.add(getMin(differences)+i*grouping_);
  		if (i>0)
  		{
  			int count = 0;
  			double xSum = 0.0;
    		double ySum = 0.0;
    		for (AnchorPointVO dataPoint : differences)
      	{
      		if (dataPoint.getxValue()>=intervals.get(i-1) && dataPoint.getxValue()<=intervals.get(i))
      		{
      			xSum += dataPoint.getxValue();
      			ySum += dataPoint.getyValue();
      			count += 1;
      		}
      	}
    		if (count > 0 && 
    				!xValues.contains(xSum/count)) //this can happen in very rare cases when the interval border coincides exactly with a value
    		{
    			clustered.add(new Pair<Double,Double>(xSum/count, ySum/count));
    			xValues.add(xSum/count);
    		}
  		}
  	}
  	return clustered;
  }
	
	private double getMin(ArrayList<AnchorPointVO> data){  
    double min = Integer.MAX_VALUE;  
    for(int i=0; i<data.size(); i++){  
      if(data.get(i).getxValue()<min)  
        min = data.get(i).getxValue();  
    }  
    return min;  
  } 
	
	private double getMax(ArrayList<AnchorPointVO> data){  
    double max = Integer.MIN_VALUE;  
    for(int i=0; i<data.size(); i++){  
      if(data.get(i).getxValue()>max)  
        max = data.get(i).getxValue();  
    }  
    return max;  
  }
	
	private double getMinClustered(ArrayList<Pair<Double,Double>> data){  
    double min = Integer.MAX_VALUE;  
    for(int i=0; i<data.size(); i++){  
      if(data.get(i).getKey()<min)  
        min = data.get(i).getKey();  
    }  
    return min;  
  } 
	
	private double getMaxClustered(ArrayList<Pair<Double,Double>> data){  
    double max = Integer.MIN_VALUE;  
    for(int i=0; i<data.size(); i++){  
      if(data.get(i).getKey()>max)  
        max = data.get(i).getKey();  
    }  
    return max;  
  }
	
	protected void removeDataPoint(AnchorPointVO dataPoint)
	{
		this.differences_.remove(dataPoint);
		initRegression();
	}
	
	public double getTargetRT(double referenceRT) throws OutOfRangeException
	{
		return referenceRT-function_.value(referenceRT);
	}
	
	protected ArrayList<AnchorPointVO> getDifferences()
	{
		return differences_;
	}
	
	/**
	 * Removing a data point can lead to either one or two less cluster points.
	 * This is relevant when checking if the minimum number of data points will be kept.
	 * @param dataPoint
	 * @return
	 */
	protected int getClusteredWithoutDataPointSize(AnchorPointVO dataPoint)
	{
		ArrayList<AnchorPointVO> differences = new ArrayList<AnchorPointVO>(differences_);
		differences.remove(dataPoint);
		ArrayList<Pair<Double,Double>> clustered = clusterKeys(differences);
		return clustered.size();
	}
	
	public ArrayList<Pair<Double,Double>> getClustered()
	{
		return clustered_;
	}
	
	protected PolynomialSplineFunction getFunction()
	{
		return function_;
	}
	
	protected double getLowerRTLimit()
	{
		return getMinClustered(clustered_);
	}
	
	protected double getUpperRTLimit()
	{
		return getMaxClustered(clustered_);
	}

	public String getDataType()
	{
		return dataType_;
	}

	public void setDataType(String dataType)
	{
		this.dataType_ = dataType;
	}

	public String getLipidClass()
	{
		return lipidClass_;
	}

	public void setLipidClass(String lipidClass)
	{
		this.lipidClass_ = lipidClass;
	}
	
	public void setGrouping(double grouping)
	{
		this.grouping_ = grouping;
		initRegression();
	}
	
}
