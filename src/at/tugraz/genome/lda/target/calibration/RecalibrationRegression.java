package at.tugraz.genome.lda.target.calibration;

import java.util.ArrayList;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;

import at.tugraz.genome.lda.utils.Pair;

public class RecalibrationRegression
{
	private ArrayList<Pair<Double,Double>> differences_;
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
	protected RecalibrationRegression(ArrayList<Pair<Double,Double>> differences, double grouping, String dataType, String lipidClass)
	{
		this(differences, grouping);
		this.dataType_ = dataType;
		this.lipidClass_ = lipidClass;
	}
	
	/**
	 * @param differences
	 * @param grouping
	 */
	private RecalibrationRegression(ArrayList<Pair<Double,Double>> differences, double grouping)
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
	
  private ArrayList<Pair<Double,Double>> clusterKeys(ArrayList<Pair<Double,Double>> differences)
  {
  	ArrayList<Pair<Double,Double>> clustered = new ArrayList<Pair<Double,Double>>();
  	ArrayList<Double> intervals = new ArrayList<Double>();
  	for (int i=0; i<=(getMax(differences)-getMin(differences))/grouping_+1;i++)
  	{
  		intervals.add(getMin(differences)+i*grouping_);
  		if (i>0)
  		{
  			int count = 0;
  			double xSum = 0.0;
    		double ySum = 0.0;
    		for (Pair<Double,Double> dataPoint : differences)
      	{
      		if (dataPoint.getKey()>=intervals.get(i-1) && dataPoint.getKey()<=intervals.get(i))
      		{
      			xSum += dataPoint.getKey();
      			ySum += dataPoint.getValue();
      			count += 1;
//      			System.out.println("Point x: "+dataPoint.getKey()+" y: "+dataPoint.getValue());
      		}
      	}
    		if (count > 0)
    		{
    			clustered.add(new Pair<Double,Double>(xSum/count, ySum/count));
//      		System.out.println("Cluster: "+new Pair<Double,Double>(xSum/count, ySum/count));
    		}
  		}
  	}
  	return clustered;
  }
	
	private double getMin(ArrayList<Pair<Double,Double>> data){  
    double min = Integer.MAX_VALUE;  
    for(int i=0; i<data.size(); i++){  
      if(data.get(i).getKey()<min)  
        min = data.get(i).getKey();  
    }  
    return min;  
  } 
	
	private double getMax(ArrayList<Pair<Double,Double>> data){  
    double max = Integer.MIN_VALUE;  
    for(int i=0; i<data.size(); i++){  
      if(data.get(i).getKey()>max)  
        max = data.get(i).getKey();  
    }  
    return max;  
  }
	
	protected void removeDataPoint(Pair<Double,Double> dataPoint)
	{
		this.differences_.remove(dataPoint);
		initRegression();
	}
	
	public double getTargetRT(double referenceRT) throws OutOfRangeException
	{
		return referenceRT-function_.value(referenceRT);
	}
	
	protected ArrayList<Pair<Double,Double>> getDifferences()
	{
		return differences_;
	}
	
	/**
	 * Removing a data point can lead to either one or two less cluster points.
	 * This is relevant when checking if the minimum number of data points will be kept.
	 * @return
	 */
	protected int getClusteredWithoutDataPointSize(Pair<Double,Double> dataPoint)
	{
		ArrayList<Pair<Double,Double>> differences = new ArrayList<Pair<Double,Double>>(differences_);
		differences.remove(dataPoint);
		ArrayList<Pair<Double,Double>> clustered = clusterKeys(differences);
		return clustered.size();
	}
	
	protected ArrayList<Pair<Double,Double>> getClustered()
	{
		return clustered_;
	}
	
	protected PolynomialSplineFunction getFunction()
	{
		return function_;
	}
	
	protected double getLowerRTLimit()
	{
		return getMin(clustered_);
	}
	
	protected double getUpperRTLimit()
	{
		return getMax(clustered_);
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
