package at.tugraz.genome.lda.target;

public abstract class AbstractRegression
{
	private double[] coefficients_;
		
	protected abstract double[] initRegression();
	
	protected double getAdjustFactor(int x)
	{
		double adjustFactor = 0;
		if (coefficients_ != null) 
		{
			for (int i=0; i<coefficients_.length; i++)
			{
				adjustFactor += coefficients_[i]*Math.pow(x,i);
			}
		}
		return adjustFactor;
	}
	
	protected double getAdjustFactor(double x)
	{
		double adjustFactor = 0;
		if (coefficients_ != null) 
		{
			for (int i=0; i<coefficients_.length; i++)
			{
				adjustFactor += coefficients_[i]*Math.pow(x,i);
			}
		}
		return adjustFactor;
	}
	
	protected void setCoefficients(double[] coefficients)
	{
		this.coefficients_ = coefficients;
	}
	
	protected double[] getCoefficients()
	{
		return coefficients_;
	}
	
}
