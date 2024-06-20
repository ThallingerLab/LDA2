package at.tugraz.genome.lda.target.calibration;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class AnchorPoint
{
	private String lipidClass_;
	private String lipidSpecies_;
	private Double xValue_;
	private Double yValue_;
	
	public AnchorPoint(String lipidClass, String lipidSpecies, Double xValue, Double yValue)
	{
		super();
		this.lipidClass_ = lipidClass;
		this.lipidSpecies_ = lipidSpecies;
		this.xValue_ = xValue;
		this.yValue_ = yValue;
	}

	public String getLipidClass()
	{
		return lipidClass_;
	}

	public String getLipidSpecies()
	{
		return lipidSpecies_;
	}

	public Double getxValue()
	{
		return xValue_;
	}

	public Double getyValue()
	{
		return yValue_;
	}
	
	
	
}
