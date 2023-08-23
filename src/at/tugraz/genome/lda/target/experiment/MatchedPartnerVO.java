package at.tugraz.genome.lda.target.experiment;

import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.target.IsotopeLabelVO;

public class MatchedPartnerVO
{
	private String lipidClass_;
	private LipidomicsMSnSet standard_;
	private LipidomicsMSnSet isotopologue_;
	private boolean useForCalibration_;
	private IsotopeLabelVO label_;
	
	protected MatchedPartnerVO(String lipidClass, LipidomicsMSnSet standard, LipidomicsMSnSet isotopologue, boolean useForCalibration, IsotopeLabelVO label)
	{
		this.lipidClass_ = lipidClass;
		this.standard_ = standard;
		this.isotopologue_ = isotopologue;
		this.useForCalibration_ = useForCalibration;
		this.label_ = label;
	}
	
	public String getLipidClass()
	{
		return lipidClass_;
	}

	public void setLipidClass(String lipidClass)
	{
		this.lipidClass_ = lipidClass;
	}

	public LipidomicsMSnSet getStandard()
	{
		return standard_;
	}

	public void setStandard(LipidomicsMSnSet standard)
	{
		this.standard_ = standard;
	}

	public LipidomicsMSnSet getIsotopologue()
	{
		return isotopologue_;
	}

	public void setIsotopologue(LipidomicsMSnSet isotopologue)
	{
		this.isotopologue_ = isotopologue;
	}

	public boolean isUseForCalibration()
	{
		return useForCalibration_;
	}

	public void setUseForCalibration(boolean useForCalibration)
	{
		this.useForCalibration_ = useForCalibration;
	}

	public IsotopeLabelVO getLabel()
	{
		return label_;
	}

	public void setLabel(IsotopeLabelVO label)
	{
		this.label_ = label;
	}
}
