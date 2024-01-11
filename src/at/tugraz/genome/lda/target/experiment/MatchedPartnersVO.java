package at.tugraz.genome.lda.target.experiment;

import java.util.Vector;

import at.tugraz.genome.lda.quantification.LipidParameterSet;

@Deprecated
public class MatchedPartnersVO
{
	private String lipidClass_;
	private Vector<LipidParameterSet> standards_;
	private Vector<MatchedPartnerVO> matchedPartners_;
	
	protected MatchedPartnersVO(String lipidClass, Vector<LipidParameterSet> standards)
	{
		this.lipidClass_ = lipidClass;
		this.standards_ = standards;
		this.matchedPartners_ = new Vector<MatchedPartnerVO>();
	}
	
	protected void addIsotopologuePartner(LipidParameterSet standard, LipidParameterSet isotopologuePartner, boolean useForCalibration)
	{
		this.matchedPartners_.add(new MatchedPartnerVO(lipidClass_, standard, isotopologuePartner, useForCalibration));
	}
	
	protected Vector<LipidParameterSet> getIsotopologuePartners(LipidParameterSet standard)
	{
		Vector<LipidParameterSet> partners = new Vector<LipidParameterSet>();
		for (MatchedPartnerVO matched : matchedPartners_)
		{
			if (matched.getStandard().equals(standard))
			{
				partners.add(matched.getIsotopologue());
			}
		}
		return partners;
	}
	
	protected String getLipidClass()
	{
		return this.lipidClass_;
	}
	
	protected Vector<LipidParameterSet> getStandards()
	{
		return this.standards_;
	}
	
	protected Vector<MatchedPartnerVO> getMatchedPartners()
	{
		return this.matchedPartners_;
	}
	
	protected class MatchedPartnerVO
	{
		private String lipidClass_;
		private LipidParameterSet standard_;
		private LipidParameterSet isotopologue_;
		private boolean useForCalibration_;
		
		protected MatchedPartnerVO(String lipidClass, LipidParameterSet standard, LipidParameterSet isotopologue, boolean useForCalibration)
		{
			this.lipidClass_ = lipidClass;
			this.standard_ = standard;
			this.isotopologue_ = isotopologue;
			this.useForCalibration_ = useForCalibration;
		}
		
		public String getLipidClass()
		{
			return lipidClass_;
		}

		public void setLipidClass(String lipidClass)
		{
			this.lipidClass_ = lipidClass;
		}

		public LipidParameterSet getStandard()
		{
			return standard_;
		}

		public void setStandard(LipidParameterSet standard)
		{
			this.standard_ = standard;
		}

		public LipidParameterSet getIsotopologue()
		{
			return isotopologue_;
		}

		public void setIsotopologue(LipidParameterSet isotopologue)
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
	}
}
