package at.tugraz.genome.lda.export;

import java.util.LinkedHashMap;

import at.tugraz.genome.lda.msn.vos.FattyAcidVO;

/**
 * Helper class for SILDA Analysis.
 * @author Leonida M. Lamp
 *
 */
public class OmegaCollector
{
	/**
	 * experiment, class, totalFA
	 * This denotes the total fatty acid content (omega inspecific), for normalization.
	 */
	LinkedHashMap<String,LinkedHashMap<String,Double>> totalFAPerClass_;
	/**
	 * experiment, class, omegapos, totalFA
	 * this is for the omega per class heatmap
	 */
	LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<Integer,Double>>> totalOmegaFAPerClass_;
	/**
	 * experiment, class, FA_VO, partner_FA_VO, totalFA
	 * this is for the fa vs fa heatmap
	 */
	LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<FattyAcidVO,LinkedHashMap<FattyAcidVO,Double>>>> faContentToPartnerContent_;
	
	public OmegaCollector()
	{
		totalFAPerClass_ = new LinkedHashMap<String,LinkedHashMap<String,Double>>();
		totalOmegaFAPerClass_ = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<Integer,Double>>>();
		faContentToPartnerContent_ = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<FattyAcidVO,LinkedHashMap<FattyAcidVO,Double>>>>();
	}
	
	public void addToTotalFAPerClass(String exp, String lClass, Double amount)
	{
		if (!totalFAPerClass_.containsKey(exp))
		{
			totalFAPerClass_.put(exp, new LinkedHashMap<String,Double>());
		}
		if (!totalFAPerClass_.get(exp).containsKey(lClass))
		{
			totalFAPerClass_.get(exp).put(lClass, amount);
		}
		else
		{
			double before = totalFAPerClass_.get(exp).get(lClass);
			totalFAPerClass_.get(exp).put(lClass, before+amount);
		}
	}
	
	public void addToTotalOmegaFAPerClass(String exp, String lClass, Integer omegaPos, Double amount)
	{
		if (!totalOmegaFAPerClass_.containsKey(exp))
		{
			totalOmegaFAPerClass_.put(exp, new LinkedHashMap<String,LinkedHashMap<Integer,Double>>());
		}
		if (!totalOmegaFAPerClass_.get(exp).containsKey(lClass))
		{
			totalOmegaFAPerClass_.get(exp).put(lClass, new LinkedHashMap<Integer,Double>());
		}
		if (!totalOmegaFAPerClass_.get(exp).get(lClass).containsKey(omegaPos))
		{
			totalOmegaFAPerClass_.get(exp).get(lClass).put(omegaPos,amount);
		}
		else
		{
			double before = totalOmegaFAPerClass_.get(exp).get(lClass).get(omegaPos);
			totalOmegaFAPerClass_.get(exp).get(lClass).put(omegaPos, before+amount);
		}
	}
	
	public void addTofaContentToPartnerContent(String exp, String lClass, FattyAcidVO vo, FattyAcidVO partnerVO, Double amount)
	{
		faContentToPartnerContent_ = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<FattyAcidVO,LinkedHashMap<FattyAcidVO,Double>>>>();
		
		
		if (!faContentToPartnerContent_.containsKey(exp))
		{
			faContentToPartnerContent_.put(exp, new LinkedHashMap<String,LinkedHashMap<FattyAcidVO,LinkedHashMap<FattyAcidVO,Double>>>());
		}
		if (!faContentToPartnerContent_.get(exp).containsKey(lClass))
		{
			faContentToPartnerContent_.get(exp).put(lClass, new LinkedHashMap<FattyAcidVO,LinkedHashMap<FattyAcidVO,Double>>());
		}
		if (!faContentToPartnerContent_.get(exp).get(lClass).containsKey(vo))
		{
			faContentToPartnerContent_.get(exp).get(lClass).put(vo, new LinkedHashMap<FattyAcidVO,Double>());
		}
		if (!faContentToPartnerContent_.get(exp).get(lClass).get(vo).containsKey(partnerVO))
		{
			faContentToPartnerContent_.get(exp).get(lClass).get(vo).put(partnerVO, amount);
		}
		else
		{
			double before = faContentToPartnerContent_.get(exp).get(lClass).get(vo).get(partnerVO);
			faContentToPartnerContent_.get(exp).get(lClass).get(vo).put(partnerVO, before+amount);
		}
	}
}
