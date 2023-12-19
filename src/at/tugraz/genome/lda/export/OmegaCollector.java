package at.tugraz.genome.lda.export;

import java.util.LinkedHashMap;
import java.util.Vector;

/**
 * Helper class for SILDA Analysis.
 * @author Leonida M. Lamp
 *
 */
public class OmegaCollector
{
	public final static String DESCRIPTION_SATURATED = "saturated";
	public final static String DESCRIPTION_ASSIGNED = "unsaturated assigned";
	public final static String DESCRIPTION_UNASSIGNED = "unsaturated unassigned (u.a.)";
	public final static String DESCRIPTION_UNASSIGNED_EVEN_CHAIN = "u.a. even chain";
	public final static String DESCRIPTION_UNASSIGNED_ODD_CHAIN = "u.a. odd chain";
	public final static String DESCRIPTION_UNASSIGNED_EVEN_CHAIN_PARTNER_KNOWN = "u.a. even chain; partner chain assigned / saturated";
	public final static String DESCRIPTION_UNASSIGNED_ODD_CHAIN_PARTNER_KNOWN = "u.a. odd chain; partner chain assigned / saturated";
	public final static String DESCRIPTION_UNASSIGNED_EVEN_CHAIN_PARTNER_UNASSIGNED = "u.a. even chain; partner chain unassigned";
	public final static String DESCRIPTION_UNASSIGNED_ODD_CHAIN_PARTNER_UNASSIGNED = "u.a. odd chain; partner chain unassigned";
	
	
	/**
	 * experiment, class, totalFA
	 * This denotes the total fatty acid content (omega inspecific), for normalization.
	 */
	LinkedHashMap<String,LinkedHashMap<String,Double>> totalFAPerClass_;
	/**
	 * experiment, description (saturated, unassigned, assigned), totalFA
	 * this is for computing the total content of FA of a description
	 */
	LinkedHashMap<String,LinkedHashMap<String,Double>> totalFAPerDescription_;
	/**
	 * experiment, class, omegapos, totalFA
	 * this is for the total fa increase validation
	 */
	LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<Integer,Double>>> totalOmegaPerClass_;
	
	
	
	/**
	 * experiment, class, faName, totalFA
	 * this is for the fa vs lClass heatmap
	 */
	LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>> totalOmegaFAPerClass_;
	/**
	 * experiment, class, faName, partner_faName, totalFA
	 * this is for the fa vs fa heatmap
	 */
	LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>>> faContentToPartnerContent_;
	
	
	
	
	/**
	 * experiment, std
	 * This denotes the total fatty acid content standard deviation (omega inspecific), for normalization.
	 */
	LinkedHashMap<String,Double> totalFASTD_;
	LinkedHashMap<String,Integer> numContributorsTotal_;
	/**
	 * experiment, std
	 * this is for the total fa increase validation STD
	 */
	LinkedHashMap<String,LinkedHashMap<Integer,Double>> omegaSTD_;
	LinkedHashMap<String,LinkedHashMap<Integer,Integer>> numContributorsOmega_;
	
	
	
	public OmegaCollector()
	{
		totalFAPerClass_ = new LinkedHashMap<String,LinkedHashMap<String,Double>>();
		totalFAPerDescription_ = new LinkedHashMap<String,LinkedHashMap<String,Double>>();
		totalOmegaPerClass_ = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<Integer,Double>>>();
		totalOmegaFAPerClass_ = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>>();
		faContentToPartnerContent_ = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>>>();
		
		totalFASTD_ = new LinkedHashMap<String,Double>();
		numContributorsTotal_ = new LinkedHashMap<String,Integer>();
		omegaSTD_ = new LinkedHashMap<String,LinkedHashMap<Integer,Double>>();
		numContributorsOmega_ = new LinkedHashMap<String,LinkedHashMap<Integer,Integer>>();
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
	
	public void addToTotalFAPerDescription(String exp, String description, Double amount)
	{
		if (!totalFAPerDescription_.containsKey(exp))
		{
			totalFAPerDescription_.put(exp, new LinkedHashMap<String,Double>());
		}
		if (!totalFAPerDescription_.get(exp).containsKey(description))
		{
			totalFAPerDescription_.get(exp).put(description, amount);
		}
		else
		{
			double before = totalFAPerDescription_.get(exp).get(description);
			totalFAPerDescription_.get(exp).put(description, before+amount);
		}
	}
	
	public void addToTotalOmegaPerClass(String exp, String lClass, Integer omegaPos, Double amount)
	{
		if (!totalOmegaPerClass_.containsKey(exp))
		{
			totalOmegaPerClass_.put(exp, new LinkedHashMap<String,LinkedHashMap<Integer,Double>>());
		}
		if (!totalOmegaPerClass_.get(exp).containsKey(lClass))
		{
			totalOmegaPerClass_.get(exp).put(lClass, new LinkedHashMap<Integer,Double>());
		}
		if (!totalOmegaPerClass_.get(exp).get(lClass).containsKey(omegaPos))
		{
			totalOmegaPerClass_.get(exp).get(lClass).put(omegaPos,amount);
		}
		else
		{
			double before = totalOmegaPerClass_.get(exp).get(lClass).get(omegaPos);
			totalOmegaPerClass_.get(exp).get(lClass).put(omegaPos, before+amount);
		}
	}
	
	public void addToTotalOmegaFAPerClass(String exp, String lClass, String faName, Double amount)
	{
		if (!totalOmegaFAPerClass_.containsKey(exp))
		{
			totalOmegaFAPerClass_.put(exp, new LinkedHashMap<String,LinkedHashMap<String,Double>>());
		}
		if (!totalOmegaFAPerClass_.get(exp).containsKey(lClass))
		{
			totalOmegaFAPerClass_.get(exp).put(lClass, new LinkedHashMap<String,Double>());
		}
		if (!totalOmegaFAPerClass_.get(exp).get(lClass).containsKey(faName))
		{
			totalOmegaFAPerClass_.get(exp).get(lClass).put(faName,amount);
		}
		else
		{
			double before = totalOmegaFAPerClass_.get(exp).get(lClass).get(faName);
			totalOmegaFAPerClass_.get(exp).get(lClass).put(faName, before+amount);
		}
	}
	
	public void addToFAContentToPartnerContent(String exp, String lClass, String vo, String partnerVO, Double amount)
	{
		if (!faContentToPartnerContent_.containsKey(exp))
		{
			faContentToPartnerContent_.put(exp, new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>>());
		}
		if (!faContentToPartnerContent_.get(exp).containsKey(lClass))
		{
			faContentToPartnerContent_.get(exp).put(lClass, new LinkedHashMap<String,LinkedHashMap<String,Double>>());
		}
		
		//just add them to both so it's symmetrical... add to first
		if (!faContentToPartnerContent_.get(exp).get(lClass).containsKey(vo))
		{
			faContentToPartnerContent_.get(exp).get(lClass).put(vo, new LinkedHashMap<String,Double>());
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
		
		//just add them to both so it's symmetrical... add to second
		if (!faContentToPartnerContent_.get(exp).get(lClass).containsKey(partnerVO))
		{
			faContentToPartnerContent_.get(exp).get(lClass).put(partnerVO, new LinkedHashMap<String,Double>());
		}
		if (!faContentToPartnerContent_.get(exp).get(lClass).get(partnerVO).containsKey(vo))
		{
			faContentToPartnerContent_.get(exp).get(lClass).get(partnerVO).put(vo, amount);
		}
		else
		{
			double before = faContentToPartnerContent_.get(exp).get(lClass).get(partnerVO).get(vo);
			faContentToPartnerContent_.get(exp).get(lClass).get(partnerVO).put(vo, before+amount);
		}
	}
	
	public Vector<String> getExperimentNames()
	{
		return new Vector<String>(totalFAPerClass_.keySet());
	}
	
	/**
	 * 
	 * @param experiment
	 * @return
	 */
	public Double getTotalFA(String experiment)
	{
		LinkedHashMap<String,Double> totalFAPerClass = totalFAPerClass_.get(experiment);
		Double totalFA = 0.0;
		for (String lClass : totalFAPerClass.keySet())
		{
			totalFA += totalFAPerClass.get(lClass);
		}
		return totalFA;
	}
	
	public LinkedHashMap<String,Double> getTotalFAPerDescription(String experiment)
	{
		return totalFAPerDescription_.get(experiment);
	}
	
	/**
	 * 
	 * @param experiment
	 * @return
	 */
	public LinkedHashMap<Integer,Double> getNormTotalOmegaFA(String experiment)
	{
		LinkedHashMap<String,LinkedHashMap<Integer,Double>> totalOmegaFAPerClass = totalOmegaPerClass_.get(experiment);
		LinkedHashMap<Integer,Double> normTotalOmegaFA = new LinkedHashMap<Integer,Double>();
		Double totalFA = getTotalFA(experiment);
		//summing it all up
		for (String lClass : totalOmegaFAPerClass.keySet())
		{
			LinkedHashMap<Integer,Double> totalOmegaFA = totalOmegaFAPerClass.get(lClass);
			for (Integer omegaPos : totalOmegaFA.keySet())
			{
				if (!normTotalOmegaFA.containsKey(omegaPos))
				{
					normTotalOmegaFA.put(omegaPos, 0.0);
				}
				Double before = normTotalOmegaFA.get(omegaPos);
				normTotalOmegaFA.put(omegaPos, before+totalOmegaFA.get(omegaPos));
			}
		}
		//normalizing
		for (Integer omegaPos : normTotalOmegaFA.keySet())
		{
			Double before = normTotalOmegaFA.get(omegaPos);
			normTotalOmegaFA.put(omegaPos, before/totalFA);
		}
		return normTotalOmegaFA;
	}
	
	public LinkedHashMap<String,LinkedHashMap<String,Double>> getTotalOmegaFAPerClass(String experiment)
	{
		return totalOmegaFAPerClass_.get(experiment);
	}
	
	public LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>> getFAContentToPartnerContent(String experiment)
	{
		return faContentToPartnerContent_.get(experiment);
	}
	
	
	
	
	
	//Starting here we got methods for the computation of standard deviations. This might not become relevant... maybe for total omega FA.
	public void addToTotalFAPerClassSTD(String exp, Double amount)
	{
		if (amount > 0)
		{
			if (!totalFASTD_.containsKey(exp))
			{
				totalFASTD_.put(exp, 0.0);
			}
			double before = totalFASTD_.get(exp);
			totalFASTD_.put(exp, before+Math.pow(amount, 2));
		
			if (!numContributorsTotal_.containsKey(exp))
			{
				numContributorsTotal_.put(exp, 0);
			}
			int prior = numContributorsTotal_.get(exp);
			numContributorsTotal_.put(exp, ++prior);
		}
	}
	
	public void addToTotalOmegaPerClassSTD(String exp, Integer omega, Double amount)
	{
		if (amount > 0)
		{
			if (!omegaSTD_.containsKey(exp))
			{
				omegaSTD_.put(exp, new LinkedHashMap<Integer,Double>());
			}
			if (!omegaSTD_.get(exp).containsKey(omega))
			{
				omegaSTD_.get(exp).put(omega, 0.0);
			}
			
			double before = omegaSTD_.get(exp).get(omega);
			omegaSTD_.get(exp).put(omega, before+Math.pow(amount, 2));
			
			if (!numContributorsOmega_.containsKey(exp))
			{
				numContributorsOmega_.put(exp, new LinkedHashMap<Integer,Integer>());
			}
			if (!numContributorsOmega_.get(exp).containsKey(omega))
			{
				numContributorsOmega_.get(exp).put(omega, 0);
			}
			
			int prior = numContributorsOmega_.get(exp).get(omega);
			numContributorsOmega_.get(exp).put(omega, ++prior);
		}
	}
	
	public Double getTotalSTD(String exp)
	{
		Double stdPow = totalFASTD_.get(exp);
		Integer numContributors = numContributorsTotal_.get(exp);
		return Math.sqrt(stdPow/numContributors);
	}
	
	public Double getOmegaSTD(String exp, Integer omega)
	{
		Double stdPow = omegaSTD_.get(exp).get(omega);
		Integer numContributors = numContributorsOmega_.get(exp).get(omega);
		return Math.sqrt(stdPow/numContributors);
	}
	
}
