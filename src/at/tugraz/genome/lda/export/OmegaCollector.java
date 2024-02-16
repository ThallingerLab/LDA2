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
	public final static String DESCRIPTION_UNASSIGNED_SUGGESTED = "u.a. suggested";
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
	 * experiment, class, faName at sn1 position, faName at sn2 position, totalFA
	 * this is for the sn1 vs sn2 heatmap
	 */
	LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>>> sn1ContentToSn2Content_;
	
	
	
	
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
		sn1ContentToSn2Content_ = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>>>();
		
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
	
	public void addToSn1ContentToSn2Content(String exp, String lClass, String sn1, String sn2, Double amount)
	{
		if (!sn1ContentToSn2Content_.containsKey(exp))
		{
			sn1ContentToSn2Content_.put(exp, new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>>());
		}
		if (!sn1ContentToSn2Content_.get(exp).containsKey(lClass))
		{
			sn1ContentToSn2Content_.get(exp).put(lClass, new LinkedHashMap<String,LinkedHashMap<String,Double>>());
		}
		
		if (!sn1ContentToSn2Content_.get(exp).get(lClass).containsKey(sn1))
		{
			sn1ContentToSn2Content_.get(exp).get(lClass).put(sn1, new LinkedHashMap<String,Double>());
		}
		if (!sn1ContentToSn2Content_.get(exp).get(lClass).get(sn1).containsKey(sn2))
		{
			sn1ContentToSn2Content_.get(exp).get(lClass).get(sn1).put(sn2, amount);
		}
		else
		{
			double before = sn1ContentToSn2Content_.get(exp).get(lClass).get(sn1).get(sn2);
			sn1ContentToSn2Content_.get(exp).get(lClass).get(sn1).put(sn2, before+amount);
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
	
	public LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>> getSn1ContentToSn2Content(String experiment)
	{
		return sn1ContentToSn2Content_.get(experiment);
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
		try
		{
			Double stdPow = omegaSTD_.get(exp).get(omega);
			Integer numContributors = numContributorsOmega_.get(exp).get(omega);
			return Math.sqrt(stdPow/numContributors);
		}
		catch (Exception ex)
		{
			System.out.println("Likely only one replicate, so no standard deviation...");
		}
		return 0.0;
	}
	
}
