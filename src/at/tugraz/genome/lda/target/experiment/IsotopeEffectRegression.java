package at.tugraz.genome.lda.target.experiment;

import java.util.Hashtable;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.target.AbstractRegression;
import at.tugraz.genome.lda.target.IsotopeLabelVO;
import at.tugraz.genome.lda.utils.StaticUtils;

public class IsotopeEffectRegression extends AbstractRegression
{
	private Vector<MatchedPartnerVO> matchedIsotopologues_;
	private final int MAX_DEUTERIUM_ALLOWED;
//	private static int MAX_DEUTERIUM_ALLOWED = 19; //TODO: this should be the highest label represented by an authentic standard, user input dependent!
	
	protected IsotopeEffectRegression(Vector<MatchedPartnerVO> matchedIsotopologues, Vector<IsotopeLabelVO> labels)
	{
		this.matchedIsotopologues_ = matchedIsotopologues;
		MAX_DEUTERIUM_ALLOWED = findMaximumDeuteriumAllowed(labels);
		super.setCoefficients(initRegression());
	}
	
	private int findMaximumDeuteriumAllowed(Vector<IsotopeLabelVO> labels)
	{
		int maxD = 0;
		for (IsotopeLabelVO label : labels)
		{
			if (label.getLabelElements().containsKey("D") && label.getLabelElements().get("D")>maxD)
			{
				maxD = label.getLabelElements().get("D");
			}
		}
		return maxD;
	}
	
	protected double[] initRegression()
	{
		double[] coefficients = new double[3];
		PolynomialCurveFitter fitter = PolynomialCurveFitter.create(2);
		WeightedObservedPoints observations = new WeightedObservedPoints();
		observations.add(10000, 0, 1.0); //with 0 isotopes the total isotope effect must be 1.0, enforcing this with a big weight
		
		for (MatchedPartnerVO matched : matchedIsotopologues_)
		{
			if (!matched.isUseForCalibration()) continue;
			
			LipidomicsMSnSet standard = matched.getStandard();
			double standardRT = standard.getPreciseRT();
			//just printing out some stuff, remove later
			System.out.println("standard: "+standard.getOmegaInformation().get(0).getDoubleBondPositionsHumanReadable()+" RT: "+standard.getPreciseRT());
			LipidomicsMSnSet isotopologue = matched.getIsotopologue();
			double isotopologueRT = isotopologue.getPreciseRT();
			//just printing out some stuff, remove later
	  	Set<String> names = isotopologue.getHumanReadableNameSet();
	  	for (String name : names)
	  	{
	  		System.out.println("partner: "+name+" RT: "+isotopologueRT);
	  	}
	  	
	  	try
			{
				Hashtable<String,Integer> chemicalFormula = StaticUtils.categorizeFormula(isotopologue.getChemicalFormula());
				int numberDeuterium = chemicalFormula.get("D");
				double totalIsotopeEffect = standardRT / isotopologueRT;
				observations.add(numberDeuterium, totalIsotopeEffect);
			}
			catch (ChemicalFormulaException ex)
			{
				ex.printStackTrace();
			}
	  	
	  	coefficients = fitter.fit(observations.toList());
			for (int i=0;i<coefficients.length-1;i++) {
				System.out.println(coefficients[i]);
			}
		}
		return coefficients;
	}
		
//		for (String lipidClass : matchedIsotopologues_.keySet())
//		{
//			MatchedPartnersVO matched = matchedIsotopologues_.get(lipidClass);
//			Vector<LipidParameterSet> standards = matched.getStandards();
//			if (standards.size()>1) //to fit a 
//			{
//				for (LipidParameterSet standard : standards) 
//				{
//					double standardRT = standard.getPreciseRT();
//					Vector<LipidParameterSet> partners = matched.getIsotopologuePartners(standard);
//					
//					//just printing out some stuff, remove later
//					System.out.println("standard: "+standard.getOmegaInformation().get(0).getDoubleBondPositionsHumanReadable()+" RT: "+standard.getPreciseRT());
//					
//					//TODO: partners size requirement is only there to avoid an issue with a wrong sn assignment of one standard, fix in the future
//					if (!partners.isEmpty() && partners.size()>1) 
//					{
//						for (LipidParameterSet partner : partners)
//						{
//							try
//							{
//								double partnerRT = partner.getPreciseRT();
//								Hashtable<String,Integer> chemicalFormula = StaticUtils.categorizeFormula(partner.getChemicalFormula());
//								if (chemicalFormula.containsKey("D"))
//								{
//									int numberDeuterium = chemicalFormula.get("D");
//									double totalIsotopeEffect = standardRT / partnerRT;
////									System.out.println(String.format("TIE: %s, standarRT: %s, deuteratedRT: %s",totalIsotopeEffect, standardRT, partnerRT));
//									observations.add(numberDeuterium, totalIsotopeEffect);								
//								}
//								else 
//								{
//									throw new ChemicalFormulaException(
//										String.format("The compound with the chemical formula '%s' should not have entered this method, as it does not contain any element 'D'!", 
//												partner.getChemicalFormula()));
//								}
//							}
//							catch (ChemicalFormulaException ex)
//							{
//								ex.printStackTrace();
//							}
//							
//							//just printing out some stuff, remove later
//							LipidomicsMSnSet partnermsn = (LipidomicsMSnSet)partner;
//					  	Set<String> names = partnermsn.getHumanReadableNameSet();
//					  	for (String name : names)
//					  	{
//					  		System.out.println("partner: "+name+" RT: "+partner.getPreciseRT());
//					  	}
//						}
//					}
//				}
//				
//				coefficients = fitter.fit(observations.toList());
//				for (int i=0;i<coefficients.length-1;i++) {
//					System.out.println(coefficients[i]);
//				}
//				
//			}
//		}
//		return coefficients;
//	}
	
	/**
	 * For plotting purposes, therefore a double for the number of deuterium is used.
	 * @param numberDeuterium
	 * @return
	 */
	public double getIsotopeEffect(double numberDeuterium)
	{
		return super.getAdjustFactor(numberDeuterium);
	}
	
	public double getIsotopeEffect(int numberDeuterium)
	{
		return super.getAdjustFactor(numberDeuterium);
	}
	
	public double getRTofUnlabeledSpecies(int numberDeuterium, double retentionTime)
	{
		double isotopeEffect = super.getAdjustFactor(numberDeuterium);
		return isotopeEffect*retentionTime;
	}
	
	public int getMaxNumDeuteriumAllowed() {
		return MAX_DEUTERIUM_ALLOWED;
	}
	
}
