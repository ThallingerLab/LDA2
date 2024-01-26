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

package at.tugraz.genome.lda.target.experiment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.target.IsotopeLabelVO;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class IsotopeEffectRegression 
{
	private Vector<MatchedPartnerVO> matchedIsotopologues_;
	private PolynomialSplineFunction function_;
	private final int MAX_DEUTERIUM_ALLOWED;
	
	protected IsotopeEffectRegression(Vector<MatchedPartnerVO> matchedIsotopologues, Vector<IsotopeLabelVO> labels)
	{
		this.matchedIsotopologues_ = matchedIsotopologues;
		MAX_DEUTERIUM_ALLOWED = findMaximumDeuteriumAllowed(labels);
		initRegression();
//		super.setCoefficients(initRegression());
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
	
	protected void initRegression()
	{
		SplineInterpolator interpolator = new SplineInterpolator();
		Hashtable<Integer,ArrayList<Double>> dataPoints = new Hashtable<Integer,ArrayList<Double>>();
		for (MatchedPartnerVO matched : matchedIsotopologues_)
		{
			if (!matched.isUseForCalibration()) continue;
			
			LipidomicsMSnSet standard = matched.getStandard();
			double standardRT = standard.getPreciseRT();
			LipidomicsMSnSet isotopologue = matched.getIsotopologue();
			double isotopologueRT = isotopologue.getPreciseRT();
	  	
	  	try
			{
				Hashtable<String,Integer> chemicalFormula = StaticUtils.categorizeFormula(isotopologue.getChemicalFormula());
				int numberDeuterium = chemicalFormula.get("D");
				double totalIsotopeEffect = standardRT / isotopologueRT;
				if (!dataPoints.containsKey(numberDeuterium))
				{
					dataPoints.put(numberDeuterium, new ArrayList<Double>());
				}
				dataPoints.get(numberDeuterium).add(totalIsotopeEffect);
			}
			catch (ChemicalFormulaException ex)
			{
				ex.printStackTrace();
			}
		}
		
		ArrayList<Integer> numD = new ArrayList<Integer>();
  	numD.addAll(dataPoints.keySet());
  	Collections.sort(numD);
  	
  	double[] xValues = new double[numD.size()+1];
		double[] yValues = new double[numD.size()+1];
  	
		xValues[0] = 0;
		yValues[0] = 1;
		
  	for (int i=1; i<numD.size()+1; i++)
  	{
  		xValues[i] = numD.get(i-1);
  		Double yValue = 0.0;
  		for (Double value : dataPoints.get(numD.get(i-1)))
  		{
  			yValue += value;
  		}
  		yValue = yValue / dataPoints.get(numD.get(i-1)).size();
  		yValues[i] = yValue;
  	}
		
		if (dataPoints.size()>1) //less than two data points cannot be interpolated
		{
			this.function_ = interpolator.interpolate(xValues, yValues);
		}
		else
		{
			this.function_ = null;
		}
	}
	
	/**
	 * For plotting purposes, therefore a double for the number of deuterium is used.
	 * @param numberDeuterium
	 * @return
	 */
	public double getIsotopeEffect(double numberDeuterium)
	{
		return function_.value(numberDeuterium);
//		return super.getAdjustFactor(numberDeuterium);
	}
	
	public double getIsotopeEffect(int numberDeuterium)
	{
		return function_.value(numberDeuterium);
//		return super.getAdjustFactor(numberDeuterium);
	}
	
	public double getRTofUnlabeledSpecies(int numberDeuterium, double retentionTime)
	{
		return getIsotopeEffect(numberDeuterium) * retentionTime;
//		double isotopeEffect = super.getAdjustFactor(numberDeuterium);
//		return isotopeEffect*retentionTime;
	}
	
	public int getMaxNumDeuteriumAllowed() {
		return MAX_DEUTERIUM_ALLOWED;
	}
	
}
