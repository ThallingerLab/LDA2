/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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
package at.tugraz.genome.lda.target.export;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import javafx.util.Pair;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class GradientAdjuster
{
	private ArrayList<Pair<Double,Double>> dataPoints_;
	private final static Integer NO_TIE = 1;
	
	protected GradientAdjuster(ArrayList<Pair<Double,Double>> dataPoints)
	{
		Collections.sort(dataPoints, new SortByKey());
		this.dataPoints_ = dataPoints;
	}
	
	protected Double getGradientAdjustedValue(Double isotopeEffect, Double preciseRT)
	{
		Double adjustedTIE = isotopeEffect;
		for (int i=0;i<dataPoints_.size()-1;i++)
		{
			Pair<Double,Double> dataPoint = dataPoints_.get(i);
			Pair<Double,Double> dataPointNext = dataPoints_.get(i+1);
			if (preciseRT > dataPoint.getKey() && preciseRT < dataPointNext.getKey())
			{
				double factor = getAdjustmentFactorForPoint(dataPoint, dataPointNext, preciseRT);
				adjustedTIE = (isotopeEffect - NO_TIE)/factor + NO_TIE;
			}
		}
		return adjustedTIE * preciseRT;
	}
	
	private Double getAdjustmentFactorForPoint(Pair<Double,Double> dataPoint, Pair<Double,Double> dataPointNext, Double preciseRT)
	{
		double m = (dataPointNext.getValue() - dataPoint.getValue()) / (dataPointNext.getKey() - dataPoint.getKey());
		double b = dataPoint.getValue()-m*dataPoint.getKey();
		return m*preciseRT+b;
	}
  
	private class SortByKey implements Comparator<Pair<Double,Double>>
	{
		@Override
		public int compare(Pair<Double,Double> arg0, Pair<Double,Double> arg1)
		{
			return Double.compare(arg0.getKey(), arg1.getKey());
		}	
	}
}
