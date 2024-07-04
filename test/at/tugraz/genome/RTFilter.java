/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2024 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome;

import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.export.QuantificationResultExporter;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class RTFilter
{
	private static final String RESULT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\labeled_n3\\28_D5-18-3(n-3)_327_30min_negative_C.xlsx";
	private static final String LABEL = "C";
	
	public static void main(String[] args)
  {
    // MainFrame frame = new MainFrame(new TestClass(), 1024, 1024);
    filterRT(100.0, 110.0);
  }
	
	/**
	 * 
	 * @param lowerBound 		lowest allowed rt
	 * @param upperBound		highest allowed rt
	 */
	private static void filterRT(double lowerBound, double upperBound)
	{
		try
		{
			QuantificationResult quantRes = LDAResultReader.readResultFile(RESULT_PATH, new Hashtable<String,Boolean>());
			Hashtable<String,Vector<LipidParameterSet>> results = quantRes.getIdentifications();
			Hashtable<String,Vector<LipidParameterSet>> newResults = new Hashtable<String,Vector<LipidParameterSet>>();
			for (String lipidClass : results.keySet())
			{
				newResults.put(lipidClass, new Vector<LipidParameterSet>());
				Vector<LipidParameterSet> params = results.get(lipidClass);
				for (LipidParameterSet param : params)
				{
					if (!param.getNameStringWithoutRt().contains(LABEL))
					{
						newResults.get(lipidClass).add(param);
					}
					else if (param.getPreciseRT() >= lowerBound && param.getPreciseRT() <= upperBound)
					{
						newResults.get(lipidClass).add(param);
					}
				}
			}
			quantRes.setIdentifications(newResults);
			QuantificationResultExporter.writeResultsToExcel(RESULT_PATH, quantRes);
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}
}
