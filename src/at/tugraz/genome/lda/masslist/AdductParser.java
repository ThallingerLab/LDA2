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

package at.tugraz.genome.lda.masslist;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.vos.AdductVO;


/**
 * 
 * @author Leonida M. Lamp
 * 
 */
public class AdductParser
{
	public final static String ADDUCT_FOLDER = "./massListCreation/adducts";
	public final static String ADDUCT_SUFFIX = ".txt";
	
	public final static String PROPERTY_NAME = "name";
	public final static String PROPERTY_FORMULA = "formula";
	public final static String PROPERTY_CHARGE = "charge";
	
	public AdductParser(){}
	
	public ArrayList<AdductVO> parse() throws FileNotFoundException, IOException, ChemicalFormulaException
	{
		ArrayList<AdductVO> allDefinedAdducts = new ArrayList<AdductVO>();
		ArrayList<String> uniqueNames = new ArrayList<String>();
		
		File folder = new File(ADDUCT_FOLDER);
		if (!folder.exists())
		{
			throw new IOException(String.format("The adduct folder '%s' does not exist!", ADDUCT_FOLDER));
		}
		File[] fileCandidates = folder.listFiles();
		for (int i=0; i<fileCandidates.length;i++)
		{
			if (fileCandidates[i].getName().endsWith(ADDUCT_SUFFIX))
			{
				try (FileInputStream in= new FileInputStream(fileCandidates[i]))
				{
					Properties properties = new Properties();
					properties.load(in);
					String name = properties.getProperty(PROPERTY_NAME, null);
					String formula = properties.getProperty(PROPERTY_FORMULA, null);
					int charge = Integer.parseInt(properties.getProperty(PROPERTY_CHARGE, "0"));
					if (name != null && formula != null && charge != 0)
					{
						if (uniqueNames.contains(name))
							throw new IOException(String.format("Duplicate name definition '%s' detected in '%s'!", name, ADDUCT_FOLDER));
						else
							uniqueNames.add(name);
						
						AdductVO adduct = new AdductVO(name, formula, charge, fileCandidates[i].getName());
						allDefinedAdducts.add(adduct);
					}
					else
					{
						System.out.format("The adduct definition file '%s' does not adhere to the required format!", fileCandidates[i].getName());
					}
				}
			}
		}
		return allDefinedAdducts;
	}
	
}
