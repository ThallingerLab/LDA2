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
public class LipidClassParser
{
	public final static String LIPID_CLASS_FOLDER = "./massListCreation/lipidClasses";
	public final static String LIPID_CLASS_SUFFIX = ".txt";
	public final static String ADDUCT_SEPARATOR = ",";
	public final static String LIPID_CLASS_NAME = "name";
	public final static String LIPID_CLASS_RT_FILTER = "adductInsensitiveRtFilter";
	public final static String LIPID_CLASS_PICK_BEST = "pickBestMatchBySpectrumCoverage";
	public final static String LIPID_CLASS_OH_RANGE_FROM = "OH_range_from";
	public final static String LIPID_CLASS_OH_RANGE_TO = "OH_range_to";
	public final static String LIPID_CLASS_RT_RANGE_FROM = "RT_range_from";
	public final static String LIPID_CLASS_RT_RANGE_TO = "RT_range_to";
	public final static String LIPID_CLASS_ADDUCTS = "adducts";
	public final static String LIPID_CLASS_FORMULA = "headgroup_formula";
	public final static String LIPID_CLASS_MIN_CHAIN_C = "min_chain_C";
	public final static String LIPID_CLASS_MAX_CHAIN_C = "max_chain_C";
	public final static String LIPID_CLASS_MIN_CHAIN_DB = "min_chain_DB";
	public final static String LIPID_CLASS_MAX_CHAIN_DB = "max_chain_DB";
	public final static String LIPID_NUMBER_FA_CHAINS = "number_of_FA_chains";
	public final static String LIPID_NUMBER_LCB_CHAINS = "number_of_LCB_chains";
	public final static String FA_LIST_NAME = "FA_chain_list_name";
	public final static String LCB_LIST_NAME = "LCB_chain_list_name";
	
	private ArrayList<AdductVO> allDefinedAdducts_;
	
	public LipidClassParser(ArrayList<AdductVO> allDefinedAdducts)
	{
		this.allDefinedAdducts_ = allDefinedAdducts;
	}
	
	public ArrayList<LipidClassVO> parse() throws IOException, ChemicalFormulaException
	{
		ArrayList<LipidClassVO> allDefinedLipidClasses = new ArrayList<LipidClassVO>();
		
		File folder = new File(LIPID_CLASS_FOLDER);
		if (!folder.exists())
		{
			throw new IOException(String.format("The adduct folder '%s' does not exist!", LIPID_CLASS_FOLDER));
		}
		File[] fileCandidates = folder.listFiles();
		for (int i=0; i<fileCandidates.length;i++)
		{
			if (fileCandidates[i].getName().endsWith(LIPID_CLASS_SUFFIX))
			{
				try (FileInputStream in= new FileInputStream(fileCandidates[i]))
				{
					Properties properties = new Properties();
					properties.load(in);
					String name = properties.getProperty(LIPID_CLASS_NAME, null);
					boolean adductInsensitiveRtFilter = Boolean.parseBoolean(properties.getProperty(LIPID_CLASS_RT_FILTER, null));
					boolean pickBestMatchBySpectrumCoverage = Boolean.parseBoolean(properties.getProperty(LIPID_CLASS_PICK_BEST, null));
					int ohRangeFrom = Integer.parseInt(properties.getProperty(LIPID_CLASS_OH_RANGE_FROM, null));
					int ohRangeTo = Integer.parseInt(properties.getProperty(LIPID_CLASS_OH_RANGE_TO, null));
					double rtRangeFrom = Double.parseDouble(properties.getProperty(LIPID_CLASS_RT_RANGE_FROM, null));
					double rtRangeTo = Double.parseDouble(properties.getProperty(LIPID_CLASS_RT_RANGE_FROM, null));
					String adducts = properties.getProperty(LIPID_CLASS_ADDUCTS, null);
					String headgroupFormula = properties.getProperty(LIPID_CLASS_FORMULA, null);
					int minChainC = Integer.parseInt(properties.getProperty(LIPID_CLASS_MIN_CHAIN_C, null));
					int maxChainC = Integer.parseInt(properties.getProperty(LIPID_CLASS_MAX_CHAIN_C, null));
					int minChainDB = Integer.parseInt(properties.getProperty(LIPID_CLASS_MIN_CHAIN_DB, null));
					int maxChainDB = Integer.parseInt(properties.getProperty(LIPID_CLASS_MAX_CHAIN_DB, null));
					int numberOfFAChains = Integer.parseInt(properties.getProperty(LIPID_NUMBER_FA_CHAINS, null));
					int numberOfLCBChains = Integer.parseInt(properties.getProperty(LIPID_NUMBER_LCB_CHAINS, null));
					String faChainList = properties.getProperty(FA_LIST_NAME, null);
					String lcbChainList = properties.getProperty(LCB_LIST_NAME, null);
					
					if (name != null && adducts != null && headgroupFormula != null && faChainList != null)
					{
						ArrayList<AdductVO> selectedAdducts = matchSelectedAdducts(adducts, fileCandidates[i].getName());
						
						LipidClassVO vo = new LipidClassVO(name,
								adductInsensitiveRtFilter, pickBestMatchBySpectrumCoverage, ohRangeFrom, ohRangeTo, 
								rtRangeFrom, rtRangeTo, selectedAdducts, headgroupFormula,
								minChainC, maxChainC, minChainDB, maxChainDB, numberOfFAChains, numberOfLCBChains, faChainList, lcbChainList);
						
						allDefinedLipidClasses.add(vo);
					}
					else
					{
						throw new IOException(String.format(
								"The lipid class definition file '%s' does not adhere to the required format!", 
								fileCandidates[i].getName()));
					}
				}
			}
		}	
		return allDefinedLipidClasses;
	}
	
	private ArrayList<AdductVO> matchSelectedAdducts(String adducts, String fileName) throws IOException
	{
		ArrayList<AdductVO> selectedAdducts = new ArrayList<AdductVO>();
		String[] adductNames = adducts.split(ADDUCT_SEPARATOR);
		for (int j=0;j<adductNames.length;j++)
		{
			boolean matched = false;
			for (AdductVO adduct : allDefinedAdducts_)
			{
				if (adduct.getAdductName().equalsIgnoreCase(adductNames[j]))
				{
					selectedAdducts.add(adduct);
					matched = true;
				}
			}
			if (!matched)
			{
				throw new IOException(String.format(
						"The provided adduct name '%s' for the lipid class definition file '%s' is not defined!", 
						adductNames[j], fileName));
			}
		}
		return selectedAdducts;
	}
	
}
