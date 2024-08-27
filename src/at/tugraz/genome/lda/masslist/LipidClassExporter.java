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

import java.io.FileOutputStream;
import java.io.IOException;

import javax.swing.JFrame;

import at.tugraz.genome.lda.WarningMessage;

/**
 * 
 * @author Leonida M. Lamp
 * 
 */
public class LipidClassExporter
{
	private LipidClassVO toExport_;
	
	public LipidClassExporter(LipidClassVO toExport)
	{
		this.toExport_ = toExport;
	}
	
	public void export()
	{
		try (FileOutputStream out= new FileOutputStream(buildLipidClassPath(toExport_.getLipidClass()));)
		{
			out.write("## The parameter 'name' defines the name of this lipid class.\n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_NAME, toExport_.getLipidClass()).getBytes());
			out.write("## When set to 'true', the parameter 'adductInsensitiveRtFilter' forces the RT filter to be calculated based on all modifications\n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_RT_FILTER, toExport_.isAdductInsensitiveRtFilter()).getBytes());
			out.write("## When set to 'true', the parameter 'pickBestMatchBySpectrumCoverage' the best matches are picked from duplicates (same lipid class and scan number) by spectrum coverage.\n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_PICK_BEST, toExport_.isPickBestMatchBySpectrumCoverage()).getBytes());
			out.write("## Number of 'OH groups' present in the fatty acyl (FA) and long chain base (LCB) chains of this lipid class, \n".getBytes());
			out.write("## and 'OH_range_from' and 'OH_range_to' define the lower and upper limit of the number of 'OH groups', respectively, which the algorithm identifies compounds for. \n".getBytes());
			out.write("## Which number is chosen for 'OH_number' from the biologically relevant range for sphingolipids does not influence the algorithm. \n".getBytes());
			out.write("## ATTENTION: this parameter should only be used for sphingolipids. \n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_OH_RANGE_FROM, toExport_.getOhRangeFrom()).getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_OH_RANGE_TO, toExport_.getOhRangeTo()).getBytes());
			out.write("## Optional definition of the retention time (RT) range (the lower limit is defined with 'RT_range_from' and the upper limit with 'RT_range_to') that this lipid class elutes at. \n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_RT_RANGE_FROM, toExport_.getRtRangeFrom()).getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_RT_RANGE_TO, toExport_.getRtRangeTo()).getBytes());
			out.write("## The parameter 'adducts' allows for the definition of adducts relevant for this lipid class. \n".getBytes());
			out.write("## The given names must correspond to names (parameter 'name') of adducts defined in .txt files in the folder ./massListCreation/adducts \n".getBytes());
			out.write("## Multiple names are separated with a comma and no spaces, e.g. adductName1,adductName2,adductName3 \n".getBytes());
			StringBuilder builder = new StringBuilder();
			builder.append(toExport_.getAdducts().get(0).getAdductName());
			for (int i=1;i<toExport_.getAdducts().size();i++)
			{
				builder.append(LipidClassParser.ADDUCT_SEPARATOR);
				builder.append(toExport_.getAdducts().get(i).getAdductName());
			}
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_ADDUCTS, builder.toString()).getBytes());
			out.write("## The parameter 'headgroup_formula' defines the chemical formula of the headgroup (without chains) \n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_FORMULA, toExport_.getHeadGroupFormulaString()).getBytes());
			out.write("## The parameters 'min_chain_C' and 'max_chain_C' describe the minimum and maximum total number of C atoms in chains to be included in the mass list \n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_MIN_CHAIN_C, toExport_.getMinChainC()).getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_MAX_CHAIN_C, toExport_.getMaxChainC()).getBytes());
			out.write("## The parameters 'min_chain_DB' and 'max_chain_DB' describe the minimum and maximum total number of double bonds in chains to be considered in the mass list \n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_MIN_CHAIN_DB, toExport_.getMinChainDB()).getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_CLASS_MAX_CHAIN_DB, toExport_.getMaxChainDB()).getBytes());
			out.write("## The parameter 'number_of_FA_chains' allows for the definition of the number of FA chains present in this lipid class. \n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_NUMBER_FA_CHAINS, toExport_.getNumberOfFAChains()).getBytes());
			out.write("## The parameter 'number_of_LCB_chains' allows for the definition of the number of LCB chains present in this lipid class. \n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LIPID_NUMBER_LCB_CHAINS, toExport_.getNumberOfLCBChains()).getBytes());
			out.write("## The parameter 'FA_chain_list_name' allows for the definition of a FA chain list name (without suffix, must be present in the folder ./fattyAcids) to base the mass list on. \n".getBytes());
			out.write("## Only compounds that are possible given the entries in the provided chain list will be included in the resulting mass list. \n".getBytes());
			out.write("## If the chain list includes stable isotope labeled chains, those will be included in all possible combinations in the generated mass list file. \n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.FA_LIST_NAME, toExport_.getFAChainList()).getBytes());
			out.write("## The parameter 'LCB_chain_list_name' allows for the definition of a LCB chain list name (without suffix, must be present in the folder ./fattyAcids) to base the mass list on. \n".getBytes());
			out.write("## Only compounds that are possible given the entries in the provided chain list will be included in the resulting mass list. \n".getBytes());
			out.write("## If the chain list includes stable isotope labeled chains, those will be included in all possible combinations in the generated mass list file. \n".getBytes());
			out.write(String.format("%s=%s\n", LipidClassParser.LCB_LIST_NAME, toExport_.getLCBChainList()).getBytes());
		}
		catch (IOException ex) 
		{
			new WarningMessage(new JFrame(), "Error", "The export of the adduct definition file failed. Error message: "+ex.getMessage());
		}
	}
	
	public static String buildLipidClassPath(String lipidClass)
	{
		return LipidClassParser.LIPID_CLASS_FOLDER+"/"+lipidClass+LipidClassParser.LIPID_CLASS_SUFFIX;
	}
	
}
