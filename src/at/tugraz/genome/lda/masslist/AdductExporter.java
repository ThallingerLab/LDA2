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
import at.tugraz.genome.lda.vos.AdductVO;

/**
 * 
 * @author Leonida M. Lamp
 * 
 */
public class AdductExporter
{
	private AdductVO toExport_;
	
	public AdductExporter(AdductVO toExport)
	{
		this.toExport_ = toExport;
	}
	
	public void export()
	{
		try (FileOutputStream out= new FileOutputStream(buildAdductPath(toExport_.getFileName()));)
		{
			out.write("## The name of the adduct.\n".getBytes());
			out.write(String.format("%s=%s\n", AdductParser.PROPERTY_NAME, toExport_.getAdductName()).getBytes());
			out.write("## The chemical formula of the adduct. Chemical formulas prefixed with a '-' and '+' are subtracted from and added to the total mass, respectively. Syntax e.g. +Na-H2\n".getBytes());
			out.write(String.format("%s=%s\n", AdductParser.PROPERTY_FORMULA, toExport_.getFormulaString()).getBytes());
			out.write("## The charge state of the molecule with this adduct.\n".getBytes());
			out.write(String.format("%s=%s\n", AdductParser.PROPERTY_CHARGE, toExport_.getCharge()).getBytes());
		}
		catch (IOException ex) 
		{
			new WarningMessage(new JFrame(), "Error", "The export of the adduct definition file failed. Error message: "+ex.getMessage());
		}
	}
	
	public static String buildAdductPath(String fileName)
	{
		return AdductParser.ADDUCT_FOLDER+"/"+fileName;
	}
	
}
