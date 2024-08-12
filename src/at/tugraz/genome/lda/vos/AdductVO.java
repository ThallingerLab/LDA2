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

package at.tugraz.genome.lda.vos;

import java.util.Hashtable;
import java.util.Objects;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.utils.StaticUtils;

public class AdductVO
{
	private final static String ADDUCT_SEPARATOR = ";";
	private String adductName_;
	private String formulaString_;
	private Hashtable<String,Integer> formula_;
	private int charge_;
	private String fileName_;
	
	public AdductVO(String name, String formula, int charge, String fileName) throws ChemicalFormulaException
	{
		this.fileName_=fileName;
		this.adductName_ = name;
		this.formulaString_ = formula;
		this.formula_ = StaticUtils.categorizeFormula(formula, true);
		this.charge_ = charge;
	}
	
	public AdductVO(String composite) throws ChemicalFormulaException
	{
		String[] split = composite.split(ADDUCT_SEPARATOR);
		this.adductName_ = split[0];
		this.formulaString_ = split[1];
		this.formula_ = StaticUtils.categorizeFormula(split[1], true);
		this.charge_ = Integer.parseInt(split[2]);
	}

	public AdductVO(AdductVO other) throws ChemicalFormulaException
	{
		this(other.getAdductName(), other.getFormulaString(), other.getCharge(), other.getFileName());
	}

	public String getAdductName()
	{
		return adductName_;
	}
	
	public String getFormulaString()
	{
		return formulaString_;
	}
	
	public Hashtable<String,Integer> getFormula()
	{
		return this.formula_;
	}

	public int getCharge()
	{
		return charge_;
	}
	
	public String getFileName()
	{
		return fileName_;
	}

	public void setAdductName(String adductName)
	{
		this.adductName_ = adductName;
	}

	public void setFormulaString(String formulaString) throws ChemicalFormulaException
	{
		this.formulaString_ = formulaString;
		this.formula_ = null;
		this.formula_ = StaticUtils.categorizeFormula(formulaString, true);
	}

	public void setCharge(int charge)
	{
		this.charge_ = charge;
	}

	public void setFileName(String filePath)
	{
		this.fileName_ = filePath;
	}

	@Override
	public int hashCode()
	{
		return Objects.hash(adductName_, charge_, formula_, fileName_);
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AdductVO other = (AdductVO) obj;
		return Objects.equals(adductName_, other.adductName_)
				&& charge_ == other.charge_
				&& Objects.equals(formula_, other.formula_)
				&& Objects.equals(fileName_, other.fileName_);
	}
	
	
}
