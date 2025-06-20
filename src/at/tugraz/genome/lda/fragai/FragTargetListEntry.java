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

package at.tugraz.genome.lda.fragai;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Objects;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.vos.AdductVO;

public class FragTargetListEntry
{
	private String lipidClass_;
	private String species_;
	private Hashtable<String,Integer> sumFormula_;
	private ArrayList<AdductVO> adducts_;
	private Double retentionTime_;
	private Double tolerance_;
	
	
	public FragTargetListEntry(String lipidClass, String species,
			Hashtable<String,Integer> sumFormula, ArrayList<AdductVO> adducts,
			Double retentionTime, Double tolerance)
	{
		this.lipidClass_ = lipidClass;
		this.species_ = species;
		this.sumFormula_ = sumFormula;
		this.adducts_ = adducts;
		this.retentionTime_ = retentionTime;
		this.tolerance_ = tolerance;
	}


	public String getLipidClass()
	{
		return lipidClass_;
	}

	public String getSpecies()
	{
		return species_;
	}


	public Hashtable<String,Integer> getSumFormula()
	{
		return sumFormula_;
	}


	public ArrayList<AdductVO> getAdducts()
	{
		return adducts_;
	}
	
	/**
	 * Computes the theoretical precursor mz value considering the adduct and the charge
	 * @param adductName 		the name of the relevant adduct
	 * @return
	 */
	public Double computeTheoreticalPrecursorMZValue(String adductName)
	{
		for (AdductVO adduct : adducts_)
		{
			if (adduct.getAdductName().equalsIgnoreCase(adductName))
			{
				Double fullMz = 0.0;
				for (String element : getSumFormula().keySet())
				{
					fullMz += Settings.getElementParser().getElementDetails(element).getMonoMass()*getSumFormula().get(element);
				}
				Double adductMz = fullMz;
				for (String element : adduct.getFormula().keySet())
				{
					adductMz += Settings.getElementParser().getElementDetails(element).getMonoMass()*adduct.getFormula().get(element);
				}
				return Math.abs(adductMz/adduct.getCharge());
			}
		}
		return null;
	}

	public Double getRetentionTime()
	{
		return retentionTime_;
	}


	public Double getTolerance()
	{
		return tolerance_;
	}


	@Override
	public int hashCode()
	{
		return Objects.hash(adducts_, lipidClass_, species_, retentionTime_, sumFormula_,
				tolerance_);
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
		FragTargetListEntry other = (FragTargetListEntry) obj;
		return Objects.equals(adducts_, other.adducts_)
				&& Objects.equals(lipidClass_, other.lipidClass_)
				&& Objects.equals(species_, other.species_)
				&& Objects.equals(retentionTime_, other.retentionTime_)
				&& Objects.equals(sumFormula_, other.sumFormula_)
				&& Objects.equals(tolerance_, other.tolerance_);
	}
	
}
