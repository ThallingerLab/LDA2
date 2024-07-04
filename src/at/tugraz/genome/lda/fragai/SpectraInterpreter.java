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

public class SpectraInterpreter
{
	private ArrayList<SpectrumContainer> spectra_;
	/** the mz tolerance in ppm */
	private static final Integer MZ_TOLERANCE = 10;
	
	public SpectraInterpreter(ArrayList<SpectrumContainer> spectra)
	{
		this.spectra_ = spectra;
	}
	
	/**
	 * Combines spectra based on lipid class and adduct. Any combined spectra with enough data to compute merged spectra are returned for export.
	 * @return
	 */
	public ArrayList<CombinedSpectrumContainer> interpretSpectra()
	{
		Hashtable<String, ArrayList<SpectrumContainer>> spectraForLipid = groupSpectraContainer();
		ArrayList<CombinedSpectrumContainer> combinedSpectra = new ArrayList<CombinedSpectrumContainer>();
		for (String lipid : spectraForLipid.keySet())
		{
			CombinedSpectrumContainer combined = new CombinedSpectrumContainer(lipid.substring(0, lipid.indexOf("_")), lipid.substring(lipid.indexOf("_")+1, lipid.length()), MZ_TOLERANCE);
			for (SpectrumContainer container : spectraForLipid.get(lipid))
			{
				combined.addContainerCopy(container);
			}
			if (combined.isViableData())
			{
				combinedSpectra.add(combined);
			}
		}
		return combinedSpectra;
	}
	
	private String getSpectrumID(SpectrumContainer spectrum)
	{
		return spectrum.getEntry().getLipidClass()+"_"+spectrum.getAdduct().getAdductName();
	}
	
	private Hashtable<String, ArrayList<SpectrumContainer>> groupSpectraContainer()
	{
		Hashtable<String, ArrayList<SpectrumContainer>> spectraForLipid = new Hashtable<String, ArrayList<SpectrumContainer>>();
		for (SpectrumContainer spectrum : spectra_)
		{
			String id = getSpectrumID(spectrum);
			if (!spectraForLipid.containsKey(id))
			{
				spectraForLipid.put(id, new ArrayList<SpectrumContainer>());
			}
			spectraForLipid.get(id).add(spectrum);
		}
		return spectraForLipid;
	}
	
}
