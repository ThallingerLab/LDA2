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

package at.tugraz.genome.lda.target.calibration;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.quantification.LipidParameterSet;


/**
 * Contains identifications of the same lipid species, 
 * within a retention time grouping parameter across several input files.
 * 
 * @author Leonida M. Lamp
 *
 */
public class IdentificationVO
{
	private ArrayList<File> files_ = new ArrayList<File>();
	private String lipidClass_;
	private String lipidSpecies_;
	private boolean isMSnAvailable_ = false;
	private ArrayList<Double> retentionTimes_ = new ArrayList<Double>();
	private ArrayList<Double> areas_ = new ArrayList<Double>();
  //TODO: this should be a setting
	private double groupingParameter_ = 0.05;
	
	/**
	 * 
	 * @param file					the file this analyte stems from for sanity check purposes
	 * @param lipidClass
	 * @param param
	 */
	protected IdentificationVO(File file, String lipidClass, LipidParameterSet param)
	{
		this.files_.add(file);
		this.lipidClass_ = lipidClass;
		this.lipidSpecies_ = param.getNameStringWithoutRt();
		this.isMSnAvailable_ = param instanceof LipidomicsMSnSet;
		this.retentionTimes_.add(param.getPreciseRT());
		this.areas_.add(new Double(param.getArea()));
	}
	
	/**
	 * Adds an analyte to the object if its retention time is within the grouping parameter
	 * Make sure the lipid class and lipid species match!
	 * @param file			the file this analyte stems from for sanity check purposes
	 * @param param
	 * @return
	 */
	protected boolean addParam(File file, LipidParameterSet param)
	{
		Double paramRT = param.getPreciseRT();
		for (int i=0; i<retentionTimes_.size(); i++)
		{
			Double rT = retentionTimes_.get(i);
			if (paramRT >= (rT - groupingParameter_) &&
					paramRT <= (rT + groupingParameter_))
			{
				this.files_.add(file);
				this.isMSnAvailable_ = param instanceof LipidomicsMSnSet ? true : isMSnAvailable_;
				this.retentionTimes_.add(param.getPreciseRT());
	  		this.areas_.add(new Double(param.getArea()));
	  		return true;
			}
		}
		return false;
	}
	
	Double getAverageRT()
	{
		Double sum = 0.0;
		for (Double rT : retentionTimes_)
		{
			sum += rT;
		}
		return sum/retentionTimes_.size();
	}
	
	Double getAverageArea()
	{
		Double sum = 0.0;
		for (Double area : areas_)
		{
			sum += area;
		}
		return sum/areas_.size();
	}
	
	boolean isMSnAvailable()
	{
		return isMSnAvailable_;
	}
	
	/**
	 * If an identification was contributed to more than once by the same result file, 
	 * it cannot be used to match identifications. If this happens in too many cases, then
	 * either the grouping parameter needs to be lowered or the experiment is not suitable.
	 * @return 		true if the identification was contributed to only once per file
	 */
	protected boolean isIdentificationUnambiguous()
	{
		Set<File> uniqueFiles = new HashSet<File>(files_);
		return uniqueFiles.size() == files_.size();
	}
	
	/**
	 * Returns the number of files this identification was found in.
	 * @return
	 */
	int getIdentificationCount()
	{
		Set<File> uniqueFiles = new HashSet<File>(files_);
		return uniqueFiles.size();
	}
	
	protected String getLipidClass()
	{
		return this.lipidClass_;
	}
	
	protected String getLipidSpecies()
	{
		return this.lipidSpecies_;
	}

	@Override
	public String toString()
	{
		return "IdentificationVO [files_=" + files_ + ", lipidClass_="
				+ lipidClass_ + ", lipidSpecies_=" + lipidSpecies_
				+ ", isMSnAvailable_=" + isMSnAvailable_ + ", retentionTimes_="
				+ retentionTimes_ + ", areas_=" + areas_ + "]";
	}
}
