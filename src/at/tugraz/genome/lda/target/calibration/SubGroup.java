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

import java.util.ArrayList;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class SubGroup
{
	private String groupName_;
	private ArrayList<String> lipidClasses_;
	
	public SubGroup(String groupName, ArrayList<String> lipidClasses)
	{
		this.groupName_ = groupName;
		this.lipidClasses_ = lipidClasses;
	}

	public String getGroupName()
	{
		return groupName_;
	}

	public ArrayList<String> getLipidClasses()
	{
		return lipidClasses_;
	}
	
	public String getId()
	{
		StringBuilder lipidClassesBuilder = new StringBuilder();
		for (int i=0;i<lipidClasses_.size();i++)
		{
			lipidClassesBuilder.append(lipidClasses_.get(i));
			if (i+1<lipidClasses_.size())
			{
				lipidClassesBuilder.append(", ");
			}
		}
		return String.format("%s; lipid classes %s", groupName_, lipidClassesBuilder.toString()); 
	}
}
