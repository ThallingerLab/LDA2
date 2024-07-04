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

package at.tugraz.genome.lda.target.calibration;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class AnchorPointVO
{
	private String lipidClass_;
	private String lipidSpecies_;
	private Double xValue_;
	private Double yValue_;
	
	public AnchorPointVO(String lipidClass, String lipidSpecies, Double xValue, Double yValue)
	{
		super();
		this.lipidClass_ = lipidClass;
		this.lipidSpecies_ = lipidSpecies;
		this.xValue_ = xValue;
		this.yValue_ = yValue;
	}

	public String getLipidClass()
	{
		return lipidClass_;
	}

	public String getLipidSpecies()
	{
		return lipidSpecies_;
	}

	public Double getxValue()
	{
		return xValue_;
	}

	public Double getyValue()
	{
		return yValue_;
	}
	
	
	
}
