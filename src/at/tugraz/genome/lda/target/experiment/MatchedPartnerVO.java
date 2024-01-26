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

package at.tugraz.genome.lda.target.experiment;

import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.target.IsotopeLabelVO;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class MatchedPartnerVO
{
	private String lipidClass_;
	private LipidomicsMSnSet standard_;
	private LipidomicsMSnSet isotopologue_;
	private boolean useForCalibration_;
	private IsotopeLabelVO label_;
	
	protected MatchedPartnerVO(String lipidClass, LipidomicsMSnSet standard, LipidomicsMSnSet isotopologue, boolean useForCalibration, IsotopeLabelVO label)
	{
		this.lipidClass_ = lipidClass;
		this.standard_ = standard;
		this.isotopologue_ = isotopologue;
		this.useForCalibration_ = useForCalibration;
		this.label_ = label;
	}
	
	public String getLipidClass()
	{
		return lipidClass_;
	}

	public void setLipidClass(String lipidClass)
	{
		this.lipidClass_ = lipidClass;
	}

	public LipidomicsMSnSet getStandard()
	{
		return standard_;
	}

	public void setStandard(LipidomicsMSnSet standard)
	{
		this.standard_ = standard;
	}

	public LipidomicsMSnSet getIsotopologue()
	{
		return isotopologue_;
	}

	public void setIsotopologue(LipidomicsMSnSet isotopologue)
	{
		this.isotopologue_ = isotopologue;
	}

	public boolean isUseForCalibration()
	{
		return useForCalibration_;
	}

	public void setUseForCalibration(boolean useForCalibration)
	{
		this.useForCalibration_ = useForCalibration;
	}

	public IsotopeLabelVO getLabel()
	{
		return label_;
	}

	public void setLabel(IsotopeLabelVO label)
	{
		this.label_ = label;
	}
}
