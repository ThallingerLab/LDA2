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

import java.util.ArrayList;
import java.util.Hashtable;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.AdductVO;

public class LipidClassVO
{
	private String lClass_;
	private boolean adductInsensitiveRtFilter_;
	private boolean pickBestMatchBySpectrumCoverage_;
	private int ohNumber_;
	private int ohRangeFrom_;
	private int ohRangeTo_;
	private double rtRangeFrom_;
	private double rtRangeTo_;
	private ArrayList<AdductVO> adducts_;
	private String headGroupFormulaString_;
	private Hashtable<String,Integer> headgroupFormula_;
	private int minChainC_;
	private int maxChainC_;
	private int minChainDB_;
	private int maxChainDB_;
	private int numberOfFAChains_;
	private int numberOfLCBChains_;
	private String faChainList_;
	private String lcbChainList_;
	
	
	public LipidClassVO(String lClass, boolean adductInsensitiveRtFilter,
			boolean pickBestMatchBySpectrumCoverage, int ohNumber, int ohRangeFrom,
			int ohRangeTo, double rtRangeFrom, double rtRangeTo, 
			ArrayList<AdductVO> adducts, String headgroupFormula, int minChainC,
			int maxChainC, int minChainDB, int maxChainDB, int numberOfFAChains,
			int numberOfLCBChains, String faChainList, String lcbChainList) throws ChemicalFormulaException
	{
		super();
		this.lClass_ = lClass;
		this.adductInsensitiveRtFilter_ = adductInsensitiveRtFilter;
		this.pickBestMatchBySpectrumCoverage_ = pickBestMatchBySpectrumCoverage;
		this.ohNumber_ = ohNumber;
		this.ohRangeFrom_ = ohRangeFrom;
		this.ohRangeTo_ = ohRangeTo;
		this.rtRangeFrom_ = rtRangeFrom;
		this.rtRangeTo_ = rtRangeTo;
		this.adducts_ = adducts;
		this.headGroupFormulaString_ = headgroupFormula;
		this.headgroupFormula_ = StaticUtils.categorizeFormula(headgroupFormula, true);
		this.minChainC_ = minChainC;
		this.maxChainC_ = maxChainC;
		this.minChainDB_ = minChainDB;
		this.maxChainDB_ = maxChainDB;
		this.numberOfFAChains_ = numberOfFAChains;
		this.numberOfLCBChains_ = numberOfLCBChains;
		this.faChainList_ = faChainList;
		this.lcbChainList_ = lcbChainList;
	}

	public LipidClassVO(LipidClassVO other) throws ChemicalFormulaException
	{
		this(other.getLipidClass(), other.isAdductInsensitiveRtFilter(), other.isPickBestMatchBySpectrumCoverage(), other.getOhNumber(), other.getOhRangeFrom(), other.getOhRangeTo(),
				other.getRtRangeFrom(), other.getRtRangeTo(), other.getAdducts(), other.getHeadGroupFormulaString(), other.getMinChainC(),
				other.getMaxChainC(), other.getMinChainDB(), other.getMaxChainDB(), other.getNumberOfFAChains(), other.getNumberOfLCBChains(), other.getFAChainList(), other.getLCBChainList());
	}

	public String getLipidClass()
	{
		return lClass_;
	}

	public boolean isAdductInsensitiveRtFilter()
	{
		return adductInsensitiveRtFilter_;
	}


	public boolean isPickBestMatchBySpectrumCoverage()
	{
		return pickBestMatchBySpectrumCoverage_;
	}


	public int getOhNumber()
	{
		return ohNumber_;
	}


	public int getOhRangeFrom()
	{
		return ohRangeFrom_;
	}


	public int getOhRangeTo()
	{
		return ohRangeTo_;
	}


	public double getRtRangeFrom()
	{
		return rtRangeFrom_;
	}


	public double getRtRangeTo()
	{
		return rtRangeTo_;
	}


	public ArrayList<AdductVO> getAdducts()
	{
		return adducts_;
	}

	public String getHeadGroupFormulaString()
	{
		return headGroupFormulaString_;
	}

	public Hashtable<String,Integer> getHeadgroupFormula()
	{
		return headgroupFormula_;
	}


	public int getMinChainC()
	{
		return minChainC_;
	}


	public int getMaxChainC()
	{
		return maxChainC_;
	}


	public int getMinChainDB()
	{
		return minChainDB_;
	}


	public int getMaxChainDB()
	{
		return maxChainDB_;
	}
	
	public int getNumberOfFAChains()
	{
		return numberOfFAChains_;
	}
	
	public int getNumberOfLCBChains()
	{
		return numberOfLCBChains_;
	}

	public String getFAChainList()
	{
		return faChainList_;
	}
	
	public String getLCBChainList()
	{
		return lcbChainList_;
	}

	public String getFAChainListPath()
	{
		return MassListCreatorPanel.CHAIN_LIST_FOLDER+faChainList_+MassListCreatorPanel.CHAIN_LIST_SUFFIX;
	}
	
	public String getLCBChainListPath()
	{
		return MassListCreatorPanel.CHAIN_LIST_FOLDER+lcbChainList_+MassListCreatorPanel.CHAIN_LIST_SUFFIX;
	}

	public void setLipidClass(String lClass)
	{
		this.lClass_ = lClass;
	}

	public void setAdductInsensitiveRtFilter(boolean adductInsensitiveRtFilter)
	{
		this.adductInsensitiveRtFilter_ = adductInsensitiveRtFilter;
	}

	public void setPickBestMatchBySpectrumCoverage(
			boolean pickBestMatchBySpectrumCoverage)
	{
		this.pickBestMatchBySpectrumCoverage_ = pickBestMatchBySpectrumCoverage;
	}

	public void setOhNumber(int ohNumber)
	{
		this.ohNumber_ = ohNumber;
	}

	public void setOhRangeFrom(int ohRangeFrom)
	{
		this.ohRangeFrom_ = ohRangeFrom;
	}

	public void setOhRangeTo(int ohRangeTo)
	{
		this.ohRangeTo_ = ohRangeTo;
	}

	public void setRtRangeFrom(double rtRangeFrom)
	{
		this.rtRangeFrom_ = rtRangeFrom;
	}

	public void setRtRangeTo(double rtRangeTo)
	{
		this.rtRangeTo_ = rtRangeTo;
	}

	public void setAdducts(ArrayList<AdductVO> adducts)
	{
		this.adducts_ = adducts;
	}

	public void setHeadGroupFormulaString(String headGroupFormulaString) throws ChemicalFormulaException
	{
		this.headGroupFormulaString_ = headGroupFormulaString;
		this.headgroupFormula_ = null;
		this.headgroupFormula_ = StaticUtils.categorizeFormula(headGroupFormulaString, true);
	}
	
	public void setMinChainC(int minChainC)
	{
		this.minChainC_ = minChainC;
	}

	public void setMaxChainC(int maxChainC)
	{
		this.maxChainC_ = maxChainC;
	}

	public void setMinChainDB(int minChainDB)
	{
		this.minChainDB_ = minChainDB;
	}

	public void setMaxChainDB(int maxChainDB)
	{
		this.maxChainDB_ = maxChainDB;
	}

	public void setNumberOfFAChains(int numberOfFAChains)
	{
		this.numberOfFAChains_ = numberOfFAChains;
	}
	
	public void setNumberOfLCBChains(int numberOfLCBChains)
	{
		this.numberOfLCBChains_ = numberOfLCBChains;
	}

	public void setFaChainList(String faChainList)
	{
		this.faChainList_ = faChainList;
	}

	public void setLCBChainList(String lcbChainList)
	{
		this.lcbChainList_ = lcbChainList;
	}
	
	
	
	
	
}
