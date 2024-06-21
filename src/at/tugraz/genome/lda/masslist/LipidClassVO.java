package at.tugraz.genome.lda.masslist;

import java.util.ArrayList;
import java.util.Hashtable;

import at.tugraz.genome.lda.vos.AdductVO;

public class LipidClassVO
{
	private String lClass_;
	private boolean adductInsensitiveRtFilter_;
	private boolean pickBestMatchBySpectrumCoverage_;
	private int ohNumber_;
	private int ohRangeFrom_;
	private int ohRangeTo_;
	private int rtRangeFrom_;
	private int rtRangeTo_;
	private int oxRangeFrom_;
	private int oxRangeTo_;
	private ArrayList<AdductVO> adducts_;
	private Hashtable<String,Integer> headgroupFormula_;
	private int minChainC_;
	private int maxChainC_;
	private int minChainDB_;
	private int maxChainDB_;
	private int numberOfChains_;
	private String faChainListPath_;
	
	
	public LipidClassVO(String lClass, boolean adductInsensitiveRtFilter,
			boolean pickBestMatchBySpectrumCoverage, int ohNumber, int ohRangeFrom,
			int ohRangeTo, int rtRangeFrom, int rtRangeTo, int oxRangeFrom,
			int oxRangeTo, ArrayList<AdductVO> adducts,
			Hashtable<String,Integer> headgroupFormula, int minChainC,
			int maxChainC, int minChainDB, int maxChainDB, int numberOfChains,
			String faChainListPath)
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
		this.oxRangeFrom_ = oxRangeFrom;
		this.oxRangeTo_ = oxRangeTo;
		this.adducts_ = adducts;
		this.headgroupFormula_ = headgroupFormula;
		this.minChainC_ = minChainC;
		this.maxChainC_ = maxChainC;
		this.minChainDB_ = minChainDB;
		this.maxChainDB_ = maxChainDB;
		this.numberOfChains_ = numberOfChains;
		this.faChainListPath_ = faChainListPath;
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


	public int getRtRangeFrom()
	{
		return rtRangeFrom_;
	}


	public int getRtRangeTo()
	{
		return rtRangeTo_;
	}


	public int getOxRangeFrom()
	{
		return oxRangeFrom_;
	}


	public int getOxRangeTo()
	{
		return oxRangeTo_;
	}


	public ArrayList<AdductVO> getAdducts()
	{
		return adducts_;
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


	public int getNumberOfChains()
	{
		return numberOfChains_;
	}


	public String getFaChainListPath()
	{
		return faChainListPath_;
	}
	
	
	
	
	
}
