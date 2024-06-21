package at.tugraz.genome.lda.masslist;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Properties;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.AdductVO;

public class LipidClassParser
{
	private final static String LIPID_CLASS_FOLDER = "./massListCreation/lipidClasses";
	private final static String LIPID_CLASS_SUFFIX = ".txt";
	private final static String ADDUCT_SEPARATOR = ",";
	private final static String CHAIN_LIST_FOLDER = "./fattyAcids";
	private final static String CHAIN_LIST_SUFFIX = ".xlsx";
	private ArrayList<AdductVO> allDefinedAdducts_;
	
	public LipidClassParser(ArrayList<AdductVO> allDefinedAdducts)
	{
		this.allDefinedAdducts_ = allDefinedAdducts;
	}
	
	public ArrayList<LipidClassVO> parse() throws IOException, ChemicalFormulaException
	{
		ArrayList<LipidClassVO> allDefinedLipidClasses = new ArrayList<LipidClassVO>();
		
		File folder = new File(LIPID_CLASS_FOLDER);
		if (!folder.exists())
		{
			throw new IOException(String.format("The adduct folder '%s' does not exist!", LIPID_CLASS_FOLDER));
		}
		File[] fileCandidates = folder.listFiles();
		for (int i=0; i<fileCandidates.length;i++)
		{
			if (fileCandidates[i].getName().endsWith(LIPID_CLASS_SUFFIX))
			{
				try (FileInputStream in= new FileInputStream(fileCandidates[i]))
				{
					Properties properties = new Properties();
					properties.load(in);
					String name = properties.getProperty("name", null);
					boolean adductInsensitiveRtFilter = Boolean.parseBoolean(properties.getProperty("adductInsensitiveRtFilter", null));
					boolean pickBestMatchBySpectrumCoverage = Boolean.parseBoolean(properties.getProperty("pickBestMatchBySpectrumCoverage", null));
					int ohNumber = Integer.parseInt(properties.getProperty("OH_number", null));
					int ohRangeFrom = Integer.parseInt(properties.getProperty("OH_range_from", null));
					int ohRangeTo = Integer.parseInt(properties.getProperty("OH_range_to", null));
					int rtRangeFrom = Integer.parseInt(properties.getProperty("RT_range_from", null));
					int rtRangeTo = Integer.parseInt(properties.getProperty("RT_range_to", null));
					int oxRangeFrom = Integer.parseInt(properties.getProperty("Ox_range_from", null));
					int oxRangeTo = Integer.parseInt(properties.getProperty("Ox_range_to", null));
					String adducts = properties.getProperty("adducts", null);
					String headgroupFormula = properties.getProperty("headgroup_formula", null);
					int minChainC = Integer.parseInt(properties.getProperty("min_chain_C", null));
					int maxChainC = Integer.parseInt(properties.getProperty("max_chain_C", null));
					int minChainDB = Integer.parseInt(properties.getProperty("min_chain_DB", null));
					int maxChainDB = Integer.parseInt(properties.getProperty("max_chain_DB", null));
					int numberOfChains = Integer.parseInt(properties.getProperty("number_of_chains", null));
					String faChainList = properties.getProperty("FA_chain_list_name", null);
					
					if (name != null && adducts != null && headgroupFormula != null && faChainList != null)
					{
						ArrayList<AdductVO> selectedAdducts = matchSelectedAdducts(adducts, fileCandidates[i].getName());
						Hashtable<String,Integer> headgroupFormulaTranslated = StaticUtils.categorizeFormula(headgroupFormula);
						String chainListPath = CHAIN_LIST_FOLDER+"/"+faChainList+CHAIN_LIST_SUFFIX;
						
						LipidClassVO vo = new LipidClassVO(name,
								adductInsensitiveRtFilter, pickBestMatchBySpectrumCoverage, ohNumber, ohRangeFrom, ohRangeTo, 
								rtRangeFrom, rtRangeTo, oxRangeFrom, oxRangeTo, selectedAdducts, 
								headgroupFormulaTranslated, minChainC, maxChainC, minChainDB, maxChainDB, numberOfChains, chainListPath);
						
						allDefinedLipidClasses.add(vo);
					}
					else
					{
						throw new IOException(String.format(
								"The lipid class definition file '%s' does not adhere to the required format!", 
								fileCandidates[i].getName()));
					}
				}
			}
		}	
		return allDefinedLipidClasses;
	}
	
	private ArrayList<AdductVO> matchSelectedAdducts(String adducts, String fileName) throws IOException
	{
		ArrayList<AdductVO> selectedAdducts = new ArrayList<AdductVO>();
		String[] adductNames = adducts.split(ADDUCT_SEPARATOR);
		for (int j=0;j<adductNames.length;j++)
		{
			boolean matched = false;
			for (AdductVO adduct : allDefinedAdducts_)
			{
				if (adduct.getAdductName().equalsIgnoreCase(adductNames[j]))
				{
					selectedAdducts.add(adduct);
					matched = true;
				}
			}
			if (!matched)
			{
				throw new IOException(String.format(
						"The provided adduct name '%s' for the lipid class definition file '%s' is not defined!", 
						adductNames[j], fileName));
			}
		}
		return selectedAdducts;
	}
	
}
