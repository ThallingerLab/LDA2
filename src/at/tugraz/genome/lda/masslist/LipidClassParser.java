package at.tugraz.genome.lda.masslist;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.vos.AdductVO;

public class LipidClassParser
{
	public final static String LIPID_CLASS_FOLDER = "./massListCreation/lipidClasses";
	public final static String LIPID_CLASS_SUFFIX = ".txt";
	public final static String ADDUCT_SEPARATOR = ",";
	public final static String LIPID_CLASS_NAME = "name";
	public final static String LIPID_CLASS_RT_FILTER = "adductInsensitiveRtFilter";
	public final static String LIPID_CLASS_PICK_BEST = "pickBestMatchBySpectrumCoverage";
	public final static String LIPID_CLASS_OH_NUMBER = "OH_number";
	public final static String LIPID_CLASS_OH_RANGE_FROM = "OH_range_from";
	public final static String LIPID_CLASS_OH_RANGE_TO = "OH_range_to";
	public final static String LIPID_CLASS_RT_RANGE_FROM = "RT_range_from";
	public final static String LIPID_CLASS_RT_RANGE_TO = "RT_range_to";
	public final static String LIPID_CLASS_OX_RANGE_FROM = "Ox_range_from";
	public final static String LIPID_CLASS_OX_RANGE_TO = "Ox_range_to";
	public final static String LIPID_CLASS_ADDUCTS = "adducts";
	public final static String LIPID_CLASS_FORMULA = "headgroup_formula";
	public final static String LIPID_CLASS_MIN_CHAIN_C = "min_chain_C";
	public final static String LIPID_CLASS_MAX_CHAIN_C = "max_chain_C";
	public final static String LIPID_CLASS_MIN_CHAIN_DB = "min_chain_DB";
	public final static String LIPID_CLASS_MAX_CHAIN_DB = "max_chain_DB";
	public final static String LIPID_NUMBER_CHAINS = "number_of_chains";
	public final static String FA_LIST_NAME = "FA_chain_list_name";
	
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
					String name = properties.getProperty(LIPID_CLASS_NAME, null);
					boolean adductInsensitiveRtFilter = Boolean.parseBoolean(properties.getProperty(LIPID_CLASS_RT_FILTER, null));
					boolean pickBestMatchBySpectrumCoverage = Boolean.parseBoolean(properties.getProperty(LIPID_CLASS_PICK_BEST, null));
					int ohNumber = Integer.parseInt(properties.getProperty(LIPID_CLASS_OH_NUMBER, null));
					int ohRangeFrom = Integer.parseInt(properties.getProperty(LIPID_CLASS_OH_RANGE_FROM, null));
					int ohRangeTo = Integer.parseInt(properties.getProperty(LIPID_CLASS_OH_RANGE_TO, null));
					double rtRangeFrom = Double.parseDouble(properties.getProperty(LIPID_CLASS_RT_RANGE_FROM, null));
					double rtRangeTo = Double.parseDouble(properties.getProperty(LIPID_CLASS_RT_RANGE_FROM, null));
					int oxRangeFrom = Integer.parseInt(properties.getProperty(LIPID_CLASS_OX_RANGE_FROM, null));
					int oxRangeTo = Integer.parseInt(properties.getProperty(LIPID_CLASS_OX_RANGE_TO, null));
					String adducts = properties.getProperty(LIPID_CLASS_ADDUCTS, null);
					String headgroupFormula = properties.getProperty(LIPID_CLASS_FORMULA, null);
					int minChainC = Integer.parseInt(properties.getProperty(LIPID_CLASS_MIN_CHAIN_C, null));
					int maxChainC = Integer.parseInt(properties.getProperty(LIPID_CLASS_MAX_CHAIN_C, null));
					int minChainDB = Integer.parseInt(properties.getProperty(LIPID_CLASS_MIN_CHAIN_DB, null));
					int maxChainDB = Integer.parseInt(properties.getProperty(LIPID_CLASS_MAX_CHAIN_DB, null));
					int numberOfChains = Integer.parseInt(properties.getProperty(LIPID_NUMBER_CHAINS, null));
					String faChainList = properties.getProperty(FA_LIST_NAME, null);
					
					if (name != null && adducts != null && headgroupFormula != null && faChainList != null)
					{
						ArrayList<AdductVO> selectedAdducts = matchSelectedAdducts(adducts, fileCandidates[i].getName());
						
						LipidClassVO vo = new LipidClassVO(name,
								adductInsensitiveRtFilter, pickBestMatchBySpectrumCoverage, ohNumber, ohRangeFrom, ohRangeTo, 
								rtRangeFrom, rtRangeTo, oxRangeFrom, oxRangeTo, selectedAdducts, headgroupFormula,
								minChainC, maxChainC, minChainDB, maxChainDB, numberOfChains, faChainList);
						
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
