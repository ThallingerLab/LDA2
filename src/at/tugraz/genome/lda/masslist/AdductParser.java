package at.tugraz.genome.lda.masslist;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.vos.AdductVO;

public class AdductParser
{
	private final static String ADDUCT_FOLDER = "./massListCreation/adducts";
	private final static String ADDUCT_SUFFIX = ".txt";
	
	public AdductParser(){}
	
	public ArrayList<AdductVO> parse() throws FileNotFoundException, IOException, ChemicalFormulaException
	{
		ArrayList<AdductVO> allDefinedAdducts = new ArrayList<AdductVO>();
		
		File folder = new File(ADDUCT_FOLDER);
		if (!folder.exists())
		{
			throw new IOException(String.format("The adduct folder '%s' does not exist!", ADDUCT_FOLDER));
		}
		File[] fileCandidates = folder.listFiles();
		for (int i=0; i<fileCandidates.length;i++)
		{
			if (fileCandidates[i].getName().endsWith(ADDUCT_SUFFIX))
			{
				try (FileInputStream in= new FileInputStream(fileCandidates[i]))
				{
					Properties properties = new Properties();
					properties.load(in);
					String name = properties.getProperty("name", null);
					String formula = properties.getProperty("formula", null);
					int charge = Integer.parseInt(properties.getProperty("charge", "0"));
					if (name != null && formula != null && charge != 0)
					{
						AdductVO adduct = new AdductVO(name, formula, charge);
						allDefinedAdducts.add(adduct);
					}
					else
					{
						System.out.format("The adduct definition file '%s' does not adhere to the required format!", fileCandidates[i].getName());
					}
				}
			}
		}	
		return allDefinedAdducts;
	}
	
}
