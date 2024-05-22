package at.tugraz.genome;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JFrame;

//import org.apache.commons.math3.util.Precision;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.dhatim.fastexcel.Worksheet;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.vos.RTCheckedVO;
import javafx.util.Pair;

/**
 * Class for analysis of the ASC output
 * 
 * @author Leonida M. Lamp
 *
 */
public class ASCAnalyzer
{
//	private final static String FILE_PATH_ASC = "D:\\Collaborator_Files\\Christer\\MS-DIAL_comparison\\LCMSdata_liver_MSDIALsearch_rtChecked.xlsx";
//	private final static String FILE_PATH_ASC = "D:\\Collaborator_Files\\Christer\\MS-DIAL_comparison\\LCMSdata_liver_targetsearch_rtChecked_allScore_Na_12.5_lessHarsh_combRemovalAll.xlsx";
	private final static String FILE_PATH_ASC = "D:\\Collaborator_Files\\Christer\\LC-MS_LDA_liver\\LCMSdata_liver_targetsearch_rtChecked_allScore_Na_12.5_lessHarsh_combRemovalAll.xlsx";
	private final static String FOLDER_PATH_NEG = "D:\\Collaborator_Files\\Christer\\LC-MS_LDA_liver\\negative\\";
	private final static String FOLDER_PATH_POS = "D:\\Collaborator_Files\\Christer\\LC-MS_LDA_liver\\positive\\";
	
	public final static String HEADER_ADDUCT = "Adduct";
	public final static String HEADER_ALEX_SCORE = "ALEX score";
	public final static String HEADER_FRAGMENT_TYPE = "Fragment type";
	public final static String HEADER_LIPID_CLASS = "Lipid class";
	public final static String HEADER_LIPID_SPECIES = "Lipid species";
	public final static String HEADER_MOLECULAR_LIPID_SPECIES = "Molecular lipid species";
	public final static String HEADER_POLARITY = "Polarity";
	public final static String HEADER_RT_GROUP = "RT group";
	public final static String HEADER_TP = "TP";
	public final static String HEADER_COMMENT = "Comment";
	
	
	private final static Double RT_GROUPING = 0.4;
	private final static Integer HEADER_ROW = 0;
	
	
	public static void main(String[] args)
  {
		new ASCAnalyzer();
  }
	
	private ASCAnalyzer()
	{
		try
		{
//			exportGroupAssignments();
			exportASCfilteredByLDA();
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
		}
	}
	
	private void exportASCfilteredByLDA()
	{
		Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>> allLDAResults = new Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>>();
		try
		{
			Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>> allASCAssignments = parsePreviousAssignements();
			File file = new File(FOLDER_PATH_NEG);
			File[] files = file.listFiles();
			parseResultFiles(files, allLDAResults);
			file = new File(FOLDER_PATH_POS);
			files = file.listFiles();
			parseResultFiles(files, allLDAResults);
			filterASC(allLDAResults, allASCAssignments);
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
		}
	}
	
	/**
	 * 
	 * @throws ExcelInputFileException
	 * @throws IOException
	 */
	private void filterASC(Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>> allLDAResults, Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>> allASCAssignments) throws ExcelInputFileException, IOException
	{
		ArrayList<RTCheckedVO> toExport = new ArrayList<RTCheckedVO>();
		for (String lipidClass : allASCAssignments.keySet())
		{
			if (allLDAResults.containsKey(lipidClass))
			{
				for (String lipidSpecies : allASCAssignments.get(lipidClass).keySet())
				{
					if (allLDAResults.get(lipidClass).containsKey(lipidSpecies))
					{
						Vector<RTCheckedVO> rtCheckedASC = allASCAssignments.get(lipidClass).get(lipidSpecies);
						Vector<RTCheckedVO> rtCheckedLDA = allLDAResults.get(lipidClass).get(lipidSpecies);
						if (rtCheckedLDA.isEmpty())
						{
//							System.out.format("LDA is empty for %s.\n", lipidSpecies);
							continue;
						}
						
						for (RTCheckedVO checkedASC : rtCheckedASC)
						{
							if (!checkedASC.getFragType().equalsIgnoreCase("MLSSF")) continue;
							for (RTCheckedVO checkedLDA : rtCheckedLDA)
							{
//								if (lipidClass.equals("TAG"))
//								{
//									System.out.format("ASC: %s, LDA: %s\n", checkedASC.getMolSpec(), checkedLDA.getMolSpec());
//									System.out.println("hi");
//								}
								if (checkedASC.getAdduct().equals(checkedLDA.getAdduct()) && checkedASC.getMolSpec().equals(checkedLDA.getMolSpec()))
								{
									Double rtASC = Double.parseDouble(checkedASC.getRtGroup());
									Double rtLDA = Double.parseDouble(checkedLDA.getRtGroup());
									if (Math.abs(rtASC - rtLDA) <= RT_GROUPING)
									{
										if (!toExport.contains(checkedASC)) toExport.add(checkedASC);
										continue;
									}
								}
							}
						}
//						System.out.format("Match for %s.\n", lipidSpecies);
					}
					else
					{
//						System.out.format("Lipid species %s is not in the LDA results. \n", lipidSpecies);
					}
				}
			}
			else
			{
				System.out.format("Lipid class %s is not in the LDA results. \n", lipidClass);
			}
		}
		Collections.sort(toExport);
		exportASCFile(toExport);
	}
	
	private void exportASCFile(ArrayList<RTCheckedVO> toExport) throws FileNotFoundException, IOException
	{
		try (BufferedOutputStream out = new BufferedOutputStream(
				new FileOutputStream(FILE_PATH_ASC.substring(0,FILE_PATH_ASC.lastIndexOf("."))+"_LDA_filtered.xlsx"));) 
  	{
			String s = Settings.VERSION;
			// the constructor can only take a version number in the format xx.yyyy
			org.dhatim.fastexcel.Workbook wb = new org.dhatim.fastexcel.Workbook(out, "Lipid Data Analyzer",
					s.substring(0, s.indexOf(".", s.indexOf(".") + 1)));

			Worksheet ws1 = wb.newWorksheet("LDA filtered");
			
			createHeaderExtended(ws1);
			writeRowsExtended(ws1, toExport);
			wb.finish();
			System.out.println("done");
		}
	}
	
	private void writeRowsExtended(Worksheet ws, ArrayList<RTCheckedVO> entries)
	{
		for (int row = 1; row < entries.size(); row++) {
			ws.value(row, 0, entries.get(row-1).getAdduct());
			ws.value(row, 1, entries.get(row-1).getScore());
			ws.value(row, 2, entries.get(row-1).getFragType());
			ws.value(row, 3, entries.get(row-1).getlClass());
			ws.value(row, 4, entries.get(row-1).getSpecies());
			ws.value(row, 5, entries.get(row-1).getMolSpec());
			ws.value(row, 6, entries.get(row-1).getPolarity());
			ws.value(row, 7, entries.get(row-1).getRtGroup());
			ws.value(row, 8, entries.get(row-1).getTruePos());
			ws.value(row, 9, entries.get(row-1).getComment());
		}
	}
	
	private void createHeaderExtended(Worksheet ws)
	{
		List<String> headerTitles = new ArrayList<String>();
		headerTitles.add(HEADER_ADDUCT);
		headerTitles.add(HEADER_ALEX_SCORE);
		headerTitles.add(HEADER_FRAGMENT_TYPE);
		headerTitles.add(HEADER_LIPID_CLASS);
		headerTitles.add(HEADER_LIPID_SPECIES);
		headerTitles.add(HEADER_MOLECULAR_LIPID_SPECIES);
		headerTitles.add(HEADER_POLARITY);
		headerTitles.add(HEADER_RT_GROUP);
		headerTitles.add(HEADER_TP);
		headerTitles.add(HEADER_COMMENT);
		for (int i = 0; i < headerTitles.size(); i++) {
			ws.value(HEADER_ROW, i, headerTitles.get(i));
		}
		ws.range(HEADER_ROW, HEADER_ROW, HEADER_ROW, headerTitles.size()).style()
				.bold().horizontalAlignment("center").fontName("Arial").fontSize(12)
				.set();
	}
	
	private void parseResultFiles(File[] files, Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>> allResults) throws ExcelInputFileException, LipidCombinameEncodingException
  {
  	for (File file : files)
  	{
  		if (!file.isDirectory() && file.getName().endsWith(".xlsx"))
  		{
  			QuantificationResult quantResOriginal = LDAResultReader.readResultFile(file.getAbsolutePath(), new Hashtable<String,Boolean>());
    		groupEntries(file, quantResOriginal, allResults);
  		}
  	}
  }
	
  /**
   * Groups all analytes of the same lipid class
   * @param file
   * @param quantRes
   * @param identifications
   * @throws LipidCombinameEncodingException 
   */
	private void groupEntries(File file, QuantificationResult quantRes, Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>> allResults) throws LipidCombinameEncodingException
  {
  	Hashtable<String,Vector<LipidParameterSet>> quantIdentifications = quantRes.getIdentifications();
  	for (String lipidClass : quantIdentifications.keySet())
  	{
  		String newLipidClass = lipidClass;
  		if (lipidClass.equals("P-PE"))
  		{
  			newLipidClass = "PE O-";
  		}
  		else if (lipidClass.equals("TG"))
			{
				newLipidClass = "TAG";
			}
			else if (lipidClass.equals("DG"))
			{
				newLipidClass = "DAG";
			}
  		if (!allResults.containsKey(newLipidClass)) allResults.put(newLipidClass, new Hashtable<String,Vector<RTCheckedVO>>());
  		for (LipidParameterSet set : quantIdentifications.get(lipidClass))
  		{
  			String sumComp = set.getNameStringWithoutRt();
  			if (lipidClass.equals("P-PE"))
    		{
    			Integer num = Integer.parseInt(sumComp.substring(sumComp.indexOf(":")+1,sumComp.length()));
    			num++;
    			sumComp = sumComp.substring(0, sumComp.indexOf(":")+1) + num.toString();
    		}
  			if (sumComp.startsWith("d") || sumComp.startsWith("t"))
  			{
  				sumComp = replaceSphingoHydroxyStrings(sumComp);
  			}
  			String space = newLipidClass.endsWith("-") ? "" : " ";
  			sumComp = newLipidClass+space+sumComp;
  			
  			String adduct = set.getModificationName();
  			if (adduct.equals("-H"))
  			{
  				adduct = "-H+";
  			}
  			else if (adduct.equals("H"))
  			{
  				adduct = "+H+";
  			}
  			else if (adduct.equals("Na"))
  			{
  				adduct = "+Na+";
  			}
  			else if (adduct.equals("HCOO"))
  			{
  				adduct = "+HCOO-";
  			}
  			else if (adduct.equals("NH4"))
  			{
  				adduct = "+NH4+";
  			}
  			else
  			{
  				continue; //not in DB
  			}
  			
  			if (!allResults.get(newLipidClass).containsKey(sumComp)) allResults.get(newLipidClass).put(sumComp, new Vector<RTCheckedVO>());
  			if (set instanceof LipidomicsMSnSet && ((LipidomicsMSnSet) set).getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED)
  			{
  				LipidomicsMSnSet msnSet = (LipidomicsMSnSet) set;
  				
  				if (lipidClass.equals("PG") || lipidClass.equals("PC") || lipidClass.equals("PE") || lipidClass.equals("PS") || 
  						lipidClass.equals("TG") || lipidClass.equals("DG") || lipidClass.equals("PI") || lipidClass.equals("CL") || 
  						lipidClass.equals("DMPE") || lipidClass.equals("MMPE")) //from fragmentguide, lipidclasses without assigned positions
  				{
  					Vector<String> combis = msnSet.getValidChainCombinations();
    				
    				for (String combi : combis)
    				{
    					Vector<FattyAcidVO> chains = StaticUtils.decodeLipidNamesFromChainCombi(combi);
    					String molSpecies = StaticUtils.getHumanReadableCombiName(chains, Settings.getFaHydroxyEncoding(), Settings.getLcbHydroxyEncoding()).replaceAll("_", "-");
    					molSpecies = newLipidClass+" "+molSpecies;
    					RTCheckedVO vo = new RTCheckedVO(adduct, null, "MLSSF", lipidClass, sumComp, molSpecies, null, set.getRt(), null, null);
    					allResults.get(newLipidClass).get(sumComp).add(vo);
    				}
  				}
  				else
  				{
  					for (String molSpecies : msnSet.getMSnIdentificationNamesWithSNPositions())
    				{
  						if (newLipidClass.equals("PE O-"))
  						{
  							Integer num = Integer.parseInt(molSpecies.substring(molSpecies.indexOf(":")+1,molSpecies.indexOf("/")));
  		    			num++;
  		    			molSpecies = "PE O"+molSpecies.substring(1, molSpecies.indexOf(":")+1)+num.toString()+molSpecies.substring(molSpecies.indexOf("/"),molSpecies.length());
  						}
  						else if (lipidClass.equals("Cer") || lipidClass.equals("SM"))
  						{
  							String first = molSpecies.substring(0, molSpecies.indexOf("/"));
  							String second = molSpecies.substring(molSpecies.indexOf("/")+1, molSpecies.length());
  							first = replaceSphingoHydroxyStrings(first);
  							molSpecies = replaceSphingoHydroxyStrings(first)+"/"+replaceSphingoHydroxyStrings(second);
  						}
  						else
  						{
  							molSpecies = newLipidClass+" "+molSpecies;
  						}					
    					RTCheckedVO vo = new RTCheckedVO(adduct, null, "MLSSF", lipidClass, sumComp, molSpecies, null, set.getRt(), null, null);
    					allResults.get(newLipidClass).get(sumComp).add(vo);
    				}
  				}
  			}
  		}
  	}
  }
	
	private String replaceSphingoHydroxyStrings(String name)
	{
		if (name.startsWith("n"))
		{
			name = name.substring(1, name.length());
		}
		else if (name.startsWith("d"))
		{
			name = name.substring(1, name.length())+";2";
		}
		else if (name.startsWith("t"))
		{
			name = name.substring(1, name.length())+";3";
		}
		return name;
	}
	
	/**
	 * Merges entries in an ASC output file according to the Adducts in one sheet and Adducts plus RTs in a second.
	 * @throws ExcelInputFileException
	 * @throws IOException
	 */
	private void exportGroupAssignments() throws ExcelInputFileException, IOException
	{
		Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>> assignemnts = parsePreviousAssignements();
		ArrayList<RTCheckedVO> sameMol = new ArrayList<RTCheckedVO>();
		ArrayList<RTCheckedVO> sameMolRt = new ArrayList<RTCheckedVO>();
		for (String lipidClass : assignemnts.keySet())
		{
			for (String lipidSpecies : assignemnts.get(lipidClass).keySet())
			{
				Vector<RTCheckedVO> rtChecked = assignemnts.get(lipidClass).get(lipidSpecies);
				Hashtable<String, Vector<RTCheckedVO>> sameMolSpecies = new Hashtable<String, Vector<RTCheckedVO>>();
				for (RTCheckedVO vo : rtChecked)
				{
					if (!sameMolSpecies.containsKey(vo.getMolSpec()))
					{
						sameMolSpecies.put(vo.getMolSpec(), new Vector<RTCheckedVO>());
					}
					sameMolSpecies.get(vo.getMolSpec()).add(vo);
				}
				
				if (rtChecked.size()>sameMolSpecies.keySet().size())
				{
					for (String molSpecies : sameMolSpecies.keySet())
					{
						Vector<RTCheckedVO> sameMolecularSpecies = sameMolSpecies.get(molSpecies);
						if (sameMolecularSpecies.get(0).getFragType().equalsIgnoreCase("MLSSF"))
						{
							if (sameMolecularSpecies.size()>1)
							{
								Vector<RTCheckedVO> consensusRT = getConsensusRT(sameMolecularSpecies);
								for (RTCheckedVO vo : consensusRT)
								{
									sameMolRt.add(vo); //"merged adducts"
								}
								RTCheckedVO consensusTP = getConsensus(sameMolecularSpecies, true);
								if (consensusTP != null)
									sameMol.add(consensusTP);
								RTCheckedVO consensusFP = getConsensus(sameMolecularSpecies, false);
								if (consensusFP != null)
									sameMol.add(consensusFP);
							}
							else
							{
								sameMolRt.add(sameMolecularSpecies.get(0));
								sameMol.add(sameMolecularSpecies.get(0));
							}
						}
					}
				}
			}
		}
		Collections.sort(sameMol);
		Collections.sort(sameMolRt);
		exportAssignments(sameMol, sameMolRt);
	}
	
	private void exportAssignments(ArrayList<RTCheckedVO> sameMol, ArrayList<RTCheckedVO> sameMolRt) throws IOException
	{
		try (BufferedOutputStream out = new BufferedOutputStream(
				new FileOutputStream(FILE_PATH_ASC.substring(0,FILE_PATH_ASC.lastIndexOf("."))+"_merged.xlsx"));) 
		{
			String s = Settings.VERSION;
			// the constructor can only take a version number in the format xx.yyyy
			org.dhatim.fastexcel.Workbook wb = new org.dhatim.fastexcel.Workbook(out, "Lipid Data Analyzer",
					s.substring(0, s.indexOf(".", s.indexOf(".") + 1)));

			Worksheet ws1 = wb.newWorksheet("Merged adducts");
			Worksheet ws2 = wb.newWorksheet("Fully merged");
			
			createHeader(ws1);
			writeRows(ws1, sameMolRt);
			createHeader(ws2);
			writeRows(ws2, sameMol);
			wb.finish();
			System.out.println("done");
		}
	}
	
	private void writeRows(Worksheet ws, ArrayList<RTCheckedVO> entries)
	{
		for (int row = 1; row < entries.size(); row++) {
			ws.value(row, 0, entries.get(row-1).getAdduct());
			ws.value(row, 1, entries.get(row-1).getScore());
			ws.value(row, 2, entries.get(row-1).getlClass());
			ws.value(row, 3, entries.get(row-1).getSpecies());
			ws.value(row, 4, entries.get(row-1).getMolSpec());
			ws.value(row, 5, entries.get(row-1).getFragType());
			ws.value(row, 6, entries.get(row-1).getRtGroup());
			ws.value(row, 7, entries.get(row-1).getTruePos());
		}
	}
	
	private void createHeader(Worksheet ws)
	{
		List<String> headerTitles = new ArrayList<String>();
		headerTitles.add("Adduct(s)");
		headerTitles.add("Score");
		headerTitles.add("Lipid class");
		headerTitles.add("Lipid species");
		headerTitles.add("Molecular lipid species");
		headerTitles.add("Fragment type");
		headerTitles.add("RT group");
		headerTitles.add("TP");
		for (int i = 0; i < headerTitles.size(); i++) {
			ws.value(HEADER_ROW, i, headerTitles.get(i));
		}
		ws.range(HEADER_ROW, HEADER_ROW, HEADER_ROW, headerTitles.size()).style()
				.bold().horizontalAlignment("center").fontName("Arial").fontSize(12)
				.set();
	}
	
	
	private Vector<RTCheckedVO> getConsensusRT(Vector<RTCheckedVO> sameMolecularSpecies)
	{
		HashSet<RTCheckedVO> unique = new HashSet<RTCheckedVO>();
		Vector<Pair<Double,RTCheckedVO>> allRts = new Vector<Pair<Double,RTCheckedVO>>();
		Hashtable<Integer,Vector<RTCheckedVO>> clusteredRTs = new Hashtable<Integer,Vector<RTCheckedVO>>();
		Vector<Pair<Double,RTCheckedVO>> currentCluster = new Vector<Pair<Double,RTCheckedVO>>();
		for (RTCheckedVO vo : sameMolecularSpecies)
		{
			allRts.add(new Pair<Double,RTCheckedVO>(Double.parseDouble(vo.getRtGroup()), vo));
		}
		
		Collections.sort(allRts,new SortByDouble()); //compare by Double
		currentCluster.add(allRts.get(0));
		
		int count = 0;
		for (int i = 1; i < allRts.size(); i++) 
		{
			Pair<Double,RTCheckedVO> current = allRts.get(i);
			Pair<Double,RTCheckedVO> lastInCluster = currentCluster.get(currentCluster.size() - 1);
			
			//if is within the grouping parameter and the 
			if ((current.getKey() - lastInCluster.getKey() <= RT_GROUPING) && (current.getValue().getTruePos().equals(lastInCluster.getValue().getTruePos()))) 
			{
				currentCluster.add(allRts.get(i));
			}
			else //start new cluster
			{
				Vector<RTCheckedVO> clustered = new Vector<RTCheckedVO>();
				for (Pair<Double,RTCheckedVO> member : currentCluster)
				{
					clustered.add(member.getValue());
				}
				clusteredRTs.put(count++, clustered);
		    currentCluster = new Vector<Pair<Double,RTCheckedVO>>();
		    currentCluster.add(allRts.get(i));
			}
		}
		
		Vector<RTCheckedVO> clustered = new Vector<RTCheckedVO>();
		for (Pair<Double,RTCheckedVO> member : currentCluster)
		{
			clustered.add(member.getValue());
		}
		clusteredRTs.put(count++, clustered);
		
		
		for (Integer num : clusteredRTs.keySet())
		{
			Vector<RTCheckedVO> clusterMembers = clusteredRTs.get(num);
			if (clusterMembers.size()>1)
			{
				unique.add(getConsensus(clusterMembers, clusterMembers.get(0).getTruePos().equals("true")));
			}
			else
			{
				unique.add(clusterMembers.get(0));
			}
		}
		return new Vector<RTCheckedVO>(unique);
	}
	
	class SortByDouble implements Comparator<Pair<Double,RTCheckedVO>> 
	{

		@Override
		public int compare(Pair<Double,RTCheckedVO> arg0,
				Pair<Double,RTCheckedVO> arg1)
		{
			return Comparator.comparing((Pair<Double,RTCheckedVO> arg) -> arg.getKey())
					.thenComparing((Pair<Double,RTCheckedVO> arg) -> arg.getValue().getTruePos())
					.compare(arg0,arg1);
		}
	}
	
	private RTCheckedVO getConsensus(Vector<RTCheckedVO> sameMolecularSpecies, boolean tp)
	{
		String highestScoreString = "";
		Double highestScore = -1.0;
		String highestScoreRT = "";
		HashSet<String> adducts = new HashSet<String>();
		RTCheckedVO first = null;
		
		for (RTCheckedVO vo : sameMolecularSpecies)
		{
			if ( (tp && vo.getTruePos().equalsIgnoreCase("true")) || (!tp && !vo.getTruePos().equalsIgnoreCase("true")) )
			{
				first = vo;
				Double score = Double.parseDouble(vo.getScore());
				if (score > highestScore)
				{
					highestScore = score;
					highestScoreString = vo.getScore();
					highestScoreRT = vo.getRtGroup();
				}
				adducts.add(vo.getAdduct());
			}
		}
		if (first != null)
		{
			RTCheckedVO vo = new RTCheckedVO(first);
			vo.setAdduct(generateAdductString(adducts));
			vo.setScore(highestScoreString);
			vo.setRtGroup_(highestScoreRT);
			vo.setTruePos_(tp ? "true" : "false");
			return vo;
		}
		return null;
	}
	
	private String generateAdductString(HashSet<String> adducts)
	{
		Vector<String> sortedAdducts = new Vector<String>(adducts);
		Collections.sort(sortedAdducts);
		StringBuilder builder = new StringBuilder();
		for (int i=0;i<sortedAdducts.size();i++)
		{
			builder.append(sortedAdducts.get(i));
			if (i<sortedAdducts.size()-1)
			{
				builder.append(",");
			}
		}
		return builder.toString();
	}
	
	private Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>> parsePreviousAssignements() throws ExcelInputFileException{
    Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>> previousAssignments = new Hashtable<String,Hashtable<String,Vector<RTCheckedVO>>>();
    InputStream myxls = null;
    Workbook workbook  = null;
    try {
      myxls = new FileInputStream(FILE_PATH_ASC);
      workbook  = new XSSFWorkbook(myxls);     
      for (int sheetNumber=0;sheetNumber!=workbook.getNumberOfSheets();sheetNumber++) {
        Sheet sheet = workbook.getSheetAt(sheetNumber);
      
        int adductColumn = -1;
        int scoreColumn = -1;
        int fragTypeColumn = -1;
        int classColumn = -1;
        int speciesColumn = -1;
        int molSpecColumn = -1;
        int rtGroupColumn = -1;
        int tpColumn = -1;
        int polarityColumn = -1;
        int commentColumn = -1;
        
        int highestColumNumber = 0;
        
        int rowCount=0;
        Row row = sheet.getRow(rowCount);
        for (int i=0;  row!=null && i!=(row.getLastCellNum()+1);i++){
          Cell cell = row.getCell(i);
          String contents = "";
//          Double numeric = null;
          int cellType = -1;
          if (cell!=null) cellType = cell.getCellType();
          if (cellType==Cell.CELL_TYPE_STRING){
            contents = cell.getStringCellValue();
          }else if (cellType==Cell.CELL_TYPE_NUMERIC || cellType==Cell.CELL_TYPE_FORMULA){
           contents = String.valueOf(cell.getNumericCellValue());
          }
          contents = contents.trim();
          if (contents==null || contents.length()==0)
            continue;
          if (contents.equalsIgnoreCase("Adduct")) {
            adductColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("ALEX score")) {
            scoreColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("MSDIAL score (avg)")) {
          	scoreColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("Fragment type")) {
            fragTypeColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("Lipid class")) {
            classColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("Lipid species")) {
            speciesColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("ALEX annotation")) //take ALEX annotation over molecular lipid species, relevant for MS-DIAL files
          {
            molSpecColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("Molecular lipid species")) {
            molSpecColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("RT group")) {
            rtGroupColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("Polarity")) {
            polarityColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("TP")) {
            tpColumn = i;
            highestColumNumber = i;
          } else if (contents.equalsIgnoreCase("Comment")) {
            commentColumn = i;
            highestColumNumber = i;
          }
        }
        if (adductColumn<0 || scoreColumn<0 || classColumn<0 || speciesColumn<0 || molSpecColumn<0 || rtGroupColumn<0 || tpColumn<0)
          continue;
        for (rowCount=1;rowCount<(sheet.getLastRowNum()+1);rowCount++){
          row = sheet.getRow(rowCount);
          
          String adduct = null;
          String score = null;
          String fragType = "MLSSF"; //default for MS-DIAL file
          String lClass = null;
          String species = null;
          String molSpec = null;
          String rtGroup = null;
          String truePos = null;
          String polarity = null;
          String comment = null;
          
          for (int i=0;  row!=null && i!=(highestColumNumber<row.getLastCellNum() ? (highestColumNumber+1) : (row.getLastCellNum()+1));i++){
            Cell cell = row.getCell(i);
            String contents = "";
            Double numeric = null;
            int cellType = -1;
            if (cell!=null) cellType = cell.getCellType();
            if (cellType==Cell.CELL_TYPE_STRING){
              contents = cell.getStringCellValue();
              try{ numeric = new Double(contents);}catch(NumberFormatException nfx){};
            }else if (cellType==Cell.CELL_TYPE_NUMERIC || cellType==Cell.CELL_TYPE_FORMULA){
             numeric = cell.getNumericCellValue();
             contents = String.valueOf(numeric);
            }else if (cellType==Cell.CELL_TYPE_BOOLEAN) {
              contents = String.valueOf(cell.getBooleanCellValue());
            }
            
            if (cellType == Cell.CELL_TYPE_STRING)
              contents = cell.getStringCellValue();
            else if (cellType == Cell.CELL_TYPE_NUMERIC){
              double cellValue = -1;
              cellValue = cell.getNumericCellValue();
              contents = String.valueOf(cellValue);
            }
            
            if (i==adductColumn)
              adduct = contents;
            else if (i==scoreColumn)
              score = contents;
            else if (i==fragTypeColumn)
              fragType = contents;
            else if (i==classColumn)
              lClass = contents;
            else if (i==speciesColumn)
              species = contents;
            else if (i==molSpecColumn)
              molSpec = contents;
            else if (i==rtGroupColumn)
              rtGroup = contents;
            else if (i==polarityColumn)
            	polarity = contents;
            else if (i==tpColumn)
              truePos = contents;
            else if (i==commentColumn)
              comment = contents;
          }

          if (adduct==null || adduct.length()==0 || score==null || score.length()==0 || fragType==null || fragType.length()==0 || lClass==null || lClass.length()==0 || species==null || species.length()==0 ||
              molSpec==null || molSpec.length()==0 || rtGroup==null || rtGroup.length()==0 || truePos==null || truePos.length()==0)
            continue;

          RTCheckedVO vo = new RTCheckedVO(adduct, score, fragType, lClass, species, molSpec, polarity, rtGroup, truePos, comment);
          Hashtable<String,Vector<RTCheckedVO>> speciesHash = new Hashtable<String,Vector<RTCheckedVO>>();
          if (previousAssignments.containsKey(lClass)) speciesHash = previousAssignments.get(lClass);
          Vector<RTCheckedVO> diffRts = new Vector<RTCheckedVO>();
          if (speciesHash.containsKey(species)) diffRts = speciesHash.get(species);
          diffRts.add(vo);
          speciesHash.put(species, diffRts);
          previousAssignments.put(lClass, speciesHash);
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
      new WarningMessage(new JFrame(), "ERROR", e.getMessage()+"; it does not seem to be Microsoft Excel");
      throw new ExcelInputFileException(e);
    } catch (Exception e) {
      e.printStackTrace();
      new WarningMessage(new JFrame(), "ERROR", e.getMessage());
      throw new ExcelInputFileException(e);
    } finally {
      if (workbook!=null) {
       try {workbook.close();}catch (IOException e) {
         e.printStackTrace();
       }
      }
      if (myxls!=null) {
        try {myxls.close();}catch (IOException e) {
         e.printStackTrace();
       }
      }
    }
    return previousAssignments;
    
  }
}
