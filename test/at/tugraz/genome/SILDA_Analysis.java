package at.tugraz.genome;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFCellStyle;
import org.apache.poi.hssf.usermodel.HSSFFont;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFFont;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.QuantificationThread;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.QuantVO;

/**
 * This class will allow for 
 * 1. Reading in masslists and printing an excel file with:
 * 		a) all found FA with count in how many different molecular species they were found to be incorporated into
 * 		b) all found FA per C=C position with count
 * 		c) all found FA per lipid class with count
 * 		d) all found FA per lipid class + C=C position (maybe)
 * @author Leonida M. Lamp
 *
 */
public class SILDA_Analysis
{
//	private static final String TARGET_LIST_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\recalibration\\SILDA_5.xlsx";
//	private static final String TARGET_LIST_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\analyzedTL\\SILDA_7_analyzed.xlsx";
//	private static final String OUT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\analyzedTL\\SILDA_7_FA_count.xlsx";
	private static final String TARGET_LIST_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_30min_final\\SILDA_30min.xlsx";
	private static final String OUT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_30min_final\\TL_analysis\\SILDA_30min_count.xlsx";
	
	private final static int TOTALS_ROW = 0;
	private final static int HEADER_ROW = 1;
	private final static String HEADER_UNIQUE_MOL = "Unique molecular species with assigned C=C positions";
	private final static String HEADER_ALL_MOL = "All molecular species with assigned C=C positions";
  private final static String HEADER_ALL_FA = "All FA with assigned C=C positions";
	
  public static void main(String[] args)
  {
  	analyseTargetList();
  }
  
	@SuppressWarnings("unchecked")
	public static void analyseTargetList()
	{
		
		try (	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(OUT_PATH));
				XSSFWorkbook workbook = new XSSFWorkbook();)
		{
			Vector<?> quantInfo = QuantificationThread.getCorrectAnalyteSequence(TARGET_LIST_PATH,false);
			LinkedHashMap<String,Integer> classSequence = (LinkedHashMap<String,Integer>)quantInfo.get(0);
			LinkedHashMap<String,Vector<String>> analyteSequence = (LinkedHashMap<String,Vector<String>>)quantInfo.get(1);
			Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo.get(4);
			
			XSSFCellStyle headerStyle = getHeaderStyle(workbook); 
			Set<Integer> allDoubleBondPositions = new HashSet<Integer>();
			Hashtable<String, Vector<DoubleBondPositionVO>> allDoubleBondPositionVOsForClass = new Hashtable<String, Vector<DoubleBondPositionVO>>();
			Vector<String> allFAwithDBPositions = new Vector<String>();
			Hashtable<String, Vector<String>> allFAwithDBPositionsForClass = new Hashtable<String, Vector<String>>();
			
			
	    for (String cName : classSequence.keySet()) 
	    {
	    	if (!allDoubleBondPositionVOsForClass.containsKey(cName))
	    	{
	    		allDoubleBondPositionVOsForClass.put(cName, new Vector<DoubleBondPositionVO>());
	    		allFAwithDBPositionsForClass.put(cName, new Vector<String>());
	    	}
	    	for (String analyte : analyteSequence.get(cName)) 
		    {
	    		Hashtable<String,QuantVO> quantAnalytes = quantObjects.get(cName).get(analyte);
		      for (String mod : quantAnalytes.keySet()) 
		      {
		        QuantVO quant = quantAnalytes.get(mod);
		        Vector<DoubleBondPositionVO> quantDB = quant.getInfoForOmegaAssignment();
		        Vector<String> fattyAcidsWithDoubleBondPositions = collectFAwithDBPositions(quantDB);
		        allDoubleBondPositions.addAll(collectDoubleBondPositions(quantDB));
		        allDoubleBondPositionVOsForClass.get(cName).addAll(quantDB);
		        allFAwithDBPositions.addAll(fattyAcidsWithDoubleBondPositions);
		        allFAwithDBPositionsForClass.get(cName).addAll(fattyAcidsWithDoubleBondPositions);
		      }
		    }
	    }
	    
	    Sheet sheet1 = workbook.createSheet("All molecular species");
	    Sheet sheet2 = workbook.createSheet("Molecular species for class");
	    writeAllMolecularSpecies(sheet1, headerStyle, allDoubleBondPositionVOsForClass);
	    writeAllMolecularSpeciesForClass(sheet2, headerStyle, allDoubleBondPositionVOsForClass);
	    for (int doubleBondPosition : allDoubleBondPositions)
	    {
	    	Sheet sheet = workbook.createSheet(String.format("n-%s Molecular Species", doubleBondPosition));
	    	writeMolecularSpecieswithDB(sheet, headerStyle, allDoubleBondPositionVOsForClass, doubleBondPosition);
	    }
	    Sheet sheet3 = workbook.createSheet("All FA");
	    Sheet sheet4 = workbook.createSheet("FA for class");
	    writeAllFA(sheet3, headerStyle, allFAwithDBPositions, allDoubleBondPositions);
	    writeAllFAForClass(sheet4, headerStyle, allFAwithDBPositionsForClass);
	    for (int doubleBondPosition : allDoubleBondPositions)
	    {
	    	Sheet sheet = workbook.createSheet(String.format("n-%s FA", doubleBondPosition));
	    	writeFAwithDB(sheet, headerStyle, allFAwithDBPositionsForClass, doubleBondPosition);
	    }
	    
	    workbook.write(out);
	    System.out.println("done");
		}
		catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void writeMolecularSpecieswithDB(Sheet sheet, XSSFCellStyle headerStyle, Hashtable<String, Vector<DoubleBondPositionVO>> allDoubleBondPositionVOsForClass, int doubleBondPosition)
	{
		Hashtable <String, Vector<DoubleBondPositionVO>> doubleBondPositionVOsForClass = new Hashtable <String, Vector<DoubleBondPositionVO>>();
		for (String cName : allDoubleBondPositionVOsForClass.keySet())
		{
			doubleBondPositionVOsForClass.put(cName, new Vector<DoubleBondPositionVO>());
			Vector<DoubleBondPositionVO> analytes = allDoubleBondPositionVOsForClass.get(cName);
			for (DoubleBondPositionVO analyte : analytes)
			{
				if (analyte.getDoubleBondPositionsHumanReadable().contains(String.format("(n-%s)", doubleBondPosition)))
				{
					doubleBondPositionVOsForClass.get(cName).add(analyte);
				}
			}
		}
		writeAllMolecularSpeciesForClass(sheet, headerStyle, doubleBondPositionVOsForClass);
	}
	
	private static void writeFAwithDB(Sheet sheet, XSSFCellStyle headerStyle, Hashtable<String, Vector<String>> allFAwithDBPositionsForClass, int doubleBondPosition)
	{
		Hashtable <String, Vector<String>> allFAwithDBPositionForClass = new Hashtable <String, Vector<String>>();
		for (String cName : allFAwithDBPositionsForClass.keySet())
		{
			allFAwithDBPositionForClass.put(cName, new Vector<String>());
			Vector<String> analytes = allFAwithDBPositionsForClass.get(cName);
			for (String analyte : analytes)
			{
				if (analyte.contains(String.format("(n-%s)", doubleBondPosition)))
				{
					allFAwithDBPositionForClass.get(cName).add(analyte);
				}
			}
		}
		writeAllFAForClass(sheet, headerStyle, allFAwithDBPositionForClass);
	}
	
	private static void writeAllFAForClass(Sheet sheet, XSSFCellStyle headerStyle, Hashtable<String, Vector<String>> allFAwithDBPositionsForClass)
	{
		List<String> headerTitles = new ArrayList<String>();
		
		ArrayList<String> classes = new ArrayList<String>(allFAwithDBPositionsForClass.keySet());
		Collections.sort(classes);
		for (String cName : classes)
		{
			Vector<String> fas = allFAwithDBPositionsForClass.get(cName);
			Collections.sort(fas);
			Set<String> uniqueFAs = new LinkedHashSet<String>(fas);
			headerTitles.add(cName);
			
			int rowCount = HEADER_ROW+1;
			for (String fa : uniqueFAs)
			{
				Row row = sheet.getRow(rowCount);
				if (row == null)
				{
					row = sheet.createRow(rowCount);
				}
				Cell cell = row.createCell(headerTitles.indexOf(cName),HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue(fa);
        rowCount++;
			}
			
		}
		createHeader(sheet, headerTitles, headerStyle);
	}
	
	private static void writeAllFA(Sheet sheet, XSSFCellStyle headerStyle, Vector<String> allFAwithDBPositions, Set<Integer> allDoubleBondPositions)
	{
		Collections.sort(allFAwithDBPositions);
		Set<String> uniqueFAs = new LinkedHashSet<String>(allFAwithDBPositions);
		List<String> headerTitles = new ArrayList<String>();
		List<Integer> countsForDB = new ArrayList<Integer>();
		headerTitles.add(HEADER_ALL_FA);
		int rowCount = HEADER_ROW+1;
		countsForDB.add(rowCount); //for the "All FA" column
		for (Integer db : allDoubleBondPositions)
		{
			headerTitles.add(String.format("n-%s", db));
			countsForDB.add(rowCount); //for each db column
		}
		
		createHeader(sheet, headerTitles, headerStyle);
		
		sheet.setColumnWidth(headerTitles.indexOf(HEADER_ALL_FA), 75 * 256);
		
		for (String fa : uniqueFAs)
		{
			Row row = sheet.createRow(rowCount);
			Cell cell = row.createCell(headerTitles.indexOf(HEADER_ALL_FA),HSSFCell.CELL_TYPE_STRING);
      cell.setCellValue(fa);
      int db = Integer.parseInt(fa.substring(fa.indexOf("(n-")+3, fa.indexOf(")")));
      int dbIndex = headerTitles.indexOf(String.format("n-%s", db));
      Row dbRow = sheet.getRow(countsForDB.get(dbIndex));
      Cell dbCell = dbRow.createCell(dbIndex,HSSFCell.CELL_TYPE_STRING);
      dbCell.setCellValue(fa);
      countsForDB.set(dbIndex, countsForDB.get(dbIndex)+1);
      rowCount++;
		}
	}
	
	private static void writeAllMolecularSpeciesForClass(Sheet sheet, XSSFCellStyle headerStyle, Hashtable<String, Vector<DoubleBondPositionVO>> allDoubleBondPositionVOsForClass)
	{
		List<String> headerTitles = new ArrayList<String>();
		
		ArrayList<String> classes = new ArrayList<String>(allDoubleBondPositionVOsForClass.keySet());
		Collections.sort(classes);
		
		Row row = sheet.createRow(TOTALS_ROW);
		for (String cName : classes)
		{
			Set<String> uniqueNames = new LinkedHashSet<String>();
			Hashtable<String, ArrayList<DoubleBondPositionVO>> toRemoveLookup = new Hashtable<String, ArrayList<DoubleBondPositionVO>>();
			Vector<DoubleBondPositionVO> vos = allDoubleBondPositionVOsForClass.get(cName);
			Collections.sort(vos);
			for (DoubleBondPositionVO vo : vos)
			{
				String id = vo.getEncodedDetailed(true, false);
				if (!toRemoveLookup.containsKey(id))
				{
					toRemoveLookup.put(id, new ArrayList<DoubleBondPositionVO>());
				}
				uniqueNames.add(vo.getDoubleBondPositionsHumanReadable());
				toRemoveLookup.get(id).add(vo);
			}
			
			ArrayList<String> uniqueNamesList = removeSpeciesWithLessOmegaInfo(toRemoveLookup, uniqueNames);
			
			headerTitles.add(cName);
			
			row = sheet.getRow(TOTALS_ROW);
			Cell cell = row.createCell(headerTitles.indexOf(cName),HSSFCell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("Total: %s", uniqueNamesList.size()));
			
			int rowCount = HEADER_ROW+1;
			for (String name : uniqueNamesList)
			{
				row = sheet.getRow(rowCount);
				if (row == null)
				{
					row = sheet.createRow(rowCount);
				}
				cell = row.createCell(headerTitles.indexOf(cName),HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue(name);
        rowCount++;
			}
			
		}
		createHeader(sheet, headerTitles, headerStyle);
	}
	
	private static void writeAllMolecularSpecies(Sheet sheet, XSSFCellStyle headerStyle, Hashtable<String, Vector<DoubleBondPositionVO>> allDoubleBondPositionVOsForClass)
	{
		List<String> headerTitles = new ArrayList<String>();
		headerTitles.add(HEADER_ALL_MOL);
		headerTitles.add(HEADER_UNIQUE_MOL);
		
		createHeader(sheet, headerTitles, headerStyle);
		
		sheet.setColumnWidth(headerTitles.indexOf(HEADER_ALL_MOL), 75 * 256);
		sheet.setColumnWidth(headerTitles.indexOf(HEADER_UNIQUE_MOL), 75 * 256);
		
		Integer countAll = 0;
		Integer countUnique = 0;
		
		int rowCount = HEADER_ROW+1;
		ArrayList<String> classes = new ArrayList<String>(allDoubleBondPositionVOsForClass.keySet());
		Collections.sort(classes);
		for (String cName : classes)
		{
			Set<String> uniqueNames = new LinkedHashSet<String>();
			Hashtable<String, ArrayList<DoubleBondPositionVO>> toRemoveLookup = new Hashtable<String, ArrayList<DoubleBondPositionVO>>();
			ArrayList<String> allNames = new ArrayList<String>();
			Vector<DoubleBondPositionVO> vos = allDoubleBondPositionVOsForClass.get(cName);
			Collections.sort(vos);
			for (DoubleBondPositionVO vo : vos)
			{
				String id = vo.getEncodedDetailed(true, false);
				if (!toRemoveLookup.containsKey(id))
				{
					toRemoveLookup.put(id, new ArrayList<DoubleBondPositionVO>());
				}
				uniqueNames.add(vo.getDoubleBondPositionsHumanReadable());
				allNames.add(vo.getDoubleBondPositionsHumanReadable());
				toRemoveLookup.get(id).add(vo);
			}
			
			ArrayList<String> uniqueNamesList = removeSpeciesWithLessOmegaInfo(toRemoveLookup, uniqueNames);
			countAll += allNames.size();
			countUnique += uniqueNamesList.size();
			
			if (cName.contains("-")) // O- or P- prefix is part of the mol species name
			{
				cName = cName.substring(cName.indexOf("-")+1);
			}
			
			for (int i=0; i<allNames.size(); i++)
			{
				Row row = sheet.createRow(rowCount);
				Cell cell = row.createCell(headerTitles.indexOf(HEADER_ALL_MOL),HSSFCell.CELL_TYPE_STRING);
	      cell.setCellValue(String.format("%s %s", cName, allNames.get(i)));
	      if (i<uniqueNamesList.size())
	      {
	      	cell = row.createCell(headerTitles.indexOf(HEADER_UNIQUE_MOL),HSSFCell.CELL_TYPE_STRING);
	        cell.setCellValue(String.format("%s %s", cName, uniqueNamesList.get(i)));
	      }
	      rowCount++;
			}
		}
		
		Row row = sheet.createRow(TOTALS_ROW);
		Cell cell = row.createCell(headerTitles.indexOf(HEADER_ALL_MOL),HSSFCell.CELL_TYPE_STRING);
		cell.setCellValue(String.format("Total: %s", countAll));
		cell = row.createCell(headerTitles.indexOf(HEADER_UNIQUE_MOL),HSSFCell.CELL_TYPE_STRING);
		cell.setCellValue(String.format("Total: %s", countUnique));
	}
	
	/** 
	 * This method removes VOs which have a matching partner with more assigned C=C positions for the count.
	 * @param toRemoveLookup
	 * @param uniqueNames
	 * @return
	 */
	private static ArrayList<String> removeSpeciesWithLessOmegaInfo(Hashtable<String, ArrayList<DoubleBondPositionVO>> toRemoveLookup, Set<String> uniqueNames)
	{
		ArrayList<String> uniqueNamesList = new ArrayList<String>(uniqueNames);
		
		for (String id : toRemoveLookup.keySet())
		{
			ArrayList<DoubleBondPositionVO> vosOfSameMolSpecies = toRemoveLookup.get(id);
			Set<Vector<Integer>> patterns = new HashSet<Vector<Integer>>();
			for (DoubleBondPositionVO vo : vosOfSameMolSpecies)
			{
				Vector<Integer> pattern = vo.getPositionAssignmentPattern();
				if (!pattern.contains(-1)) //add only fully assigned species
				{
					patterns.add(vo.getPositionAssignmentPattern());
				}
			}
			for (DoubleBondPositionVO vo : vosOfSameMolSpecies)
			{
				Vector<Integer> pattern = vo.getPositionAssignmentPattern();
				if (pattern.contains(-1))
				{
					for (Vector<Integer> fullyAssigned : patterns)
					{
						for (int i=0; i<fullyAssigned.size();i++)
						{
							if (!pattern.get(i).equals(-1) && pattern.get(i).equals(fullyAssigned.get(i))) //I rely on the fact, that the FA in the string are sorted.
							{
								uniqueNamesList.remove(vo.getDoubleBondPositionsHumanReadable());
							}
						}
					}
				}
			}
		}
		return uniqueNamesList;
	}
	
	private static Vector<Integer> collectDoubleBondPositions(Vector<DoubleBondPositionVO> vos)
	{
		Vector<Integer> dbs = new Vector<Integer>();
		for (DoubleBondPositionVO vo : vos)
		{
			for (FattyAcidVO fa : vo.getChainCombination())
			{
				if (fa.getOmegaPosition() > 0)
				{
					dbs.add(fa.getOmegaPosition());
				}
			}
		}
		return dbs;
	}
	
	private static Vector<String> collectFAwithDBPositions(Vector<DoubleBondPositionVO> vos) throws LipidCombinameEncodingException
	{
		Vector<String> fattyAcidsWithDoubleBondPositions = new Vector<String>();
		for (DoubleBondPositionVO vo : vos)
		{
			Vector<FattyAcidVO> chainCombination = vo.getChainCombination();
			for (FattyAcidVO fa : chainCombination)
			{
				if (fa.getOmegaPosition()>0)
				{
					String name = StaticUtils.getHumanReadableChainName(fa, Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), fa.getOhNumber()>0);
					fattyAcidsWithDoubleBondPositions.add(name);
				}
			}
		}
		return fattyAcidsWithDoubleBondPositions;
	}
	
	/**
   * Creates a formatted header row with the given header titles in the row number given by HEADER_ROW
   * @param ws Excel worksheet to write the header to
   * @param headerTitles List of header titles
   */
  private static void createHeader(Sheet sheet, List<String> headerTitles, XSSFCellStyle headerStyle) 
  {
  	Row row = sheet.createRow(HEADER_ROW);
  	Cell cell;
  	
  	for (int i=0; i<headerTitles.size(); i++) 
  	{
  		cell = row.createCell(i, HSSFCell.CELL_TYPE_STRING);
  		cell.setCellValue(headerTitles.get(i));
  		cell.setCellStyle(headerStyle);
  		sheet.setColumnWidth(i, 20 * 256);
  	}
  }
	
	public static XSSFCellStyle getHeaderStyle(XSSFWorkbook wb)
	{
    XSSFCellStyle arial12style = wb.createCellStyle();
    XSSFFont arial12font = wb.createFont();
    arial12font.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);
    arial12font.setFontName("Arial");
    arial12font.setFontHeightInPoints((short)12);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(HSSFCellStyle.ALIGN_CENTER);
    return arial12style;
  }
}
