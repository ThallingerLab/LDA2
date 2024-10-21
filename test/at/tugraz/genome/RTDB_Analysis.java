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

import org.apache.commons.math3.util.Precision;
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
public class RTDB_Analysis
{
//	private static final String TARGET_LIST_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\recalibration\\SILDA_5.xlsx";
//	private static final String TARGET_LIST_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\analyzedTL\\SILDA_7_analyzed.xlsx";
//	private static final String OUT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\analyzedTL\\SILDA_7_FA_count.xlsx";
//	private static final String TARGET_LIST_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_30min_final\\final_TL\\SILDA_30min.xlsx";
//	private static final String OUT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_30min_final\\TL_analysis\\SILDA_30min_count.xlsx";
	
	private static final String TARGET_LIST_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\final_TL\\SILDA_60min.xlsx";
	private static final String OUT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\final_TL\\TL_Analysis\\SILDA_60min_analysis.xlsx";
	
	private final static int TOTALS_ROW = 0;
	private final static int HEADER_ROW = 1;
	private final static String HEADER_UNIQUE_MOL_NO_CC = "Unique molecular species";
	private final static String HEADER_UNIQUE_MOL = "Unique molecular species with assigned C=C positions";
	private final static String HEADER_ALL_MOL = "All RT-DB entries";
  private final static String HEADER_ALL_FA = "All FA with assigned C=C positions";
  
  private static ArrayList<Double> diffPerOmega_ = new ArrayList<Double>();
	
  public static void main(String[] args)
  {
  	analyseTargetList();
  	printDiffPerOmega();
  }
  
  private static void printDiffPerOmega()
  {
  	Double sum = 0.0;
  	for (Double rt : diffPerOmega_) sum += rt;
  	System.out.println(sum/diffPerOmega_.size());
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
	    Sheet sheet2 = workbook.createSheet("Per lipid (sub)class => All RT-DB entries");
	    Sheet sheet3 = workbook.createSheet("Per lipid (sub)class => Unique molecular species with assigned C=C positions");
	    writeAllMolecularSpecies(sheet1, headerStyle, allDoubleBondPositionVOsForClass);
	    writeAllMolecularSpeciesForClass(sheet2, headerStyle, allDoubleBondPositionVOsForClass, false);
	    writeAllMolecularSpeciesForClass(sheet3, headerStyle, allDoubleBondPositionVOsForClass, true);
	    
	    for (int doubleBondPosition : allDoubleBondPositions)
	    {
	    	Sheet sheet = workbook.createSheet(String.format("Unique n-%s molecular species", doubleBondPosition));
	    	writeMolecularSpecieswithDB(sheet, headerStyle, allDoubleBondPositionVOsForClass, doubleBondPosition);
	    }
	    Sheet sheet4 = workbook.createSheet("All FA");
	    Sheet sheet5 = workbook.createSheet("All FA => bond-type independent");
	    Sheet sheet6 = workbook.createSheet("FA per lipid (sub)class");
	    writeAllFA(sheet4, headerStyle, allFAwithDBPositions, allDoubleBondPositions, false);
	    writeAllFA(sheet5, headerStyle, allFAwithDBPositions, allDoubleBondPositions, true);
	    writeAllFAForClass(sheet6, headerStyle, allFAwithDBPositionsForClass);
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
		writeAllMolecularSpeciesForClass(sheet, headerStyle, doubleBondPositionVOsForClass, true);
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
	
	private static void writeAllFA(Sheet sheet, XSSFCellStyle headerStyle, Vector<String> allFAwithDBPositions, Set<Integer> allDoubleBondPositions, boolean bondTypeIndependent)
	{
		Collections.sort(allFAwithDBPositions);
		Set<String> uniqueFAs = new LinkedHashSet<String>(allFAwithDBPositions);
		if (bondTypeIndependent)
		{
			Vector<String> unique = new Vector<String>();
			for (String fa : uniqueFAs)
			{
				String name = fa;
				if (name.startsWith("P-") || name.startsWith("O-"))
					name = fa.substring(2, fa.length());
				unique.add(name);
			}
			Collections.sort(unique);
			uniqueFAs = new LinkedHashSet<String>(unique);
		}
		
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
	
	private static void writeAllMolecularSpeciesForClass(Sheet sheet, XSSFCellStyle headerStyle, Hashtable<String, Vector<DoubleBondPositionVO>> allDoubleBondPositionVOsForClass, boolean unique)
	{
		List<String> headerTitles = new ArrayList<String>();
		
		ArrayList<String> classes = new ArrayList<String>(allDoubleBondPositionVOsForClass.keySet());
		Collections.sort(classes);
		
		Row row = sheet.createRow(TOTALS_ROW);
		for (String cName : classes)
		{
			Set<String> uniqueNames = new LinkedHashSet<String>();
			ArrayList<String> allNames = new ArrayList<String>();
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
				allNames.add(String.format("%s; RT=%s min", vo.getDoubleBondPositionsHumanReadable(), Precision.round(vo.getExpectedRetentionTime(), 2)));
			}
			
			ArrayList<String> namesList = unique ? removeSpeciesWithLessOmegaInfo(toRemoveLookup, uniqueNames) : allNames;
			
			headerTitles.add(cName);
			
			row = sheet.getRow(TOTALS_ROW);
			Cell cell = row.createCell(headerTitles.indexOf(cName),HSSFCell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("Total: %s", namesList.size()));
			
			int rowCount = HEADER_ROW+1;
			for (String name : namesList)
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
		headerTitles.add(HEADER_UNIQUE_MOL_NO_CC);
		
		createHeader(sheet, headerTitles, headerStyle);
		
		sheet.setColumnWidth(headerTitles.indexOf(HEADER_ALL_MOL), 75 * 256);
		sheet.setColumnWidth(headerTitles.indexOf(HEADER_UNIQUE_MOL), 75 * 256);
		sheet.setColumnWidth(headerTitles.indexOf(HEADER_UNIQUE_MOL_NO_CC), 75 * 256);
		
		Integer countAll = 0;
		Integer countUniqueCC = 0;
		Integer countUniqueMol = 0;
		
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
				allNames.add(String.format("%s; RT=%s min", vo.getDoubleBondPositionsHumanReadable(), Precision.round(vo.getExpectedRetentionTime(), 2)));
				toRemoveLookup.get(id).add(vo);
			}
			
			ArrayList<String> uniqueNamesList = removeSpeciesWithLessOmegaInfo(toRemoveLookup, uniqueNames);
			ArrayList<String> uniqueMolList = removeOmegaIDs(uniqueNamesList);
			
			countAll += allNames.size();
			countUniqueCC += uniqueNamesList.size();
			countUniqueMol += uniqueMolList.size();
			
			analyzeRTDiffs(allNames, uniqueMolList);
			
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
	      if (i<uniqueMolList.size())
	      {
	      	cell = row.createCell(headerTitles.indexOf(HEADER_UNIQUE_MOL_NO_CC),HSSFCell.CELL_TYPE_STRING);
	        cell.setCellValue(String.format("%s %s", cName, uniqueMolList.get(i)));
	      }
	      rowCount++;
			}
		}
		
		Row row = sheet.createRow(TOTALS_ROW);
		Cell cell = row.createCell(headerTitles.indexOf(HEADER_ALL_MOL),HSSFCell.CELL_TYPE_STRING);
		cell.setCellValue(String.format("Total: %s", countAll));
		cell = row.createCell(headerTitles.indexOf(HEADER_UNIQUE_MOL),HSSFCell.CELL_TYPE_STRING);
		cell.setCellValue(String.format("Total: %s", countUniqueCC));
		cell = row.createCell(headerTitles.indexOf(HEADER_UNIQUE_MOL_NO_CC),HSSFCell.CELL_TYPE_STRING);
		cell.setCellValue(String.format("Total: %s", countUniqueMol));
	}
	
	private static void analyzeRTDiffs(ArrayList<String> allNames, ArrayList<String> uniqueMolList)
	{
		Hashtable<String,ArrayList<String>> grouped = new Hashtable<String,ArrayList<String>>();
		for (String name : uniqueMolList)
		{
			grouped.put(name, new ArrayList<String>());
		}
		for (String name : allNames)
		{
			grouped.get(name.substring(0, name.indexOf(";")).replaceAll("\\([^)]*\\)", "")).add(name);
		}
		for (String name : grouped.keySet())
		{
			ArrayList<String> allMolsOfName = grouped.get(name);
			if (allMolsOfName.size() > 1)
			{
				Hashtable<Integer,ArrayList<Double>> groupedMols = new Hashtable<Integer,ArrayList<Double>>();
				ArrayList<Integer> omegas = new ArrayList<Integer>();
				for (String molName : allMolsOfName)
				{
					Double rt = Double.parseDouble(molName.substring(molName.indexOf("RT=")+3, molName.indexOf(" min")));
					Integer sumOmega = Integer.parseInt(molName.substring(molName.indexOf("(n-")+3, molName.indexOf(")")));
					String lessOmega = molName.replaceFirst("\\([^)]*\\)", "");
					String temp = molName.replaceAll(":", "");
					if (temp.length()<molName.length()-1 && !molName.contains(":0") && !lessOmega.contains("(n-")) 
					{
						continue;
					}
					if (lessOmega.contains("(n-"))
					{
						sumOmega += Integer.parseInt(lessOmega.substring(lessOmega.indexOf("(n-")+3, lessOmega.indexOf(")")));
					}
					if (!groupedMols.containsKey(sumOmega)) 
					{
						groupedMols.put(sumOmega, new ArrayList<Double>());
						omegas.add(sumOmega);
					}
					groupedMols.get(sumOmega).add(rt);
				}
				if (omegas.size()>1)
				{
					Collections.sort(omegas);
					Integer previousOmega = omegas.get(0);
					Double sumRts = 0.0;
					for (Double rt : groupedMols.get(previousOmega)) sumRts += rt;
					Double previousRTs = sumRts / groupedMols.get(previousOmega).size();
					
					for (int i=1; i<omegas.size();i++)
					{
						Integer omegadiff = omegas.get(i)-previousOmega;
						Double sumRtsNew = 0.0;
						for (Double rt : groupedMols.get(omegas.get(i))) sumRtsNew += rt;
						Double rtsNew = sumRtsNew / groupedMols.get(omegas.get(i)).size();
						Double rtDiff = Math.abs(previousRTs-rtsNew);
						diffPerOmega_.add(rtDiff/omegadiff);
					}
				}
			}
		}
	}
	
	
	private static ArrayList<String> removeOmegaIDs(ArrayList<String> uniqueNamesList)
	{
		ArrayList<String> unique = new ArrayList<String>();
		for (String name : uniqueNamesList)
		{
			String uniqueName = name.replaceAll("\\([^)]*\\)", "");
			if (!unique.contains(uniqueName))
			{
				unique.add(uniqueName);
			}	
		}
		return unique;
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
