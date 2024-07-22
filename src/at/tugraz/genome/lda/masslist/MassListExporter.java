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

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import javax.swing.JFrame;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.exception.SheetNotPresentException;
import at.tugraz.genome.lda.msn.parser.FALibParser;
import at.tugraz.genome.lda.msn.parser.SPBLibParser;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.AdductVO;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;

public class MassListExporter
{
	public final static String OPTION_ADDUCT_INSENSITIVE_RT_FILTER = "adductInsensitiveRtFilter";
	public final static String OPTION_PICK_BEST_MATCH_BY_SPECTRUM_COVERAGE = "pickBestMatchBySpectrumCoverage";
	public final static String OPTION_START_RT = "Start-RT:";
	public final static String OPTION_STOP_RT = "Stop-RT:";
	public final static String OPTION_OH_NUMBER = "OH-Number:";
	public final static String OPTION_OH_RANGE = "OH-Range:";
	
	public final static String HEADER_NAME = "Name";
	public final static String HEADER_COLON = "";
	public final static String HEADER_DBS = "dbs";
	public final static String HEADER_MASS_NEUTRAL = "neutral mass";
  public final static String HEADER_PSM = "PSM";
  public final static String HEADER_RETENTION_TIME = "tR (min)";
  
  public final static String HEADER_CD_CLASS = "Lipid Class";
  public final static String HEADER_CD_SPECIES = "Lipid Species";
  public final static String HEADER_CD_FORMULA = "Chemical Formula";
  public final static String HEADER_CD_MASS_NEUTRAL = "Neutral Mass";
  public final static String HEADER_CD_ADDUCT_NAME = "Adduct Name";
  public final static String HEADER_CD_ADDUCT_FORMULA = "Adduct Formula";
  public final static String HEADER_CD_ADDUCT_MASS = "Adduct m/z";
  
  private final static int OPTIONS_ROW_LDA = 0;
  private final static int HEADER_ROW_LDA = 2;
  private final static int FIRST_VALUE_ROW_LDA = 3;
	
	private String outPath_;
	private ArrayList<LipidClassVO> lipidClasses_;
	private String exportIonMode_;
	private String exportFormat_;
	
	public MassListExporter(String outPath, ArrayList<LipidClassVO> lipidClasses, String exportIonMode, String exportFormat)
	{
		this.outPath_ = outPath;
		this.lipidClasses_ = lipidClasses;
		this.exportIonMode_ = exportIonMode;
		this.exportFormat_ = exportFormat;
	}
	
	public void export()
	{
		try (	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outPath_));
				XSSFWorkbook workbook = new XSSFWorkbook();)
		{
			if (exportFormat_.equals(MassListCreatorPanel.EXPORT_FORMAT_LDA))
			{
				exportLDAFormat(workbook);
			}
			else if (exportFormat_.equals(MassListCreatorPanel.EXPORT_FORMAT_COMPOUND_DISCOVERER))
			{
				exportCDFormat(workbook);
			}
			
			workbook.write(out);
	    System.out.println("workbook written!");
		}
		catch (IOException | RulesException | SheetNotPresentException | ChemicalFormulaException ex) {
			new WarningMessage(new JFrame(), "Error", "The following error occurred during the export: "+ex.getMessage());
		}
	}
	
	
	private void exportCDFormat(XSSFWorkbook workbook) throws IOException, RulesException, SheetNotPresentException, ChemicalFormulaException
	{
		Sheet sheet = workbook.createSheet("Mass List");
		XSSFCellStyle headerStyle = ExcelUtils.getMassListHeaderStyle(workbook);
		List<String> headerTitles = createCDHeaderTitles();
		createHeader(0, sheet, headerTitles, headerStyle);
		int rowCount = 0;
		Row row;
    Cell cell;
		for (LipidClassVO lClassVO : lipidClasses_)
		{
			ArrayList<FattyAcidVO> chainsFA = new ArrayList<FattyAcidVO>();
	  	ArrayList<FattyAcidVO> chainsLCB = new ArrayList<FattyAcidVO>();
	    if (lClassVO.getNumberOfFAChains() > 0)
	    {
	    	FALibParser faParser = new FALibParser(lClassVO.getFAChainListPath());
		    faParser.parseFile();
	    	chainsFA = faParser.getFattyAcidSet(false); //do not include the oxstate, as that is currently buggy
	    }
	    if (lClassVO.getNumberOfLCBChains() > 0)
	    {
	    	SPBLibParser faParser = new SPBLibParser(new File(lClassVO.getLCBChainListPath()));
		    faParser.parseFile();
	    	chainsLCB = faParser.getFattyAcidSet(false); //do not include the oxstate, as that is currently buggy
	    }
			
	    for (int i=lClassVO.getMinChainC(); i<=lClassVO.getMaxChainC(); i++)
			{
				for (int j=lClassVO.getMinChainDB(); j<=lClassVO.getMaxChainDB(); j++)
				{
					for (int k=lClassVO.getOhRangeFrom(); k<=lClassVO.getOhRangeTo(); k++)
					{
						ArrayList<ArrayList<FattyAcidVO>> filtered = getPossCombis(lClassVO.getNumberOfFAChains(),lClassVO.getNumberOfLCBChains(),i,j,k,
								chainsFA,chainsLCB, new ArrayList<FattyAcidVO>()); 
						if (filtered.isEmpty()) continue; //no combinations possible
						
						Hashtable<String,Hashtable<String,Integer>> prefixElements = new Hashtable<String,Hashtable<String,Integer>>();
						prefixElements = computeElementsFromChains(prefixElements,lClassVO,filtered);
						for (String label : prefixElements.keySet())
	          {
							Hashtable<String,Integer> elements = prefixElements.get(label);
							double massNeutral = computeNeutralMass(elements);
							for (AdductVO adduct : lClassVO.getAdducts())
	            {
	            	if ( isAdductExport(adduct.getCharge()) )
								{
	            		double massAdduct = computeAdductMass(massNeutral, adduct);
	            		row = sheet.createRow(++rowCount);
	            		cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_CD_CLASS),HSSFCell.CELL_TYPE_STRING);
	                cell.setCellValue(lClassVO.getLipidClass());
	                cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_CD_SPECIES),HSSFCell.CELL_TYPE_STRING);
	                cell.setCellValue(buildLipidSpeciesString(lClassVO, i,j,k));
	                cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_CD_FORMULA),HSSFCell.CELL_TYPE_STRING);
	                cell.setCellValue(StaticUtils.getFormulaInHillNotation(elements, false));
	                cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_CD_MASS_NEUTRAL),HSSFCell.CELL_TYPE_NUMERIC);
	                cell.setCellValue(massNeutral);
	                cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_CD_ADDUCT_NAME),HSSFCell.CELL_TYPE_STRING);
	                cell.setCellValue(adduct.getAdductName());
	                cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_CD_ADDUCT_FORMULA),HSSFCell.CELL_TYPE_STRING);
	                cell.setCellValue(StaticUtils.getFormulaInHillNotation(computeAdductFormula(elements, adduct), false));
	                cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_CD_ADDUCT_MASS),HSSFCell.CELL_TYPE_NUMERIC);
	                cell.setCellValue(massAdduct);
								}
	            }
	          }
					}
				}
			}
		}
	}
	
	private String buildLipidSpeciesString(LipidClassVO lClassVO, int cAtoms, int doubleBonds, int ohNumber)
	{
		String ohString = "";
		if (ohNumber == 1)
		{
			ohString = ";O";
		}
		else if (ohNumber > 1)
		{
			ohString = ";O"+ohNumber;
		}
		
		return String.format("%s %s:%s%s", lClassVO.getLipidClass(), cAtoms, doubleBonds, ohString);
	}
	
	private void exportLDAFormat(XSSFWorkbook workbook) throws IOException, RulesException, SheetNotPresentException, ChemicalFormulaException
	{
		for (LipidClassVO lClassVO : lipidClasses_)
		{
			Sheet sheet = workbook.createSheet(lClassVO.getLipidClass());
			
			XSSFCellStyle headerStyle = ExcelUtils.getMassListHeaderStyle(workbook);
			createHeader(OPTIONS_ROW_LDA, sheet, createLDAOptionsTitles(lClassVO), headerStyle);
			List<String> headerTitles = createLDAHeaderTitles(lClassVO);
			createHeader(HEADER_ROW_LDA, sheet, headerTitles, headerStyle);
			
	  	ArrayList<FattyAcidVO> chainsFA = new ArrayList<FattyAcidVO>();
	  	ArrayList<FattyAcidVO> chainsLCB = new ArrayList<FattyAcidVO>();
	    if (lClassVO.getNumberOfFAChains() > 0)
	    {
	    	FALibParser faParser = new FALibParser(lClassVO.getFAChainListPath());
		    faParser.parseFile();
	    	chainsFA = faParser.getFattyAcidSet(false); //do not include the oxstate, as that is currently buggy
	    }
	    if (lClassVO.getNumberOfLCBChains() > 0)
	    {
	    	SPBLibParser faParser = new SPBLibParser(new File(lClassVO.getLCBChainListPath()));
		    faParser.parseFile();
	    	chainsLCB = faParser.getFattyAcidSet(false); //do not include the oxstate, as that is currently buggy
	    }
	    
	    String psmString = getPSMString(lClassVO);
	    int rowCount = FIRST_VALUE_ROW_LDA;
			Row row;
	    Cell cell;
			for (int i=lClassVO.getMinChainC(); i<=lClassVO.getMaxChainC(); i++)
			{
				for (int j=lClassVO.getMinChainDB(); j<=lClassVO.getMaxChainDB(); j++)
				{
					ArrayList<ArrayList<FattyAcidVO>> filtered = getPossCombis(lClassVO.getNumberOfFAChains(),lClassVO.getNumberOfLCBChains(),i,j,lClassVO.getOhNumber(),
							chainsFA,chainsLCB, new ArrayList<FattyAcidVO>()); 
					if (filtered.isEmpty()) continue; //no combinations possible
					Hashtable<String,Hashtable<String,Integer>> prefixElements = new Hashtable<String,Hashtable<String,Integer>>();
					prefixElements = computeElementsFromChains(prefixElements,lClassVO,filtered);
					for (String label : prefixElements.keySet())
          {
						row = sheet.createRow(rowCount);
          	Hashtable<String,Integer> elements = prefixElements.get(label);
          	
            cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_NAME),HSSFCell.CELL_TYPE_STRING);
	          cell.setCellValue(label+i);
	          cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_COLON),HSSFCell.CELL_TYPE_STRING);
	          cell.setCellValue(":");
	          cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_DBS),HSSFCell.CELL_TYPE_NUMERIC);
	          cell.setCellValue(j);
	          for (String element : determineLDAHeaderElements())
	          {
	          	cell = row.createCell(headerTitles.indexOf(element),HSSFCell.CELL_TYPE_NUMERIC);
	          	if (!elements.containsKey(element))
	          	{
	          		cell.setCellValue(0);
	          	}
	          	else
	          	{
	          		cell.setCellValue(elements.get(element));
	          	}
	          }
	          double massNeutral = computeNeutralMass(elements);
	          cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_MASS_NEUTRAL),HSSFCell.CELL_TYPE_NUMERIC);
	          cell.setCellValue(massNeutral);
            for (AdductVO adduct : lClassVO.getAdducts())
            {
            	if ( isAdductExport(adduct.getCharge()) )
							{
            		double massAdduct = computeAdductMass(massNeutral, adduct);
              	cell = row.createCell(headerTitles.indexOf(getAdductHeader(adduct)),HSSFCell.CELL_TYPE_NUMERIC);
                cell.setCellValue(massAdduct);
							}
            }
            cell = row.createCell(headerTitles.indexOf(MassListExporter.HEADER_PSM),HSSFCell.CELL_TYPE_STRING);
	          cell.setCellValue(psmString);
            rowCount++;
          }
				}
			}
		}
	}
	
	
	private boolean isAdductExport(int charge)
	{
		if ( (exportIonMode_.equals(MassListCreatorPanel.EXPORT_OPTION_NEG) && charge < 0)
    		|| (exportIonMode_.equals(MassListCreatorPanel.EXPORT_OPTION_POS) && charge > 0)
    		|| (exportIonMode_.equals(MassListCreatorPanel.EXPORT_OPTION_BOTH) && charge != 0) )
		{
			return true;
		}
		return false;
	}
	
	/**
	 * TODO: ensure the syntax is fully correct in edge cases.
	 * @param lClassVO
	 * @return
	 */
	private String getPSMString(LipidClassVO lClassVO)
	{
		StringBuilder builder = new StringBuilder();
		if (lClassVO.getNumberOfLCBChains() > 0)
		{
			return "";
		}
		for (int i=lClassVO.getOhRangeFrom(); i<=lClassVO.getOhRangeTo();i++)
		{
			if (i<1) continue;
			builder.append(";O"+i);
		}
		return builder.toString();
	}
	
	private ArrayList<ArrayList<FattyAcidVO>> getPossCombis(int numberOfFAChains, int numberOfLCBChains, int cAtoms, int dbs, int oxNum, 
			ArrayList<FattyAcidVO> fas, ArrayList<FattyAcidVO> lcbs, ArrayList<FattyAcidVO> added){
		ArrayList<ArrayList<FattyAcidVO>> combis = new ArrayList<ArrayList<FattyAcidVO>>();
		int totalChains = numberOfFAChains + numberOfLCBChains;
		
		ArrayList<FattyAcidVO> toIterate = fas;
		int faMinus = 1;
		int lcbMinus = 0;
		if (numberOfFAChains == 0 && totalChains>0)
		{
			toIterate = lcbs;
			faMinus = 0;
			lcbMinus = 1;
		}
		
		int maxC = checkMaxCAvailable(added, cAtoms);
		int maxD = checkMaxDBAvailable(added, dbs);
		int maxO = checkMaxOxAvailable(added, oxNum);
		
    for (FattyAcidVO fa : toIterate) {
    	if (fa.getcAtoms()<=maxC && fa.getDoubleBonds()<=maxD && fa.getOhNumber()<=maxO)
    	{
    		ArrayList<FattyAcidVO> toAdd = new ArrayList<FattyAcidVO>(added);
        toAdd.add(fa);
        if (totalChains>1) {
        	ArrayList<ArrayList<FattyAcidVO>> combis2 = getPossCombis(numberOfFAChains-faMinus, numberOfLCBChains-lcbMinus, cAtoms, dbs, oxNum, fas, lcbs, toAdd);
          combis.addAll(combis2);
        }else{
          if (checkTotalCDbsOxValid(toAdd,cAtoms,dbs,oxNum)) {
            combis.add(toAdd);
          }
        }
    	}
    }
    return combis;
  }
	
	private int checkMaxCAvailable(ArrayList<FattyAcidVO> added, int cAtoms)
	{
		int tot = cAtoms;
    for (FattyAcidVO fa : added) {
    	tot -= fa.getcAtoms();
    }
    return tot;
	}
	
	private int checkMaxDBAvailable(ArrayList<FattyAcidVO> added, int dbs)
	{
		int tot = dbs;
    for (FattyAcidVO fa : added) {
    	tot -= fa.getDoubleBonds();
    }
    return tot;
	}
	
	private int checkMaxOxAvailable(ArrayList<FattyAcidVO> added, int ox)
	{
		int tot = ox;
    for (FattyAcidVO fa : added) {
    	tot -= fa.getOhNumber();
    }
    return tot;
	}
	
	private boolean checkTotalCDbsOxValid(ArrayList<FattyAcidVO> toAdd, int cAtoms, int dbs, int oxNum) {
		int totC = 0;
    int totD = 0;
    int totO = 0;
    for (FattyAcidVO fa : toAdd) {
      totC += fa.getcAtoms();
      totD += fa.getDoubleBonds();
      totO += fa.getOhNumber();
    }
    return (totC==cAtoms && totD==dbs && totO==oxNum);
  }
	
	/**
	 * 
	 * @param labelElements
	 * @param fattyAcidC
	 * @param dbs
	 * @param lClassVO
	 * @param filtered
	 * @return
	 * @throws ChemicalFormulaException
	 */
	private Hashtable<String,Hashtable<String,Integer>> computeElementsFromChains(Hashtable<String,Hashtable<String,Integer>>labelElements, 
  		LipidClassVO lClassVO, ArrayList<ArrayList<FattyAcidVO>> filtered) throws ChemicalFormulaException
  {
		for (ArrayList<FattyAcidVO> comb : filtered)
    {
    	Vector<String> labels = new Vector<String>();
    	Hashtable<String,Integer> combElements = new Hashtable<String,Integer>(lClassVO.getHeadgroupFormula());
    	for (FattyAcidVO fa : comb) 
    	{
    		Hashtable<String,Integer> formula = StaticUtils.categorizeFormula(fa.getFormula(), true);
    		for (String element : formula.keySet())
    		{
    			if (!combElements.containsKey(element))
    			{
    				combElements.put(element, 0);
    			}
    			combElements.put(element, combElements.get(element) + formula.get(element));
    		}
    		combElements.put("H", combElements.get("H")-1); //removing a hydrogen from each FA
    		labels.add(fa.getPrefix());
    	}
    	Collections.sort(labels);
    	StringBuilder builder = new StringBuilder();
    	for (String label : labels) 
  		{
  			builder.append(label);
  		}
    	String labelID = builder.toString();
    	//removing a hydroxy group for each FA chain that is esterified to an LCB chain (a hydrogen is already removed earlier)
    	if (lClassVO.getNumberOfLCBChains()>0)
    	{
    		combElements.put("H", combElements.get("H")-lClassVO.getNumberOfFAChains());
    		combElements.put("O", combElements.get("O")-lClassVO.getNumberOfFAChains());
    	}
    	labelElements.put(labelID, combElements);
    }
  	return labelElements;
  }
	
	private Hashtable<String,Integer> computeAdductFormula(Hashtable<String,Integer> neutralElements, AdductVO adduct)
	{
		Hashtable<String,Integer> elements = new Hashtable<String,Integer>();
		Hashtable<String,Integer> elementsAdduct = new Hashtable<String,Integer>(adduct.getFormula());
		Set<String> allElements = new HashSet<String>(neutralElements.keySet());
		allElements.addAll(elementsAdduct.keySet());
		for (String element : allElements)
		{
			int neutralEle = neutralElements.get(element) != null ? neutralElements.get(element) : 0;
			int adductEle = elementsAdduct.get(element) != null ? elementsAdduct.get(element) : 0;
			elements.put(element, neutralEle+adductEle);
		}
		return elements;
	}
	
	private double computeNeutralMass(Hashtable<String,Integer> elements)
	{
		double neutralMass = 0.0;
		ElementConfigParser parser = Settings.getElementParser();
		
		for (String element : elements.keySet())
		{
			neutralMass += parser.getElementDetails(element).getMonoMass()*elements.get(element);
		}
		return neutralMass;
	}
	
	private double computeAdductMass(double neutralMass, AdductVO adduct)
	{
		double massAdduct = neutralMass;
  	for (String element : adduct.getFormula().keySet())
		{
  		massAdduct += Settings.getElementParser().getElementDetails(element).getMonoMass()*adduct.getFormula().get(element);
		}
		return Math.abs(massAdduct/adduct.getCharge());
	}
	
	private List<String> createLDAHeaderTitles(LipidClassVO lClassVO) 
  {
    List<String> headerTitles = new ArrayList<String>();
    headerTitles.add(MassListExporter.HEADER_NAME);
    headerTitles.add(MassListExporter.HEADER_COLON);
    headerTitles.add(MassListExporter.HEADER_DBS);
    for (String element : determineLDAHeaderElements()) 
    {
    	headerTitles.add(element);
    }
    headerTitles.add(MassListExporter.HEADER_MASS_NEUTRAL);
    for (AdductVO adduct : lClassVO.getAdducts())
    {
    	if ( isAdductExport(adduct.getCharge()) )
    	{
    		headerTitles.add(getAdductHeader(adduct));
    	}
    }
    
    headerTitles.add(MassListExporter.HEADER_PSM);
    headerTitles.add(MassListExporter.HEADER_RETENTION_TIME);
    
    return headerTitles;
  }
	
	private List<String> createCDHeaderTitles() 
  {
    List<String> headerTitles = new ArrayList<String>();
    headerTitles.add(MassListExporter.HEADER_CD_CLASS);
    headerTitles.add(MassListExporter.HEADER_CD_SPECIES);
    headerTitles.add(MassListExporter.HEADER_CD_FORMULA);
    headerTitles.add(MassListExporter.HEADER_CD_MASS_NEUTRAL);
    headerTitles.add(MassListExporter.HEADER_CD_ADDUCT_NAME);
    headerTitles.add(MassListExporter.HEADER_CD_ADDUCT_FORMULA);
    headerTitles.add(MassListExporter.HEADER_CD_ADDUCT_MASS);
    
    return headerTitles;
  }

	
	private String getAdductHeader(AdductVO adduct)
	{
		//TODO: absolute value for charge state for backward compatibility (12.07.2024)
		return String.format("mass(form[%s] name[%s] charge=%s)", adduct.getFormulaString(), adduct.getAdductName(), Math.abs(adduct.getCharge()));
	}
	
	private List<String> createLDAOptionsTitles(LipidClassVO lClassVO) 
  {
    List<String> titles = new ArrayList<String>();
    if (lClassVO.isAdductInsensitiveRtFilter())
    {
    	titles.add(OPTION_ADDUCT_INSENSITIVE_RT_FILTER);
    }
    if (lClassVO.isPickBestMatchBySpectrumCoverage())
    {
    	titles.add(OPTION_PICK_BEST_MATCH_BY_SPECTRUM_COVERAGE);
    }
    if (lClassVO.getRtRangeFrom() > -1 && lClassVO.getRtRangeTo() > -1 && lClassVO.getRtRangeFrom() < lClassVO.getRtRangeTo())
    {
    	titles.add(OPTION_START_RT+lClassVO.getRtRangeFrom());
    	titles.add(OPTION_STOP_RT+lClassVO.getRtRangeTo());
    }
    if (lClassVO.getOhNumber()>0)
    {
    	titles.add(OPTION_OH_NUMBER+lClassVO.getOhNumber());
    }
    if (lClassVO.getOhRangeFrom() > -1 && lClassVO.getOhRangeTo() > -1 && lClassVO.getOhRangeFrom() < lClassVO.getOhRangeTo())
    {
    	titles.add(OPTION_OH_RANGE+lClassVO.getOhRangeFrom()+"-"+lClassVO.getOhRangeTo());
    }
    return titles;
  }
	
	//TODO: write out formula for each compound instead, this will be more flexible
	private ArrayList<String> determineLDAHeaderElements()
	{
		ArrayList<String> elements = new ArrayList<String>();
		elements.add("C");
		elements.add("Cc");
		elements.add("H");
		elements.add("D");
		elements.add("N");
		elements.add("O");
		elements.add("P");
		return elements;
	}
	
	/**
	 * Creates a formatted header row with the given header titles in the row number given by @param rowNum
	 * @param rowNum
	 * @param sheet
	 * @param headerTitles
	 * @param headerStyle
	 */
  private void createHeader(int rowNum, Sheet sheet, List<String> headerTitles, XSSFCellStyle headerStyle) 
  {
  	Row row = sheet.createRow(rowNum);
  	Cell cell;
  	for (int i=0; i<headerTitles.size(); i++) 
  	{
  		cell = row.createCell(i, HSSFCell.CELL_TYPE_STRING);
  		cell.setCellValue(headerTitles.get(i));
  		cell.setCellStyle(headerStyle);
  		sheet.setColumnWidth(i, 10 * 256);
  	}
  }
	
	
	
}
