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
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.exception.SheetNotPresentException;
import at.tugraz.genome.lda.msn.parser.FALibParser;
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
	
	private final static int OPTIONS_ROW = 0;
	
	public final static String HEADER_NAME = "Name";
	public final static String HEADER_COLON = "";
	public final static String HEADER_DBS = "dbs";
	public final static String HEADER_MASS_NEUTRAL = "M";
  public final static String HEADER_PSM = "PSM";
  public final static String HEADER_RETENTION_TIME = "tR (min)";
  
  private final static int HEADER_ROW = 2;
  private final static int FIRST_VALUE_ROW = 3;
	
	private String outPath_;
	private ArrayList<LipidClassVO> lipidClasses_;
	private String exportOption_;
	
	public MassListExporter(String outPath, ArrayList<LipidClassVO> lipidClasses, String exportOption)
	{
		this.outPath_ = outPath;
		this.lipidClasses_ = lipidClasses;
		this.exportOption_ = exportOption;
	}
	
	public void export()
	{
		try (	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outPath_));
				XSSFWorkbook workbook = new XSSFWorkbook();)
		{
			for (LipidClassVO lClassVO : lipidClasses_)
			{
				Sheet sheet = workbook.createSheet(lClassVO.getLipidClass());
				
				XSSFCellStyle headerStyle = ExcelUtils.getMassListHeaderStyle(workbook);
				createHeader(OPTIONS_ROW, sheet, createOptionsTitles(lClassVO), headerStyle);
				List<String> headerTitles = createHeaderTitles(lClassVO);
				createHeader(HEADER_ROW, sheet, headerTitles, headerStyle);
				FALibParser faParser = new FALibParser(lClassVO.getFaChainListPath());
		    faParser.parseFile();
		    ArrayList<FattyAcidVO> chains = faParser.getFattyAcidSet(Settings.getFaHydroxyEncoding().getEncodedPrefix((short) 0));
		    String psmString = getPSMString(lClassVO);
		    int rowCount = FIRST_VALUE_ROW;
				Row row;
		    Cell cell;
				for (int i=lClassVO.getMinChainC(); i<=lClassVO.getMaxChainC();i++)
				{
					for (int j=lClassVO.getMinChainDB(); j<=lClassVO.getMaxChainDB(); j++)
					{
						ArrayList<ArrayList<FattyAcidVO>> filtered = getPossCombis(lClassVO.getNumberOfChains(),i,j,chains,new ArrayList<FattyAcidVO>()); 
						if (filtered.isEmpty()) continue; //no combinations possible
						Hashtable<String,Hashtable<String,Integer>> prefixElements = new Hashtable<String,Hashtable<String,Integer>>();
						prefixElements = computeElementsFromFAChains(prefixElements,i,j,lClassVO,filtered);
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
		          for (String element : elements.keySet())
		          {
		          	cell = row.createCell(headerTitles.indexOf(element),HSSFCell.CELL_TYPE_NUMERIC);
	              cell.setCellValue(elements.get(element));
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
			workbook.write(out);
	    System.out.println("workbook written!");
		}
		catch (IOException | RulesException | SheetNotPresentException | HydroxylationEncodingException | ChemicalFormulaException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private boolean isAdductExport(int charge)
	{
		if ( (exportOption_.equals(MassListCreatorPanel.EXPORT_OPTION_NEG) && charge < 0)
    		|| (exportOption_.equals(MassListCreatorPanel.EXPORT_OPTION_POS) && charge > 0)
    		|| (exportOption_.equals(MassListCreatorPanel.EXPORT_OPTION_BOTH) && charge != 0) )
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
		if (lClassVO.getOxRangeFrom() < 1 && lClassVO.getOxRangeTo() < 1)
		{
			return "";
		}
		for (int i=lClassVO.getOxRangeFrom(); i<=lClassVO.getOxRangeTo();i++)
		{
			if (i<1) continue;
			builder.append(";O"+i);
		}
		return builder.toString();
	}
	
	private ArrayList<ArrayList<FattyAcidVO>> getPossCombis(int numberOfChains, int cAtoms, int dbs, ArrayList<FattyAcidVO> fas, ArrayList<FattyAcidVO> added){
		ArrayList<ArrayList<FattyAcidVO>> combis = new ArrayList<ArrayList<FattyAcidVO>>();
    for (FattyAcidVO fa : fas) {
    	ArrayList<FattyAcidVO> toAdd = new ArrayList<FattyAcidVO>(added);
      toAdd.add(fa);
      if (numberOfChains>1) {
      	ArrayList<ArrayList<FattyAcidVO>> combis2 = getPossCombis(numberOfChains-1, cAtoms, dbs, fas, toAdd);
        if (combis2.size()>0) {
          combis.addAll(combis2);
        }
      }else{
        if (checkTotalCDbsValid(toAdd,cAtoms,dbs)) {
          combis.add(toAdd);
        }
      }
    }
    return combis;
  }
	
	private boolean checkTotalCDbsValid(ArrayList<FattyAcidVO> toAdd, int cAtoms, int dbs) {
    int totC = 0;
    int totD = 0;
    for (FattyAcidVO fa : toAdd) {
      totC += fa.getcAtoms();
      totD += fa.getDoubleBonds();
    }
    return (totC==cAtoms && totD==dbs);
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
	private Hashtable<String,Hashtable<String,Integer>> computeElementsFromFAChains(Hashtable<String,Hashtable<String,Integer>>labelElements, 
  		int fattyAcidC, int dbs, LipidClassVO lClassVO, ArrayList<ArrayList<FattyAcidVO>> filtered) throws ChemicalFormulaException
  {
  	labelElements = computeElementsFromUnlabeledChains(labelElements,fattyAcidC,dbs,lClassVO);
  	for (ArrayList<FattyAcidVO> combs : filtered) 
    {
    	Vector<String> labels = new Vector<String>();
    	Hashtable<String,Integer> classElements = new Hashtable<String,Integer>(labelElements.get(""));
    	for (FattyAcidVO fa : combs) 
    	{
    		if (fa.getPrefix().length()>0) 
    		{
    			labels.add(fa.getPrefix());
    			Hashtable<String,Integer> formula = StaticUtils.categorizeFormula(fa.getFormula(), true);
    			if (formula.containsKey("D"))
    			{
    				classElements.put("D", classElements.get("D")+formula.get("D"));
    				classElements.put("H", classElements.get("H")-formula.get("D"));
    			}
    			if (formula.containsKey("Cc"))
    			{
    				classElements.put("Cc", classElements.get("Cc")+formula.get("Cc"));
    				classElements.put("C", classElements.get("C")-formula.get("Cc"));
    			}
    		}
    	}
    	if (!labels.isEmpty()) 
    	{
    		Collections.sort(labels);
      	StringBuilder builder = new StringBuilder();
      	for (String label : labels) 
    		{
    			builder.append(label);
    		}
      	String labelID = builder.toString();
      	labelElements.put(labelID, classElements);
    	}
    }
  	return labelElements;
  }
	
	private Hashtable<String,Hashtable<String,Integer>> computeElementsFromUnlabeledChains(Hashtable<String,Hashtable<String,Integer>> labelElements, 
  		int chainC, int dbs, LipidClassVO lClassVO)
  {
  	Hashtable<String,Integer> elements = new Hashtable<String,Integer>(lClassVO.getHeadgroupFormula());
  	Integer adjustedC = elements.get("C")+chainC;
  	elements.put("C", adjustedC);
  	Integer adjustedH = elements.get("H") + chainC*2-dbs*2;
  	elements.put("H", adjustedH);
  	elements.put("Cc", 0);
  	elements.put("D", 0);
  	labelElements.put("", elements);
  	return labelElements;
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
	
	private List<String> createHeaderTitles(LipidClassVO lClassVO) 
  {
    List<String> headerTitles = new ArrayList<String>();
    headerTitles.add(MassListExporter.HEADER_NAME);
    headerTitles.add(MassListExporter.HEADER_COLON);
    headerTitles.add(MassListExporter.HEADER_DBS);
    for (String element : determineHeaderElements(lClassVO)) 
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
	
	private String getAdductHeader(AdductVO adduct)
	{
		return String.format("mass(form[%s] name[%s] charge=%s)", adduct.getFormulaString(), adduct.getAdductName(), adduct.getCharge());
	}
	
	private List<String> createOptionsTitles(LipidClassVO lClassVO) 
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
	
	private ArrayList<String> determineHeaderElements(LipidClassVO lClassVO)
	{
		ArrayList<String> elements = new ArrayList<String>(lClassVO.getHeadgroupFormula().keySet());
		elements.add("Cc");
		elements.add("D");
		Collections.sort(elements);
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
