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

package at.tugraz.genome.lda.fragai;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.dhatim.fastexcel.reader.Cell;
import org.dhatim.fastexcel.reader.CellType;
import org.dhatim.fastexcel.reader.ReadableWorkbook;
import org.dhatim.fastexcel.reader.Row;
import org.dhatim.fastexcel.reader.Sheet;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.AdductVO;

public class ExcelTargetListParser
{
	private final static String HEADER_CLASS = "Class";
	private final static String HEADER_SPECIES = "Species";
	private final static String HEADER_SUM_FORMULA = "Sum Formula";
	private final static String HEADER_ADDUCTS= "Adducts";
	private final static String HEADER_RETENTION_TIME = "Retention Time";
	private final static String HEADER_TOLERANCE= "Tolerance (+-)";
	private final static String ADDUCT_START= "{";
	private final static String ADDUCT_END= "}";
	private final static String ADDUCT_REGEX= "\\}\\{";
	private final static int HEADER_ROW = 0;
	private File directory_;
	private List<String> headerTitles_;
  private ArrayList<FragTargetListEntry> targetListEntries_;
  
  /**
	 * Constructor specifying the directory containing the target list in excel format.
	 * @param directory
	 */
	public ExcelTargetListParser(File directory)
	{
		this.directory_ = directory;
		this.targetListEntries_ = new ArrayList<FragTargetListEntry>();
		
	}
	
  /**
   * Parses the excel file.
   * @throws IOException		if something went wrong parsing the file
   */
  public void parse() throws IOException
  {
  	String supportedFormat = ".xlsx";
    if (!directory_.getAbsolutePath().endsWith(supportedFormat))
    {
    	throw new IOException(String.format("Only the file format '%s' is supported!", supportedFormat));
    }
    try (InputStream is = new FileInputStream(directory_.getAbsolutePath());
        ReadableWorkbook wb = new ReadableWorkbook(is);
        Stream<Sheet> sheets = wb.getSheets();) {
      			sheets.forEach((s) -> { 
            				try { readSheet(s); } 
            				catch (IOException | ChemicalFormulaException ex) 
            				{
            					headerTitles_ = null; //this indicates that the parsing was not successful
            				} });     
        }
    if (headerTitles_ == null)
    {
    	throw new IOException(String.format("The file '%s' was not parsed successfully. Ensure the required worksheet(s) and headers are named according to the template and entries are in the correct format.\n", directory_.toString()));
    }
  }
  
  /**
	 * Read the excel sheet.
	 * @param sheet					the excel sheet to be read
	 * @throws IOException	if there is something wrong with the excel sheet.
   * @throws ChemicalFormulaException if there is something wrong with formulas in the excel sheet.
	 */
  private void readSheet(Sheet sheet) throws IOException, ChemicalFormulaException
	{
		List<Row> rows = null;
    rows = sheet.read();
    Row headerRow = rows.get(HEADER_ROW);
    this.headerTitles_ = readSheetHeaderTitles(headerRow);
    List<Row> contentRows = rows.subList(HEADER_ROW+1, rows.size());
    readContentRows(contentRows);
	}
  
  /**
   * Reads all rows except the header.
   * @param contentRows
   * @throws IOException
   * @throws ChemicalFormulaException
   */
  private void readContentRows(List<Row> contentRows) throws IOException, ChemicalFormulaException
  {
  	if (headerTitles_.contains(HEADER_CLASS) && 
  			headerTitles_.contains(HEADER_SPECIES) && 
  			headerTitles_.contains(HEADER_SUM_FORMULA) &&
  			headerTitles_.contains(HEADER_ADDUCTS) && 
  			headerTitles_.contains(HEADER_RETENTION_TIME) && 
  			headerTitles_.contains(HEADER_TOLERANCE)
  			)
		{
			for (Row row : contentRows) 
			{
	      List<Cell> cells = collectRowCells(row);
	      String lipidClass = null;
	      String species = null;
	    	Hashtable<String,Integer> sumFormula = null;
	    	ArrayList<AdductVO> adducts = null;
	    	Double retentionTime = null;
	    	Double tolerance = null;
	    	
	    	int index;
	      String rawValue;
	      
	      for (Cell cell : cells) {
	        index  = cell.getColumnIndex();
	        rawValue = cell.getRawValue();
	        
	        if (index == headerTitles_.indexOf(HEADER_CLASS)) {
	        	lipidClass = rawValue;
	        } else if (index == headerTitles_.indexOf(HEADER_SPECIES)) {
	        	species = rawValue;
	        } else if (index == headerTitles_.indexOf(HEADER_SUM_FORMULA)) {
						try {sumFormula = StaticUtils.categorizeFormula(rawValue);} catch (Exception ex) {}
	        } else if (index == headerTitles_.indexOf(HEADER_ADDUCTS)) {
	        	adducts = new ArrayList<AdductVO>();
	        	String[] split = rawValue.split(ADDUCT_REGEX);
	        	for (int i=0;i<split.length;i++)
	        	{
	        		adducts.add(new AdductVO(split[i].replace(ADDUCT_START, "").replace(ADDUCT_END, "")));
	        	}
	        } else if (index == headerTitles_.indexOf(HEADER_RETENTION_TIME)) {
	        	retentionTime = Double.parseDouble(rawValue);
	        } else if (index == headerTitles_.indexOf(HEADER_TOLERANCE)) {
	        	tolerance = Double.parseDouble(rawValue);
	        }
	      }
	      if (lipidClass != null && sumFormula != null && adducts != null && retentionTime != null && tolerance != null)
	      {
	      	targetListEntries_.add(new FragTargetListEntry(lipidClass, species, sumFormula, adducts, retentionTime, tolerance));
	      }
	      else
	      {
	      	throw new IOException(String.format("The content row number %s does not contain all required entries! \n", row.getRowNum()));
	      }
	    }
		}
		else
		{
			throw new IOException("The target list does not contain all required headers.");
		}
  }
  
  /**
   * Parses the header of an Excel sheet
   * @param headerRow 	row number of the header
   * @return List of the header titles
   */
  private List<String> readSheetHeaderTitles(Row headerRow) {
    try (Stream<Cell> cells = headerRow.stream();) {
      return cells.map((c) -> (!(c==null || c.getType().equals(CellType.ERROR)) ? c.getText() : "null")).collect(Collectors.toList());
    } 
  }
  
  /**
   * Collects all cells of a row that are neither null nor contain erroneous values.
   * @param row
   * @return
   */
  private List<Cell> collectRowCells(Row row)
  {
  	return row.stream().filter((c) -> !(c==null || c.getType().equals(CellType.ERROR))).collect(Collectors.toList());
  }
  
  public ArrayList<FragTargetListEntry> getTargetListEntries()
  {
  	return this.targetListEntries_;
  }
  
  
}