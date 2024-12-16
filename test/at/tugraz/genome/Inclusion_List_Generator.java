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
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.Vector;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.util.Precision;
import org.dhatim.fastexcel.Workbook;
import org.dhatim.fastexcel.Worksheet;
import org.dhatim.fastexcel.reader.Cell;
import org.dhatim.fastexcel.reader.CellType;
import org.dhatim.fastexcel.reader.ReadableWorkbook;
import org.dhatim.fastexcel.reader.Row;
import org.dhatim.fastexcel.reader.Sheet;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.utils.StaticUtils;

public class Inclusion_List_Generator
{
//	private final static String MASS_LIST_PATH = "D:\\Collaborator_Files\\SILDA\\massLists\\masslists_all_labels";
//	private final static String MZ_VALUES_FILE_PATH = "D:\\Collaborator_Files\\SILDA\\massLists\\masslists_all_labels\\all_mz_values_for_inclusion_list.xlsx";
	private final static String MASS_LIST_PATH = "D:\\Collaborator_Files\\SILDA\\massLists\\masslists_all_labels\\less_compounds";
	private final static String MZ_VALUES_FILE_PATH = "D:\\Collaborator_Files\\SILDA\\massLists\\masslists_all_labels\\less_compounds\\less_mz_values_for_inclusion_list.xlsx";
	private final static int HEADER_ROW_IN = 0;
	private final static int HEADER_ROW_OUT = 0;
	private final static int DECIMALS = 5;
	private final static String FORMATE_MZ = "mass(form[+HCOO] name[HCOO])";
	private final static String DEPROTONATED_MZ = "mass(form[-H] name[-H])";
	
	private Set<Double> mzValues_ = new HashSet<Double>();
	
	public void loadMassLists()
	{
		Vector<File> resultFiles = new Vector<File>();
		File folder = new File(MASS_LIST_PATH);
		if (folder.exists() && folder.isDirectory())
		{
			File[] fileCandidates = folder.listFiles();
			for (int i=0; i!=fileCandidates.length;i++)
			{
				if (fileCandidates[i].isFile())
				{
					String fileName = StaticUtils.extractFileName(fileCandidates[i].getAbsolutePath());
					String suffix = fileName.substring(fileName.lastIndexOf(".")+1);
					if (suffix.equalsIgnoreCase("xlsx"))
					{
						resultFiles.add(fileCandidates[i]);
					}
				}
			}
		}
		else
		{
			System.out.println("This is not a folder!");
		}
		parseMassLists(resultFiles);
	}
	
	private void parseMassLists(Vector<File> resultFiles)
	{
		for (File resultFile : resultFiles)
		{
			try (InputStream is = new FileInputStream(resultFile);
	        ReadableWorkbook wb = new ReadableWorkbook(is);
	        Stream<Sheet> sheets = wb.getSheets();) {
				sheets.forEach((s) -> {
					mzValues_.addAll(readSheet(s));
				});
	    } catch (IOException ex){
	      ex.printStackTrace();
	    }
		}
	}
	
	private Vector<Double> readSheet(Sheet sheet)
	{
		Vector<Double> mzValues = new Vector<Double>();
		List<Row> rows = null;
    try {
      rows = sheet.read();
    } catch (IOException ex) {}
    int headerRowNum = HEADER_ROW_IN;
    Row headerRow = rows.get(headerRowNum);
    List<String> headerTitles = readSheetHeaderTitles(headerRow);
    while (!(headerTitles.contains(FORMATE_MZ) || headerTitles.contains(DEPROTONATED_MZ)))
    {
    	headerRowNum++;
    	headerRow = rows.get(headerRowNum);
    	headerTitles = readSheetHeaderTitles(headerRow);
    }
    List<Row> contentRows = rows.subList(headerRowNum+1, rows.size());
    int index;
    String rawValue;
    for (Row row : contentRows) {
      List<Cell> cells = row.stream().filter((c) -> !(c==null || c.getType().equals(CellType.ERROR))).collect(Collectors.toList());
      for (Cell cell : cells) {
        index = cell.getColumnIndex();
        rawValue = cell.getRawValue();
        if (rawValue == null) continue;
        
        if (index == headerTitles.indexOf(FORMATE_MZ)) {
        	mzValues.add(Precision.round(Double.parseDouble(rawValue), DECIMALS));
        }
        else if (index == headerTitles.indexOf(DEPROTONATED_MZ)) {
        	mzValues.add(Precision.round(Double.parseDouble(rawValue), DECIMALS));
        }
      }
    }
    return mzValues;   
	}
	
	/**
   * Parses the header of an Excel sheet
   * @param headerRow row number of the header
   * @return List of the header titles
   */
  private List<String> readSheetHeaderTitles(Row headerRow) {
    try (Stream<Cell> cells = headerRow.stream();) {
      return cells.map((c) -> (!(c==null || c.getType().equals(CellType.ERROR)) ? c.getText() : "null")).collect(Collectors.toList());
    } 
  } 
  
  public void writeMzValuesToExcel()
  {
  	try (BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(MZ_VALUES_FILE_PATH));) {
      String s = Settings.VERSION;
      //the constructor can only take a version number in the format xx.yyyy
      Workbook wb = new Workbook(out, "Lipid Data Analyzer", s.substring(0, s.indexOf(".", s.indexOf(".") + 1)));
      Worksheet ws = wb.newWorksheet("mz_values");
      List<String> headerTitles = new ArrayList<String>();
      headerTitles.add("m/z values");
      createHeader(ws, headerTitles);
      int row = HEADER_ROW_OUT+1;
      for (Double value : mzValues_)
      {
      	ws.value(row, 0, value);
      	row++;
      }
      wb.finish();
  	}
  	catch (Exception ex)
  	{
  		ex.printStackTrace();
  	}
  }
  
  /**
   * Creates a formatted header row with the given header titles in the row number given by HEADER_ROW_OUT
   * @param ws Excel worksheet to write the header to
   * @param headerTitles List of header titles
   */
  private static void createHeader(Worksheet ws, List<String> headerTitles) {
    for (int i=0; i<headerTitles.size(); i++) {
      ws.value(HEADER_ROW_OUT, i, headerTitles.get(i));
    }
    ws.range(HEADER_ROW_OUT, HEADER_ROW_OUT, HEADER_ROW_OUT, headerTitles.size()).style()
    .bold().horizontalAlignment("center").fontName("Arial").fontSize(12).set();
  }
}
