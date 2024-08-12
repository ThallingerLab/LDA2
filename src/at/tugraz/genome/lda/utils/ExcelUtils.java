/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.utils;

import java.util.Hashtable;

import org.apache.poi.hssf.usermodel.HSSFCellStyle;
import org.apache.poi.hssf.usermodel.HSSFFont;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFFont;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
public class ExcelUtils
{
  /** bolder letters are broader - this factor accounts for this fact*/
  public final static double BOLD_MULT = 1.4d;
  /** the Excel width is declared by (# of letters)/256 - this value can be used to calculate an ideal cell width*/
  public final static int CHAR_MULT = 256;
  
  public final static String EXCEL_TEMP_PREFIX = "~$";
  
  /**
   * parses an Excel row and returns the entries in an hash, where the key is the column index - if numeric: Double is returned, else: String
   * @param row Excel row to be parsed
   * @stringOnly shall only the string values be returned
   * @return hash, where the key is the column index - if numeric: Double is returned, else: String
   */
  public static Hashtable<Integer,Object> getEntriesOfOneRow(Row row, boolean stringOnly){
    Hashtable<Integer,Object> rowEntries = new Hashtable<Integer,Object>();
    for (int i=0;  row!=null && i!=(row.getLastCellNum()+1);i++){
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
        if (contents.endsWith(".0"))contents = contents.substring(0,contents.length()-".0".length());
      }
      if (contents.length()>0){
        if (numeric!=null && !stringOnly) rowEntries.put(i, numeric);
        else rowEntries.put(i, contents);
      }
    }
    return rowEntries;
  }
  
  public static XSSFCellStyle getMassListHeaderStyle(XSSFWorkbook wb)
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
  
  public static XSSFCellStyle getMassListNumberStyle(XSSFWorkbook wb)
  {
  	XSSFCellStyle numberStyle = wb.createCellStyle();
  	numberStyle.setDataFormat(2);
  	return numberStyle;
  }


}
