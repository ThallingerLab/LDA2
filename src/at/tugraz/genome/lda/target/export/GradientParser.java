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
package at.tugraz.genome.lda.target.export;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.dhatim.fastexcel.reader.Cell;
import org.dhatim.fastexcel.reader.CellType;
import org.dhatim.fastexcel.reader.ReadableWorkbook;
import org.dhatim.fastexcel.reader.Row;
import org.dhatim.fastexcel.reader.Sheet;

import at.tugraz.genome.lda.parser.LDAResultReader;
import javafx.util.Pair;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class GradientParser
{
	private final static String HEADER_TIME = "Adjustment time /min";
	private final static String HEADER_FACTOR = "Adjustment Factor";
	private final static Integer HEADER_ROW = 0;
	
	protected static GradientAdjuster parseGradient(File gradientFile) throws Exception
	{
		try (InputStream is = new FileInputStream(gradientFile.getAbsolutePath());
        ReadableWorkbook wb = new ReadableWorkbook(is);) 
		{
			return new GradientAdjuster(parseGradientSheet(wb.getFirstSheet()));
		}
	}
	
	/**
	 * Parses an excel sheet specifying details about a gradient. 
	 * HEADER_TIME and HEADER_FACTOR must be present.
	 * @param sheet
	 * @return
	 * @throws IOException
	 */
	private static ArrayList<Pair<Double,Double>> parseGradientSheet(Sheet sheet) throws IOException
	{
		ArrayList<Pair<Double,Double>> dataPoints = new ArrayList<Pair<Double,Double>>();
		List<Row> rows = null;
    rows = sheet.read();
    Row headerRow = rows.get(HEADER_ROW);
    List<String> headerTitles = LDAResultReader.readSheetHeaderTitles(headerRow);
    List<Row> contentRows = rows.subList(HEADER_ROW+1, rows.size());
    int index;
    String rawValue;
    Double time;
  	Double factor;
    
    for (Row row : contentRows) 
    {
    	time = null;
    	factor = null;
    	List<Cell> cells = row.stream().filter((c) -> !(c==null || c.getType().equals(CellType.ERROR))).collect(Collectors.toList());
      for (Cell cell : cells) {
        index = cell.getColumnIndex();
        rawValue = cell.getRawValue();
        
        if (index == headerTitles.indexOf(HEADER_TIME))
        {
        	time = Double.parseDouble(rawValue);
        }
        else if (index == headerTitles.indexOf(HEADER_FACTOR))
        {
        	factor = Double.parseDouble(rawValue);
        }
      }
      if (time != null && factor != null)
      {
      	dataPoints.add(new Pair<Double,Double>(time,factor));
      }
    }
    return dataPoints;
	}
  
}
