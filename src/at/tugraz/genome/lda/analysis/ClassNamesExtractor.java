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

package at.tugraz.genome.lda.analysis;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Collection;
import java.util.Hashtable;
import java.util.List;
import java.util.Optional;
import java.util.TreeMap;
import java.util.Vector;
import java.util.stream.Stream;

import javax.swing.JFrame;

import org.dhatim.fastexcel.reader.ReadableWorkbook;
import org.dhatim.fastexcel.reader.Row;
import org.dhatim.fastexcel.reader.Sheet;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.export.QuantificationResultExporter;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
public class ClassNamesExtractor
{

  protected Vector<File> resultFiles_;
  protected Vector<String> lipidClasses_;
  private Hashtable<String,String> classesHash_;

  
  public ClassNamesExtractor(Vector<File> resultFiles){
    this.resultFiles_ = resultFiles;
  }

  public void parseInput(int statisticsViewMode, boolean combineOxWithNonOx) throws ExcelInputFileException, LipidCombinameEncodingException{
    extractInformation(statisticsViewMode, combineOxWithNonOx);
  }

  protected void extractInformation(int statisticsViewMode, boolean combineOxWithNonOx) throws ExcelInputFileException, LipidCombinameEncodingException{
    lipidClasses_ = new Vector<String>();
    classesHash_ = new Hashtable<String,String>();
    for (int i=0; i!=resultFiles_.size();i++){
      File resultFile = resultFiles_.get(i);
      parseResultFile(resultFile);
    }
  }
  
  protected void parseResultFile(File resultFile) throws ExcelInputFileException, LipidCombinameEncodingException{
  	if (resultFile.getAbsolutePath().endsWith(".xlsx")){
    	parseResultFileFastExcel(resultFile);
    } 
  	else 
  	{
    	new WarningMessage(new JFrame(), "ERROR", "The specified file format is not supported!");
      throw new ExcelInputFileException("The specified file format is not supported!");
    }
  }
  
  protected void parseResultFileFastExcel(File resultFile) throws ExcelInputFileException, LipidCombinameEncodingException{
  	try (InputStream is = new FileInputStream(resultFile);
        ReadableWorkbook wb = new ReadableWorkbook(is);) 
  	{
  		Optional<Sheet> lipidClassSheet = wb.getSheets().filter((s) -> (s.getName().equals(QuantificationResultExporter.SHEET_LIPID_CLASS_LOOKUP))).findFirst();
  		if (lipidClassSheet.isPresent())
  		{
  			TreeMap<String,String> map = parseLipidClassSheet(lipidClassSheet.get());
  			Collection<String> lipidClasses = map.values();
  			lipidClasses_.addAll(lipidClasses);
  			for (String lipidClass : lipidClasses)
  			{
  				classesHash_.put(lipidClass, lipidClass);
  			}
  		}
  		else
  		{
  			wb.getSheets().filter((s) -> (!s.getName().equals(QuantificationResultExporter.SHEET_CONSTANTS)&&
										          !s.getName().endsWith(QuantificationResultExporter.ADDUCT_OMEGA_SHEET)&&
										          !s.getName().endsWith(QuantificationResultExporter.ADDUCT_OVERVIEW_SHEET)&&
										          !s.getName().endsWith(QuantificationResultExporter.ADDUCT_MSN_SHEET)))
							.forEach((s) -> {
								if (!classesHash_.containsKey(s.getName())) 
								{
									lipidClasses_.add(s.getName());
									classesHash_.put(s.getName(), s.getName());
								}
							});  
  		} 
  	}
  	catch (IOException ex)
  	{
      ex.printStackTrace();
      new WarningMessage(new JFrame(), "ERROR", ex.getMessage());
      throw new ExcelInputFileException(ex);
    }
  }
  
  public static TreeMap<String,String> parseLipidClassSheet(Sheet sheet) throws IOException
  {
  	TreeMap<String,String> map = new TreeMap<String,String>();
  	List<Row> rows = sheet.read();
    List<Row> contentRows = rows.subList(1, rows.size());

    for (Row row : contentRows) {
      
      String key = row.getCellText(0);
      String lipidClass = row.getCellText(1);
      map.put(key, lipidClass);
    }
    return map;
  }

  public Vector<String> getLipidClasses(){
    return this.lipidClasses_;
  }
}
