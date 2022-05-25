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
import java.util.Hashtable;
import java.util.Vector;
import java.util.stream.Stream;

import javax.swing.JFrame;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.dhatim.fastexcel.reader.ReadableWorkbook;

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

  public void parseInput() throws ExcelInputFileException, LipidCombinameEncodingException{
    extractInformation();
  }

  protected void extractInformation() throws ExcelInputFileException, LipidCombinameEncodingException{
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
  	//for backwards compatibility, in case there are files in the old excel format
  	else if (resultFile.getAbsolutePath().endsWith(".xls")) {
      parseResultFileApachePOI(resultFile);
    } else {
    	new WarningMessage(new JFrame(), "ERROR", "The specified file format is not supported!");
      throw new ExcelInputFileException("The specified file format is not supported!");
    }
  }
  
  protected void parseResultFileFastExcel(File resultFile) throws ExcelInputFileException, LipidCombinameEncodingException{
  	try (InputStream is = new FileInputStream(resultFile);
        ReadableWorkbook wb = new ReadableWorkbook(is);
        Stream<org.dhatim.fastexcel.reader.Sheet> sheets = wb.getSheets();) 
  	{
  		sheets.filter((s) -> (!s.getName().equals(QuantificationResultExporter.SHEET_CONSTANTS)&&
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
  	catch (IOException ex)
  	{
      ex.printStackTrace();
      new WarningMessage(new JFrame(), "ERROR", ex.getMessage());
      throw new ExcelInputFileException(ex);
    }
  }
  
  @Deprecated
  //slower Apache POI reader is only used for old xls excel files
  //TODO: remove completely after an adequate transition period *note written: 25.05.2022*
  //(1 year should suffice as xls files are not written since a few years already and newer LDA versions do not support such old file formats anymore anyway) 
  protected void parseResultFileApachePOI(File resultFile) throws ExcelInputFileException, LipidCombinameEncodingException{
    try (InputStream myxls = new FileInputStream(resultFile);
        Workbook workbook = new HSSFWorkbook(myxls);)
    {
      for (int sheetNumber=0;sheetNumber!=workbook.getNumberOfSheets();sheetNumber++){
        Sheet sheet = workbook.getSheetAt(sheetNumber);
        String className = sheet.getSheetName();
        if (!className.equals(QuantificationResultExporter.SHEET_CONSTANTS)&&
            !className.endsWith(QuantificationResultExporter.ADDUCT_OMEGA_SHEET)&&
            !className.endsWith(QuantificationResultExporter.ADDUCT_OVERVIEW_SHEET)&&
            !className.endsWith(QuantificationResultExporter.ADDUCT_MSN_SHEET)){
          if (!classesHash_.containsKey(className)){
            lipidClasses_.add(className);
            classesHash_.put(className, className);
          }
        }
      }
    }
    catch (IOException e) {
      e.printStackTrace();
    }

  }

  public Vector<String> getLipidClasses(){
    return this.lipidClasses_;
  }
}
