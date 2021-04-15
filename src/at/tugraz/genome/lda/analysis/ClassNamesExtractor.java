/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.QuantificationThread;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;

/**
 * 
 * @author Juergen Hartler
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
    try {
      InputStream myxls = new FileInputStream(resultFile);
      Workbook workbook = null;
      if (resultFile.getAbsolutePath().endsWith(".xlsx")) workbook = new XSSFWorkbook(myxls);
      else if (resultFile.getAbsolutePath().endsWith(".xls")) workbook     = new HSSFWorkbook(myxls);
      for (int sheetNumber=0;sheetNumber!=workbook.getNumberOfSheets();sheetNumber++){
        Sheet sheet = workbook.getSheetAt(sheetNumber);
        String className = sheet.getSheetName();
        if (!className.contains("Overview")&&!className.endsWith(QuantificationThread.OVERVIEW_SHEET_ADDUCT)&&!className.equalsIgnoreCase(QuantificationThread.CONSTANTS_SHEET)&&!className.endsWith(QuantificationThread.MSN_SHEET_ADDUCT)){
          if (!classesHash_.containsKey(className)){
            lipidClasses_.add(className);
            classesHash_.put(className, className);
          }
        }
      }
      myxls.close();
    }
    catch (IOException e) {
      e.printStackTrace();
    }

  }

  public Vector<String> getLipidClasses(){
    return this.lipidClasses_;
  }
}
