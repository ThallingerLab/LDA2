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
package at.tugraz.genome;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Vector;

import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.LDAResultReader;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LDAToFaConverter extends ConverterBase
{
  
  private String fileName_;
  
  public LDAToFaConverter(String fileName){
    fileName_ = fileName;
  }
  
  public void convert() throws ExcelInputFileException, FileNotFoundException, LipidCombinameEncodingException{
    QuantificationResult returnParam = LDAResultReader.readResultFile(fileName_, new Hashtable<String,Boolean>());
    if (returnParam==null){
      System.out.println("It was not possible to translate the file: "+fileName_);
      return;
    }
    Hashtable<String,Vector<LipidParameterSet>> results = returnParam.getIdentifications();
    if (results == null || results.keySet().size()==0){
      System.out.println("The file does not contain any data: "+fileName_);
      return;      
    }
    boolean allZero = true;
    for (Vector<LipidParameterSet> sets : results.values()){
      if (sets.size()!=0){
        allZero = false;
        break;
      }
    }
    if (allZero){
      System.out.println("The file does not contain any data: "+fileName_);
      return;      
    }
    String outFilename = fileName_.substring(0,fileName_.lastIndexOf("."))+"_FA"+fileName_.substring(fileName_.lastIndexOf("."));
    BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outFilename));
    try{
      Workbook resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      for (String className : results.keySet()){
        Sheet sheet = resultWorkbook.createSheet(className);
        int rowCount = 0;
        Row row = sheet.createRow(rowCount);
        rowCount++;

        int longestSpecies = 0;
        int longestAdduct = 0;
        int longestRt = 0;
        int longestSpeciesIntensity = 0;
        int longestChain = 0;
        int longestChainPercent = 0;
        int longestChainIntensity = 0;
        int longestMolecularSpecies = 0;
        
        this.createCell(row, headerStyle, COLUMN_LIPID_SPECIES, TEXT_LIPID_SPECIES);
        this.createCell(row, headerStyle, COLUMN_ADDUCT, TEXT_ADDUCT);
        this.createCell(row, headerStyle, COLUMN_RT, TEXT_RT);
        this.createCell(row, headerStyle, COLUMN_LIPID_SPECIES_INTENSITY, TEXT_LIPID_SPECIES_INTENSITY);
        this.createCell(row, headerStyle, COLUMN_CHAIN, TEXT_CHAIN);
        this.createCell(row, headerStyle, COLUMN_CHAIN_PERCENT, TEXT_CHAIN_PERCENT);
        this.createCell(row, headerStyle, COLUMN_CHAIN_INTENSITY, TEXT_CHAIN_INTENSITY);
        this.createCell(row, headerStyle, COLUMN_MOLECULAR_SPECIES, TEXT_MOLECULAR_SPECIES);

        for (LipidParameterSet set : results.get(className)){
          if (set.getIsotopicProbes().size()==0) continue;
          String displayName = set.getNameStringWithoutRt();
          String rt = set.getRt();
          float ms1Area = getMS1Area(set);
          String ms1AreaString = String.valueOf(ms1Area);
          String mod = String.valueOf(set.getModificationName());
          if (StaticUtils.isThereChainInformationAvailable(displayName, set)){
            LipidomicsMSnSet msn = (LipidomicsMSnSet) set;
            LinkedHashMap<String,String[]> faDetails = extractFaIntensityDetails(msn, ms1AreaString, ms1Area);
            for (String fa : faDetails.keySet()){
              String[] areaResults = faDetails.get(fa);
              String percent = areaResults[1];
              String faAreaString = areaResults[2];
              String combiName = areaResults[3];
              row = sheet.createRow(rowCount);
              rowCount++;
              createCell(row, null, COLUMN_LIPID_SPECIES, displayName);
              if (displayName.length()>longestSpecies) longestSpecies = displayName.length();
              createCell(row, null, COLUMN_ADDUCT, mod);
              if (mod.length()>longestAdduct) longestAdduct = mod.length();
              createNumericCell(row, null, COLUMN_RT, rt);
              if (rt.length()>longestRt) longestRt = rt.length();
              createNumericCell(row, null, COLUMN_LIPID_SPECIES_INTENSITY, ms1AreaString);
              if (ms1AreaString.length()>longestSpeciesIntensity) longestSpeciesIntensity = ms1AreaString.length();
              createCell(row, null, COLUMN_CHAIN, fa);
              if (fa.length()>longestChain) longestChain = fa.length();
              createNumericCell(row, null, COLUMN_CHAIN_PERCENT, percent);
              if (percent.length()>longestChainPercent) longestChainPercent = percent.length();
              createNumericCell(row, null, COLUMN_CHAIN_INTENSITY, faAreaString);
              if (faAreaString.length()>longestChainIntensity) longestChainIntensity = faAreaString.length();
              createCell(row, null, COLUMN_MOLECULAR_SPECIES, combiName);
              if (combiName.length()>longestMolecularSpecies) longestMolecularSpecies = combiName.length();
            }            
          } else {
            row = sheet.createRow(rowCount);
            rowCount++;
            createCell(row, null, COLUMN_LIPID_SPECIES, displayName);
            if (displayName.length()>longestSpecies) longestSpecies = displayName.length();
            createCell(row, null, COLUMN_ADDUCT, mod);
            if (mod.length()>longestAdduct) longestAdduct = mod.length();
            createNumericCell(row, null, COLUMN_RT, rt);
            if (rt.length()>longestRt) longestRt = rt.length();
            createNumericCell(row, null, COLUMN_LIPID_SPECIES_INTENSITY, ms1AreaString);
            if (ms1AreaString.length()>longestSpeciesIntensity) longestSpeciesIntensity = ms1AreaString.length();

          }
        }
        
        setColumnWidth(sheet, COLUMN_LIPID_SPECIES, TEXT_LIPID_SPECIES, longestSpecies);
        setColumnWidth(sheet, COLUMN_ADDUCT, TEXT_ADDUCT, longestAdduct);
        setColumnWidth(sheet, COLUMN_RT, TEXT_RT, longestRt);
        setColumnWidth(sheet, COLUMN_LIPID_SPECIES_INTENSITY, TEXT_LIPID_SPECIES_INTENSITY, longestSpeciesIntensity);
        setColumnWidth(sheet, COLUMN_CHAIN, TEXT_CHAIN, longestChain);
        setColumnWidth(sheet, COLUMN_CHAIN_PERCENT, TEXT_CHAIN_PERCENT, longestChainPercent);
        setColumnWidth(sheet, COLUMN_CHAIN_INTENSITY, TEXT_CHAIN_INTENSITY, longestChainIntensity);
        setColumnWidth(sheet, COLUMN_MOLECULAR_SPECIES, TEXT_MOLECULAR_SPECIES, longestMolecularSpecies);
      }
      resultWorkbook.write(out);
      resultWorkbook.close();
      out.close();
    } catch (IOException iox) {
      try {
        out.close();
      } catch (IOException e) {
        e.printStackTrace();
      }      
    } finally{
      try {
        out.close();
      }
      catch (IOException e) {
        e.printStackTrace();
      }
    }
  }
  
}
