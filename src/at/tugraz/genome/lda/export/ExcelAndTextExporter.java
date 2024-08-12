/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.export;

import java.awt.Color;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.util.Precision;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFColor;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.analysis.ComparativeResultsLookup;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RetentionTimeGroupingException;
import at.tugraz.genome.lda.export.vos.SpeciesExportVO;
import at.tugraz.genome.lda.export.vos.SummaryVO;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.swing.LipidomicsTableCellRenderer;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.ExportOptionsVO;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;

/**
 * Class for exporting and writing data from the heat map to Excel or tab-delimited TXT format
 * 
 * @author Juergen Hartler
 *
 */
public class ExcelAndTextExporter extends LDAExporter
{
  
  /**
   * central static method for executing the export
   * @param includeResultFiles true when the original result files shall be read (in some cases not possible, e.g. in the classes overview)
   * @param speciesType structural level of data export (lipid species, chain level, position level - for details see LipidomicsConstants.EXPORT_ANALYTE_TYPE)
   * @param exportDoubleBondPositionsForClass true when double bond positions shall be exported for the lipid class of this analyte
   * @param sheetName analyte class
   * @param out the output stream for writing the results
   * @param excelFile true in case of an Excel file; false for a tab-delimited file
   * @param maxIsotope the highest number of allowed isotopes
   * @param molNames the names of the exported molecules
   * @param isRtGrouped true when the peaks were grouped by a certain retention time
   * @param isGrouped true when sample groups are exported
   * @param expIdNames the names of the experiments
   * @param expNames a lookup from the original names to the displayed ones in the heat map
   * @param expFullPaths a lookup from the original names to the full paths of the Excel files
   * @param expsOfGroup a sorted hash map containing the sample groups as key, and the abbreviated experiment names belonging to this group as values
   * @param results the values according to the heat map/bar chart selection
   * @param preferredUnit a description of the physical unit
   * @param headerText the text written in the header of the exported file
   * @param expVO value object specifying which values shall be exported
   * @param compLookup lookup to the comparative values over several LDA results (used in Statistical Analysis section)
   * @param modifications the modifications for this analyte class
   * @throws ExportException when there is something wrong
   * @throws SpectrummillParserException when there are elements missing in the elementconfig.xml
   * @throws ExcelInputFileException when there is something wrong with the LDA results
   * @throws IOException when there is a general IO problem
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static void exportToFile(boolean includeResultFiles, short speciesType, boolean exportDoubleBondPositionsForClass, String sheetName, OutputStream out, boolean excelFile, int maxIsotope, Vector<String> molNames, 
      boolean isRtGrouped, boolean isGrouped, Vector<String> expIdNames, Hashtable<String,String> expNames, LinkedHashMap<String,String> expFullPaths,
      LinkedHashMap<String,Vector<String>> expsOfGroup, Hashtable<String,Hashtable<String,Vector<Double>>> results,
      String preferredUnit, String headerText, ExportOptionsVO expVO, ComparativeResultsLookup compLookup, ArrayList<String> modifications)
          throws ExportException,
      SpectrummillParserException, ExcelInputFileException, IOException, LipidCombinameEncodingException, RetentionTimeGroupingException{
    
    //first read the original Excel results of this analyte class
    Hashtable<String,QuantificationResult> originalExcelResults = new Hashtable<String,QuantificationResult>();
    HydroxyEncoding faEncoding = null;
    HydroxyEncoding lcbEncoding = null;
    for (String expId : expFullPaths.keySet()){
      QuantificationResult result = LDAResultReader.readResultFile(expFullPaths.get(expId),  new Hashtable<String,Boolean>(), sheetName);
      if (result.getFaHydroxyEncoding()!=null) faEncoding = result.getFaHydroxyEncoding();
      if (result.getLcbHydroxyEncoding()!=null) lcbEncoding = result.getLcbHydroxyEncoding();
      originalExcelResults.put(expId, result);
    }
    LinkedHashMap<String,Boolean> adductsSorted = extractAdductsSortedByAbundance(sheetName,originalExcelResults);
    
    //extract the values according to the molecular species
    boolean rtGroupingUsed = isRtGrouped;
    LinkedHashMap<String,SummaryVO> molSpeciesDetails = new LinkedHashMap<String,SummaryVO>();
    //Quantification
    for (String molName : molNames){
      Hashtable<String,Hashtable<String,Vector<LipidParameterSet>>> relevantOriginals = null;
      if (includeResultFiles){
        relevantOriginals = new Hashtable<String,Hashtable<String,Vector<LipidParameterSet>>>();
        for (String expId : expFullPaths.keySet()){
          ResultAreaVO areaVO = compLookup.getResultAreaVO(sheetName,molName,expId);
          if (!isGrouped && !expNames.containsKey(expId))
            continue;
          if (areaVO!=null){
            if (areaVO.isAStandard())
              rtGroupingUsed = false;
            Hashtable<String,Vector<LipidParameterSet>> sets = getRelevantOriginalResults(originalExcelResults.get(expId).getIdentifications().get(sheetName),areaVO,false);
            relevantOriginals.put(expId, sets);
          }
        }
      }
      SpeciesExportVO exportVO = LDAExporter.extractExportableSummaryInformation(speciesType, exportDoubleBondPositionsForClass, false, 0, 0, false, 0, 0, rtGroupingUsed, adductsSorted, new Vector<String>(expFullPaths.keySet()),
        expsOfGroup,  molName, results.get(molName), relevantOriginals, maxIsotope, faEncoding, lcbEncoding);
      for (SummaryVO sumVO : exportVO.getSummaries()){
        String id = "";
        //only for the classes overview, the name of the analyte class is not added
        if (includeResultFiles)
          id += sheetName+" ";
        id += sumVO.getSpeciesId();
        if (sumVO.getMolecularId()!=null)
          id += " | "+sumVO.getMolecularId();
        molSpeciesDetails.put(id, sumVO);
      }
    }
    try{
      writeToFile(sheetName,out,excelFile,molSpeciesDetails, isGrouped, expIdNames, expNames, expsOfGroup, preferredUnit, expVO,modifications);
    }finally{
      if (out!=null)
      out.close();
    }
  }
  
  /**
   * writes the prepared results finally to an Excel or tab-delimited format
   * @param sheetName the name of the sheet (i.e. the lipid class)
   * @param out the output stream to write to
   * @param excelFile true in case of an Excel file; false for a tab-delimited file
   * @param molSpeciesDetails the details for every line
   * @param isGrouped true when sample groups are exported
   * @param expIdNames the names of the experiments
   * @param expNames a lookup from the original names to the displayed ones in the heat map
   * @param expsOfGroup a sorted hash map containing the sample groups as key, and the abbreviated experiment names belonging to this group as values
   * @param preferredUnit a description of the physical unit
   * @param expVO value object specifying which values shall be exported
   * @param modifications the modifications for this analyte class
   * @throws IOException when there is a general IO problem
   */
  private static void writeToFile(String sheetName, OutputStream out, boolean excelFile, LinkedHashMap<String,SummaryVO> molSpeciesDetails,
      boolean isGrouped, Vector<String> expIdNames, Hashtable<String,String> expNames, LinkedHashMap<String,Vector<String>> expsOfGroup,
      String preferredUnit, ExportOptionsVO expVO, ArrayList<String> modifications) throws IOException{
    Workbook workbook = null;
    Sheet sheet = null;
    Row row = null;
    
    CellStyle headerStyle = null;
    CellStyle normalStyle = null;
    CellStyle zeroStyle = null;
    CellStyle ms2Unambiguous = null;
    CellStyle ms2Split = null;
    CellStyle ms2NoSplitPossible = null;
    
    if (excelFile){
      workbook = new XSSFWorkbook();
      sheet = workbook.createSheet(sheetName);
      headerStyle = getHeaderStyle(workbook);
      normalStyle = getNormalStyle(workbook);
      zeroStyle = getZeroStyle(workbook);
      ms2Unambiguous = getMs2UnambiguousStyle((XSSFWorkbook)workbook);
      ms2Split = getMs2SplitStyle((XSSFWorkbook)workbook);
      ms2NoSplitPossible = getMs2NoSplitPossibleStyle((XSSFWorkbook)workbook);
      row = sheet.createRow(0);
    }
    boolean expInColumn = false;
    Vector<String> columnNames = new Vector<String>(molSpeciesDetails.keySet());
    Vector<String> rowNames = new Vector<String>(expIdNames);
    if (!expVO.isAnalyteInColumn()){
      columnNames = new Vector<String>(expIdNames);
      rowNames = new Vector<String>(molSpeciesDetails.keySet());
      expInColumn = true;
    } 
    StaticUtils.printOutLabelToFile(sheetName, out, row, headerStyle, (columnNames.size()+2)/2, excelFile, false);
    if (!excelFile)
      out.write("\n".getBytes());   
    int rowCount = 1; 
    int columnCount = 0;
    Hashtable<Integer,String> headerValues = new Hashtable<Integer,String>();
    Hashtable<Integer,Integer> longestValues = new Hashtable<Integer,Integer>();
    String unit = "Arb.Units / AU";
    if (preferredUnit!=null && preferredUnit.length()>0)
      unit = preferredUnit;
    String unitString = unit;
    if (excelFile)
      row = sheet.createRow(rowCount);
    StaticUtils.printOutLabelToFile(unitString, out, row, headerStyle, columnCount, excelFile, false);
    headerValues.put(columnCount, unitString);
    longestValues.put(columnCount, 0);
    columnCount++;
    StaticUtils.printOutLabelToFile("", out, row, headerStyle, columnCount, excelFile, true);
    headerValues.put(columnCount, "");
    longestValues.put(columnCount, 0);
    columnCount++;
    for (String columnName : columnNames){
      String valueToPrint = columnName;
      if (expInColumn)
        valueToPrint = expNames.get(columnName);
      StaticUtils.printOutLabelToFile(valueToPrint, out, row, headerStyle, columnCount, excelFile, true);
      headerValues.put(columnCount, valueToPrint);
      longestValues.put(columnCount, 0);
      columnCount++;
    }
    rowCount++;
    if (excelFile)
      row = sheet.createRow(rowCount);
    if (!excelFile)
      out.write("\n".getBytes());
    double area;
    for (String rowName : rowNames){
      columnCount = 0;
      String valueToPrint = rowName;
      if (!expInColumn)
        valueToPrint = expNames.get(rowName);
      StaticUtils.printOutLabelToFile(valueToPrint, out, row, headerStyle, columnCount, excelFile, false);
      if (valueToPrint.length()>headerValues.get(columnCount).length())
        headerValues.put(columnCount,valueToPrint);
      columnCount++;
      StaticUtils.printOutLabelToFile("Value", out, row, headerStyle, columnCount, excelFile, true);
      if ("Value".length()>headerValues.get(columnCount).length())
        headerValues.put(columnCount,"Value");
      columnCount++;
      // write out the quantitative value
      for (String columnName : columnNames){
        String toWrite = "0";
        String molName = columnName;
        String expName = rowName;
        if (!expVO.isAnalyteInColumn()){
          molName = rowName;
          expName = columnName;
        }
        SummaryVO sumVO = molSpeciesDetails.get(molName);
        if (isGrouped)
          area = sumVO.getMeanArea(expName);
        else
          area = sumVO.getArea(expName);
        toWrite = StaticUtils.extractDisplayValue(area,expVO.getCommaPositions());        
        CellStyle style = null;
        if (toWrite.equalsIgnoreCase("0")||toWrite.equalsIgnoreCase("0.0"))
          style = zeroStyle;
        else
          style = getCorrespondingCellStyle(sumVO.getEvidenceReliabilty(expName),normalStyle,
        		  ms2Unambiguous, ms2Split, ms2NoSplitPossible);
        StaticUtils.printOutNumberToFile(toWrite, out, row, style, columnCount, excelFile, true);
        if (toWrite.length()>longestValues.get(columnCount))
          longestValues.put(columnCount,toWrite.length());
        columnCount++;
      }    
      rowCount++;
      if (excelFile)
        row = sheet.createRow(rowCount);
      if (!excelFile)
        out.write("\n".getBytes());
      // write out the standard deviation
      if (expVO!=null && expVO.getExportType() != ExportOptionsVO.EXPORT_NO_DEVIATION){
        columnCount = 0;
        StaticUtils.printOutLabelToFile("", out, row, headerStyle, columnCount, excelFile, false);
        columnCount++;
        String title = "";
        if (expVO.getExportType() == ExportOptionsVO.EXPORT_SD_DEVIATION){
          title+="SD-"+expVO.getSdValue();
        } else if (expVO.getExportType() == ExportOptionsVO.EXPORT_SD_ERROR){
          title+="SE";
        }
        StaticUtils.printOutLabelToFile(title, out, row, headerStyle, columnCount, excelFile, true);
        if (title.length()>headerValues.get(columnCount).length())
          headerValues.put(columnCount,title);

        columnCount++;
        for (String columnName : columnNames){
          String toWrite = "-1.0";
          String molName = columnName;
          String expName = rowName;
          if (!expVO.isAnalyteInColumn()){
            molName = rowName;
            expName = columnName;
          }
          SummaryVO sumVO = molSpeciesDetails.get(molName);
          double value = sumVO.calculateDeviationValue(expVO, expsOfGroup.get(expName));
          if (!Double.isInfinite(value) && !Double.isNaN(value))
            toWrite = StaticUtils.extractDisplayValue(value);
          CellStyle style = null;
          if (toWrite.equalsIgnoreCase("-1")||toWrite.equalsIgnoreCase("-1.0"))
            style = zeroStyle;
          else
            style = getCorrespondingCellStyle(sumVO.getEvidenceReliabilty(expName),normalStyle,
          		  ms2Unambiguous, ms2Split, ms2NoSplitPossible);
          StaticUtils.printOutNumberToFile(toWrite, out, row, style, columnCount, excelFile, true);
          if (toWrite.length()>longestValues.get(columnCount))
            longestValues.put(columnCount,toWrite.length());
          columnCount++;
        }
        rowCount++;       
        if (!excelFile)
          out.write("\n".getBytes());
        else
          row = sheet.createRow(rowCount);
      }
      //export of retention times
      for (String modName : modifications){
        if (expVO!=null && expVO.isExportRT()){
          columnCount = 0;
          StaticUtils.printOutLabelToFile("", out, row, headerStyle, columnCount, excelFile, false);
          columnCount++;
          String rtName = "RT";
          if (modifications.size()>1)
            rtName+="-"+modName;
          StaticUtils.printOutLabelToFile(rtName, out, row, headerStyle, columnCount, excelFile, true);
          if (rtName.length()>headerValues.get(columnCount).length())
            headerValues.put(columnCount,rtName);

          columnCount++;
          for (String columnName : columnNames){
            String toWrite = "-1.0";
            String molName = columnName;
            String expName = rowName;
            if (!expVO.isAnalyteInColumn()){
              molName = rowName;
              expName = columnName;
            }
            SummaryVO sumVO = molSpeciesDetails.get(molName);
            Double rt;
            if (isGrouped)
              rt = sumVO.getMeanRetentionTime(modName, expName);
            else
              rt = sumVO.getRetentionTime(modName, expName);
            if (rt!=null && !Double.isInfinite(rt) && !Double.isNaN(rt))
                toWrite = StaticUtils.extractDisplayValue(rt);
            CellStyle style = null;
            if (toWrite.equalsIgnoreCase("-1")||toWrite.equalsIgnoreCase("-1.0"))
              style = zeroStyle;
            else
              style = getCorrespondingCellStyle(sumVO.getEvidenceReliabilty(expName),normalStyle,
            		  ms2Unambiguous, ms2Split, ms2NoSplitPossible);
            StaticUtils.printOutNumberToFile(toWrite, out, row, style, columnCount, excelFile, true);
            if (toWrite.length()>longestValues.get(columnCount))
              longestValues.put(columnCount,toWrite.length());
            columnCount++;
          }  
          rowCount++;
          
          if (excelFile)
            row = sheet.createRow(rowCount);
          else
            out.write("\n".getBytes());
        }
        //export of standard deviations of the retention times
        if (expVO!=null && expVO.isExportRTDev()){
          columnCount = 0;
          StaticUtils.printOutLabelToFile("", out, row, headerStyle, columnCount, excelFile, false);
          columnCount++;
          StaticUtils.printOutLabelToFile("RT-StDev", out, row, headerStyle, columnCount, excelFile, true);
          if ("RT-StDev".length()>headerValues.get(columnCount).length())
            headerValues.put(columnCount,"RT-StDev");

          columnCount++;
          for (String columnName : columnNames){
            String toWrite = "-1.0";
            String molName = columnName;
            String expName = rowName;
            if (!expVO.isAnalyteInColumn()){
              molName = rowName;
              expName = columnName;
            }
            SummaryVO sumVO = molSpeciesDetails.get(molName);
            Double rtStdev = sumVO.getStdevRetentionTime(modName, expName);
            if (rtStdev!=null && !Double.isInfinite(rtStdev) && !Double.isNaN(rtStdev))
              toWrite = StaticUtils.extractDisplayValue(rtStdev);
            CellStyle style = null;
            if (toWrite.equalsIgnoreCase("-1")||toWrite.equalsIgnoreCase("-1.0"))
              style = zeroStyle;
            else
              style = getCorrespondingCellStyle(sumVO.getEvidenceReliabilty(expName),normalStyle,
            	ms2Unambiguous, ms2Split, ms2NoSplitPossible);
            StaticUtils.printOutNumberToFile(toWrite, out, row, style, columnCount, excelFile, true);
            columnCount++;
          }  
          rowCount++;         
          if (excelFile)
            row = sheet.createRow(rowCount);           
          else
            out.write("\n".getBytes());
        }
      }      
    }  
    
    if (excelFile){
      for (int column : headerValues.keySet()){
        setColumnWidth(sheet, column, headerValues.get(column), longestValues.get(column));
      }
      workbook.write(out);
    }
  }
  
  /**
   * returns the style for conventional Excel cells
   * @param wb the Excel workbook
   * @return the style for conventional Excel cells
   */
  private static CellStyle getNormalStyle(Workbook wb){
    CellStyle arialStyle = wb.createCellStyle();
    org.apache.poi.ss.usermodel.Font arialfont = getArialFont(wb);
    arialStyle.setFont(arialfont);
    return arialStyle;
  }

  /**
   * returns the style for an Excel cell containing a zero
   * @param wb the Excel workbook
   * @return style for an Excel cell containing a zero
   */
  private static CellStyle getZeroStyle(Workbook wb){
    CellStyle arialStyle = wb.createCellStyle();
    org.apache.poi.ss.usermodel.Font arialfont = getArialRedFont(wb);
    arialStyle.setFont(arialfont);
    return arialStyle;
  }
  
  /**
   * returns the style for a header Excel cell
   * @param wb the Excel workbook
   * @return style for a header Excel cell
   */
  private static CellStyle getHeaderStyle(Workbook wb){
    CellStyle arial12style = wb.createCellStyle();
    org.apache.poi.ss.usermodel.Font arial12font = getArial12BoldFont(wb);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(CellStyle.ALIGN_CENTER);
    return arial12style;
  }
  
  /**
   * returns the style for an MS2 unambiguously identified Excel cell
   * @param wb the Excel workbook
   * @return style for a header Excel cell
   */
  private static CellStyle getMs2UnambiguousStyle(XSSFWorkbook wb){
    XSSFCellStyle arialStyleGreen = wb.createCellStyle();
    org.apache.poi.ss.usermodel.Font arialfont = getArialFont(wb);
    arialStyleGreen.setFont(arialfont);
    XSSFColor color = new XSSFColor(LipidomicsTableCellRenderer.BRIGHT_GREEN);
    arialStyleGreen.setFillForegroundColor(color);
    arialStyleGreen.setFillPattern(CellStyle.SOLID_FOREGROUND);
    return arialStyleGreen;
  }

  /**
   * returns the style for an MS2-split identified Excel cell
   * @param wb the Excel workbook
   * @return style for a header Excel cell
   */
  private static CellStyle getMs2SplitStyle(XSSFWorkbook wb) {
    XSSFCellStyle arialStyleYellow = wb.createCellStyle();
    org.apache.poi.ss.usermodel.Font arialfont = getArialFont(wb);
    arialStyleYellow.setFont(arialfont);
    XSSFColor color = new XSSFColor(Color.YELLOW);
    arialStyleYellow.setFillForegroundColor(color);
    arialStyleYellow.setFillPattern(CellStyle.SOLID_FOREGROUND);
    return arialStyleYellow;
  }

  /**
   * returns the style for an MS2 no-MS1-split-possible identified Excel cell
   * @param wb the Excel workbook
   * @return style for a header Excel cell
   */
  private static CellStyle getMs2NoSplitPossibleStyle(XSSFWorkbook wb){
    XSSFCellStyle arialStyleOrange = wb.createCellStyle();
    org.apache.poi.ss.usermodel.Font arialfont = getArialFont(wb);
    arialStyleOrange.setFont(arialfont);
    XSSFColor color = new XSSFColor(Color.ORANGE);
    arialStyleOrange.setFillForegroundColor(color);
    arialStyleOrange.setFillPattern(CellStyle.SOLID_FOREGROUND);
    return arialStyleOrange;	  
  }
  
  /**
   * returns an Arial 12 font in bold
   * @param wb the Excel workbook
   * @return an Arial 12 font in bold
   */
  private static org.apache.poi.ss.usermodel.Font getArial12BoldFont(Workbook wb){
    org.apache.poi.ss.usermodel.Font arial12font = getArialFont(wb);
    arial12font.setFontHeightInPoints((short)12);
    arial12font.setBoldweight(org.apache.poi.ss.usermodel.Font.BOLDWEIGHT_BOLD);
    return arial12font;
  }
  
  /**
   * returns the standard Arial font in red
   * @param wb the Excel workbook
   * @return Arial font in red
   */
  private static org.apache.poi.ss.usermodel.Font getArialRedFont(Workbook wb){
    org.apache.poi.ss.usermodel.Font arialfont = getArialFont(wb);
    arialfont.setColor(org.apache.poi.ss.usermodel.Font.COLOR_RED);
    return arialfont;
  }

  /**
   * returns the standard Arial font
   * @param wb the Excel workbook
   * @return standard Arial font
   */
  private static org.apache.poi.ss.usermodel.Font getArialFont(Workbook wb){
    org.apache.poi.ss.usermodel.Font arial12font = wb.createFont();
    arial12font.setFontName("Arial");
    return arial12font;
  }

  /**
   * sets the column width of an Excel cell according to the entries
   * @param sheet the Excel tab where the width should be set
   * @param column the column number
   * @param headerValue the longest header value entry
   * @param longestValue length of the longest cell entry in this column
   */
  private static void setColumnWidth(Sheet sheet, int column, String headerValue, int longestValue){
    int columnWidth = (int)((headerValue.length()*256)*ExcelUtils.BOLD_MULT);
    if ((longestValue+1)*256>columnWidth) columnWidth =  (longestValue+1)*256;
    sheet.setColumnWidth(column,columnWidth); 
  }
  
  /**
   * returns the corresponding evidence depending on the value of the found evidence
   * @param evidence the evidence encoding specified in the SummaryVO
   * @param normalStyle the standard cell style without any color
   * @param ms2Unambiguous style when ms2 found and no overlap
   * @param ms2Split style when ms2 is found, but there is an overlap with an isobar from another species, but the two hits can be separated
   * @param ms2NoSplitPossible style when ms2 is found, but there is an overlap with an isobar from another species, and the two hits cannot be separated*
   * @return
   */
  private static CellStyle getCorrespondingCellStyle(short evidence, CellStyle normalStyle, CellStyle ms2Unambiguous,
      CellStyle ms2Split, CellStyle ms2NoSplitPossible) {
    if (evidence==SummaryVO.EVIDENCE_MS2_UNAMBIGUOUS)
      return ms2Unambiguous;
    else if (evidence==SummaryVO.EVIDENCE_MS2_SPLIT)
      return ms2Split;
    else if (evidence==SummaryVO.EVIDENCE_MS2_NO_SPLIT_POSSIBLE)
      return ms2NoSplitPossible;
    else	
      return normalStyle;
  }
  
  
  
  
  
  
  //TODO: Starting here, this is just for SILDA analysis, remove when done
  
  
  
  public static void writeExcelSheetOmegaSummary(OmegaCollector omegaCollector, Sheet sheet, XSSFWorkbook workbook, OutputStream out, boolean includeResultFiles, short speciesType, boolean exportDoubleBondPositionsForClass, 
  		int maxIsotope, Vector<String> molNames, boolean isRtGrouped, boolean isGrouped, Vector<String> expIdNames, Hashtable<String,String> expNames, 
  		LinkedHashMap<String,String> expFullPaths, LinkedHashMap<String,Vector<String>> expsOfGroup, Hashtable<String,Hashtable<String,Vector<Double>>> results,
      String preferredUnit, String headerText, ExportOptionsVO expVO, ComparativeResultsLookup compLookup, ArrayList<String> modifications)
          throws ExportException,
      SpectrummillParserException, ExcelInputFileException, IOException, LipidCombinameEncodingException, RetentionTimeGroupingException{

    //first read the original Excel results of this analyte class
    Hashtable<String,QuantificationResult> originalExcelResults = new Hashtable<String,QuantificationResult>();
    HydroxyEncoding faEncoding = null;
    HydroxyEncoding lcbEncoding = null;
    for (String expId : expFullPaths.keySet()){
      QuantificationResult result = LDAResultReader.readResultFile(expFullPaths.get(expId),  new Hashtable<String,Boolean>(), sheet.getSheetName());
      if (result.getFaHydroxyEncoding()!=null) faEncoding = result.getFaHydroxyEncoding();
      if (result.getLcbHydroxyEncoding()!=null) lcbEncoding = result.getLcbHydroxyEncoding();
      originalExcelResults.put(expId, result);
    }
    LinkedHashMap<String,Boolean> adductsSorted = extractAdductsSortedByAbundance(sheet.getSheetName(),originalExcelResults);

    //extract the values according to the molecular species
    boolean rtGroupingUsed = isRtGrouped;
    LinkedHashMap<String,SummaryVO> molSpeciesDetails = new LinkedHashMap<String,SummaryVO>();
    LinkedHashMap<String,Vector<DoubleBondPositionVO>> omegaSuggestionAvailable = new LinkedHashMap<String,Vector<DoubleBondPositionVO>>();
    //Quantification
    for (String molName : molNames){
    	if (!omegaSuggestionAvailable.containsKey(molName)) omegaSuggestionAvailable.put(molName, new Vector<DoubleBondPositionVO>());
      Hashtable<String,Hashtable<String,Vector<LipidParameterSet>>> relevantOriginals = null;
      if (includeResultFiles){
        relevantOriginals = new Hashtable<String,Hashtable<String,Vector<LipidParameterSet>>>();
        for (String expId : expFullPaths.keySet()){
          ResultAreaVO areaVO = compLookup.getResultAreaVO(sheet.getSheetName(),molName,expId);
          for (String mod : modifications)
          {
          	if (areaVO == null) continue;
          	for (LipidParameterSet set : areaVO.getLipidParameterSets(mod))
          	{
          		if (!set.getOmegaInformation().isEmpty())
          		{
          			omegaSuggestionAvailable.get(molName).addAll(set.getOmegaInformation());
          		}
          	}
          }
          if (!isGrouped && !expNames.containsKey(expId))
            continue;
          if (areaVO!=null){
            if (areaVO.isAStandard())
              rtGroupingUsed = false;
            Hashtable<String,Vector<LipidParameterSet>> sets = getRelevantOriginalResults(originalExcelResults.get(expId).getIdentifications().get(sheet.getSheetName()),areaVO,false);
            relevantOriginals.put(expId, sets);
          }
        }
      }
      SpeciesExportVO exportVO = LDAExporter.extractExportableSummaryInformation(speciesType, exportDoubleBondPositionsForClass, false, 0, 0, false, 0, 0, rtGroupingUsed, adductsSorted, new Vector<String>(expFullPaths.keySet()),
        expsOfGroup,  molName, results.get(molName), relevantOriginals, maxIsotope, faEncoding, lcbEncoding);
      for (SummaryVO sumVO : exportVO.getSummaries()){
        String id = "";
        //only for the classes overview, the name of the analyte class is not added
        if (includeResultFiles)
          id += sheet.getSheetName()+" ";
        id += sumVO.getSpeciesId();
        if (sumVO.getMolecularId()!=null)
          id += " | "+sumVO.getMolecularId();
        molSpeciesDetails.put(id, sumVO);
        omegaSuggestionAvailable.put(id, omegaSuggestionAvailable.get(molName));
      }
    }
    writeExcelSheetToFile(omegaCollector, sheet,workbook,out,molSpeciesDetails, isGrouped, expIdNames, expNames, expsOfGroup, preferredUnit, expVO, modifications, omegaSuggestionAvailable);
  }


  //TODO: maybe just write out the summary part and not the stuff that's already done in the normal excel sheet.
  private static void writeExcelSheetToFile(OmegaCollector omegaCollector, Sheet sheet, XSSFWorkbook workbook, OutputStream out, LinkedHashMap<String,SummaryVO> molSpeciesDetails,
      boolean isGrouped, Vector<String> expIdNames, Hashtable<String,String> expNames, LinkedHashMap<String,Vector<String>> expsOfGroup,
      String preferredUnit, ExportOptionsVO expVO, ArrayList<String> modifications, LinkedHashMap<String,Vector<DoubleBondPositionVO>> omegaSuggestionAvailable) throws IOException{  
    CellStyle headerStyle = getHeaderStyle(workbook);
    CellStyle numberStyle = ExcelUtils.getMassListNumberStyle(workbook);
    CellStyle normalStyle = getNormalStyle(workbook);
    CellStyle zeroStyle = getZeroStyle(workbook);
    CellStyle ms2Unambiguous = getMs2UnambiguousStyle((XSSFWorkbook)workbook);
    CellStyle ms2Split = getMs2SplitStyle((XSSFWorkbook)workbook);
    CellStyle ms2NoSplitPossible = getMs2NoSplitPossibleStyle((XSSFWorkbook)workbook);
    Row row = sheet.createRow(0);

    boolean excelFile = true; //only allowing text export
    Vector<String> columnNames = new Vector<String>(expIdNames);
    Collections.sort(columnNames);
    
    Vector<String> rowNames = new Vector<String>(molSpeciesDetails.keySet());
    boolean expInColumn = true; //not allowing this setting, as it messes with omega pos summary 

    StaticUtils.printOutLabelToFile(sheet.getSheetName(), out, row, headerStyle, (columnNames.size()+2)/2, excelFile, false);
    int rowCount = 1; 
    int columnCount = 0;
    Hashtable<Integer,String> headerValues = new Hashtable<Integer,String>();
    Hashtable<Integer,Integer> longestValues = new Hashtable<Integer,Integer>();
    String unit = "Arb.Units / AU";
    if (preferredUnit!=null && preferredUnit.length()>0)
      unit = preferredUnit;
    String unitString = unit;
    row = sheet.createRow(rowCount);
    StaticUtils.printOutLabelToFile(unitString, out, row, headerStyle, columnCount, excelFile, false);
    headerValues.put(columnCount, unitString);
    longestValues.put(columnCount, 0);
    columnCount++;
    StaticUtils.printOutLabelToFile("", out, row, headerStyle, columnCount, excelFile, true);
    headerValues.put(columnCount, "");
    longestValues.put(columnCount, 0);
    columnCount++;
    Hashtable<String,Hashtable<String,OmegaSummary>> valuesForOmega = new Hashtable<String,Hashtable<String,OmegaSummary>>();
    for (String columnName : columnNames){
      String valueToPrint = columnName;
      if (expInColumn)
        valueToPrint = expNames.get(columnName);
      StaticUtils.printOutLabelToFile(valueToPrint, out, row, headerStyle, columnCount, excelFile, true);
      headerValues.put(columnCount, valueToPrint);
      longestValues.put(columnCount, 0);
      valuesForOmega.put(valueToPrint, new Hashtable<String,OmegaSummary>());
      columnCount++;
    }
    row = sheet.createRow(++rowCount);
    double area;
    LinkedHashMap<String,ArrayList<Pair<Double,Double>>> weightedRTCollector = new LinkedHashMap<String,ArrayList<Pair<Double,Double>>>();
    for (String rowName : rowNames){
    	double areaSum = 0.0;
    	String molNameForWeightedCollector = rowName;
    	
    	if ((rowName.contains("|") || rowName.startsWith("L")) && !rowName.contains("IntS") && expVO.isExportRT()) //this means there's a molecular species and we export RTs
      {
      	if (rowName.contains("|"))
      	{
      		molNameForWeightedCollector = rowName.substring(rowName.indexOf("|")+2);
      	}
      	else if (rowName.startsWith("L"))
      	{
      		molNameForWeightedCollector = rowName.substring(rowName.indexOf(" ")+1, rowName.indexOf("_"));
      	}
      	
      	ArrayList<String> insensitive = new ArrayList<String>();
      	insensitive.add(molNameForWeightedCollector);
      	if (molNameForWeightedCollector.contains("/"))
      	{
      		insensitive = new ArrayList<String>(Arrays.asList(molNameForWeightedCollector.split("/")));
      	}
      	else if (molNameForWeightedCollector.contains("_"))
      	{
      		insensitive = new ArrayList<String>(Arrays.asList(molNameForWeightedCollector.split("_")));
      	}
      	molNameForWeightedCollector = insensitive.toString();
      	
      	if (!weightedRTCollector.containsKey(molNameForWeightedCollector))
      	{
      		weightedRTCollector.put(molNameForWeightedCollector, new ArrayList<Pair<Double,Double>>());
      	}
      }
    	
    	//extract omega collector info
    	fillOmegaCollector(omegaCollector, sheet.getSheetName(), molSpeciesDetails,
    				columnNames, expNames, expsOfGroup,
    				expVO, modifications, rowName, omegaSuggestionAvailable.get(rowName));
    	
    	String omegaId = extractOmegaPos(rowName);
      columnCount = 0;
      String valueToPrint = rowName;
      if (!expInColumn)
        valueToPrint = expNames.get(rowName);
      StaticUtils.printOutLabelToFile(valueToPrint, out, row, headerStyle, columnCount, excelFile, false);
      if (valueToPrint.length()>headerValues.get(columnCount).length())
        headerValues.put(columnCount,valueToPrint);
      columnCount++;
      StaticUtils.printOutLabelToFile("Value", out, row, headerStyle, columnCount, excelFile, true);
      if ("Value".length()>headerValues.get(columnCount).length())
        headerValues.put(columnCount,"Value");
      columnCount++;
      
      if (!isAreaBelowDetectorLimit(columnNames, rowName, expVO, isGrouped, molSpeciesDetails)) continue;
      
      for (String columnName : columnNames){
      	if (!valuesForOmega.get(columnName).containsKey(omegaId))
      	{
      		valuesForOmega.get(columnName).put(omegaId, new OmegaSummary());
      	}
        String toWrite = "0";
        String molName = columnName;
        String expName = rowName;
        if (!expVO.isAnalyteInColumn()){
          molName = rowName;
          expName = columnName;
        }
        SummaryVO sumVO = molSpeciesDetails.get(molName);
        
      	if (isGrouped)
          area = sumVO.getMeanArea(expName);
        else
          area = sumVO.getArea(expName);
        
      	areaSum += area;
        valuesForOmega.get(columnName).get(omegaId).addArea(area);
        toWrite = StaticUtils.extractDisplayValue(area,expVO.getCommaPositions());   
        
        CellStyle style = null;
        if (toWrite.equalsIgnoreCase("0")||toWrite.equalsIgnoreCase("0.0"))
          style = zeroStyle;
        else
          style = numberStyle;
//          getCorrespondingCellStyle(sumVO.getEvidenceReliabilty(expName),normalStyle,
//        		  ms2Unambiguous, ms2Split, ms2NoSplitPossible);
        StaticUtils.printOutNumberToFile(toWrite, out, row, style, columnCount, excelFile, true);
        if (toWrite.length()>longestValues.get(columnCount))
          longestValues.put(columnCount,toWrite.length());
        columnCount++;
      }    
      row = sheet.createRow(++rowCount);
      // write out the standard deviation
      if (expVO!=null && expVO.getExportType() != ExportOptionsVO.EXPORT_NO_DEVIATION){
        columnCount = 0;
        StaticUtils.printOutLabelToFile("", out, row, headerStyle, columnCount, excelFile, false);
        columnCount++;
        String title = "";
        if (expVO.getExportType() == ExportOptionsVO.EXPORT_SD_DEVIATION){
          title+="SD-"+expVO.getSdValue();
        } else if (expVO.getExportType() == ExportOptionsVO.EXPORT_SD_ERROR){
          title+="SE";
        }
        StaticUtils.printOutLabelToFile(title, out, row, headerStyle, columnCount, excelFile, true);
        if (title.length()>headerValues.get(columnCount).length())
          headerValues.put(columnCount,title);

        columnCount++;
        for (String columnName : columnNames){
          String toWrite = "-1.0";
          String molName = columnName;
          String expName = rowName;
          if (!expVO.isAnalyteInColumn()){
            molName = rowName;
            expName = columnName;
          }
          SummaryVO sumVO = molSpeciesDetails.get(molName);
          double value = sumVO.calculateDeviationValue(expVO, expsOfGroup.get(expName));
          valuesForOmega.get(columnName).get(omegaId).addStdev(value);
          if (!Double.isInfinite(value) && !Double.isNaN(value))
            toWrite = StaticUtils.extractDisplayValue(value);
          CellStyle style = null;
          if (toWrite.equalsIgnoreCase("-1")||toWrite.equalsIgnoreCase("-1.0"))
            style = zeroStyle;
          else
            style = numberStyle;
//            getCorrespondingCellStyle(sumVO.getEvidenceReliabilty(expName),normalStyle,
//          		  ms2Unambiguous, ms2Split, ms2NoSplitPossible);
          StaticUtils.printOutNumberToFile(toWrite, out, row, style, columnCount, excelFile, true);
          if (toWrite.length()>longestValues.get(columnCount))
            longestValues.put(columnCount,toWrite.length());
          columnCount++;
        }      
        row = sheet.createRow(++rowCount);
      }
      //export of retention times
      for (String modName : modifications){
        if (expVO!=null && expVO.isExportRT()){
          columnCount = 0;
          StaticUtils.printOutLabelToFile("", out, row, headerStyle, columnCount, excelFile, false);
          columnCount++;
          String rtName = "RT";
          if (modifications.size()>1)
            rtName+="-"+modName;
          StaticUtils.printOutLabelToFile(rtName, out, row, headerStyle, columnCount, excelFile, true);
          if (rtName.length()>headerValues.get(columnCount).length())
            headerValues.put(columnCount,rtName);

          columnCount++;
          for (String columnName : columnNames){
            String toWrite = "-1.0";
            String molName = columnName;
            String expName = rowName;
            if (!expVO.isAnalyteInColumn()){
              molName = rowName;
              expName = columnName;
            }
            SummaryVO sumVO = molSpeciesDetails.get(molName);
            Double rt;
            if (isGrouped)
              rt = sumVO.getMeanRetentionTime(modName, expName);
            else
              rt = sumVO.getRetentionTime(modName, expName);
            
            Double rtForWeighted = rt!=null && !Double.isInfinite(rt) && !Double.isNaN(rt) ? rt : 0.0;
            if (weightedRTCollector.containsKey(molNameForWeightedCollector))
            {
            	weightedRTCollector.get(molNameForWeightedCollector).add(new Pair<Double,Double>(areaSum, rtForWeighted));
            }
            
            if (rt!=null && !Double.isInfinite(rt) && !Double.isNaN(rt))
                toWrite = StaticUtils.extractDisplayValue(rt);
            CellStyle style = null;
            if (toWrite.equalsIgnoreCase("-1")||toWrite.equalsIgnoreCase("-1.0"))
              style = zeroStyle;
            else
              style = getCorrespondingCellStyle(sumVO.getEvidenceReliabilty(expName),normalStyle,
            		  ms2Unambiguous, ms2Split, ms2NoSplitPossible);
            StaticUtils.printOutNumberToFile(toWrite, out, row, style, columnCount, excelFile, true);
            if (toWrite.length()>longestValues.get(columnCount))
              longestValues.put(columnCount,toWrite.length());
            columnCount++;
          }  

          row = sheet.createRow(++rowCount);
        }
        //export of standard deviations of the retention times
        if (expVO!=null && expVO.isExportRTDev()){
          columnCount = 0;
          StaticUtils.printOutLabelToFile("", out, row, headerStyle, columnCount, excelFile, false);
          columnCount++;
          StaticUtils.printOutLabelToFile("RT-StDev", out, row, headerStyle, columnCount, excelFile, true);
          if ("RT-StDev".length()>headerValues.get(columnCount).length())
            headerValues.put(columnCount,"RT-StDev");

          columnCount++;
          for (String columnName : columnNames){
            String toWrite = "-1.0";
            String molName = columnName;
            String expName = rowName;
            if (!expVO.isAnalyteInColumn()){
              molName = rowName;
              expName = columnName;
            }
            SummaryVO sumVO = molSpeciesDetails.get(molName);
            Double rtStdev = sumVO.getStdevRetentionTime(modName, expName);
            if (rtStdev!=null && !Double.isInfinite(rtStdev) && !Double.isNaN(rtStdev))
              toWrite = StaticUtils.extractDisplayValue(rtStdev);
            CellStyle style = null;
            if (toWrite.equalsIgnoreCase("-1")||toWrite.equalsIgnoreCase("-1.0"))
              style = zeroStyle;
            else
              style = getCorrespondingCellStyle(sumVO.getEvidenceReliabilty(expName),normalStyle,
            	ms2Unambiguous, ms2Split, ms2NoSplitPossible);
            StaticUtils.printOutNumberToFile(toWrite, out, row, style, columnCount, excelFile, true);
            columnCount++;
          }          
          row = sheet.createRow(++rowCount);           
        }
      }      
    }
    
    rowCount += 3;
    for (String molName : weightedRTCollector.keySet())
    {
    	row = sheet.createRow(++rowCount); 
    	Cell cell = row.createCell(0,Cell.CELL_TYPE_STRING);
    	cell.setCellValue(molName);
	    cell = row.createCell(1,Cell.CELL_TYPE_NUMERIC);
	    Double weightedRT = 0.0;
	    Double sum = 0.0;
	    Double sumWeights = 0.0;
	    Integer num = 0;
	    for (Pair<Double,Double> pair : weightedRTCollector.get(molName))
	    {
	    	if (pair.getFirst()>0 && pair.getSecond()>0)
	    	{
	    		sum += pair.getFirst()*pair.getSecond();
	    		sumWeights += pair.getFirst();
	    		num++;
	    	}
	    }
	    weightedRT = sum/sumWeights;
	    cell.setCellValue(Precision.round(weightedRT, 2));
    }

    writeValuesForOmega(valuesForOmega, sheet, workbook, out, rowCount+3, headerStyle, numberStyle);

    if (excelFile){
      for (int column : headerValues.keySet()){
        setColumnWidth(sheet, column, headerValues.get(column), longestValues.get(column));
      }
    }
  }
  
  private static void fillOmegaCollector(OmegaCollector omegaCollector, String lClass, LinkedHashMap<String,SummaryVO> molSpeciesDetails,
			Vector<String> columnNames, Hashtable<String,String> expNames, LinkedHashMap<String,Vector<String>> expsOfGroup,
			ExportOptionsVO expVO, ArrayList<String> modifications, String rowName, Vector<DoubleBondPositionVO> vector)
  {
  	if (lClass.equals("Cer") || lClass.equals("SM") || !rowName.contains("_")) return; //TODO: we ignore these for now, if rowName does not contain "_" it is an internal standard
  	
  	int numFA = 2;
  	if (lClass.startsWith("L")) //lysos have 1 FA
  	{
  		numFA = 1;
  	}
  	
  	//summing up the total for all FA (assumes this is set to MSn Only)
  	for (String experimentName : columnNames){
  		
      SummaryVO sumVO = molSpeciesDetails.get(rowName);
      Double area = sumVO.getMeanArea(experimentName);	
      Double stdev = sumVO.calculateDeviationValue(expVO, expsOfGroup.get(experimentName));
      
      //summing up the total for all FA. mol dependent on num of FA
      omegaCollector.addToTotalFAPerClass(experimentName, lClass, area*numFA);
      
      if (rowName.contains("|") || rowName.startsWith("L")) //this means there's a molecular species
      {
      	String molName = rowName;
      	if (rowName.contains("|"))
      	{
      		molName = rowName.substring(rowName.indexOf("|")+2);
      	}
      	else if (rowName.startsWith("L"))
      	{
      		molName = rowName.substring(rowName.indexOf(" ")+1, rowName.indexOf("_"));
      	}
      	try
      	{
      		Vector<FattyAcidVO> chains = StaticUtils.decodeFAsFromHumanReadableName(molName,Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(),false,null);
      		Vector<FattyAcidVO> sameChainType = new Vector<FattyAcidVO>();
      		for (FattyAcidVO chain : chains)
      		{
      			//mapping all chain types to acyl chains
//      			FattyAcidVO sameChain = new FattyAcidVO(
//      					LipidomicsConstants.CHAIN_TYPE_FA_ACYL, "", chain.getcAtoms(), chain.getDoubleBonds(), chain.getOhNumber(), chain.getMass(), chain.getFormula(), chain.getOxState());
      			
      			FattyAcidVO sameChain = new FattyAcidVO(
      					chain.getChainType(), "", chain.getcAtoms(), chain.getDoubleBonds(), chain.getOhNumber(), chain.getMass(), chain.getFormula(), chain.getOxState());
      			
      			sameChain.setOmegaPosition(chain.getOmegaPosition());
      			
      			sameChainType.add(sameChain);
      			
      			if (sameChain.getOmegaPosition()>0)
      			{
      				//filling info for the per description count
      				omegaCollector.addToTotalFAPerDescription(experimentName, OmegaCollector.DESCRIPTION_ASSIGNED, area);
      				
      				//filling info of the omega total validation, one mol for each
      				omegaCollector.addToTotalOmegaPerClass(experimentName, lClass, sameChain.getOmegaPosition(), area);
      				omegaCollector.addToTotalOmegaPerClassSTD(experimentName, sameChain.getOmegaPosition(), stdev);
      				
      				//filling info for the fa to lClass heatmap, one mol for each
      				omegaCollector.addToTotalOmegaFAPerClass(experimentName, lClass, sameChain.getCarbonDbsId(), area);
      			}
      			else
      			{
      				String description = sameChain.getDoubleBonds()>0 ? OmegaCollector.DESCRIPTION_UNASSIGNED : OmegaCollector.DESCRIPTION_SATURATED;
      				//filling info for the per description count
      				omegaCollector.addToTotalFAPerDescription(experimentName, description, area);
      				
      				// go into further detail
      				if (description.equals(OmegaCollector.DESCRIPTION_UNASSIGNED))
      				{
      					if (!vector.isEmpty())
      					{
      						omegaCollector.addToTotalFAPerDescription(experimentName, OmegaCollector.DESCRIPTION_UNASSIGNED_SUGGESTED, area);
      					}
      					
      					//first odd or even chain, then partner
      					if (sameChain.getcAtoms()%2 == 1) //if uneven
      					{
      						omegaCollector.addToTotalFAPerDescription(experimentName, OmegaCollector.DESCRIPTION_UNASSIGNED_ODD_CHAIN, area);
      						if (chains.size()>1)
      						{
      							FattyAcidVO partner = chains.indexOf(chain) == 1 ? chains.get(0) : chains.get(1); 
      							if (partner.getOmegaPosition() > 0 || partner.getDoubleBonds() == 0)
      							{
      								omegaCollector.addToTotalFAPerDescription(experimentName, OmegaCollector.DESCRIPTION_UNASSIGNED_ODD_CHAIN_PARTNER_KNOWN, area);
      							}
      							else
      							{
      								omegaCollector.addToTotalFAPerDescription(experimentName, OmegaCollector.DESCRIPTION_UNASSIGNED_ODD_CHAIN_PARTNER_UNASSIGNED, area);
      							}
      						}
      						else
      						{
      							omegaCollector.addToTotalFAPerDescription(experimentName, OmegaCollector.DESCRIPTION_UNASSIGNED_ODD_CHAIN_PARTNER_KNOWN, area);
      						}
      					}
      					else
      					{
      						omegaCollector.addToTotalFAPerDescription(experimentName, OmegaCollector.DESCRIPTION_UNASSIGNED_EVEN_CHAIN, area);
      						if (chains.size()>1)
      						{
      							FattyAcidVO partner = chains.indexOf(chain) == 1 ? chains.get(0) : chains.get(1); 
      							if (partner.getOmegaPosition() > 0 || partner.getDoubleBonds() == 0)
      							{
      								omegaCollector.addToTotalFAPerDescription(experimentName, OmegaCollector.DESCRIPTION_UNASSIGNED_EVEN_CHAIN_PARTNER_KNOWN, area);
      							}
      							else
      							{
      								omegaCollector.addToTotalFAPerDescription(experimentName, OmegaCollector.DESCRIPTION_UNASSIGNED_EVEN_CHAIN_PARTNER_UNASSIGNED, area);
      							}
      						}
      						else
      						{
      							omegaCollector.addToTotalFAPerDescription(experimentName, OmegaCollector.DESCRIPTION_UNASSIGNED_EVEN_CHAIN_PARTNER_KNOWN, area);
      						}
      					}
      				}
      			}
      		}
      		
      		
      		//ensuring that the first entry 
      		if (sameChainType.size() > 1)
      		{
      			FattyAcidVO first = sameChainType.get(0);
      			FattyAcidVO second = sameChainType.get(1);
      			if (first.getDoubleBonds() < second.getDoubleBonds())
      			{
      				sameChainType.clear();
      				sameChainType.add(second);
      				sameChainType.add(first);
      			}
      		}
      		
      		boolean bothAssigned = true;
      		for (FattyAcidVO chain : sameChainType)
      		{
      			if (chain.getDoubleBonds() > 0 && chain.getOmegaPosition()<1) bothAssigned = false;
      		}
      		
      		//filling partner info. mol dependent on num of FA, as that's also part of the total calculation.
      		if (bothAssigned)
      		{
      			omegaCollector.addToFAContentToPartnerContent(experimentName, lClass, 
            		sameChainType.size()>1 ? StaticUtils.getHumanReadableChainName(sameChainType.get(1), Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), sameChainType.get(0).getOhNumber()>0) : "",
            		StaticUtils.getHumanReadableChainName(sameChainType.get(0), Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), sameChainType.get(0).getOhNumber()>0),
            		area*numFA);
      		}
          if (molName.contains("/"))
          {
          	omegaCollector.addToSn1ContentToSn2Content(experimentName, lClass, sameChainType.get(0).getCarbonDbsId(), sameChainType.get(1).getCarbonDbsId(), area*numFA);
          }
      	}
      	
      	catch (Exception ex)
      	{
      		ex.printStackTrace();
      	}
      }
  	}
  }

  public static void writeHeatMapData(OmegaCollector omegaCollector, XSSFWorkbook workbook, OutputStream out)
  {	
  	Sheet sheetA = workbook.createSheet("Validation");
  	Row headerRowA = sheetA.createRow(0);
  	Sheet sheetD = workbook.createSheet("Description Count");
  	Row headerRowD = sheetD.createRow(0);
  	Sheet sheetE = workbook.createSheet("FAvsFA Combined");
  	CellStyle headerStyle = getHeaderStyle(workbook);
    CellStyle numberStyle = ExcelUtils.getMassListNumberStyle(workbook);
  	
  	ArrayList<LinkedHashMap<String,LinkedHashMap<String,Double>>> fafaData = new ArrayList<LinkedHashMap<String,LinkedHashMap<String,Double>>>();
  	Vector<String> experimentNames = omegaCollector.getExperimentNames();
  	String experimentNameUnsuppl = "F) RAW264.7 unsupplemented"; //ensure this is named this way
  	for (int i=0;i<experimentNames.size();i++)
  	{
  		int columnNr = i+1;
  		String experimentName = experimentNames.get(i);
  		System.out.println(experimentName);
  		Double totalFAContent = omegaCollector.getTotalFA(experimentName);
  		
//  		if (!experimentName.startsWith("F")) continue; //testing with just unsupplemented
	  	Sheet sheetB = workbook.createSheet("HeatMap LC "+experimentName.substring(0,2)); //TODO: change for other data types.
	  	Sheet sheetC = workbook.createSheet("HeatMap FA "+experimentName.substring(0,2));
	  	
	    Row headerRowB = sheetB.createRow(0);
	    Row headerRowC = sheetC.createRow(0);
  		
  		
  		//for validation
	    writeValidation(omegaCollector, experimentName, experimentNameUnsuppl, columnNr, i, sheetA, headerRowA, headerStyle, numberStyle, totalFAContent);
  		
  		//for fa to lclass heatmap
  		writeFAvsLClassHeatmap(omegaCollector, experimentName, sheetB, headerRowB, headerStyle, numberStyle, totalFAContent);
  		
  		//for fa to fa heatmap
  		fafaData.add(writeFAvsFAHeatmap(omegaCollector, experimentName, sheetC, headerRowC, headerStyle, numberStyle, totalFAContent));
  		
  		//for description
  		writeDescription(omegaCollector, experimentName, columnNr, i, sheetD, headerRowD, headerStyle, numberStyle, totalFAContent);
  		
  	}
  	
  	writeFAvsFABarChartData(fafaData, sheetE, headerStyle, numberStyle);
  	
  }
  
  private static void writeFAvsFABarChartData(ArrayList<LinkedHashMap<String,LinkedHashMap<String,Double>>> fafaData,
  		Sheet sheet, CellStyle headerStyle, CellStyle numberStyle)
  {
  	LinkedHashSet<String> set = new LinkedHashSet<String>();
  	for (LinkedHashMap<String,LinkedHashMap<String,Double>> data : fafaData)
  	{
  		set.addAll(data.keySet());
  	}
  	Hashtable<String,Hashtable<String,ArrayList<Double>>> valueStd = new Hashtable<String,Hashtable<String,ArrayList<Double>>>();
  	for (String fa : set)
  	{
  		valueStd.put(fa, new Hashtable<String,ArrayList<Double>>());
  		for (LinkedHashMap<String,LinkedHashMap<String,Double>> data : fafaData)
    	{
    		if (data.containsKey(fa))
    		{
    			for (String partner : data.get(fa).keySet())
    			{
    				if (!valueStd.get(fa).containsKey(partner))
    				{
    					valueStd.get(fa).put(partner, new ArrayList<Double>());
    				}
    				valueStd.get(fa).get(partner).add(data.get(fa).get(partner));
    			}
    		}
    	}
  	}
  	
  	ArrayList<String> fas = new ArrayList<String>(set);
  	int rowNr = 0;
  	Row headerRow = sheet.createRow(rowNr);
  	for (int i=0; i<fas.size(); i++)
  	{
  		String fa = fas.get(i);
  		Hashtable<String,ArrayList<Double>> values = valueStd.get(fa);
  		int columnNr = i+1;
  		Cell cell = headerRow.createCell(columnNr,Cell.CELL_TYPE_STRING);
  		cell.setCellValue(fa);
  		cell.setCellStyle(headerStyle);
  		
  		for (int j=0; j<fas.size(); j++)
    	{
  			String partner = fas.get(j);
  			rowNr = j+1;
  			Row row;
  			if (i==0) 
  			{
  				row = sheet.createRow(rowNr);
  				cell = row.createCell(0,Cell.CELL_TYPE_STRING);
    			cell.setCellValue(partner);
    			cell.setCellStyle(headerStyle);
  			}
  			else 
  			{
  				row = sheet.getRow(rowNr);
  			}
  			
  			cell = row.createCell(columnNr, Cell.CELL_TYPE_NUMERIC);
  			double number = 0.0;
  			
  			if (values.containsKey(partner))
  			{
  				ArrayList<Double> nonZero = values.get(partner);
  				nonZero.removeAll(Collections.singleton(0.0));
  				if (!nonZero.isEmpty())
  				{
  					for (double x : nonZero)
      			{
      				number += x;
      			}
      			number = number / nonZero.size();
  				}
  			} 			
  			cell.setCellValue(number);
  			cell.setCellStyle(numberStyle);
    		
    	}
  	}
  	
  	rowNr += 3;
  	headerRow = sheet.createRow(rowNr);
  	int rowNrLocal = rowNr;
  	for (int i=0; i<fas.size(); i++)
  	{
  		String fa = fas.get(i);
  		Hashtable<String,ArrayList<Double>> values = valueStd.get(fa);
  		int columnNr = i+1;
  		Cell cell = headerRow.createCell(columnNr,Cell.CELL_TYPE_STRING);
  		cell.setCellValue(fa);
  		cell.setCellStyle(headerStyle);
  		
  		for (int j=0; j<fas.size(); j++)
    	{
  			String partner = fas.get(j);
  			rowNr = rowNrLocal+j+1;
  			Row row;
  			if (i==0) 
  			{
  				row = sheet.createRow(rowNr);
  				cell = row.createCell(0,Cell.CELL_TYPE_STRING);
    			cell.setCellValue(partner);
    			cell.setCellStyle(headerStyle);
  			}
  			else 
  			{
  				row = sheet.getRow(rowNr);
  			}
  			
  			cell = row.createCell(columnNr, Cell.CELL_TYPE_NUMERIC);
  			
  			double standardDeviation = 0.0;
  			if (values.containsKey(partner) && values.get(partner).size() > 2)
  			{
  				ArrayList<Double> nonZero = values.get(partner);
  				nonZero.removeAll(Collections.singleton(0.0));
  				if (nonZero.size() > 2)
  				{
  					double sum = 0.0;
    				for (double x : nonZero)
      			{
      				sum += x;
      			}
    				
    				if (sum > 0)
    				{
    					double mean = sum / nonZero.size();
    					sum = 0.0;
      				for (double x : nonZero)
        			{
      					sum += Math.pow(x-mean, 2);
        			}
      				standardDeviation = Math.sqrt(sum / (nonZero.size()-1));
    				}
    				cell.setCellValue(standardDeviation);
      			cell.setCellStyle(numberStyle);
  				}
  			}
    	}
  	}
  	
  }
  
  
  private static void writeValuesForOmega(Hashtable<String,Hashtable<String,OmegaSummary>> valuesForOmega, 
  		Sheet sheet, XSSFWorkbook workbook, OutputStream out, int rowCount, CellStyle headerStyle, CellStyle normalStyle)
  {
  	Row row = sheet.createRow(rowCount);
  	ArrayList<String> groupNames = new ArrayList<String>(valuesForOmega.keySet());
		Collections.sort(groupNames);

  	ArrayList<String> labelList = generateLabelListForOmega();
  	int columnCount = 1;
  	for (String groupName : groupNames)
  	{
  		StaticUtils.addExcelLabel(headerStyle, row, ++columnCount, groupName);
  	}
  	
  //write areas
  	row = sheet.createRow(++rowCount);
  	StaticUtils.addExcelLabel(headerStyle, row, 0, "area [AU]");
  	for (String label : labelList)
  	{
  		columnCount = 1;
  		row = sheet.createRow(++rowCount);
  		StaticUtils.addExcelLabel(headerStyle, row, columnCount, label);
  		for (String groupName : groupNames)
    	{
  			double area = 0.0;
  			if (valuesForOmega.get(groupName).containsKey(label))
  			{
  				area = valuesForOmega.get(groupName).get(label).getAverageArea();
  				if ((Double.isInfinite(area)||Double.isNaN(area)))
    			{
    				area = 0.0;
    			}
  			}
  			Cell number = row.createCell(++columnCount,Cell.CELL_TYPE_NUMERIC);
  	    number.setCellValue(area);
  	    number.setCellStyle(normalStyle);
    	}
  	}

  	//write standard deviations
  	rowCount+=2;
  	row = sheet.createRow(++rowCount);
  	StaticUtils.addExcelLabel(headerStyle, row, 0, "SD-1.0");
  	for (String label : labelList)
  	{
  		columnCount = 1;
  		row = sheet.createRow(++rowCount);
  		StaticUtils.addExcelLabel(headerStyle, row, columnCount, label);
  		for (String groupName : groupNames)
    	{
  			double stdev = 0.0;
  			if (valuesForOmega.get(groupName).containsKey(label))
  			{
  				stdev = valuesForOmega.get(groupName).get(label).getStdev();
  				if ((Double.isInfinite(stdev)||Double.isNaN(stdev)))
    			{
  					stdev = 0.0;
    			}
  			}
  			Cell number = row.createCell(++columnCount,Cell.CELL_TYPE_NUMERIC);
  	    number.setCellValue(stdev);
  	    number.setCellStyle(normalStyle);
    	}
  	}
  }

  private static String extractOmegaPos(String name)
  {
  	if (!name.contains(LipidomicsConstants.OMEGA_POSITION_START)) return "none";
  	String omegaId = "";
  	ArrayList<Integer> omegaPos = new ArrayList<Integer>();
  	omegaPos.add(Integer.parseInt(name.substring(
  			name.indexOf(LipidomicsConstants.OMEGA_POSITION_START)+LipidomicsConstants.OMEGA_POSITION_START.length(),
  			name.indexOf(LipidomicsConstants.OMEGA_POSITION_END))));
  	if (name.indexOf(LipidomicsConstants.OMEGA_POSITION_END)+LipidomicsConstants.OMEGA_POSITION_END.length() < name.length())
  	{
  		String remainder = name.substring(
  				name.indexOf(LipidomicsConstants.OMEGA_POSITION_END)+LipidomicsConstants.OMEGA_POSITION_END.length(),
  				name.length());
  		if (remainder.contains(LipidomicsConstants.OMEGA_POSITION_START))
  		{
  			omegaPos.add(Integer.parseInt(remainder.substring(
  					remainder.indexOf(LipidomicsConstants.OMEGA_POSITION_START)+LipidomicsConstants.OMEGA_POSITION_START.length(),
  					remainder.indexOf(LipidomicsConstants.OMEGA_POSITION_END))));
  			Collections.sort(omegaPos);
  		}
  	}
  	if (omegaPos.size() == 1 && name.contains(":0")) //this will not work for lyso species
  	{
  		omegaId = String.format("(n-%s)", omegaPos.get(0));
  	}
  	else if (omegaPos.size() == 1) //we must have a second unknown FA
  	{
  		omegaId = String.format("(n-%s)+(n-x)", omegaPos.get(0));
  	}
  	else if (omegaPos.size() == 2)
  	{
  		omegaId = String.format("(n-%s)+(n-%s)", omegaPos.get(0),omegaPos.get(1));
  	}
  	return omegaId;
  }
  
  private static ArrayList<String> generateLabelListForOmega()
  {
  	ArrayList<String> labelList = new ArrayList<String>();
  	labelList.add("(n-3)");
  	labelList.add("(n-3)+(n-3)");
  	labelList.add("(n-3)+(n-6)");
  	labelList.add("(n-3)+(n-7)");
  	labelList.add("(n-3)+(n-9)");
  	labelList.add("(n-3)+(n-10)");
  	labelList.add("(n-6)");
  	labelList.add("(n-6)+(n-6)");
  	labelList.add("(n-6)+(n-7)");
  	labelList.add("(n-6)+(n-9)");
  	labelList.add("(n-6)+(n-10)");
  	labelList.add("(n-7)");
  	labelList.add("(n-7)+(n-7)");
  	labelList.add("(n-7)+(n-9)");
  	labelList.add("(n-7)+(n-10)");
  	labelList.add("(n-9)");
  	labelList.add("(n-9)+(n-9)");
  	labelList.add("(n-9)+(n-10)");
  	labelList.add("(n-10)");
  	labelList.add("(n-10)+(n-10)");
  	return labelList;
  }
  
  private static boolean isAreaBelowDetectorLimit(
  		Vector<String> columnNames, String rowName, ExportOptionsVO expVO, boolean isGrouped, LinkedHashMap<String,SummaryVO> molSpeciesDetails)
  {
//  	double area = 0.0;
//  	for (String columnName : columnNames)
//    {
//    	String molName = columnName;
//      String expName = rowName;
//    	if (!expVO.isAnalyteInColumn()){
//        molName = rowName;
//        expName = columnName;
//      }
//      SummaryVO sumVO = molSpeciesDetails.get(molName);
//      if (isGrouped)
//        area = sumVO.getMeanArea(expName);
//      else
//        area = sumVO.getArea(expName);
//      if (area > 9.5E9)
//      {
//      	return false;
//      }
//    }
  	return true;
  }
  
  private static void writeValidation(
  		OmegaCollector omegaCollector, String experimentName, String experimentNameUnsuppl, int columnNr, int i,
  		Sheet sheetA, Row headerRowA, CellStyle headerStyle, CellStyle numberStyle, Double totalFAContent)
  {
  	LinkedHashMap<Integer,Double> normTotalOmegaFA = omegaCollector.getNormTotalOmegaFA(experimentName);
//		LinkedHashMap<Integer,Double> normTotalOmegaFAUnsuppl = omegaCollector.getNormTotalOmegaFA(experimentNameUnsuppl);
		Vector<Integer> omegaPos = new Vector<Integer>(normTotalOmegaFA.keySet());
		Collections.sort(omegaPos);
		
		Cell cell = headerRowA.createCell(columnNr,Cell.CELL_TYPE_STRING);
		cell.setCellValue(experimentName);
		cell.setCellStyle(headerStyle);
		
		//print values
		for (int j=0; j<omegaPos.size();j++)
		{
			Integer omega = omegaPos.get(j);
			Double normalized = normTotalOmegaFA.get(omega);
//			Double refNorm = normTotalOmegaFAUnsuppl.get(omega);
			
			Double value = normalized;
			
			int rowNr = j+1;
			Row row;
			if (i==0) 
			{
				row = sheetA.createRow(rowNr);
				cell = row.createCell(0,Cell.CELL_TYPE_STRING);
  			cell.setCellValue(String.format("omega %s", omega));
  			cell.setCellStyle(headerStyle);
			}
			else 
			{
				row = sheetA.getRow(rowNr);
			}
			
			cell = row.createCell(columnNr, Cell.CELL_TYPE_NUMERIC);
			cell.setCellValue(value);
			cell.setCellStyle(numberStyle);
		}
		
		//print std
		for (int j=0; j<omegaPos.size();j++)
		{
			Integer omega = omegaPos.get(j);
//			Double refNorm = normTotalOmegaFAUnsuppl.get(omega);
			Double stdev = omegaCollector.getOmegaSTD(experimentName, omega);
			
			Double value = stdev/totalFAContent;
			
			int rowNr = j+(omegaPos.size()+2);
			Row row;
			if (i==0) 
			{
				row = sheetA.createRow(rowNr);
				cell = row.createCell(0,Cell.CELL_TYPE_STRING);
  			cell.setCellValue(String.format("omega %s", omega));
  			cell.setCellStyle(headerStyle);
			}
			else 
			{
				row = sheetA.getRow(rowNr);
			}
			
			cell = row.createCell(columnNr, Cell.CELL_TYPE_NUMERIC);
			cell.setCellValue(value);
			cell.setCellStyle(numberStyle);
		}
  }
  
  private static void writeDescription(OmegaCollector omegaCollector, String experimentName, int columnNr, int i,
  		Sheet sheetD, Row headerRowD, CellStyle headerStyle, CellStyle numberStyle, Double totalFAContent)
  {
  	LinkedHashMap<String,Double> descriptionCount = omegaCollector.getTotalFAPerDescription(experimentName);
		
		Cell cell = headerRowD.createCell(columnNr,Cell.CELL_TYPE_STRING);
		cell.setCellValue(experimentName);
		cell.setCellStyle(headerStyle);
		
		Row row;
		if (i==0) 
		{
			row = sheetD.createRow(1);
			cell = row.createCell(0,Cell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("%s", OmegaCollector.DESCRIPTION_SATURATED));
			cell.setCellStyle(headerStyle);
			
			row = sheetD.createRow(2);
			cell = row.createCell(0,Cell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("%s", OmegaCollector.DESCRIPTION_ASSIGNED));
			cell.setCellStyle(headerStyle);
			
			row = sheetD.createRow(3);
			cell = row.createCell(0,Cell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("%s", OmegaCollector.DESCRIPTION_UNASSIGNED));
			cell.setCellStyle(headerStyle);
			
			row = sheetD.createRow(4);
			cell = row.createCell(0,Cell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("%s", OmegaCollector.DESCRIPTION_UNASSIGNED_SUGGESTED));
			cell.setCellStyle(headerStyle);
			
			row = sheetD.createRow(5);
			cell = row.createCell(0,Cell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("%s", OmegaCollector.DESCRIPTION_UNASSIGNED_EVEN_CHAIN));
			cell.setCellStyle(headerStyle);
			
			row = sheetD.createRow(6);
			cell = row.createCell(0,Cell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("%s", OmegaCollector.DESCRIPTION_UNASSIGNED_ODD_CHAIN));
			cell.setCellStyle(headerStyle);
			
			row = sheetD.createRow(7);
			cell = row.createCell(0,Cell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("%s", OmegaCollector.DESCRIPTION_UNASSIGNED_EVEN_CHAIN_PARTNER_KNOWN));
			cell.setCellStyle(headerStyle);
			
			row = sheetD.createRow(8);
			cell = row.createCell(0,Cell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("%s", OmegaCollector.DESCRIPTION_UNASSIGNED_ODD_CHAIN_PARTNER_KNOWN));
			cell.setCellStyle(headerStyle);
			
			row = sheetD.createRow(9);
			cell = row.createCell(0,Cell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("%s", OmegaCollector.DESCRIPTION_UNASSIGNED_EVEN_CHAIN_PARTNER_UNASSIGNED));
			cell.setCellStyle(headerStyle);
			
			row = sheetD.createRow(10);
			cell = row.createCell(0,Cell.CELL_TYPE_STRING);
			cell.setCellValue(String.format("%s", OmegaCollector.DESCRIPTION_UNASSIGNED_ODD_CHAIN_PARTNER_UNASSIGNED));
			cell.setCellStyle(headerStyle);
		}
		
		row = sheetD.getRow(1);
		cell = row.createCell(columnNr, Cell.CELL_TYPE_NUMERIC);
		cell.setCellValue(descriptionCount.get(OmegaCollector.DESCRIPTION_SATURATED) /totalFAContent);
		cell.setCellStyle(numberStyle);
		
		row = sheetD.getRow(2);
		cell = row.createCell(columnNr, Cell.CELL_TYPE_NUMERIC);
		cell.setCellValue(descriptionCount.get(OmegaCollector.DESCRIPTION_ASSIGNED) /totalFAContent);
		cell.setCellStyle(numberStyle);
		
		row = sheetD.getRow(3);
		cell = row.createCell(columnNr, Cell.CELL_TYPE_NUMERIC);
		cell.setCellValue(descriptionCount.get(OmegaCollector.DESCRIPTION_UNASSIGNED) /totalFAContent);
		cell.setCellStyle(numberStyle);
		
		row = sheetD.getRow(4);
		cell = row.createCell(columnNr, Cell.CELL_TYPE_NUMERIC);
		cell.setCellValue(descriptionCount.get(OmegaCollector.DESCRIPTION_UNASSIGNED_SUGGESTED) /totalFAContent);
		cell.setCellStyle(numberStyle);
		
		row = sheetD.getRow(5);
		cell = row.createCell(columnNr,Cell.CELL_TYPE_NUMERIC);
		cell.setCellValue(descriptionCount.get(OmegaCollector.DESCRIPTION_UNASSIGNED_EVEN_CHAIN) /totalFAContent);
		cell.setCellStyle(numberStyle);
		
		row = sheetD.getRow(6);
		cell = row.createCell(columnNr,Cell.CELL_TYPE_NUMERIC);
		cell.setCellValue(descriptionCount.get(OmegaCollector.DESCRIPTION_UNASSIGNED_ODD_CHAIN) /totalFAContent);
		cell.setCellStyle(numberStyle);
		
		row = sheetD.getRow(7);
		cell = row.createCell(columnNr,Cell.CELL_TYPE_NUMERIC);
		cell.setCellValue(descriptionCount.get(OmegaCollector.DESCRIPTION_UNASSIGNED_EVEN_CHAIN_PARTNER_KNOWN) /totalFAContent);
		cell.setCellStyle(numberStyle);
		
		row = sheetD.getRow(8);
		cell = row.createCell(columnNr,Cell.CELL_TYPE_NUMERIC);
		cell.setCellValue(descriptionCount.get(OmegaCollector.DESCRIPTION_UNASSIGNED_ODD_CHAIN_PARTNER_KNOWN) /totalFAContent);
		cell.setCellStyle(numberStyle);
		
		row = sheetD.getRow(9);
		cell = row.createCell(columnNr,Cell.CELL_TYPE_NUMERIC);
		cell.setCellValue(descriptionCount.get(OmegaCollector.DESCRIPTION_UNASSIGNED_EVEN_CHAIN_PARTNER_UNASSIGNED) /totalFAContent);
		cell.setCellStyle(numberStyle);
		
		row = sheetD.getRow(10);
		cell = row.createCell(columnNr,Cell.CELL_TYPE_NUMERIC);
		cell.setCellValue(descriptionCount.get(OmegaCollector.DESCRIPTION_UNASSIGNED_ODD_CHAIN_PARTNER_UNASSIGNED) /totalFAContent);
		cell.setCellStyle(numberStyle);
  }
  
  private static void writeFAvsLClassHeatmap(
  		OmegaCollector omegaCollector, String experimentName, Sheet sheetB, Row headerRowB, CellStyle headerStyle, CellStyle numberStyle, Double totalFAContent)
  {
		LinkedHashMap<String,LinkedHashMap<String,Double>> totalOmegaFAPerClass = omegaCollector.getTotalOmegaFAPerClass(experimentName);
		Vector<String> lClasses = new Vector<String>(totalOmegaFAPerClass.keySet());
		Collections.sort(lClasses);
		
		Set<String> allFASet = new HashSet<String>();
		for (String lClass : lClasses)
		{
			LinkedHashMap<String,Double> entries = totalOmegaFAPerClass.get(lClass);
			for (String faName : entries.keySet())
  		{
				Double area = entries.get(faName);
				if (area > 0)
				{
					allFASet.add(faName);
				}
  		}
		}
		
		Vector<FattyAcidVO> allFAs = new Vector<FattyAcidVO>();
		for (String faName : allFASet)
		{
			try
			{
				FattyAcidVO vo = StaticUtils.decodeHumanReadableChain(faName, Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(),false,null);
  			allFAs.add(vo);
			} catch (Exception ex){ex.printStackTrace();}
		}
		sortFAs(allFAs);
		
		//writing it to sheet
		for (int i=0; i<lClasses.size();i++)
		{
			int columnNr = i+1;
			Cell cell = headerRowB.createCell(columnNr,Cell.CELL_TYPE_STRING);
			cell.setCellValue(lClasses.get(i));
			cell.setCellStyle(headerStyle);
			
			for (int j=0; j<allFAs.size();j++)
			{
				int rowNr = j+1;
				Row row;
				if (i==0) 
				{
					row = sheetB.createRow(rowNr);
					cell = row.createCell(0,Cell.CELL_TYPE_STRING);
    			cell.setCellValue(allFAs.get(j).getCarbonDbsId());
    			cell.setCellStyle(headerStyle);
				}
				else 
				{
					row = sheetB.getRow(rowNr);
				}
				
				cell = row.createCell(columnNr, Cell.CELL_TYPE_NUMERIC);
				double number = 0.0;
				if (totalOmegaFAPerClass.get(lClasses.get(i)).containsKey(allFAs.get(j).getCarbonDbsId()))
				{
					number = totalOmegaFAPerClass.get(lClasses.get(i)).get(allFAs.get(j).getCarbonDbsId()) /totalFAContent;
				}
				cell.setCellValue(number);
  			cell.setCellStyle(numberStyle);
			}
		}
  }
  
  /**
   * 
   * @param omegaCollector
   * @param experimentName
   * @param sheetC
   * @param headerRowC
   * @param headerStyle
   * @param numberStyle
   * @param totalFAContent
   * @return a map with: key1: faName, key2: partnerFAName; list of intensities.
   */
  private static LinkedHashMap<String,LinkedHashMap<String,Double>> writeFAvsFAHeatmap(
  		OmegaCollector omegaCollector, String experimentName, Sheet sheetC, Row headerRowC, CellStyle headerStyle, CellStyle numberStyle, Double totalFAContent)
  {
  	LinkedHashMap<String,LinkedHashMap<String,Double>> collection = new LinkedHashMap<String,LinkedHashMap<String,Double>>();
  	
  	LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>> faContentToPartnerContent = omegaCollector.getFAContentToPartnerContent(experimentName);
//  	LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,Double>>> faContentToPartnerContent = omegaCollector.getSn1ContentToSn2Content(experimentName);
		//just want the total fa to fa, so summing up over classes
		LinkedHashMap<String,LinkedHashMap<String,Double>> totalFAFA = new LinkedHashMap<String,LinkedHashMap<String,Double>>();
		Set<String> allFASet1 = new HashSet<String>(); //to figure out rows and columns, this are the species with less double bonds if applicable
		Set<String> allFASet2 = new HashSet<String>(); //to figure out rows and columns, this are the species with more double bonds if applicable
		for (String lClass : faContentToPartnerContent.keySet())
		{
			LinkedHashMap<String,LinkedHashMap<String,Double>> fafaForClass = faContentToPartnerContent.get(lClass);
			for (String fa : fafaForClass.keySet())
			{
				LinkedHashMap<String,Double> partnerFAs = fafaForClass.get(fa);
				if (!totalFAFA.containsKey(fa))
				{
					totalFAFA.put(fa, new LinkedHashMap<String,Double>());
				}
				for (String partnerFA : partnerFAs.keySet())
				{
					Double value = partnerFAs.get(partnerFA);
					
					if (value>0 && 
							(fa.contains(LipidomicsConstants.OMEGA_POSITION_START) || partnerFA.contains(LipidomicsConstants.OMEGA_POSITION_START) || fa.length() <=1 ))  //one of the FA has to have an omega pos, or lyso
					{
						allFASet1.add(fa);
						allFASet2.add(partnerFA);
						allFASet2.add(fa);
						allFASet1.add(partnerFA);
					}
					if (!totalFAFA.get(fa).containsKey(partnerFA))
  				{
  					totalFAFA.get(fa).put(partnerFA, 0.0);
  				}
					Double before = totalFAFA.get(fa).get(partnerFA);
					totalFAFA.get(fa).put(partnerFA, before+value);
				}
			}
		}
		
		Double maxFAContent = 0.0;
		
		for (String fa : totalFAFA.keySet())
		{
			for (String partnerFA : totalFAFA.get(fa).keySet())
			{
				if (totalFAFA.get(fa).get(partnerFA) > maxFAContent)
				{
					maxFAContent = totalFAFA.get(fa).get(partnerFA);
				}
			}
		}
		
		//sorting fas
		Vector<FattyAcidVO> allFAs1 = new Vector<FattyAcidVO>();
		for (String faName : allFASet1)
		{
			if (faName.equals("")) 
			{
				if (!allFAs1.contains(null))
					allFAs1.add(null);
				continue;
			}
			try
			{
				FattyAcidVO vo = StaticUtils.decodeHumanReadableChain(faName, Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(),false,null);
  			allFAs1.add(vo);
			} catch (Exception ex){ex.printStackTrace();}
		}
		sortFAs(allFAs1);
		
		//sorting fas
		Vector<FattyAcidVO> allFAs2 = new Vector<FattyAcidVO>();
		for (String faName : allFASet2)
		{
			if (faName.equals(""))
			{
				if (!allFAs2.contains(null))
					allFAs2.add(null);
				continue;
			}
			try
			{
				FattyAcidVO vo = StaticUtils.decodeHumanReadableChain(faName, Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(),false,null);
  			allFAs2.add(vo);
			} catch (Exception ex){ex.printStackTrace();}
		}
		sortFAs(allFAs2);
		
		//writing fafa to sheet
		try 
		{
			for (int i=0; i<allFAs1.size();i++)
			{
				int columnNr = i+1;
				String key1 = allFAs1.get(i) == null ? "" : StaticUtils.getHumanReadableChainName(allFAs1.get(i), Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), allFAs1.get(i).getOhNumber()>0);
				String faName = key1.length() <= 1 ? "lyso" : key1;
				collection.put(faName, new LinkedHashMap<String,Double>());
				
				Cell cell = headerRowC.createCell(columnNr,Cell.CELL_TYPE_STRING);
				cell.setCellValue(faName);
				cell.setCellStyle(headerStyle);
				
				for (int j=0; j<allFAs2.size();j++)
				{
					String key2 = allFAs2.get(j) == null ? "" : StaticUtils.getHumanReadableChainName(allFAs2.get(j), Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), allFAs2.get(j).getOhNumber()>0);
					String faNamePartner = key2.length() <= 1 ? "lyso" : key2;
					
					int rowNr = j+1;
					Row row;
					if (i==0) 
					{
						row = sheetC.createRow(rowNr);
						cell = row.createCell(0,Cell.CELL_TYPE_STRING);
	    			cell.setCellValue(faNamePartner);
	    			cell.setCellStyle(headerStyle);
					}
					else 
					{
						row = sheetC.getRow(rowNr);
					}
					
					cell = row.createCell(columnNr, Cell.CELL_TYPE_NUMERIC);
					double number = 0.0;
					if (
//							columnNr <= rowNr && 
							totalFAFA.containsKey(key1) && totalFAFA.get(key1).containsKey(key2))
					{
						number = totalFAFA.get(key1).get(key2) /totalFAContent;
//						number = totalFAFA.get(key1).get(key2) /maxFAContent; //normalizing to highest occurrence
					}
					collection.get(faName).put(faNamePartner, number);
					
					cell.setCellValue(number);
	  			cell.setCellStyle(numberStyle);
				}
			}
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
		System.out.println("hihi");
		return collection;
  }
  
  private static void sortFAs(Vector<FattyAcidVO> fas)
  {
  	Collections.sort(fas, new Comparator<FattyAcidVO>() {
      @Override
      public int compare(FattyAcidVO fa1, FattyAcidVO fa2) {
      	if (fa1 == null && fa2 == null) return 0;
      	else if (fa1 == null) return 1;
      	else if (fa2 == null) return -1;
      	
      	return Comparator.comparing(FattyAcidVO::getcAtoms)
      			.thenComparing(FattyAcidVO::getDoubleBonds)
      			.thenComparing(FattyAcidVO::getOmegaPosition)
      			.thenComparing(FattyAcidVO::getChainType)
						.compare(fa1, fa2);
      	
//      	return Comparator.comparing(FattyAcidVO::getOmegaPosition)
//      			.thenComparing(FattyAcidVO::getcAtoms)
//      			.thenComparing(FattyAcidVO::getDoubleBonds)
//						.compare(fa1, fa2);
      }
		});
  }

  private static class OmegaSummary
  {
  	double sumArea_ = 0.0;
  	double sumVariance_ = 0.0;
  	int numContributions_ = 0;

  	private OmegaSummary() {}

  	void addArea(double area)
  	{
  		if (!(Double.isInfinite(area)||Double.isNaN(area)))
  		{
  			sumArea_ += area;
  			numContributions_++;
  		}
  	}

  	void addStdev(double stdev)
  	{
  		if (!(Double.isInfinite(stdev)||Double.isNaN(stdev)))
  		{
  			sumVariance_ += Math.pow(stdev, 2);
  		}
  	}

  	double getAverageArea()
  	{
  		return sumArea_/numContributions_;
  	}

  	double getStdev()
  	{
  		return Math.sqrt(sumVariance_/numContributions_);
  	}
  }  
}
