/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Vector;

import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFColor;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.analysis.ComparativeResultsLookup;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RetentionTimeGroupingException;
import at.tugraz.genome.lda.export.vos.SpeciesExportVO;
import at.tugraz.genome.lda.export.vos.SummaryVO;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.swing.LipidomicsTableCellRenderer;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.StaticUtils;
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
}
