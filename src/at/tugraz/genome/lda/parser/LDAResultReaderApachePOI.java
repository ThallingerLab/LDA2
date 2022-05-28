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
package at.tugraz.genome.lda.parser;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JFrame;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.analysis.ComparativeNameExtractor;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.export.QuantificationResultExporter;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.msn.vos.IntensityChainVO;
import at.tugraz.genome.lda.msn.vos.IntensityPositionVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.StaticUtils;
//import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.Probe3D;

/**
 * 
 * @author Juergen Hartler
 *
 */
//TODO: features supporting omega identifications are commented out - since this class is deprecated only minimal maintenance required

@Deprecated
//slower Apache POI reader is only used for old xls excel files
//TODO: remove completely after an adequate transition period *note written: 25.05.2022*
//(1 year should suffice as xls files are not written since a few years already and newer LDA versions do not support such old file formats anymore anyway) 
public class LDAResultReaderApachePOI
{

  /**
   * reads an LDA results file in Excel format
   * @param filePath the absolute path to the Excel file
   * @param showModifications this hash is filled by the method and gives information whether there are more than one modifications present; key: lipid class
   * @return the contents of the Excel file stored in the corresponding value object
   * @throws ExcelInputFileException when there is something wrong with the Excel file
   */
  public static QuantificationResult readResultFile(String filePath, Hashtable<String,Boolean> showModifications)
      throws ExcelInputFileException{
    return LDAResultReaderApachePOI.readResultFile(filePath, showModifications, null);
  }
  
  /**
   * reads an LDA results file in Excel format
   * @param filePath filePath the absolute path to the Excel file
   * @param showModifications this hash is filled by the method and gives information whether there are more than one modifications present; key: lipid class
   * @param specificClass filter for parsing only the results of one analyte class; enter null when no filter is required
   * @return the contents of the Excel file stored in the corresponding value object
   * @throws ExcelInputFileException when there is something wrong with the Excel file
   */
  public static QuantificationResult readResultFile(String filePath, Hashtable<String,Boolean> showModifications,
      String specificClass) throws ExcelInputFileException{
    Hashtable<String,Vector<LipidParameterSet>> resultParams = new Hashtable<String,Vector<LipidParameterSet>>();
    Hashtable<String,Integer> msLevels = new Hashtable<String,Integer>();
    LipidomicsConstants readConstants = null;
    HydroxyEncoding faOhEncodings = null;
    HydroxyEncoding lcbOhEncodings = null;
    String suffix = "";
    if (filePath!=null && filePath.length()>3)
      suffix = filePath.substring(filePath.lastIndexOf("."));
    if (!(suffix.equalsIgnoreCase(".xls")||(suffix.equalsIgnoreCase(".xlsx")))){
      new WarningMessage(new JFrame(), "ERROR", "The specified file is not Microsoft Excel!");
      throw new ExcelInputFileException("The specified file is not Microsoft Excel!");
    } 
    try {
      InputStream myxls = new FileInputStream(filePath);
      Workbook workbook  = null;
      if (suffix.equalsIgnoreCase(".xlsx")) workbook = new XSSFWorkbook(myxls);
      else if (suffix.equalsIgnoreCase(".xls")) workbook  = new HSSFWorkbook(myxls);
      //Workbook workbook = Workbook.getWorkbook(new File(filePath));
      
      for (int sheetNumber=0;sheetNumber!=workbook.getNumberOfSheets();sheetNumber++){       
        Sheet sheet = workbook.getSheetAt(sheetNumber);
        if (!sheet.getSheetName().endsWith(QuantificationResultExporter.ADDUCT_OMEGA_SHEET)&&!sheet.getSheetName().endsWith(QuantificationResultExporter.ADDUCT_OVERVIEW_SHEET)&&
            !sheet.getSheetName().endsWith(QuantificationResultExporter.ADDUCT_MSN_SHEET) && !sheet.getSheetName().equalsIgnoreCase(QuantificationResultExporter.SHEET_CONSTANTS)){
          if (specificClass!=null && !sheet.getSheetName().equalsIgnoreCase(specificClass))
            continue;
         
          Vector<LipidParameterSet> resultPrms = new Vector<LipidParameterSet>();
          int nameColumn = -1;
          int dbsColumn = -1;
          int ohColumn = -1;
          int modificationColumn = -1;
          int formulaColumn = -1;
          int modFormulaColumn = -1;
          int rtColumn = -1;
          int isotopeColumn = -1;
          int areaColumn = -1;
          int areaErrorColumn = -1;
          int backgroundColumn = -1;
          int chargeColumn = -1;
          int mzColumn = -1;
          int mzToleranceColumn = -1;
          int peakColumn = -1;
          int lowerValleyColumn = -1;
          int upperValleyColumn = -1;
          int lowMzColumn = -1;
          int upMzColumn = -1;
          int ellipseTimePosColumn = -1;
          int ellipseMzPosColumn = -1;
          int ellipseTimeStretchColumn = -1;
          int ellipseMzStretchColumn = -1;
          int lowerHardLimitColumn = -1;
          int upperHardLimitColumn = -1;
          int percentalSplitColumn = -1;
          int apexIntensityColumn = -1;
          int lowerValley10PcColumn = -1;
          int lowerValley50PcColumn = -1;
          int upperValley10PcColumn = -1;
          int upperValley50PcColumn = -1;
          int lowerMz10PcColumn = -1;
          int lowerMz50PcColumn = -1;
          int upperMz10PcColumn = -1;
          int upperMz50PcColumn = -1;
          
          
          int msLevel=1;
          LipidParameterSet params = null;
          boolean showModification = false;
          Hashtable<String,String> analyteNames = new Hashtable<String,String>();
          for (int rowCount=0;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
            Row row = sheet.getRow(rowCount);
            String name = null;
            int dbs = -1;
            int oh = LipidomicsConstants.EXCEL_NO_OH_INFO;
            int paramCharge = 1;
            String modification = null;
            String formula = null;
            String modFormula = null;
            String rtString = "";
            float area = 0f;
            float areaError = 0f;
            float background = 0f; 
            int charge = -1;
            float mz = 0f;
            float mzTolerance = 0f;
            float peak = 0f;
            float lowerValley = 0f;
            float upperValley = 0f;
            float apexIntensity = 0f;
            float lowerValley10Pc = 0f;
            float upperValley10Pc = 0f;
            float lowerValley50Pc = 0f;
            float upperValley50Pc = 0f;
            float lowerMz10Pc = -1f;
            float upperMz10Pc = -1f;
            float lowerMz50Pc = -1f;
            float upperMz50Pc = -1f;

            int isotope = -1;
            float lowMz = -1;
            float upMz = -1;
            float ellipseTimePosition = -1f;
            float ellipseMzPosition = -1f;
            float ellipseTimeStretch = -1f;
            float ellipseMzStretch = -1f;
            float lowerRtHardLimit = -1f;
            float upperRtHardLimit = -1f;
            float percentalSplit = -1f;
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
              }
              
              if (cellType == Cell.CELL_TYPE_STRING)
                contents = cell.getStringCellValue();
              else if (cellType == Cell.CELL_TYPE_NUMERIC){
                double cellValue = -1;
                cellValue = cell.getNumericCellValue();
                contents = String.valueOf(cellValue);
              }
              if (rowCount==0){
                if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_NAME))
                  nameColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_DBS))
                  dbsColumn = i;
                else if (contents.equalsIgnoreCase(LipidomicsConstants.EXCEL_MS_OH))
                  ohColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_MODIFICATION))
                  modificationColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_FORMULA))
                  formulaColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_MOD_FORMULA))
                  modFormulaColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_RT))
                  rtColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_ISOTOPE))
                  isotopeColumn = i;            
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_AREA))
                  areaColumn = i;            
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_AREA_ERROR))
                  areaErrorColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_BACKGROUND))
                  backgroundColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_CHARGE))
                  chargeColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_MZ_MS1))
                  mzColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_MZ_TOLERANCE))
                  mzToleranceColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_PEAK))
                  peakColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_LOWER_VALLEY))
                  lowerValleyColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_UPPER_VALLEY))
                  upperValleyColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_LOWER_MZ))
                  lowMzColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_UPPER_MZ))
                  upMzColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_ELL_CENT_TIME))
                  ellipseTimePosColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_ELL_CENT_MZ))
                  ellipseMzPosColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_ELL_STRETCH_TIME))
                  ellipseTimeStretchColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_ELL_STRETCH_MZ))
                  ellipseMzStretchColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_LOWER_RT_HARD_LIMIT))
                  lowerHardLimitColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_UPPER_RT_HARD_LIMIT))
                  upperHardLimitColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_PERCENTAL_SPLIT))
                  percentalSplitColumn = i;
                else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_RAW_APEX)){
                  apexIntensityColumn = i;
                }else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_LOWER_VALLEY10PC)){
                  lowerValley10PcColumn = i;
                }else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_LOWER_VALLEY50PC)){
                  lowerValley50PcColumn = i;
                }else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_UPPER_VALLEY10PC)){
                  upperValley10PcColumn = i;
                }else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_UPPER_VALLEY50PC)){
                  upperValley50PcColumn = i;
                }else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_LOWER_MZ10PC)){
                  lowerMz10PcColumn = i;
                }else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_LOWER_MZ50PC)){
                  lowerMz50PcColumn = i;
                }else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_UPPER_MZ10PC)){
                  upperMz10PcColumn = i;
                }else if (contents.equalsIgnoreCase(QuantificationResultExporter.HEADER_UPPER_MZ50PC)){
                  upperMz50PcColumn = i;
                }else if (contents.startsWith(QuantificationResultExporter.HEADER_MS_LEVEL)){
                  String levelString = contents.substring(QuantificationResultExporter.HEADER_MS_LEVEL.length()).trim();
                  msLevel = Integer.valueOf(levelString);
                }
                
              }else{
                if (i==nameColumn)
                  name = contents;
                if (i==dbsColumn&&contents!=null&&contents.length()>0){
                  dbs = numeric.intValue();
                } else if (i==ohColumn&&contents!=null&&contents.length()>0){
                  oh = numeric.intValue();
                }if (i==modificationColumn)
                  modification = contents;
                if (i==formulaColumn)
                  formula = contents;
                if (i==modFormulaColumn)
                  modFormula = contents;
                if (i==rtColumn)
                  rtString = contents;
                if (i==isotopeColumn && contents!=null&&contents.length()>0){
                  isotope = numeric.intValue();
                }  
                if (i==chargeColumn&&contents!=null&&contents.length()>0)
                  paramCharge = numeric.intValue();
                if (i==areaColumn && contents!=null && contents.length()>0)
                  area = numeric.floatValue();
                if (i==areaErrorColumn && contents!=null && contents.length()>0)
                  areaError = numeric.floatValue();
                if (i==backgroundColumn && contents!=null && contents.length()>0)
                  background = numeric.floatValue();
                if (i==chargeColumn && contents!=null && contents.length()>0){
                  charge = numeric.intValue();
                }  
                if (i==mzColumn && contents!=null && contents.length()>0)
                  mz = numeric.floatValue();
                if (i==mzToleranceColumn && contents!=null && contents.length()>0)
                  mzTolerance = numeric.floatValue();
                if (i==peakColumn && contents!=null && contents.length()>0)
                  peak = numeric.floatValue();
                if (i==lowerValleyColumn && contents!=null && contents.length()>0)
                  lowerValley = numeric.floatValue();
                if (i==apexIntensityColumn && contents!=null && contents.length()>0)
                  apexIntensity = numeric.floatValue();
                if (i==lowerValley10PcColumn && contents!=null && contents.length()>0)
                  lowerValley10Pc = numeric.floatValue();
                if (i==lowerValley50PcColumn && contents!=null && contents.length()>0)
                  lowerValley50Pc = numeric.floatValue();
                if (i==upperValley10PcColumn && contents!=null && contents.length()>0)
                  upperValley10Pc = numeric.floatValue();
                if (i==upperValley50PcColumn && contents!=null && contents.length()>0)
                  upperValley50Pc = numeric.floatValue();               
                if (i==lowerMz10PcColumn && contents!=null && contents.length()>0)
                  lowerMz10Pc = numeric.floatValue();
                if (i==lowerMz50PcColumn && contents!=null && contents.length()>0)
                  lowerMz50Pc = numeric.floatValue();
                if (i==upperMz10PcColumn && contents!=null && contents.length()>0)
                  upperMz10Pc = numeric.floatValue();
                if (i==upperMz50PcColumn && contents!=null && contents.length()>0)
                  upperMz50Pc = numeric.floatValue();               
                
                if (i==upperValleyColumn && contents!=null && contents.length()>0)
                  upperValley = numeric.floatValue();
                if (i==lowMzColumn && contents!=null && contents.length()>0)
                  lowMz = numeric.floatValue();
                if (i==upMzColumn && contents!=null && contents.length()>0)
                  upMz = numeric.floatValue();
                if (i==ellipseTimePosColumn && contents!=null && contents.length()>0)
                  ellipseTimePosition = numeric.floatValue();
                if (i==ellipseMzPosColumn && contents!=null && contents.length()>0)
                  ellipseMzPosition = numeric.floatValue();
                if (i==ellipseTimeStretchColumn && contents!=null && contents.length()>0)
                  ellipseTimeStretch = numeric.floatValue();
                if (i==ellipseMzStretchColumn && contents!=null && contents.length()>0)
                  ellipseMzStretch = numeric.floatValue();
                
                if (i==lowerHardLimitColumn && contents!=null && contents.length()>0)
                  lowerRtHardLimit = numeric.floatValue();
                if (i==upperHardLimitColumn && contents!=null && contents.length()>0)
                  upperRtHardLimit = numeric.floatValue();
                if (i==percentalSplitColumn && contents!=null && contents.length()>0)
                  percentalSplit = numeric.floatValue();

              }
            }
            if (name!=null&&name.length()>0){
              if (params!=null){
                // this is for backward compatibility
                if (params.ProbeCount()>0)
                  params.setCharge(params.Probe(0).Charge);
                resultPrms.add(params);
                if (analyteNames.containsKey(params.getNameString())) showModification = true;
                analyteNames.put(params.getNameString(), params.getNameString());
              }
              //this is for backward compatibility
              if (modificationColumn==-1 || formulaColumn==-1 || modFormulaColumn==-1){
                Object[] components = ComparativeNameExtractor.splitOldNameStringToComponents(name);
                name = (String)components[0];
                dbs = (Integer)components[1];
                formula = (String)components[2];
                modification = "";
                modFormula = "";
              }
              params = new LipidParameterSet(mz, name, dbs, oh, modification, rtString, formula, modFormula,paramCharge);
              if (lowerRtHardLimit>=0) params.setLowerRtHardLimit(lowerRtHardLimit);
              if (upperRtHardLimit>=0) params.setUpperRtHardLimit(upperRtHardLimit);
              if (percentalSplit>=0) params.setPercentalSplit(percentalSplit);
              params.Area = area;
              params.LowerMzBand = mzTolerance;
              params.UpperMzBand = mzTolerance;
            }else{
              if (params!=null){
                if (charge!=-1){
                  CgProbe probe = new CgProbe(0,charge);
                  if (area>0 || params.getLowerRtHardLimit()>=0 || params.getUpperRtHardLimit()>=0){
                    probe.AreaStatus = CgAreaStatus.OK;
                    probe.Area = area;
                    probe.AreaError = areaError;
                    probe.Background = background;
                    probe.Peak = peak;
                    probe.LowerValley = lowerValley;
                    probe.UpperValley = upperValley;
                    if (upperValley10Pc>0f){
                      probe.setApexIntensity(apexIntensity);
                      probe.setLowerValley10(lowerValley10Pc);
                      probe.setLowerValley50(lowerValley50Pc);
                      probe.setUpperValley10(upperValley10Pc);
                      probe.setUpperValley50(upperValley50Pc);
                    }
                  }else{
                    probe.AreaStatus = CgAreaStatus.TooSmall;
                  }            
                  probe.Mz = mz;
                  probe.LowerMzBand = mzTolerance;
                  probe.UpperMzBand = mzTolerance;
                  probe.isotopeNumber = isotope;
                  if (ellipseTimePosition>0&&ellipseMzPosition>0&&ellipseTimeStretch>0&&ellipseMzStretch>0&&
                      lowMz>0&&upMz>0){
                    Probe3D probe3D = new Probe3D(probe,ellipseTimePosition,ellipseMzPosition,
                        ellipseTimeStretch,ellipseMzStretch,-1f,-1f,lowerMz10Pc,upperMz10Pc,lowerMz50Pc,upperMz50Pc);
                    probe3D.LowerMzBand = lowMz;
                    probe3D.UpperMzBand = upMz;
                    if (params.getLowerRtHardLimit()>=0) probe3D.setLowerHardRtLimit(params.getLowerRtHardLimit());
                    if (params.getUpperRtHardLimit()>=0) probe3D.setUpperHardRtLimit(params.getUpperRtHardLimit());
                    params.AddProbe(probe3D);
                  }else
                    params.AddProbe(probe);
                }
              }
            }
        
          }
          if (params!=null){
            // this is for backward compatibility
            if (params.ProbeCount()>0)
              params.setCharge(params.Probe(0).Charge);
            resultPrms.add(params);
            if (analyteNames.containsKey(params.getNameString())) showModification = true;
            analyteNames.put(params.getNameString(), params.getNameString());
          }
          if (!showModification){
            String modificationString = null;
            for (LipidParameterSet set : resultPrms){
              if (modificationString==null) modificationString = set.getModificationName();
              if (!modificationString.equalsIgnoreCase(set.getModificationName())){
                showModification = true;
                break;
              }
            }
          }
          resultParams.put(sheet.getSheetName(), resultPrms);
          showModifications.put(sheet.getSheetName(), showModification);
          msLevels.put(sheet.getSheetName(), msLevel);
        } else if (sheet.getSheetName().equalsIgnoreCase(QuantificationResultExporter.SHEET_CONSTANTS)){
          Object[] settings = LipidomicsConstants.readSettingsFromExcelApachePOI(sheet);
          readConstants = (LipidomicsConstants)settings[0];
          faOhEncodings = (HydroxyEncoding)settings[1];
          lcbOhEncodings = (HydroxyEncoding)settings[2];
        }
      }
      for (int sheetNumber=0;sheetNumber!=workbook.getNumberOfSheets();sheetNumber++){  
        Sheet sheet = workbook.getSheetAt(sheetNumber);
        if (sheet.getSheetName().endsWith(QuantificationResultExporter.ADDUCT_MSN_SHEET)){
          if (specificClass!=null && !sheet.getSheetName().equalsIgnoreCase((specificClass+QuantificationResultExporter.ADDUCT_MSN_SHEET)))
            continue;
          String lipidClass = sheet.getSheetName().substring(0,sheet.getSheetName().lastIndexOf(QuantificationResultExporter.ADDUCT_MSN_SHEET));
          Vector<LipidParameterSet> resultPrms = resultParams.get(lipidClass);
          resultPrms = readMSnEvidence(sheet,resultPrms,readConstants,faOhEncodings,lcbOhEncodings);
          resultParams.put(lipidClass,resultPrms);
        } 
//        else if (sheet.getSheetName().endsWith(QuantificationResultExporter.ADDUCT_OMEGA_SHEET)) {
//          if (specificClass!=null && !sheet.getSheetName().equalsIgnoreCase((specificClass+QuantificationResultExporter.ADDUCT_OMEGA_SHEET)))
//            continue;
//          String lipidClass = sheet.getSheetName().replace(QuantificationResultExporter.ADDUCT_OMEGA_SHEET, "");
//          Vector<LipidParameterSet> resultPrms = resultParams.get(lipidClass);
//          readOmegaAssignments(sheet,resultPrms,faOhEncodings,lcbOhEncodings);
//          LipidParameterSet.setOmegaInformationAvailable(true);
//        }
      }
      myxls.close();
    } catch (IOException e) {
      e.printStackTrace();
      new WarningMessage(new JFrame(), "ERROR", e.getMessage()+"; it does not seem to be Microsoft Excel");
      throw new ExcelInputFileException(e);
    } catch (Exception e) {
      e.printStackTrace();
      new WarningMessage(new JFrame(), "ERROR", e.getMessage());
      throw new ExcelInputFileException(e);
    }
    
    return new QuantificationResult(resultParams,readConstants,msLevels,faOhEncodings,lcbOhEncodings);
  }
  
  
  /**
   * Reads double bond position evidence from an Excel sheet. The results are stored in a Vector of DoubleBondPositionVOs for each LipidParameterSet.
   * @param sheet Omega Excel sheet
   * @param params the LipidParameterSets for this analyte class
   * @param faHydroxyEncoding the character encoding of the number of hydroxylation sites for the
   * @param lcbHydroxyEncoding the character encoding of the number of hydroxylation sites for the LCB
   */
//  public static void readOmegaAssignments(
//      Sheet sheet, Vector<LipidParameterSet> params, HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding) {
//    
//    Hashtable<String,LipidParameterSet> msHash = new Hashtable<String,LipidParameterSet>();
//    for (LipidParameterSet param : params){
//      msHash.put(param.getNamePlusModHumanReadable(), param);
//    }
//    
//    Map<String, Integer> headerMap = new HashMap<String,Integer>();
//    Row headerRow = sheet.getRow(0);
//    for (Cell cell : headerRow) {
//      headerMap.put(cell.getStringCellValue(),cell.getColumnIndex());
//    }
//    
//    String identifier = null;
//    
//    for (Row row : sheet) {
//      if (row.getRowNum() == 0) continue;
//      
//      int lastCellNum = row.getLastCellNum();
//      String molecularSpecies = null;
//      String doubleBondPosition = null;
//      int accuracy = 0;
//      boolean isAssigned = false;
//      float expectedRetentionTime = 0f;
//      
//      for (int i=0;  row!=null && i<lastCellNum; i++){
//        Cell cell = row.getCell(i);
//        String contentsString = "";
//        Double contentsNumeric = null;
//        boolean contentsBoolean = false;
//        if (cell!=null) {
//          switch (cell.getCellType()) {
//            case Cell.CELL_TYPE_STRING:
//              contentsString = cell.getStringCellValue();
//              break;
//            case Cell.CELL_TYPE_NUMERIC:
//              contentsNumeric = cell.getNumericCellValue();
//              break;
//            case Cell.CELL_TYPE_BOOLEAN:
//              contentsBoolean = cell.getBooleanCellValue();
//              break;
//          }
//        }
//        
//        if (i == headerMap.get(QuantificationResultExporter.HEADER_IDENTIFIER)) {
//          if (!contentsString.equals("")) {
//            identifier = contentsString;
//          }
//        } else if (i == headerMap.get(QuantificationResultExporter.HEADER_MOLECULAR_SPECIES)) {
//          molecularSpecies = contentsString;
//        } else if (i == headerMap.get(QuantificationResultExporter.HEADER_DOUBLE_BOND_POSITION_LEVEL)) {
//          doubleBondPosition = contentsString;
//        } else if (i == headerMap.get(QuantificationResultExporter.HEADER_EXPECTED_RT)) {
//          expectedRetentionTime = contentsNumeric.floatValue();
//        } else if (i == headerMap.get(QuantificationResultExporter.HEADER_ACCURACY)) {
//          accuracy = contentsNumeric.intValue();
//        } else if (i == headerMap.get(QuantificationResultExporter.HEADER_ASSIGNED)) {
//          isAssigned = contentsBoolean;
//        }
//
//      }
//      
//      try {
//        Vector<FattyAcidVO> chainCombination = StaticUtils.decodeFAsFromHumanReadableName(
//            doubleBondPosition, Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), false);
//        
//        DoubleBondPositionVO doubleBondPositionVO = new DoubleBondPositionVO(
//            chainCombination, expectedRetentionTime, accuracy, molecularSpecies, isAssigned);
//        
//        msHash.get(identifier).addOmegaInformation(doubleBondPositionVO);
//      } catch (LipidCombinameEncodingException ex) {
//        System.out.println(ex.getMessage());
//      }
//      
//    }
//  }
  
  /**
   * reads MSn evidence from Excel sheet - where applicable, the results are stored in an LipidomicsMSnSet - LipidParameterSets and LipidomicsMSnSets are returned in the vector
   * @param sheet the Excel sheet to be read
   * @param ms1Results the results from MS1 Excel reading
   * @param readConstants the lipidomics constants used for this file (for checking if this file was quantified using an Alex123 target list)
   * @param faHydroxyEncoding the character encoding of the number of hydroxylation sites for the
   * @param lcbHydroxyEncoding the character encoding of the number of hydroxylation sites for the LCB
   * @return Vector containing MS1 (LipidParameterSet) and  MS2 (LipidomicsMSnSe) results
   * @throws RulesException
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static Vector<LipidParameterSet> readMSnEvidence(Sheet sheet, Vector<LipidParameterSet> ms1Results, LipidomicsConstants readConstants,
      HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding) throws RulesException, LipidCombinameEncodingException {
    Hashtable<String,LipidParameterSet> msHash = new Hashtable<String,LipidParameterSet>();
    for (LipidParameterSet ms1 : ms1Results){
      msHash.put(ms1.getNamePlusModHumanReadable(), ms1);
    }
    
    LipidParameterSet addingMSnEvidence = null;
    Hashtable<Integer,String> columnToIdentification = new Hashtable<Integer,String>();
    Hashtable<String,Double> relativeAreas = new Hashtable<String,Double>();
    String regex = "MS(\\d+) scan RTs";
    Pattern msLevelPattern =  Pattern.compile(regex);
    Hashtable<Integer,LinkedHashMap<Integer,Float>> msnRetentionTimes = new Hashtable<Integer,LinkedHashMap<Integer,Float>>();
    boolean checkMSnAreas = false;
    boolean headGroupFragmentActive = false;
    boolean headGroupRules = false;
    boolean chainFragmentActive = false;
    boolean chainRules = false;
    boolean positionRules = false;
    String combiKey = "";

    // is the header row read and the column indices initialized 
    boolean headerRowRead = false;
    // the header row has to be read
    boolean readFragmentHeaderRow = false;

    // is the intensity header read and the column indices initialized 
    boolean intensityHeaderRead = false;
    // the intensity header row has to be read
    boolean readIntensityHeaderRow = false;
    
    
    // column indices for the fragment information rows
    int nameColumn = -1;
    int ohColumn = -1;
    int chainTypeColumn = -1;
    int formulaColumn = -1;
    int msLevelColumn = -1;
    int chargeColumn = -1;
    int mzColumn = -1;
    int mzTolColumn = -1;
    int areaColumn = -1;
    int peakColumn = -1;
    int startTimeColumn = -1;
    int stopTimeColumn = -1;
    int startMzColumn = -1;
    int stopMzColumn = -1;
    int ellTimeColumn = -1;
    int ellMzColumn = -1;
    int ellTimeRangeColumn = -1;
    int ellMzRangeColumn = -1;
    
    // column indices for the intensity information
    int ruleColumn = -1;
    int originalRuleColumn = -1;
    int ruleValuesColumn = -1;
    int ruleMissedColumn = -1;

        
    // the information to be written into the returning VO
    int status = LipidomicsMSnSet.NO_MSN_PRESENT;
    float mzTolerance = -1f;
    Hashtable<String,CgProbe> headGroupFragments = new Hashtable<String,CgProbe>();
    Hashtable<String,IntensityRuleVO> headIntensityRules = new Hashtable<String,IntensityRuleVO>();
    Hashtable<String,Hashtable<String,CgProbe>> chainFragments = new Hashtable<String,Hashtable<String,CgProbe>>();
    Hashtable<String,Hashtable<String,IntensityChainVO>> chainIntensityRules = new Hashtable<String,Hashtable<String,IntensityChainVO>>();
    Vector<String> validChainCombinations = new Vector<String>();
    Hashtable<String,Hashtable<Integer,Integer>> positionDefinition = new Hashtable<String,Hashtable<Integer,Integer>>();
    Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> positionEvidence = new Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>();
    Hashtable<Integer,Float> basePeakValues = new Hashtable<Integer,Float>();
    
    Hashtable<String,String> uniqueRules = new Hashtable<String,String>();
    int numberOfPositions = -1;
    boolean containsAlkyl = false;
    boolean containsAlkenyl = false;
    boolean usedAlexMsnTargets = false;
    
    char knownPosSep = LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS.toCharArray()[0];
    char unknownPosSep = LipidomicsConstants.CHAIN_SEPARATOR_NO_POS.toCharArray()[0];
    
    for (int rowCount=0;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
      Hashtable<Integer,Object> cellEntries =  ExcelUtils.getEntriesOfOneRow(sheet.getRow(rowCount),false);
      // check if an Alex123 target list was used for the MSn fragments of this class
      if (addingMSnEvidence==null && cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) && (cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) instanceof String) &&
          ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().startsWith(QuantificationResultExporter.HEADER_ALEX123_MSN_TARGETS_USED)){
        StringTokenizer tokenizer = new StringTokenizer(((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)),"=");
        if (tokenizer.countTokens()!=2) continue;
        tokenizer.nextToken();
        String value = tokenizer.nextToken().trim();
        if (value.equalsIgnoreCase("true") || value.equalsIgnoreCase("yes"))
          usedAlexMsnTargets = true;
      // reading the first row, containing the sum formula, and the individual identifications
      } else if (addingMSnEvidence==null && cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) && cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_FORMULA)){
        // for every new lipid MS1 species, the parameters holding the information have to be initialized
        status = LipidomicsMSnSet.NO_MSN_PRESENT;
        mzTolerance = -1f;
        headGroupFragments = new Hashtable<String,CgProbe>();
        headIntensityRules = new Hashtable<String,IntensityRuleVO>();
        chainFragments = new Hashtable<String,Hashtable<String,CgProbe>>();
        chainIntensityRules = new Hashtable<String,Hashtable<String,IntensityChainVO>>();
        validChainCombinations = new Vector<String>();
        positionDefinition = new Hashtable<String,Hashtable<Integer,Integer>>();
        positionEvidence = new Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>();
        basePeakValues = new Hashtable<Integer,Float>();
        
        columnToIdentification = new Hashtable<Integer,String>();
        String speciesName = (String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME);
        addingMSnEvidence = msHash.get(speciesName);
        int count = QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME+1;
        numberOfPositions = -1;
        containsAlkyl = false;
        containsAlkenyl = false;
        while (cellEntries.containsKey(count)){
          String lipidIdentification = (String)cellEntries.get(count);
          if (lipidIdentification!=null && lipidIdentification.length()>0){
            //this put has to be here to store the MSn scan numbers
            columnToIdentification.put(count, lipidIdentification);
            // this if was extended by "lipidIdentification.indexOf(":")!=-1" to support single chains - I hope this has no negative side effects
            if (lipidIdentification.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS)!=-1||lipidIdentification.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)!=-1||
                lipidIdentification.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_DBS)!=-1){
              boolean isAlexOhEncoding = false;
              //if there is a ";", the next character after the following numbers must be a ":"; if there is a "/" or a "_" it is the OH index of an Alex123 notation
              if (usedAlexMsnTargets){
                if (lipidIdentification.indexOf(LipidomicsConstants.CHAIN_COMBI_SEPARATOR_AMBIG_POS)!=-1){
                  boolean makeSubstring = false;
                  char[] chars = lipidIdentification.toCharArray();
                  for (int i=lipidIdentification.indexOf(LipidomicsConstants.CHAIN_COMBI_SEPARATOR_AMBIG_POS)+1; i!=chars.length; i++){
                    //it is the OH index of an Alex123 notation
                    if (chars[i]==knownPosSep || chars[i]==unknownPosSep) {
                      isAlexOhEncoding = true;
                      break;
                    //it is an LDA identification where the position is unknown -> cut
                    }else if (chars[i]==':'){
                      makeSubstring = true;
                      break;
                    }else
                      continue;
                  }
                  if (makeSubstring)
                    lipidIdentification = lipidIdentification.substring(0,lipidIdentification.indexOf(LipidomicsConstants.CHAIN_COMBI_SEPARATOR_AMBIG_POS));
                }
              }else{
                if (lipidIdentification.indexOf(LipidomicsConstants.CHAIN_COMBI_SEPARATOR_AMBIG_POS)!=-1)lipidIdentification = lipidIdentification.substring(0,lipidIdentification.indexOf(LipidomicsConstants.CHAIN_COMBI_SEPARATOR_AMBIG_POS));
              }
              /** remove omega-DB position annotations if present */
              lipidIdentification = StaticUtils.getHumanReadableWODoubleBondPositions(lipidIdentification);
              Object[] nameAndNumberOfPos = StaticUtils.cleanEmptyFAPositionsAndEncodeToLDACombiName(lipidIdentification,faHydroxyEncoding,
                  lcbHydroxyEncoding, isAlexOhEncoding);
              lipidIdentification = (String)nameAndNumberOfPos[0];
              int posNr = (Integer)nameAndNumberOfPos[1];
              if (posNr>numberOfPositions) numberOfPositions = posNr;
              containsAlkyl = (Boolean)nameAndNumberOfPos[2];
              containsAlkenyl = (Boolean)nameAndNumberOfPos[3];
              //when the correct MSn LDA decoded MSn identification name is known, put the correct one to the column
              columnToIdentification.put(count, lipidIdentification);
              validChainCombinations.add(lipidIdentification);
              //System.out.println("0. lipidIdentification: "+sheet.getSheetName()+" "+lipidIdentification);
            }
          }
          count++;
        }
        checkMSnAreas = true;
        headGroupFragmentActive = false;
        headGroupRules = false;
        chainFragmentActive = false;
        chainRules = false;
        positionRules = false;

      }
      // reading the second row, containing the relative areas of the analytes and the retention times
      else if (checkMSnAreas && addingMSnEvidence!=null){
        relativeAreas = new Hashtable<String,Double>();
        msnRetentionTimes = new Hashtable<Integer,LinkedHashMap<Integer,Float>>();
        double totalArea = (Double)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME);
        for (Integer column : columnToIdentification.keySet()){
          String lipidIdentification = columnToIdentification.get(column);
          //this is for the retention times
          Matcher msLevelMatcher = msLevelPattern.matcher(lipidIdentification);
          if (msLevelMatcher.matches()){
            int msLevel =  Integer.parseInt(msLevelMatcher.group(1));
            LinkedHashMap<Integer,Float> rts = new LinkedHashMap<Integer,Float>();
            StringTokenizer rtTokenizer = null;
            if (cellEntries.get(column) instanceof Double)
              rtTokenizer = new StringTokenizer(((Double)cellEntries.get(column)).toString(),";");
            else 
              rtTokenizer = new StringTokenizer((String)cellEntries.get(column),";");
            //this is for backward compatibility
            int count = 1;
            while (rtTokenizer.hasMoreTokens()){
              String kvPairString = rtTokenizer.nextToken();
              int scanNr = -1;
              float rt = -1;
              //this is for backward compatibility
              if (kvPairString.indexOf("=")==-1){
                scanNr = count;
                rt = new Float(kvPairString);
                count++;
              }else{
                String[] kvPair = kvPairString.split("=");
                scanNr = Integer.parseInt(kvPair[0]);
                rt = Float.parseFloat(kvPair[1]);
              }
              rts.put(scanNr, rt);
            }
            msnRetentionTimes.put(msLevel, rts);
          //this is for the relative areas
          }else
            relativeAreas.put(lipidIdentification, ((Double)cellEntries.get(column))/totalArea);
        }
        checkMSnAreas = false;
      }
      // the following if is for activating the head group fragments section
      else if (!headGroupFragmentActive && addingMSnEvidence!=null && cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)
          && cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) instanceof String 
          && ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_SECTION_HEAD_FRAGMENTS)){
        headGroupFragmentActive = true;
        headGroupRules = false;
        chainFragmentActive = false;
        chainRules = false;
        positionRules = false;
        readFragmentHeaderRow = true;
        headerRowRead = false;
        intensityHeaderRead = false;
      } 
      // the following if is for activating the head group rules section
      else if (!headGroupRules && addingMSnEvidence!=null && cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)
          && cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) instanceof String 
          && ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_SECTION_HEAD_INTENSITIES)){
        headGroupFragmentActive = false;
        headGroupRules = true;
        chainFragmentActive = false;
        chainRules = false;
        positionRules = false;
        headerRowRead = false;
        readIntensityHeaderRow = true;
        intensityHeaderRead = false;
        uniqueRules = new Hashtable<String,String>();
      }

      // the following if is for activating the chain fragments section
      else if (!chainFragmentActive && addingMSnEvidence!=null && cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)
          && cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) instanceof String 
          && ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_SECTION_CHAIN_FRAGMENTS)){
        headGroupFragmentActive = false;
        headGroupRules = false;
        chainFragmentActive = true;
        chainRules = false;
        positionRules = false;
        readFragmentHeaderRow = true;
        headerRowRead = false; 
        intensityHeaderRead = false;
      }
      // the following if is for activating the chain fragments section
      else if (!chainRules && addingMSnEvidence!=null && cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)
          && cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) instanceof String 
          && ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_SECTION_CHAIN_INTENSITIES)){
        headGroupFragmentActive = false;
        headGroupRules = false;
        chainFragmentActive = false;
        chainRules = true;
        positionRules = false;
        headerRowRead = false;
        readIntensityHeaderRow = true;
        intensityHeaderRead = false;  
        uniqueRules = new Hashtable<String,String>();
      }
      // the following ifs are for activating the position information section
      else if (addingMSnEvidence!=null && cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)
          && cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) instanceof String 
          && ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().startsWith(LipidomicsConstants.EXCEL_MSN_SECTION_POSITION_INTENSITIES)){
        headGroupFragmentActive = false;
        headGroupRules = false;
        chainFragmentActive = false;
        chainRules = false;
        positionRules = true;
        headerRowRead = false;
        readIntensityHeaderRow = true;
        intensityHeaderRead = false;
        combiKey = ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().substring(LipidomicsConstants.EXCEL_MSN_SECTION_POSITION_INTENSITIES.length());
        combiKey = combiKey.substring(combiKey.indexOf("(")+1,combiKey.indexOf(")"));
        if (combiKey.indexOf(";")!=-1) combiKey = combiKey.substring(0,combiKey.indexOf(";"));
        validChainCombinations = correctUndefinedChainCombinations(combiKey,validChainCombinations,relativeAreas,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
        Hashtable<String,Integer> chainOccurenceInCombi = new Hashtable<String,Integer>();
        combiKey = (String)StaticUtils.cleanEmptyFAPositionsAndEncodeToLDACombiName(combiKey,faHydroxyEncoding,
            lcbHydroxyEncoding,  (usedAlexMsnTargets && combiKey.indexOf(LipidomicsConstants.CHAIN_COMBI_SEPARATOR_AMBIG_POS)!=-1))[0];
        combiKey = getPermutationThatIsInValidChainCombinations(combiKey, validChainCombinations);
        
        Vector<String> fas = StaticUtils.splitChainCombiToEncodedStrings(combiKey,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
        for (String fa: fas){
          int count = 0;
          if (chainOccurenceInCombi.containsKey(fa)) count = chainOccurenceInCombi.get(fa);
          count++;
          chainOccurenceInCombi.put(fa, count);
        }
        // there is only one type of fatty acid chain
        if (chainOccurenceInCombi.size()==1 && fas.size()==numberOfPositions){
          Hashtable<Integer,Integer> definitions = new Hashtable<Integer,Integer>();
          for (int i=0; i!=fas.size();i++){
            definitions.put(i, i);
          }
          positionDefinition.put(combiKey, definitions);
          positionEvidence.put(combiKey, new Hashtable<Integer,Vector<IntensityPositionVO>>());
          uniqueRules = new Hashtable<String,String>();
          continue;
        }
      }
      // the final procedure for creating a LipidomicsMSnSet after the information was read
      else if (addingMSnEvidence!=null && !cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)){
        String speciesName = addingMSnEvidence.getNamePlusModHumanReadable();
//        for (positionDefinition)
        positionDefinition = cleanPositionDefinition(positionDefinition);
        if (isDefinitionPresent(positionDefinition)) status = LipidomicsMSnSet.POSITION_DETECTED;
        addingMSnEvidence = new LipidomicsMSnSet(addingMSnEvidence, status, mzTolerance,headGroupFragments, headIntensityRules,
        chainFragments, chainIntensityRules, validChainCombinations, relativeAreas,  positionDefinition,positionEvidence, numberOfPositions, basePeakValues,
        msnRetentionTimes, faHydroxyEncoding, lcbHydroxyEncoding);
        msHash.put(speciesName,addingMSnEvidence);
        addingMSnEvidence=null;
      }
      // reading the fragment header row for assigning the columns
      else if (readFragmentHeaderRow){
        for (Integer columnId : cellEntries.keySet()){
          if (cellEntries.get(columnId)==null || ((String)cellEntries.get(columnId)).length()==0) continue;
          String entry = (String)cellEntries.get(columnId);
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_NAME)) nameColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_OH)) ohColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_CHAIN_TYPE)) chainTypeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_FORMULA)) formulaColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_MSLEVEL)) msLevelColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_CHARGE)) chargeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ)) mzColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_TOLERANCE)) mzTolColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_AREA)) areaColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_PEAK)) peakColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_TIME_LOWER)) startTimeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_TIME_UPPER)) stopTimeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_LOWER)) startMzColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_UPPER)) stopMzColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_TIME)) ellTimeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_MZ)) ellMzColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_TIME_RANGE)) ellTimeRangeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_MZ_RANGE)) ellMzRangeColumn = columnId;
        }
        readFragmentHeaderRow = false;
        headerRowRead = true;
      }
      // reading the fragment header row for assigning the columns
      else if (readIntensityHeaderRow){
        for (Integer columnId : cellEntries.keySet()){
          if (cellEntries.get(columnId)==null || ((String)cellEntries.get(columnId)).length()==0) continue;
          String entry = (String)cellEntries.get(columnId);
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_INTENSITY_RULE)) ruleColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_INTENSITY_ORIGINAL)) originalRuleColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_INTENSITY_VALUES)) ruleValuesColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_INTENSITY_MISSED)) ruleMissedColumn = columnId;
        }
        readIntensityHeaderRow = false;
        intensityHeaderRead = true;
      } 
      // reading information about the fragment
      else if (addingMSnEvidence!=null && headerRowRead && (headGroupFragmentActive||chainFragmentActive)){
        String fragmentName = null;
        if (nameColumn>-1 && cellEntries.containsKey(nameColumn)) fragmentName = (String)cellEntries.get(nameColumn);
        int oh = 0;
        if (ohColumn>-1 && cellEntries.containsKey(ohColumn)) oh = (int)Math.rint((Double)cellEntries.get(ohColumn));
        // I have to put this to LipidomicsConstants.CHAIN_TYPE_NO_CHAIN for backward compatibility - by this I know that these results were
        // produced by an older LDA version
        short chainType = LipidomicsConstants.CHAIN_TYPE_NO_CHAIN;
        if (chainTypeColumn>-1 && cellEntries.containsKey(chainTypeColumn)) {
          String chainTypeString =  (String)cellEntries.get(chainTypeColumn);
          if (chainTypeString.equalsIgnoreCase(FragmentRuleVO.CHAIN_NAME.substring(1)))
            chainType = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
          else if (chainTypeString.equalsIgnoreCase(FragmentRuleVO.ALKYL_CHAIN_NAME.substring(1)))
            chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
          else if (chainTypeString.equalsIgnoreCase(FragmentRuleVO.ALKENYL_CHAIN_NAME.substring(1)))
            chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
          else if (chainTypeString.equalsIgnoreCase(FragmentRuleVO.LCB_NAME.substring(1)))
            chainType = LipidomicsConstants.CHAIN_TYPE_LCB;
          else
            throw new LipidCombinameEncodingException("The entry \""+chainTypeString+"\" is not allowed in the column \""+LipidomicsConstants.EXCEL_MSN_FRAGMENT_CHAIN_TYPE+"\"!");
        //this is required for backward compatibility
        } else if (chainFragmentActive) {
          chainType = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
          if (containsAlkyl && fragmentName.indexOf("("+LipidomicsConstants.ALKYL_PREFIX)!=-1)
            chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
          if (containsAlkenyl && fragmentName.indexOf("("+LipidomicsConstants.ALKENYL_PREFIX)!=-1)
            chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
        }
        String formula = null;
        if (formulaColumn>-1 && cellEntries.containsKey(formulaColumn)) formula = "+"+((String)cellEntries.get(formulaColumn)).trim().replaceAll(" ", " \\+");
        int msLevel = -1;
        if (msLevelColumn>-1 && cellEntries.containsKey(msLevelColumn)) msLevel = (int)Math.rint((Double)cellEntries.get(msLevelColumn));
        int charge = -1;
        if (chargeColumn>-1 && cellEntries.containsKey(chargeColumn)) charge = (int)Math.rint((Double)cellEntries.get(chargeColumn));
        float mz = -1f;
        if (mzColumn>-1 && cellEntries.containsKey(mzColumn)) mz = ((Double)cellEntries.get(mzColumn)).floatValue();
        if (mzTolerance<0 && mzTolColumn>-1 && cellEntries.containsKey(mzTolColumn)) mzTolerance = ((Double)cellEntries.get(mzTolColumn)).floatValue();
        float area = -1f;
        if (areaColumn>-1 && cellEntries.containsKey(areaColumn)) area = ((Double)cellEntries.get(areaColumn)).floatValue();
        float peak = -1f;
        if (peakColumn>-1 && cellEntries.containsKey(peakColumn)) peak = ((Double)cellEntries.get(peakColumn)).floatValue();
        float startTime = -1f;
        if (startTimeColumn>-1 && cellEntries.containsKey(startTimeColumn)) startTime = ((Double)cellEntries.get(startTimeColumn)).floatValue();
        float stopTime = -1f;
        if (stopTimeColumn>-1 && cellEntries.containsKey(stopTimeColumn)) stopTime = ((Double)cellEntries.get(stopTimeColumn)).floatValue();
        float startMz = -1f;
        if (startMzColumn>-1 && cellEntries.containsKey(startMzColumn)) startMz = ((Double)cellEntries.get(startMzColumn)).floatValue();
        float stopMz = -1f;
        if (stopMzColumn>-1 && cellEntries.containsKey(stopMzColumn)) stopMz = ((Double)cellEntries.get(stopMzColumn)).floatValue();
        float ellipseTimePosition = -1f;
        if (ellTimeColumn>-1 && cellEntries.containsKey(ellTimeColumn)) ellipseTimePosition = ((Double)cellEntries.get(ellTimeColumn)).floatValue();
        float ellipseMzPosition = -1f;
        if (ellMzColumn>-1 && cellEntries.containsKey(ellMzColumn)) ellipseMzPosition = ((Double)cellEntries.get(ellMzColumn)).floatValue();
        float ellipseTimeStretch = -1f;
        if (ellTimeRangeColumn>-1 && cellEntries.containsKey(ellTimeRangeColumn)) ellipseTimeStretch = ((Double)cellEntries.get(ellTimeRangeColumn)).floatValue();
        float ellipseMzStretch = -1f;
        if (ellMzRangeColumn>-1 && cellEntries.containsKey(ellMzRangeColumn)) ellipseMzStretch = ((Double)cellEntries.get(ellMzRangeColumn)).floatValue();
        //TODO: formula is only excluded for the damaged Alex123 target list
        if (fragmentName!=null && fragmentName.length()>0 && /*formula!=null && formula.length()>0 &&*/ msLevel>0 && charge>0 && mz>0f && mzTolerance>0f
            && area>0f && peak>0f && startTime>-1f && stopTime>0f){
          CgProbe probe = new CgProbe(-1, charge, msLevel, formula);
          probe.AreaStatus = CgAreaStatus.OK;
          probe.Mz = mz;
          probe.Area = area;
          probe.AreaError = 0f;
          probe.Background = 0f;
          probe.Peak = peak;
          probe.LowerValley = startTime;
          probe.UpperValley = stopTime;
          probe.LowerMzBand = mzTolerance;
          probe.UpperMzBand = mzTolerance;
          probe.isotopeNumber = 0;
          if (startMz>0f && stopMz>0f && ellipseTimePosition>0f && ellipseMzPosition>0f && ellipseTimeStretch>0f && ellipseMzStretch>0f){
            probe = new Probe3D(probe,ellipseTimePosition,ellipseMzPosition,
                ellipseTimeStretch,ellipseMzStretch,-1f,-1f,-1f,-1f,-1f,-1f);
            probe.LowerMzBand = startMz;
            probe.UpperMzBand = stopMz;
          }
          if (headGroupFragmentActive){
            headGroupFragments.put(fragmentName, probe);
            status = LipidomicsMSnSet.HEAD_GROUP_DETECTED;
          }else if (chainFragmentActive){
            status = LipidomicsMSnSet.FRAGMENTS_DETECTED;
            Object[] faAndFragment = StaticUtils.parseChainFaAndFragmentNameFromExcel(fragmentName,chainType,oh);
            FattyAcidVO faName = (FattyAcidVO)faAndFragment[0];
            fragmentName =  (String)faAndFragment[1];
            Hashtable<String,CgProbe> fragments = new Hashtable<String,CgProbe>();
            if (chainFragments.containsKey(faName.getChainId())) fragments = chainFragments.get(faName.getChainId());
            fragments.put(fragmentName, probe);
            chainFragments.put(faName.getChainId(),fragments);
          }
        }
      }
      // reading information about intensity rules
      else if (addingMSnEvidence!=null && intensityHeaderRead && (headGroupRules||chainRules||positionRules)){
        String readableRuleInterpretation = null;
        if (ruleColumn>-1 && cellEntries.containsKey(ruleColumn)) readableRuleInterpretation = (String)cellEntries.get(ruleColumn);
        String rule = null;
        if (originalRuleColumn>-1 && cellEntries.containsKey(originalRuleColumn)) rule = (String)cellEntries.get(originalRuleColumn);
        String ruleValueInterpretation = null;
        if (ruleValuesColumn>-1 && cellEntries.containsKey(ruleValuesColumn)) ruleValueInterpretation = (String)cellEntries.get(ruleValuesColumn);
        if (readableRuleInterpretation!=null && readableRuleInterpretation.length()>0 && rule!=null && rule.length()>0 &&
            ruleValueInterpretation!=null && ruleValueInterpretation.length()>0){
          String uniqueRuleString = readableRuleInterpretation+";"+rule+";"+ruleValueInterpretation;
          String missedString = "";
          if (ruleMissedColumn>-1 && cellEntries.containsKey(ruleMissedColumn)) missedString = (String)cellEntries.get(ruleMissedColumn);
          Hashtable<String,Short> missed = new Hashtable<String,Short>();
          StringTokenizer tokenizer = new StringTokenizer(missedString,";");
          while (tokenizer.hasMoreTokens()){
            String token = tokenizer.nextToken();
            //this is necessary for position rules
            if (token.indexOf("[")!=-1)
              token = token.substring(0, token.indexOf("["));
            short type = LipidomicsConstants.CHAIN_TYPE_NO_CHAIN;
            if (chainRules||positionRules){
              boolean isAChain = false;
              //is it a chain fragment in Alex notation
              if (token.startsWith("FA ") || token.startsWith("-FA ") || token.startsWith("LCB ") || token.startsWith("-LCB ") || token.startsWith("O-") || token.startsWith("-O-")){
                isAChain = true;
              }else {
                String remnant = new String(readableRuleInterpretation);
                while (remnant.indexOf(token)!=-1) {
                  remnant = remnant.substring(remnant.indexOf(token)+token.length());
                  if (remnant.startsWith("(")) {
                    isAChain = true;
                    break;
                  }
                }
              }
              if (isAChain)
                type = LipidomicsConstants.CHAIN_TYPE_MISSED;
            }           
            missed.put(token, type);
          }
          if (uniqueRules.containsKey(uniqueRuleString)) continue;
          uniqueRules.put(uniqueRuleString, uniqueRuleString);
          int currentSection = -1;
          if (headGroupRules) currentSection = FragRuleParser.HEAD_SECTION;
          else if (chainRules) currentSection = FragRuleParser.CHAINS_SECTION;
          else if (positionRules) currentSection = FragRuleParser.POSITION_SECTION;
          Hashtable<String,Short> head = new Hashtable<String,Short>();
          for (String key : headGroupFragments.keySet()) head.put(key,LipidomicsConstants.CHAIN_TYPE_NO_CHAIN);
          Hashtable<String,Short> chain = new Hashtable<String,Short>();
          for (String chainId : chainFragments.keySet()){
            FattyAcidVO fa = StaticUtils.decodeLipidNameForCreatingCombis(chainId);
            Hashtable<String,CgProbe> fragments = chainFragments.get(chainId);
            for (String name : fragments.keySet()) chain.put(name, fa.getChainType());
          }
          IntensityRuleVO ruleVO = FragRuleParser.extractIntensityVOFromEquation(rule, -1, currentSection, head, chain, null, missed);
          if (chainRules || positionRules){
            IntensityRuleVO ruleInst = null;
            if (chainRules) ruleInst = IntensityChainVO.getFattyAcidsFromReadableRule(readableRuleInterpretation,ruleVO,chainFragments,missed,
                faHydroxyEncoding, lcbHydroxyEncoding, addingMSnEvidence.getOhNumber()>0);
            else if (positionRules) ruleInst = IntensityPositionVO.getFattyAcidsFromReadableRule(readableRuleInterpretation,ruleVO,chainFragments,missed,
                faHydroxyEncoding, lcbHydroxyEncoding, addingMSnEvidence.getOhNumber()>0);
            if (ruleInst!=null) ruleVO = ruleInst;
          }
          if (ruleVO.containsBasePeak()){
            Hashtable<String,Float> values = LipidomicsMSnSet.getFragmentAreas(ruleVO,headGroupFragments,chainFragments);
            for (String name : values.keySet()){
              if (values.get(name)<0) values.put(name, 0f);
            }
            float basePeak = ruleVO.extractBasePeakValue(ruleValueInterpretation,values);
            int msLevel = ruleVO.getMSLevel(headGroupFragments, chainFragments);
            if (msLevel>1) basePeakValues.put(msLevel, basePeak);
          }
          if (ruleVO!=null && headGroupRules) {
            headIntensityRules.put(ruleVO.getRuleIdentifier(), ruleVO);
          }
          else if (ruleVO!=null && chainRules){
            //this line was replaced by the next one since this step was already performed before - I simply casted now the object
            //IntensityChainVO chainVO =  IntensityChainVO.getFattyAcidsFromReadableRule(readableRuleInterpretation,ruleVO,chainFragments,missed);
            IntensityChainVO chainVO = (IntensityChainVO)ruleVO;
            Hashtable<String,IntensityChainVO> rules = new Hashtable<String,IntensityChainVO>();
            Vector<String> ids = new Vector<String>();
            if (chainVO.getChainType()==IntensityChainVO.DIFF_CHAIN_TYPES){
              Vector<String>  chains = new Vector<String>();
              for (FattyAcidVO faVO : chainVO.getParticipatingChains())
                chains.add(faVO.getChainId());
              ids = StaticUtils.getAllAffectedChainCombinations(validChainCombinations,chains);
            } else ids.add(chainVO.getParticipatingChains().get(0).getChainId());
            for (String id : ids) {
              if (chainIntensityRules.containsKey(id)) rules = chainIntensityRules.get(id);
              rules.put(ruleVO.getReadableRuleInterpretation(faHydroxyEncoding,lcbHydroxyEncoding), chainVO);
              chainIntensityRules.put(id, rules);
            }
          } else if (ruleVO!=null && positionRules){
            //this line was replaced by the next one since this step was already perforemed before - I simply casted now the object
            //IntensityPositionVO posVO = IntensityPositionVO.getFattyAcidsFromReadableRule(readableRuleInterpretation,ruleVO,chainFragments,missed);
            IntensityPositionVO posVO = (IntensityPositionVO)ruleVO;
            Hashtable<Integer,Integer> defs = new Hashtable<Integer,Integer>();
            Hashtable<Integer,Vector<IntensityPositionVO>> evidence = new Hashtable<Integer,Vector<IntensityPositionVO>>();
            if (positionDefinition.containsKey(combiKey)) defs = positionDefinition.get(combiKey);
            if (positionEvidence.containsKey(combiKey)) evidence = positionEvidence.get(combiKey);
            Vector<String> fas = StaticUtils.splitChainCombiToEncodedStrings(combiKey,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
            Hashtable<String,String> usedFAs = new Hashtable<String,String>();
            for (int i=0;i!=fas.size();i++){
              String fa = fas.get(i);
              int position = posVO.getPositionByFA(fa);
              if ( position<0 || usedFAs.containsKey(fa)) continue;
              if (defs.containsKey(i)){
                // if the stored assignment is not the same -> no position can be assigned
                if (defs.get(i)!=(position-1)){
                  boolean otherSameFAAtPositionOrUndefined = false;
                  //check if there is the same FA more than once, and this position is assigned or undefined
                  for (int j=0;j!=fas.size();j++){
                    if (i==j || !fa.equalsIgnoreCase(fas.get(j)))continue;
                    if (!defs.containsKey(j)||defs.get(i)==(position-1)||defs.get(j)==(position-1)) otherSameFAAtPositionOrUndefined=true;
                  }
                  if (!otherSameFAAtPositionOrUndefined) defs.remove(i);
                }
              }else{
                //if several rules are fulfilled for the same FA, and this FA exists more than once, the
                //defs would store the same position twice and would not allow another rule with another
                //position for the correct assignment
                boolean sameFaAlreadySamePosition = false;
                for (int j=0; j!=i; j++){
                  if (!fa.equalsIgnoreCase(fas.get(j))) continue;
                  if (defs.get(j)==(position-1)) sameFaAlreadySamePosition = true;
                }
                if (!sameFaAlreadySamePosition){
                  defs.put(i, (position-1));
                  usedFAs.put(fa, fa);
                }
              }
              Vector<IntensityPositionVO> rules = new Vector<IntensityPositionVO>();
              if (evidence.containsKey(position)) rules = evidence.get(position);
              boolean ruleIsThere = false;
              for (IntensityPositionVO other : rules){
                if (other.getReadableRuleInterpretation(faHydroxyEncoding, lcbHydroxyEncoding).equalsIgnoreCase(posVO.getReadableRuleInterpretation(faHydroxyEncoding, lcbHydroxyEncoding))){
                  ruleIsThere = true;
                  break;
                }
              }
              if (!ruleIsThere)rules.add(posVO);
              
              evidence.put(position, rules);
            }
            positionDefinition.put(combiKey, defs);
            positionEvidence.put(combiKey, evidence);
          }
        }
      }

      // if in the next row the fragment head row has to be read - the column indices have to be initialized
      if (readFragmentHeaderRow){
        nameColumn = -1;
        ohColumn = -1;
        chainTypeColumn = -1;
        formulaColumn = -1;
        msLevelColumn = -1;
        chargeColumn = -1;
        mzColumn = -1;
        mzTolColumn = -1;
        areaColumn = -1;
        peakColumn = -1;
        startTimeColumn = -1;
        stopTimeColumn = -1;
        startMzColumn = -1;
        stopMzColumn = -1;
        ellTimeColumn = -1;
        ellMzColumn = -1;
        ellTimeRangeColumn = -1;
        ellMzRangeColumn = -1;
      }
      // if in the next row the intensity head row has to be read - the column indices have to be initialized
      if (readIntensityHeaderRow){
        ruleColumn = -1;
        originalRuleColumn = -1;
        ruleValuesColumn = -1;
      }
    }
    
    if (addingMSnEvidence!=null && (addingMSnEvidence instanceof LipidParameterSet)){
      String speciesName = addingMSnEvidence.getNamePlusModHumanReadable();
      positionDefinition = cleanPositionDefinition(positionDefinition);
      if (isDefinitionPresent(positionDefinition)) status = LipidomicsMSnSet.POSITION_DETECTED;
      addingMSnEvidence = new LipidomicsMSnSet(addingMSnEvidence, status, mzTolerance,headGroupFragments, headIntensityRules,
      chainFragments, chainIntensityRules, validChainCombinations, relativeAreas, positionDefinition,positionEvidence, numberOfPositions, basePeakValues,
      msnRetentionTimes, faHydroxyEncoding, lcbHydroxyEncoding);
      msHash.put(speciesName,addingMSnEvidence);
    }
    
    String lipidClass = sheet.getSheetName().substring(0,sheet.getSheetName().lastIndexOf(QuantificationResultExporter.ADDUCT_MSN_SHEET));
    if (usedAlexMsnTargets) readConstants.getAlexTargetlistUsed().put(lipidClass,true);
    Vector<LipidParameterSet> msnResults = new Vector<LipidParameterSet>();
    for (LipidParameterSet ms1 : ms1Results){
      msnResults.add(msHash.get(ms1.getNamePlusModHumanReadable()));
    }
    return msnResults;
  }

  /**
   * cleans unidentified position definitions from the position definition hash (-1 corresponds to unidentified)
   * @param posDefs the uncleaned position definition hash
   * @return cleaned position definition hash
   */
  private static Hashtable<String,Hashtable<Integer,Integer>> cleanPositionDefinition (Hashtable<String,Hashtable<Integer,Integer>> posDefs){
    Hashtable<String,Hashtable<Integer,Integer>> cleaned = new Hashtable<String,Hashtable<Integer,Integer>>();
    for (String combi : posDefs.keySet()){
      Hashtable<Integer,Integer> posss = posDefs.get(combi);
      Hashtable<Integer,Integer> possCleaned = new Hashtable<Integer,Integer>();
      Hashtable<Integer,Integer> usedPositions = new Hashtable<Integer,Integer>();
      for (Integer stringPos : posss.keySet()){
        if (posss.get(stringPos)>-1){
          if (usedPositions.containsKey(posss.get(stringPos))){
            //the same position is referenced more than once -> no evidence remove it
            if (possCleaned.containsKey(usedPositions.get(posss.get(stringPos))))possCleaned.remove(usedPositions.get(posss.get(stringPos)));
          }else{
            possCleaned.put(stringPos, posss.get(stringPos));
            usedPositions.put(posss.get(stringPos), stringPos);
          }
        }
      }
      cleaned.put(combi, possCleaned);
    }
    return cleaned;
  }
  
  
  /**
   * the parser adds found chain combinations by replacing the "/" by "_", where the sequence of the fatty acids may
   * be different than in the stored values - this method corrects the sequence to be exactly as in the stored Excel
   * @param combiKey fa combination String as it should be
   * @param validChainCombinations the currently available chain combinations
   * @param relAreas the relative share of each chain combination; first key: chain combination; second key: relative share on MS1 area
   * @param sep the separator for the combination
   * @return the corrected available chain combinations
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private static Vector<String> correctUndefinedChainCombinations(String combiKey, Vector<String> validChainCombinations, Hashtable<String,Double> relAreas, String sep) throws LipidCombinameEncodingException{
    Vector<String> corrected = new Vector<String>();
    for (String combi : validChainCombinations){
      if (StaticUtils.isAPermutedVersion(combiKey,combi,sep)) {
        corrected.add(combiKey);
        if (!combiKey.equalsIgnoreCase(combi)) {
          relAreas.put(combiKey, relAreas.get(combi));
          relAreas.remove(combi);
        }
      }
      else corrected.add(combi);
    }
    return corrected;
  }
  
  /**
   * checks if the position definition hash contains any defined position
   * @param posDefs position definition hash
   * @return true if definitions are present
   */
  private static boolean isDefinitionPresent(Hashtable<String,Hashtable<Integer,Integer>> posDefs){
    for (Hashtable<Integer,Integer> posss : posDefs.values()){
      if (posss.size()>0) return true;
    }
    return false;
  }
  
  
  /**
   * selects the permutation that is used in the validChainCombinations
   * @param combiKey an encoded LDA chain combination
   * @param validChainCombinations the available chain combinations in valid validChainCombinations
   * @return the one in validChainCombinations
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private static String getPermutationThatIsInValidChainCombinations(String combiKey, Vector<String> validChainCombinations) throws LipidCombinameEncodingException {
    Vector<String> chains = StaticUtils.splitChainCombiToEncodedStrings(combiKey,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
    Vector<String> permuts = StaticUtils.getPermutedChainNames(chains, LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
    for (String permut : permuts) {
      if (validChainCombinations.contains(permut))
        return permut;
    }
    return null;
  }
  
}
