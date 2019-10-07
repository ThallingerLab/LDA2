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

package at.tugraz.genome.lda;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Timer;
import java.util.TimerTask;
import java.util.Vector;

import javax.swing.JFrame;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.alex123.RdbOutputWriter;
import at.tugraz.genome.lda.alex123.TargetlistDirParser;
import at.tugraz.genome.lda.alex123.TargetlistParser;
import at.tugraz.genome.lda.alex123.vos.TargetlistEntry;
import at.tugraz.genome.lda.exception.AlexTargetlistParserException;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.MSnAnalyzer;
import at.tugraz.genome.lda.msn.OtherAdductChecker;
import at.tugraz.genome.lda.msn.PostQuantificationProcessor;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.IntensityChainVO;
import at.tugraz.genome.lda.msn.vos.IntensityPositionVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.msn.vos.RtPredictVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.RangeInteger;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleStringVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.Probe3D;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.maspectras.utils.StringUtils;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class QuantificationThread extends Thread
{
  private String chromFile_;
  private String chromFileName_;
  private String quantFile_;
  private String resultFile_;
//  private float mzTolerance_;
  private float minusTime_;
  private float plusTime_;
  private boolean finished_ = false;
  private int amountOfIsotopes_;
  private int isotopesMustMatch_;
  private String errorString_ = null;
  private boolean searchUnknownTime_ = false;
  private float basePeakCutoff_;
  private float rtShift_;
  
  private int totalAmountOfLipids_;
  private int currentLipidCount_;
  private String currentLipid_;
  private int numberOfProcessors_;
  /** the ion mode of the search: true for positive, and false for negative; required only for ALEX123*/
  private boolean ionMode_;
  /** in the case of MSnFirst: in the first round analytes not containing MSn spectra are added to an hash - in the second round normal quantitation*/
  private boolean msnRoundFinished_;
  /** was the task started by the command line interface*/
  private boolean cli_;
  
  private Hashtable<Integer,Boolean> availableThreads_;
  private Hashtable<Integer,LipidomicsAnalyzer> analyzers_;
  private Hashtable<Integer,SingleQuantThread> threads_;
  private Hashtable<String,Hashtable<String,Hashtable<String,Integer>>> quantStatus_;
  /** in the case of MSnFirst: this hash contains all spectra where no MSn spectra are present*/
  private Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> onlyMS1DataPresent_;
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> results_;
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> ms2Removed_;
  /** if a peak split has to be removed because a split partner has a wrong retention time, the unsplit peak version is stored*/
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> unsplittedPeaks_;
  
  private Hashtable<Integer,String> threadToClass_;
  private Hashtable<Integer,String> threadToAnalyte_;
  private Hashtable<Integer,String> threadToMod_;
  
  public final static String OVERVIEW_SHEET_ADDUCT = " - Overview";
  public final static String MSN_SHEET_ADDUCT = " - MSn";
  public final static String CONSTANTS_SHEET = "About";
  private long startCalcTime_;
  
  private final static String OVERVIEW_SHEET_HEADER_SPECIES = "Species";
  
  private Timer timer_;
  
  /** in the case of MSnFirst: contains LM-Models and suggestions for the next range for quantitation*/
  private Hashtable<String,Hashtable<String,RtPredictVO>> latestRtPredictions_; 
  
  private final static int STATUS_WAITING = 0;
  private final static int STATUS_CALCULATING = 1;
  private final static int STATUS_FINISHED = 2;

  private final static int MSN_ROW_FRAGMENT_OH = 1;
  private final static int MSN_ROW_FRAGMENT_CHAIN_TYPE = 2;
  
  private final static int MSN_ROW_FRAGMENT_MSLEVEL = 2;
  private final static int MSN_ROW_FRAGMENT_CHARGE = 3;
  private final static int MSN_ROW_FRAGMENT_MZ = 4;
  private final static int MSN_ROW_FRAGMENT_MZ_TOLERANCE = 5;
  private final static int MSN_ROW_FRAGMENT_AREA = 6;
  private final static int MSN_ROW_FRAGMENT_PEAK = 7;
  private final static int MSN_ROW_FRAGMENT_TIME_LOWER = 8;
  private final static int MSN_ROW_FRAGMENT_TIME_UPPER = 9;
  private final static int MSN_ROW_FRAGMENT_MZ_LOWER = 10;
  private final static int MSN_ROW_FRAGMENT_MZ_UPPER = 11;
  private final static int MSN_ROW_FRAGMENT_ELLIPSE_TIME = 12;
  private final static int MSN_ROW_FRAGMENT_ELLIPSE_MZ = 13;
  private final static int MSN_ROW_FRAGMENT_ELLIPSE_TIME_RANGE = 14;
  private final static int MSN_ROW_FRAGMENT_ELLIPSE_MZ_RANGE = 15;

  private final static int MSN_ROW_INTENSITY_RULE = 0;
  private final static int MSN_ROW_INTENSITY_ORIGINAL = 1;
  private final static int MSN_ROW_INTENSITY_VALUES = 2;
  private final static int MSN_ROW_INTENSITY_MISSED = 3;
  
    
  public QuantificationThread(String chromFile,String quantFile,String resultFile,//float mzTolerance, 
      float minusTime, float plusTime, int amountOfIsotopes, int isotopesMustMatch, boolean searchUnknownTime,
      float basePeakCutoff, float rtShift, int numberOfProcessors, boolean ionMode, boolean cli){
    super();
    this.chromFile_ = chromFile;
    this.chromFileName_ = StringUtils.getJustFileName(chromFile);
    this.quantFile_ = quantFile;
    this.resultFile_ = resultFile;
    this.minusTime_ = minusTime;
    this.plusTime_ = plusTime;
    this.amountOfIsotopes_ = amountOfIsotopes;
    this.isotopesMustMatch_ = isotopesMustMatch;
    errorString_ = null;
    this.searchUnknownTime_ = searchUnknownTime;
    this.basePeakCutoff_ = basePeakCutoff;
    rtShift_ = rtShift;
    finished_ = false;
    numberOfProcessors_ = numberOfProcessors;
    this.ionMode_ = ionMode;
    this.cli_ = cli;
  }
  
  public void run(){
    try{
      startAutomatedLipidomicsQuantification(chromFile_, quantFile_,resultFile_, //mzTolerance_,
          minusTime_,plusTime_,amountOfIsotopes_,isotopesMustMatch_, searchUnknownTime_,basePeakCutoff_,rtShift_,
          numberOfProcessors_,ionMode_);
    } catch (Exception ex){
      ex.printStackTrace();
      errorString_ = ex.toString();
      this.finished_ = true;
      
    }
  }
  
  public boolean finished(){
    if (finished_ && this.timer_!=null)
      this.timer_.cancel();
    return this.finished_;
  }
  
  public String getErrorString(){
    return this.errorString_;
  }
  
  @SuppressWarnings("unchecked")
  private void startAutomatedLipidomicsQuantification(String chromFile,String quantFile, String resultFile, //float mzTolerance,
      float minusTime, float plusTime, int amountOfIsotopes, int isotopesMustMatch, boolean searchUnknownTime, float basePeakCutoff,
      float rtShift, int numberOfProcessors, boolean ionMode) throws Exception{
    startCalcTime_ = System.currentTimeMillis();
    totalAmountOfLipids_ = -1;
    currentLipidCount_ = -1;
    currentLipid_ = "";
    msnRoundFinished_ = false;
    latestRtPredictions_ = new Hashtable<String,Hashtable<String,RtPredictVO>>();
    String pureFile = chromFile.substring(0,chromFile.lastIndexOf("."));
    String errorMessage = StaticUtils.existChromNecessaryFiles(pureFile);
    if (errorMessage!=null && errorMessage.length()>0) throw new Exception(errorMessage);
    String[] chromPaths = StringUtils.getChromFilePaths(pureFile+".chrom");
    float[] maxRetTimes = initThreadMonitors(chromPaths, numberOfProcessors, basePeakCutoff);
    float highestRetTime = maxRetTimes[1];
    float lowestRetTime = maxRetTimes[0];
    @SuppressWarnings("rawtypes")
    Vector quantContent = null;
    File quant = new File(quantFile);
    if (quant.isFile() && (quant.getName().endsWith(".xls") || quant.getName().endsWith(".xlsx"))){
      quantContent = QuantificationThread.parseQuantExcelFile(quantFile, minusTime, plusTime, amountOfIsotopes, isotopesMustMatch, searchUnknownTime, basePeakCutoff, rtShift, lowestRetTime, highestRetTime);
    } else if ((quant.isFile() && quant.getName().endsWith(".txt")) || quant.isDirectory()){
      quantContent = QuantificationThread.parseAlex123TargetList(quantFile, minusTime, plusTime, amountOfIsotopes,
          isotopesMustMatch, searchUnknownTime, basePeakCutoff, rtShift, lowestRetTime, highestRetTime, ionMode);
    }
    if (quantContent!=null){
      LinkedHashMap<String,Integer> classSequence = (LinkedHashMap<String,Integer>)quantContent.get(0);
      Hashtable<String,Vector<String>> analyteSequence = (Hashtable<String,Vector<String>>)quantContent.get(1);
    
      totalAmountOfLipids_ = 0;
      for (String className : classSequence.keySet()) totalAmountOfLipids_ += analyteSequence.get(className).size();
      timer_ = new java.util.Timer();
      timer_.schedule(new ThreadSupervisor(quantContent,basePeakCutoff,resultFile), 10, 100);
    } else {
      this.errorString_ = "The quantification file/folder does not contain any usable files";
      this.finished_ = true;
    }
  }
  
  public int getTotalAmountOfLipids()
  {
    return totalAmountOfLipids_;
  }

  public int getCurrentLipidCount()
  {
    return currentLipidCount_;
  }

  public String getCurrentLipid()
  {
    return currentLipid_;
  }

  private static boolean isWithinRemovedBoundaries(CgProbe probe, Vector<CgProbe> removedProbes){
    boolean isWithinBoundaries = false;
    
    for (CgProbe removedProbe : removedProbes){
      if ((probe.LowerValley>removedProbe.LowerValley-20 && probe.UpperValley<removedProbe.UpperValley+20)||
          (removedProbe.LowerValley<probe.Peak && removedProbe.LowerValley<removedProbe.UpperValley))
        isWithinBoundaries = true;
    }
    return isWithinBoundaries;
  }
  
  @SuppressWarnings("resource")
  public static void writeResultsToExcel(String resultFile,QuantificationResult quantRes) throws Exception{
    BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(resultFile));
    Workbook resultWorkbook = new XSSFWorkbook();
    //this line is for backward compatibility to the old Excel 2003 saved files
    if (resultFile.endsWith(".xls")) resultWorkbook = new HSSFWorkbook();
    CellStyle headerStyle = getHeaderStyle(resultWorkbook);
    boolean hasRtInfo = hasRtInfo(quantRes.getIdentifications());
    int rtPlus = 0;
    if (hasRtInfo) rtPlus=1;
    LipidomicsConstants constants = quantRes.getConstants();
    if (constants!=null){
      Sheet constantsSheet = resultWorkbook.createSheet(CONSTANTS_SHEET );
      constants.writeSettingsToExcel(constantsSheet,headerStyle,quantRes.getFaHydroxyEncoding(),quantRes.getLcbHydroxyEncoding());
    }
    for (String sheetName: quantRes.getIdentifications().keySet()){
      Vector<LipidParameterSet> params = quantRes.getIdentifications().get(sheetName);
      boolean hasMSnInformation = false;
//      boolean hasChainInfo = false;
      Sheet resultSheet = resultWorkbook.createSheet(sheetName);
      Sheet resultMSnSheet  = null;
      boolean hasOhInfo = false;
      for (LipidParameterSet param : params){
        if (param instanceof LipidomicsMSnSet) {
          hasMSnInformation = true;
//          if (((LipidomicsMSnSet)param).getStatus()>=LipidomicsMSnSet.FRAGMENTS_DETECTED)
//            hasChainInfo = true;
        }
        if (param.getOhNumber()>LipidomicsConstants.EXCEL_NO_OH_INFO)
          hasOhInfo = true;
        if (hasMSnInformation && hasOhInfo)
          break;
      }
      int ohPlus = 0;
      if (hasOhInfo) ohPlus = 1;
      Hashtable<Integer,String> msnLongestHeaders = null;
      Hashtable<Integer,Integer> msnLongestEntries = null;
      if (hasMSnInformation){
        resultMSnSheet = resultWorkbook.createSheet(sheetName+MSN_SHEET_ADDUCT);
        msnLongestHeaders = new Hashtable<Integer,String>();
        msnLongestEntries = new Hashtable<Integer,Integer>();
      }
      Sheet resultSheetOV = null;
//      HSSFRow headerRowOV = null;
//      HSSFRow valuesRowOV = null;
      Cell label = null;
      
      int amountOfIsotopes = 0;
      int massAdd = 0;
      if (Settings.isOverviewInExcelDesired()/* && params.size()<256*/){
        resultSheetOV = resultWorkbook.createSheet(sheetName+OVERVIEW_SHEET_ADDUCT);
        Row headerRowOV = resultSheetOV.createRow(0);
        int beginIndex = 0;
        int lastIndexSlash = resultFile.lastIndexOf("/");
        int lastIndexBackSlash = resultFile.lastIndexOf("\\");
        if (lastIndexSlash>lastIndexBackSlash)
          beginIndex = lastIndexSlash;
        else
          beginIndex = lastIndexBackSlash;
        beginIndex++;
        label = headerRowOV.createCell(0,Cell.CELL_TYPE_STRING);
        label.setCellValue(OVERVIEW_SHEET_HEADER_SPECIES);
        label.setCellStyle(headerStyle);
        label = headerRowOV.createCell(1,Cell.CELL_TYPE_STRING);
        label.setCellValue(resultFile.substring(beginIndex,resultFile.lastIndexOf(".")));
        label.setCellStyle(headerStyle);
        if (Settings.isMassInOverviewExcelDesired()){
          label = headerRowOV.createCell(2,Cell.CELL_TYPE_STRING);
          label.setCellValue("m/z");
          label.setCellStyle(headerStyle);
          massAdd++;
        }
        if (Settings.isIsotopeInOverviewExcelDesired()){
          for (LipidParameterSet param : params){
            int maxIso = param.getIsotopicProbes().size();
            if (maxIso>amountOfIsotopes) amountOfIsotopes = maxIso;
          }
          for (int i=0;i!=amountOfIsotopes;i++){
            label = headerRowOV.createCell(2+massAdd+i,Cell.CELL_TYPE_STRING);
            label.setCellValue("Isotope "+i);
            label.setCellStyle(headerStyle);            
          }
        }
//        headerRowOV = resultSheetOV.createRow(0);
//        valuesRowOV = resultSheetOV.createRow(1);
      }
      Row row = resultSheet.createRow(0);
      int msLevel = 1;
      if (quantRes.getMsLevels()!=null&&quantRes.getMsLevels().containsKey(sheetName)) msLevel = quantRes.getMsLevels().get(sheetName);
      
      // Create the label, specifying content and format 
      label = row.createCell(1,Cell.CELL_TYPE_STRING);
      label.setCellValue("Name");
      label.setCellStyle(headerStyle);
      label = row.createCell(2,Cell.CELL_TYPE_STRING);
      label.setCellValue("Dbs");
      label.setCellStyle(headerStyle);
      if (hasOhInfo) {
        label = row.createCell(3,Cell.CELL_TYPE_STRING);
        label.setCellValue(LipidomicsConstants.EXCEL_MS_OH);
        label.setCellStyle(headerStyle);
      }
      label = row.createCell(3+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("Modification");
      label.setCellStyle(headerStyle);
      label = row.createCell(4+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("Formula");
      label.setCellStyle(headerStyle);
      label = row.createCell(5+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("Mod-Formula");
      label.setCellStyle(headerStyle);
      if (hasRtInfo){
        label = row.createCell(6+ohPlus,Cell.CELL_TYPE_STRING);
        label.setCellValue("RT");
        label.setCellStyle(headerStyle);
      }
      label = row.createCell(6+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("Isotope");
      label.setCellStyle(headerStyle);
      label = row.createCell(7+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("Area");
      label.setCellStyle(headerStyle);
      label = row.createCell(8+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("AreaError");
      label.setCellStyle(headerStyle);
      label = row.createCell(9+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("Background");
      label.setCellStyle(headerStyle);
      label = row.createCell(10+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("Charge");
      label.setCellStyle(headerStyle);
      label = row.createCell(11+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("Mz");
      label.setCellStyle(headerStyle);
      label = row.createCell(12+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("MzTolerance");
      label.setCellStyle(headerStyle);

      label = row.createCell(13+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("Peak");
      label.setCellStyle(headerStyle);
      label = row.createCell(14+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("LowerValley");
      label.setCellStyle(headerStyle);
      label = row.createCell(15+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("UpperValley");
      label.setCellStyle(headerStyle);
      label = row.createCell(16+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("LowMz");
      label.setCellStyle(headerStyle);
      label = row.createCell(17+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("UpMz");
      label.setCellStyle(headerStyle);
      label = row.createCell(18+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("EllCentTime");
      label.setCellStyle(headerStyle);
      label = row.createCell(19+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("EllCentMz");
      label.setCellStyle(headerStyle);
      label = row.createCell(20+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("EllStretchTime");
      label.setCellStyle(headerStyle);
      label = row.createCell(21+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("EllStretchMz");
      label.setCellStyle(headerStyle);
      label = row.createCell(22+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("LowerRtHardLimit");
      label.setCellStyle(headerStyle);
      label = row.createCell(23+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("UpperRtHardLimit");
      label.setCellStyle(headerStyle);
      label = row.createCell(24+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("PercentalSplit");
      label.setCellStyle(headerStyle);
      label = row.createCell(25+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue("level="+String.valueOf(msLevel));
      label.setCellStyle(headerStyle);
      label = row.createCell(26+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue(LDAResultReader.COLUMN_APEX_INTENSITY);
      label.setCellStyle(headerStyle);
      label = row.createCell(27+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue(LDAResultReader.COLUMN_LOWER_VALLEY10PC);
      label.setCellStyle(headerStyle);
      label = row.createCell(28+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue(LDAResultReader.COLUMN_LOWER_VALLEY50PC);
      label.setCellStyle(headerStyle);
      label = row.createCell(29+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue(LDAResultReader.COLUMN_UPPER_VALLEY50PC);
      label.setCellStyle(headerStyle);
      label = row.createCell(30+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue(LDAResultReader.COLUMN_UPPER_VALLEY10PC);
      label.setCellStyle(headerStyle);
      label = row.createCell(31+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue(LDAResultReader.COLUMN_LOWER_MZ10PC);
      label.setCellStyle(headerStyle);
      label = row.createCell(32+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue(LDAResultReader.COLUMN_LOWER_MZ50PC);
      label.setCellStyle(headerStyle);
      label = row.createCell(33+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue(LDAResultReader.COLUMN_UPPER_MZ50PC);
      label.setCellStyle(headerStyle);
      label = row.createCell(34+rtPlus+ohPlus,Cell.CELL_TYPE_STRING);
      label.setCellValue(LDAResultReader.COLUMN_UPPER_MZ10PC);
      label.setCellStyle(headerStyle);


//      if (Settings.isOverviewInExcelDesired() /*&& params.size()<256*/){
//        int beginIndex = 0;
//        int lastIndexSlash = resultFile.lastIndexOf("/");
//        int lastIndexBackSlash = resultFile.lastIndexOf("\\");
//        if (lastIndexSlash>lastIndexBackSlash)
//          beginIndex = lastIndexSlash;
//        else
//          beginIndex = lastIndexBackSlash;
//        beginIndex++;
//        label = valuesRowOV.createCell(0,HSSFCell.CELL_TYPE_STRING);
//        label.setCellValue(resultFile.substring(beginIndex,resultFile.lastIndexOf(".")));
//        label.setCellStyle(headerStyle);
//      }
      int resultCount = 0;
      int resultRowCount = 0;
      int msnRowCount = 0;
      if (hasMSnInformation && constants.getAlexTargetlistUsed().containsKey(sheetName) && constants.getAlexTargetlistUsed().get(sheetName)){
        Row msnRow = resultMSnSheet.createRow(msnRowCount);
        msnRowCount++;
        Cell cell = msnRow.createCell(0);
        cell.setCellStyle(headerStyle);
        cell.setCellValue(LDAResultReader.ALEX123_MSN_TARGETS_USED+"=true");
      }
      for (LipidParameterSet param : params){
        if (param instanceof LipidomicsMSnSet){
          try {
            msnRowCount = writeMSnEvidence(msnRowCount,resultMSnSheet,(LipidomicsMSnSet)param, headerStyle, msnLongestHeaders, msnLongestEntries,
                quantRes.getFaHydroxyEncoding(),quantRes.getLcbHydroxyEncoding());
          }catch(Exception ex){
            ex.printStackTrace();
            throw ex;
          }
        }
        resultCount++;
        double totalArea = 0;
        Cell number;
        resultRowCount++;
        row = resultSheet.createRow(resultRowCount);
        number = row.createCell(0,Cell.CELL_TYPE_NUMERIC);
        number.setCellValue(new Double(resultCount));
        label = row.createCell(1,Cell.CELL_TYPE_STRING);
        label.setCellValue(param.Peptide);
        label = row.createCell(2,Cell.CELL_TYPE_STRING);
        label.setCellValue(String.valueOf(param.getDoubleBonds()));
        if (hasOhInfo) {
          label = row.createCell(3,Cell.CELL_TYPE_STRING);
          label.setCellValue(String.valueOf(param.getOhNumber()));
        }
        label = row.createCell(3+ohPlus,Cell.CELL_TYPE_STRING);
        label.setCellValue(param.getModificationName());
        label = row.createCell(4+ohPlus,Cell.CELL_TYPE_STRING);
        label.setCellValue(param.getAnalyteFormula());
        label = row.createCell(5+ohPlus,Cell.CELL_TYPE_STRING);
        label.setCellValue(param.getModificationFormula());
        if (hasRtInfo){
          label = row.createCell(6+ohPlus,Cell.CELL_TYPE_STRING);
          label.setCellValue(param.getRt());
        }  
        
        Cell totalAreaCell = row.createCell(7+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
        number = row.createCell(10+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
        number.setCellValue(param.getCharge());        
        number = row.createCell(11+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
        number.setCellValue(param.Mz[0]);
        number = row.createCell(12+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
        number.setCellValue(param.LowerMzBand);
        if (param.getLowerRtHardLimit()>=0){
          number = row.createCell(22+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
          number.setCellValue(param.getLowerRtHardLimit());
        }
        if (param.getUpperRtHardLimit()>=0){
          number = row.createCell(23+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
          number.setCellValue(param.getUpperRtHardLimit());
        }
        if (param.getPercentalSplit()>=0){
          number = row.createCell(24+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
          number.setCellValue(param.getPercentalSplit());
        }
        Row valuesRowOV = null;
        if (Settings.isOverviewInExcelDesired() /*&& params.size()<256*/){
//          label = headerRowOV.createCell(resultCount,HSSFCell.CELL_TYPE_STRING);
          valuesRowOV = resultSheetOV.createRow(resultCount);
          label = valuesRowOV.createCell(0,Cell.CELL_TYPE_STRING);
          label.setCellValue(param.getNamePlusModHumanReadable());
          label.setCellStyle(headerStyle);
        }
        Vector<Vector<CgProbe>> isotopicProbes = param.getIsotopicProbes();
        for (int j=0;j!=isotopicProbes.size();j++){
          resultRowCount++;
          row = resultSheet.createRow(resultRowCount);
          number = row.createCell(6+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
          int chargeState = j;
          if (param.getMinIsotope()<0) chargeState*=-1;
          number.setCellValue(chargeState);          
          Vector<CgProbe> probes = isotopicProbes.get(j);
          float totalIsoArea = 0f;
          for (int i=0; i!=probes.size(); i++){
            totalIsoArea+=probes.get(i).Area;
          }
          number = row.createCell(7+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
          number.setCellValue(totalIsoArea);
          totalArea+=totalIsoArea;
          for (int i=0; i!=probes.size(); i++){
            CgProbe probe = probes.get(i);
            resultRowCount++;
            row = resultSheet.createRow(resultRowCount);
            if (probe!=null&&probe.AreaStatus==CgAreaStatus.OK){
              number = row.createCell(6+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
              number.setCellValue(chargeState);
              number = row.createCell(7+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
              number.setCellValue(new Double(probe.Area));
              number = row.createCell(8+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
              number.setCellValue(new Double(probe.AreaError));
              number = row.createCell(9+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
              number.setCellValue(new Double(probe.Background));
              number = row.createCell(10+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
              number.setCellValue(probe.Charge);
              number = row.createCell(13+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
              number.setCellValue(new Double(probe.Peak));
              number = row.createCell(14+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
              number.setCellValue(new Double(probe.LowerValley));
              number = row.createCell(15+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
              number.setCellValue(new Double(probe.UpperValley));
              if (probe.getLowerValley10()!=null){
                number = row.createCell(26+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe.getApexIntensity()));
                number = row.createCell(27+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe.getLowerValley10()));
                number = row.createCell(28+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe.getLowerValley50()));
                number = row.createCell(29+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe.getUpperValley50()));
                number = row.createCell(30+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe.getUpperValley10()));
              }
              if (probe instanceof Probe3D){
                Probe3D probe3D = (Probe3D)probe;
                number = row.createCell(16+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe3D.LowerMzBand));
                number = row.createCell(17+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe3D.UpperMzBand));
                number = row.createCell(18+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe3D.getEllipseTimePosition()));
                number = row.createCell(19+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe3D.getEllipseMzPosition()));
                number = row.createCell(20+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe3D.getEllipseTimeStretch()));
                number = row.createCell(21+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                number.setCellValue(new Double(probe3D.getEllipseMzStretch()));
                if (((Probe3D) probe).getLowMz10()>-1){
                  number = row.createCell(31+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                  number.setCellValue(new Double(probe3D.getLowMz10()));
                  number = row.createCell(32+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                  number.setCellValue(new Double(probe3D.getLowMz50()));
                  number = row.createCell(33+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                  number.setCellValue(new Double(probe3D.getUpMz50()));
                  number = row.createCell(34+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
                  number.setCellValue(new Double(((Probe3D) probe).getUpMz10()));
                }
              }
            } else{
              number = row.createCell(7+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
              number.setCellValue(new Double(0));
            }
            number = row.createCell(11+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
            number.setCellValue(probe.Mz);
            number = row.createCell(12+rtPlus+ohPlus,Cell.CELL_TYPE_NUMERIC);
            number.setCellValue(param.LowerMzBand);
          }
        }
        totalAreaCell.setCellValue(totalArea);
        if (Settings.isOverviewInExcelDesired()  /*&& params.size()<256*/){
//          label = valuesRowOV.createCell(resultCount,HSSFCell.CELL_TYPE_NUMERIC);
          label = valuesRowOV.createCell(1,Cell.CELL_TYPE_NUMERIC);
          label.setCellValue(totalArea);
          if (Settings.isMassInOverviewExcelDesired()){
            label = valuesRowOV.createCell(2,Cell.CELL_TYPE_STRING);
            label.setCellValue(param.Mz[0]);
          }
          if (Settings.isIsotopeInOverviewExcelDesired()){
            for (int i=0;i!=amountOfIsotopes;i++){
              if (param.getIsotopicProbes().size()>i){
                Vector<CgProbe> probes = isotopicProbes.get(i);
                float totalIsoArea = 0f;
                for (int j=0; j!=probes.size(); j++){
                  totalIsoArea+=probes.get(j).Area;
                }
                label = valuesRowOV.createCell(2+massAdd+i,Cell.CELL_TYPE_STRING);
                label.setCellValue(totalIsoArea);
              }
            }  
          }
        }
      }
      if (hasMSnInformation){
        for (Integer column : msnLongestHeaders.keySet()) {
          setColumnWidth(resultMSnSheet, column, msnLongestHeaders.get(column), (msnLongestEntries.containsKey(column) ? msnLongestEntries.get(column) : 0));
        }
      }
    }
    resultWorkbook.write(out);
    out.close();
  }

  /**
   * writes the MSn evidence for one identified MSn species
   * @param msnRowCount the current Excel row where subsequent information can be written
   * @param sheet Excel sheet for MSn evidence for this class
   * @param param the lipid identification containing MSn evidence
   * @param headerStyle style for the header row
   * @param longestHeaders the names of the longest header entries for each column - for setting the column width
   * @param longestEntries the length of the longest entries for each column - for setting the column width
   * @param faHydroxyEncoding the character encoding of the number of hydroxylation sites for the
   * @param lcbHydroxyEncoding the character encoding of the number of hydroxylation sites for the LCB
   * @return the next row where subsequent information can be written
   * @throws RulesException
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded 
   * @throws ChemicalFormulaException thrown when something is wrong with the chemical formula
   */
  @SuppressWarnings("unchecked")
  private static int writeMSnEvidence(int msnRowCount, Sheet sheet, LipidomicsMSnSet param, CellStyle headerStyle,
      Hashtable<Integer,String> longestHeaders, Hashtable<Integer,Integer> longestEntries, HydroxyEncoding faHydroxyEncoding,
      HydroxyEncoding lcbHydroxyEncoding) throws RulesException, LipidCombinameEncodingException, ChemicalFormulaException {
    int count = msnRowCount;
    Row row = sheet.createRow(count);
    count++;
    Row areaRow = sheet.createRow(count);
    count++;
    createHeaderCell(row, 0, headerStyle, param.getNamePlusModHumanReadable(), longestHeaders);
    createValueCell(areaRow, Cell.CELL_TYPE_NUMERIC, 0, String.valueOf(param.Area), null, (double)param.Area, longestEntries);
    
    // writing of fragments and rules concerning the head group
    Hashtable<String,CgProbe> headGroupFragments = param.getHeadGroupFragments();
    if (headGroupFragments.size()>0){
     Row headGroupRow =  sheet.createRow(count);
     count++;
     createHeaderCell(headGroupRow, 0, headerStyle, LipidomicsConstants.EXCEL_MSN_SECTION_HEAD_FRAGMENTS, longestHeaders);
     headGroupRow = sheet.createRow(count);
     count++;
     writeMSnFragmentHeader(headGroupRow,headerStyle,LipidomicsConstants.CHAIN_TYPE_NO_CHAIN,longestHeaders);
     for (String name : headGroupFragments.keySet()){
       headGroupRow = sheet.createRow(count);
       count++;
       writeMSnFragment(headGroupRow,name,param.getMSnMzTolerance(),headGroupFragments.get(name),LipidomicsConstants.CHAIN_TYPE_NO_CHAIN,-1,
           longestEntries);
     }
     Hashtable<String,IntensityRuleVO> headIntRules = param.getHeadIntensityRules();
     if (headIntRules.size()>0){
       headGroupRow =  sheet.createRow(count);
       count++;
       createHeaderCell(headGroupRow, 0, headerStyle, LipidomicsConstants.EXCEL_MSN_SECTION_HEAD_INTENSITIES, longestHeaders);
       headGroupRow = sheet.createRow(count);
       count++;
       writeMSnIntensityHeader(headGroupRow,headerStyle,  longestHeaders);
       Hashtable<String,String> uniqueRules = new Hashtable<String,String>();
       for (IntensityRuleVO ruleVO : headIntRules.values()){
         if (writeMSnIntensity(sheet,count,ruleVO,param,uniqueRules,longestEntries,faHydroxyEncoding,lcbHydroxyEncoding)){
           count++;
         }
       }
     }
    }
    
    // writing of fragments and rules concerning the fragments
    Hashtable<String,Hashtable<String,CgProbe>> chainFragments = param.getChainFragments();
    if (chainFragments.size()>0){
      Row chainRow =  sheet.createRow(count);
      count++;
      createHeaderCell(chainRow, 0, headerStyle, LipidomicsConstants.EXCEL_MSN_SECTION_CHAIN_FRAGMENTS, longestHeaders);
      chainRow = sheet.createRow(count);
      count++;
      writeMSnFragmentHeader(chainRow,headerStyle,LipidomicsConstants.CHAIN_TYPE_FA_ACYL,longestHeaders);
      for (String faName : chainFragments.keySet()){
        Hashtable<String,CgProbe> fragments = chainFragments.get(faName);
        for (String name : fragments.keySet()){
          chainRow = sheet.createRow(count);
          count++;
          FattyAcidVO fa = StaticUtils.decodeLipidNameForCreatingCombis(faName);
          String displayName = StaticUtils.getChainFragmentDisplayName(name,fa.getCarbonDbsId());
          writeMSnFragment(chainRow,displayName,param.getMSnMzTolerance(),fragments.get(name),fa.getChainType(),fa.getOhNumber(),longestEntries);
        }
      }
      Hashtable<String,Hashtable<String,IntensityChainVO>> chainRules = param.getChainIntensityRules();
      if (chainRules.size()>0){
        chainRow =  sheet.createRow(count);
        count++;
        createHeaderCell(chainRow, 0, headerStyle, LipidomicsConstants.EXCEL_MSN_SECTION_CHAIN_INTENSITIES, longestHeaders);
        chainRow = sheet.createRow(count);
        count++;
        writeMSnIntensityHeader(chainRow,headerStyle,longestHeaders);
        Hashtable<String,String> uniqueRules = new Hashtable<String,String>();
        for (String faName : chainRules.keySet()){
          Hashtable<String,IntensityChainVO> rules = chainRules.get(faName);
          for (IntensityRuleVO rule : rules.values()){
            if (writeMSnIntensity(sheet,count,rule,param,uniqueRules,longestEntries,faHydroxyEncoding,lcbHydroxyEncoding)){
              count++;
            }
          }
        }
      }
    }
    
    // writing of position rules
    Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> posRules = param.getPositionEvidence();
    if (posRules.size()>0){
      for (String combiName : posRules.keySet()){
        Row positionRow =  sheet.createRow(count);
        count++;
        createHeaderCell(positionRow, 0, headerStyle, LipidomicsConstants.EXCEL_MSN_SECTION_POSITION_INTENSITIES+" ("+param.getPositionInsensitiveHumanReadableCombiName(combiName)+")", longestHeaders);
        positionRow = sheet.createRow(count);
        count++;
        writeMSnIntensityHeader(positionRow,headerStyle,longestHeaders);
        Hashtable<String,String> uniqueRules = new Hashtable<String,String>();
        Vector<IntensityRuleVO> rules = param.getFAsInSequenceAsInRule(combiName); 
        for (IntensityRuleVO rule : rules){
          if (writeMSnIntensity(sheet,count,rule,param,uniqueRules,longestEntries,faHydroxyEncoding,lcbHydroxyEncoding)){
            count++;
          }  
        }
      }
    }
    Hashtable<String,String> identifiedLipids = new Hashtable<String,String>();
    String identificationString;
    int cellCount = 1;
    Vector<Object> detected = null;
    try {detected = param.getMSnIdentificationNames();
    }catch (LipidCombinameEncodingException lcx) {
      detected = new Vector<Object>();
      lcx.printStackTrace();
    }
    for (Object nameObject : detected){
      identificationString = "";
      double area = 0d;
      if (nameObject instanceof Vector){
        area = param.getRelativeIntensity(((Vector<String>)nameObject).get(0))*((double)param.Area);
        for (String name : (Vector<String>)nameObject){
          identificationString+=name+LipidomicsConstants.CHAIN_COMBI_SEPARATOR_AMBIG_POS;
        }
        identificationString = identificationString.substring(0,identificationString.length()-1);
      }else{
        String name = (String) nameObject;
        identificationString = name;
        if (param.getStatus()==LipidomicsMSnSet.HEAD_GROUP_DETECTED)area = param.Area;
        else area = param.getRelativeIntensity(name)*((double)param.Area);
      }
      identifiedLipids.put(identificationString, identificationString);
      createHeaderCell(row, cellCount, headerStyle, identificationString, longestHeaders);
      createValueCell(areaRow, Cell.CELL_TYPE_NUMERIC, cellCount, String.valueOf(area), null, area, longestEntries);
      cellCount++;
    }
    
    //writing the retention times of the uses spectra
    List<Integer> msLevels = new ArrayList<Integer>(param.getMsnRetentionTimes().keySet());
    Collections.sort(msLevels);
    for (int msLevel : msLevels){
      LinkedHashMap<Integer,Float> rts = param.getMsnRetentionTimes().get(msLevel);
      String rtString = "";
      for (Integer scanNr : rts.keySet()) rtString += scanNr+"="+rts.get(scanNr)+";";
      if (rtString.length()>0) rtString = rtString.substring(0,rtString.length()-1);
      createHeaderCell(row, cellCount, headerStyle, "MS"+msLevel+" scan RTs", longestHeaders);
      createValueCell(areaRow, Cell.CELL_TYPE_STRING, cellCount, rtString, null, null, longestEntries);
      cellCount++;
    }
    count++;
    return count;
  }
    
  /**
   * writes the header row for subsequent fragment evidence
   * @param row Excel row that shall be used for writing
   * @param headerStyle style for the header row
   * @param chainType the chain type - allowed: LipidomicsConstants.CHAIN_TYPE_NO_CHAIN|CHAIN_TYPE_FA_ACYL|CHAIN_TYPE_FA_ALKYL|CHAIN_TYPE_FA_ALKENYL|CHAIN_TYPE_LCB
   * @param longestHeaders the names of the longest header entries for each column - for setting the column width
   */
  private static void writeMSnFragmentHeader(Row row, CellStyle headerStyle, short chainType, Hashtable<Integer,String> longestHeaders){
    createHeaderCell(row, LDAResultReader.MSN_ROW_FRAGMENT_NAME, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_NAME, longestHeaders);
    int add = 0;
    if (chainType!=LipidomicsConstants.CHAIN_TYPE_NO_CHAIN) {
      createHeaderCell(row, MSN_ROW_FRAGMENT_OH, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_OH, longestHeaders);
      createHeaderCell(row, MSN_ROW_FRAGMENT_CHAIN_TYPE, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_CHAIN_TYPE, longestHeaders);      
      add = 2;
    }
    createHeaderCell(row, LDAResultReader.MSN_ROW_FRAGMENT_FORMULA+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_FORMULA, longestHeaders);      
    createHeaderCell(row, MSN_ROW_FRAGMENT_MSLEVEL+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_MSLEVEL, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_CHARGE+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_CHARGE, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_MZ+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_MZ_TOLERANCE+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_TOLERANCE, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_AREA+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_AREA, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_PEAK+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_PEAK, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_TIME_LOWER+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_TIME_LOWER, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_TIME_UPPER+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_TIME_UPPER, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_MZ_LOWER+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_LOWER, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_MZ_UPPER+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_UPPER, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_ELLIPSE_TIME+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_TIME, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_ELLIPSE_MZ+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_MZ, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_ELLIPSE_TIME_RANGE+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_TIME_RANGE, longestHeaders);
    createHeaderCell(row, MSN_ROW_FRAGMENT_ELLIPSE_MZ_RANGE+add, headerStyle, LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_MZ_RANGE, longestHeaders);
  }
  
  
  /**
   * creates a header cell
   * @param row the row where it has to be created
   * @param column the column number
   * @param headerStyle the header style
   * @param value the value to be written in the header cell
   * @param msnLongestHeaders the names of the longest header entries for each column - for setting the column width
   * @return created cell
   */
  private static Cell createHeaderCell(Row row, int column, CellStyle headerStyle, String value, Hashtable<Integer,String> longestHeaders) {
    Cell cell = row.createCell(column,Cell.CELL_TYPE_STRING);
    cell.setCellStyle(headerStyle);
    cell.setCellValue(value);
    if (!longestHeaders.containsKey(column) || value.length()>longestHeaders.get(column).length())
      longestHeaders.put(column, value);
    return cell;
  }

  /**
   * writes a line containing found fragments
   * @param row Excel row that shall be used for writing
   * @param name display name of the fragment
   * @param mzTolerance m/z tolerance for the identification
   * @param probe VO containing information about the identified fragment
   * @param chainType the chain type - allowed: LipidomicsConstants.CHAIN_TYPE_NO_CHAIN|CHAIN_TYPE_FA_ACYL|CHAIN_TYPE_FA_ALKYL|CHAIN_TYPE_FA_ALKENYL|CHAIN_TYPE_LCB
   * @param oh the number of hydroxylations on a chain
   * @throws ChemicalFormulaException thrown when something is wrong with the chemical formula
   */
  private static void writeMSnFragment(Row row, String name, float mzTolerance, CgProbe probe,short chainType, int oh, Hashtable<Integer,Integer> longestEntries) throws ChemicalFormulaException{
    createValueCell(row, Cell.CELL_TYPE_STRING, LDAResultReader.MSN_ROW_FRAGMENT_NAME, name, null, null, longestEntries);
    int add = 0;
    if (chainType!=LipidomicsConstants.CHAIN_TYPE_NO_CHAIN) {
      createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_OH, String.valueOf(oh), oh, null, longestEntries);
      createValueCell(row, Cell.CELL_TYPE_STRING, MSN_ROW_FRAGMENT_CHAIN_TYPE, StaticUtils.getHumanReadableChainType(chainType), null, null, longestEntries);
      add = 2;
    }    
    //TODO: this is only here because of a damaged Alex123 file - delete in future version!
    String formula = "";
    if (probe.getFormula()!=null) formula = StaticUtils.getFormulaInHillNotation(StaticUtils.categorizeFormula(probe.getFormula().trim()),true);
    createValueCell(row, Cell.CELL_TYPE_STRING, LDAResultReader.MSN_ROW_FRAGMENT_FORMULA+add, formula, null, null, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_MSLEVEL+add, String.valueOf(probe.getMsLevel()), probe.getMsLevel(), null, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_CHARGE+add, String.valueOf(probe.Charge), probe.Charge, null, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_MZ+add, String.valueOf(probe.Mz), null, (double)probe.Mz, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_MZ_TOLERANCE+add, String.valueOf(mzTolerance), null, (double)mzTolerance, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_AREA+add, String.valueOf(probe.Area), null, (double)probe.Area, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_PEAK+add, String.valueOf(probe.Peak), null, (double)probe.Peak, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_TIME_LOWER+add, String.valueOf(probe.LowerValley), null, (double)probe.LowerValley, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_TIME_UPPER+add, String.valueOf(probe.UpperValley), null, (double)probe.UpperValley, longestEntries);
    if (probe instanceof Probe3D){
      Probe3D probe3D = (Probe3D)probe;
      createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_MZ_LOWER+add, String.valueOf(probe3D.LowerMzBand), null, (double)probe3D.LowerMzBand, longestEntries);
      createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_MZ_UPPER+add, String.valueOf(probe3D.UpperMzBand), null, (double)probe3D.UpperMzBand, longestEntries);
      createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_ELLIPSE_TIME+add, String.valueOf(probe3D.getEllipseTimePosition()), null, (double)probe3D.getEllipseTimePosition(), longestEntries);
      createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_ELLIPSE_MZ+add, String.valueOf(probe3D.getEllipseMzPosition()), null, (double)probe3D.getEllipseMzPosition(), longestEntries);
      createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_ELLIPSE_TIME_RANGE+add, String.valueOf(probe3D.getEllipseTimeStretch()), null, (double)probe3D.getEllipseTimeStretch(), longestEntries);
      createValueCell(row, Cell.CELL_TYPE_NUMERIC, MSN_ROW_FRAGMENT_ELLIPSE_MZ_RANGE+add, String.valueOf(probe3D.getEllipseMzStretch()), null, (double)probe3D.getEllipseMzStretch(), longestEntries);
    }
  }
  
  
  /**
   * creates a value cell
   * @param row the row where it has to be created
   * @param cellType the cell type according to the types stored in the class org.apache.poi.ss.usermodel.Cell
   * @param column the column number
   * @param strValue the string value representation
   * @param intValue the int value representation (can be null)
   * @param doubleValue the double value representation (can be null)
   * @param longestEntries the length of the longest entries for each column - for setting the column width
   * @return created cell
   */
  private static Cell createValueCell(Row row, int cellType, int column, String strValue, Integer intValue, Double doubleValue, Hashtable<Integer,Integer> longestEntries) {
    Cell cell = row.createCell(column,cellType);
    int length = strValue.length();
    if (cellType==Cell.CELL_TYPE_NUMERIC) {
      if (intValue!=null)
        cell.setCellValue(intValue);
      else if (doubleValue!=null)
        cell.setCellValue(doubleValue);
    }else {
      cell.setCellValue(strValue);
    }
    if (!longestEntries.containsKey(column) || length>longestEntries.get(column))
      longestEntries.put(column, length);
    return cell;
  }
  
  
  /**
   * writes the header row for subsequent intensity evidence
   * @param row Excel row that shall be used for writing
   * @param headerStyle style for the header row
   * @param longestHeaders the names of the longest header entries for each column - for setting the column width
   */
  private static void writeMSnIntensityHeader(Row row, CellStyle headerStyle, Hashtable<Integer,String> longestHeaders){
    createHeaderCell(row, MSN_ROW_INTENSITY_RULE, headerStyle, LipidomicsConstants.EXCEL_MSN_INTENSITY_RULE, longestHeaders);
    createHeaderCell(row, MSN_ROW_INTENSITY_ORIGINAL, headerStyle, LipidomicsConstants.EXCEL_MSN_INTENSITY_ORIGINAL, longestHeaders);
    createHeaderCell(row, MSN_ROW_INTENSITY_VALUES, headerStyle, LipidomicsConstants.EXCEL_MSN_INTENSITY_VALUES, longestHeaders);
    createHeaderCell(row, MSN_ROW_INTENSITY_MISSED, headerStyle, LipidomicsConstants.EXCEL_MSN_INTENSITY_MISSED, longestHeaders);    
  }
  
  /**
   * writes a line containing the intensity evidence
   * @param row Excel row that shall be used for writing
   * @param ruleVO IntensityRuleVO to be written
   * @param param the lipid identification containing MSn evidence
   * @param faHydroxyEncoding the character encoding of the number of hydroxylation sites for the
   * @param lcbHydroxyEncoding the character encoding of the number of hydroxylation sites for the LCB
   * @return lengths for adaption of cell width: int[0] longest amount of characters for MSN_ROW_INTENSITY_RULE cell; int[1] longest amount of characters for MSN_ROW_INTENSITY_ORIGINAL cell; int[2] longest amount of characters for MSN_ROW_INTENSITY_VALUES cell
   * @throws RulesException if something is not possible
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private static boolean writeMSnIntensity(Sheet sheet, int count, IntensityRuleVO ruleVO, LipidomicsMSnSet param, Hashtable<String,String> uniqueRules,
     Hashtable<Integer,Integer> longestEntries, HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding) throws RulesException, LipidCombinameEncodingException{
    String ruleInterpretation = ruleVO.getReadableRuleInterpretation(faHydroxyEncoding, lcbHydroxyEncoding);
    String rule = ruleVO.getRuleIdentifier();
    Hashtable<String,Float> fragmentAreas = param.getFragmentAreas(ruleVO);
    Hashtable<String,String> missedFragments = new Hashtable<String,String>();
    for (String frag : fragmentAreas.keySet()){
      if (fragmentAreas.get(frag)<0f){
        missedFragments.put(frag, frag);
        fragmentAreas.put(frag,0f);
      }
    }
    String missed = "";
    for (String name : missedFragments.keySet()) missed += name+";";
    if (missed.length()>0) missed = missed.substring(0,missed.length()-1);
    String valueInterpretation = ruleVO.getRuleValueInterpretation(fragmentAreas, param.getBasePeak(ruleVO));
    String uniqueId = ruleInterpretation+";"+rule+";"+valueInterpretation;
    if (uniqueRules.containsKey(uniqueId)) return false;
    uniqueRules.put(uniqueId, uniqueId);
    Row row = sheet.createRow(count);
    createValueCell(row, Cell.CELL_TYPE_STRING, MSN_ROW_INTENSITY_RULE, ruleInterpretation, null, null, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_STRING, MSN_ROW_INTENSITY_ORIGINAL, rule, null, null, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_STRING, MSN_ROW_INTENSITY_VALUES, valueInterpretation, null, null, longestEntries);
    createValueCell(row, Cell.CELL_TYPE_STRING, MSN_ROW_INTENSITY_MISSED, missed, null, null, longestEntries);
    return true;
  }
  
  
  public static boolean hasRtInfo(Hashtable<String,Vector<LipidParameterSet>> sheetParams){
    boolean hasRt = true;
    for (Vector<LipidParameterSet> params : sheetParams.values()){
      boolean hasEntry = false;
      for (LipidParameterSet param : params){
        if (param.getRt()!=null && param.getRt().length()>0) hasRt = true;
        else hasRt = false;
        hasEntry = true;
        break;
      }
      if (hasEntry) break;
    }
    return hasRt;
  }
  
  public static void setAnalyzerProperties(LipidomicsAnalyzer analyzer){
    if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE){
      analyzer.setShotgunParameters(LipidomicsConstants.getShogunProcessing());
    }else if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_PRM){
      analyzer.setPrmParameters(LipidomicsConstants.getChromSmoothRange(),LipidomicsConstants.getChromSmoothRepeats(),
          LipidomicsConstants.getPeakDiscardingAreaFactor(),LipidomicsConstants.getRelativeAreaCutoff(),
          LipidomicsConstants.getRelativeFarAreaCutoff(),LipidomicsConstants.getRelativeFarAreaTimeSpace(),
          LipidomicsConstants.getMs2MzTolerance());
    }else{
      analyzer.set3DParameters(LipidomicsConstants.getChromSmoothRange(),LipidomicsConstants.getChromSmoothRepeats(),
        LipidomicsConstants.removeIfOtherIsotopePresent(),LipidomicsConstants.useNoiseCutoff(), LipidomicsConstants.getNoiseCutoffDeviationValue(),
        LipidomicsConstants.getMinimumRelativeIntensity(), LipidomicsConstants.getScanStep(),
        LipidomicsConstants.getProfileMzRange(), LipidomicsConstants.getProfileTimeTolerance_(),
        LipidomicsConstants.getProfileIntThreshold_(), LipidomicsConstants.getBroaderProfileTimeTolerance_(),
        LipidomicsConstants.getProfileSmoothRange(),LipidomicsConstants.getProfileSmoothRepeats(),
        LipidomicsConstants.getProfileMeanSmoothRepeats(), LipidomicsConstants.getProfileMzMinRange(),
        LipidomicsConstants.getProfileSteepnessChange1(),LipidomicsConstants.getProfileSteepnessChange2(), 
        LipidomicsConstants.getProfileIntensityCutoff1(), LipidomicsConstants.getProfileIntensityCutoff2(), 
        LipidomicsConstants.getProfileGeneralIntCutoff(),LipidomicsConstants.getProfilePeakAcceptanceRange(), 
        LipidomicsConstants.getProfileSmoothingCorrection(), LipidomicsConstants.getProfileMaxRange(),
        LipidomicsConstants.getSmallChromMzRange(),LipidomicsConstants.getSmallChromSmoothRepeats(),
        LipidomicsConstants.getSmallChromMeanSmoothRepeats(), LipidomicsConstants.getSmallChromSmoothRange(),
        LipidomicsConstants.getSmallChromIntensityCutoff(),LipidomicsConstants.getBroadChromSmoothRepeats(),
        LipidomicsConstants.getBroadChromMeanSmoothRepeats(), LipidomicsConstants.getBroadChromSmoothRange(),
        LipidomicsConstants.getBroadChromIntensityCutoff(), LipidomicsConstants.getBroadChromSteepnessChangeNoSmall(),
        LipidomicsConstants.getBroadIntensityCutoffNoSmall(),
        LipidomicsConstants.getFinalProbeTimeCompTolerance(), LipidomicsConstants.getFinalProbeMzCompTolerance(),
        LipidomicsConstants.getOverlapDistanceDeviationFactor(), LipidomicsConstants.getOverlapPossibleIntensityThreshold(), 
        LipidomicsConstants.getOverlapSureIntensityThreshold(), LipidomicsConstants.getOverlapPeakDistanceDivisor(),
        LipidomicsConstants.getOverlapFullDistanceDivisor(), LipidomicsConstants.getPeakDiscardingAreaFactor(),
        LipidomicsConstants.getIsotopeInBetweenTime(),LipidomicsConstants.getIsoInBetweenAreaFactor(),LipidomicsConstants.getIsoInBetweenMaxTimeDistance(),
        LipidomicsConstants.getIsoNearNormalProbeTime(),LipidomicsConstants.getRelativeAreaCutoff(),
        LipidomicsConstants.getRelativeFarAreaCutoff(),LipidomicsConstants.getRelativeFarAreaTimeSpace(),
        LipidomicsConstants.getRelativeIsoInBetweenCutoff(),LipidomicsConstants.getTwinPeakMzTolerance(),
        LipidomicsConstants.getClosePeakTimeTolerance(),LipidomicsConstants.getTwinInBetweenCutoff(), LipidomicsConstants.getUnionInBetweenCutoff(),
        LipidomicsConstants.getMs2MzTolerance());
    }
  }
  
  public static String[] extractFormulaAndAdductName(String contents) throws ExcelInputFileException{
    String[] formulaAndName = new String[4];
    String formula = "";
    String adductName = "";
    String charge = "1";
    String multi = "1";
    String formNameString = contents.substring(contents.indexOf("(")+1,contents.indexOf(")"));
    if (formNameString.indexOf("form")!=-1&&formNameString.indexOf("name")!=-1){
      String nearFormString = formNameString.substring(formNameString.indexOf("form"));
      formula = nearFormString.substring(nearFormString.indexOf("[")+1,nearFormString.indexOf("]")).trim();
      String nearNameString = formNameString.substring(formNameString.indexOf("name"));
      adductName = nearNameString.substring(nearNameString.indexOf("[")+1,nearNameString.indexOf("]")).trim();
      if (formNameString.indexOf("charge=")!=-1){
        char[] afterChargeString = formNameString.substring(formNameString.indexOf("charge=")+"charge=".length()).toCharArray();
        try{
          charge = String.valueOf(extractDigitsAfterEqualSign(afterChargeString));
          if (Integer.parseInt(charge)<1){
            charge = "1";
            System.out.println("Warning: The charge entry in the column header \""+contents+"\" is smaller than 1! Setting it to \"1\"");
          }
        }catch (NumberFormatException nfx){
          System.out.println("Warning: The charge entry in the column header \""+contents+"\" is not integer format! Setting it to \"1\"");
        }
      }
      if (formNameString.indexOf("mult=")!=-1){
        char[] afterMultiString = formNameString.substring(formNameString.indexOf("mult=")+"mult=".length()).toCharArray();
        try{
          multi = String.valueOf(extractDigitsAfterEqualSign(afterMultiString));
          if (Integer.parseInt(charge)<1){
            multi = "1";
            System.out.println("Warning: The mult entry in the column header \""+contents+"\" is smaller than 1! Setting it to \"1\"");
          }
        }catch (NumberFormatException nfx){
          System.out.println("Warning: The mult entry in the column header \""+contents+"\" is not integer format! Setting it to \"1\"");
        }

      }
    } else{
      adductName = formNameString;
    }
    formulaAndName[0] = formula;
    formulaAndName[1] = adductName;
    formulaAndName[2] = charge;
    formulaAndName[3] = multi;
    return formulaAndName;
  }
    
  private static int extractDigitsAfterEqualSign(char[] afterString){
    int i=0;
    String digitsString = "";
    while (i<afterString.length&&Character.isDigit(afterString[i])){
      digitsString+=String.valueOf(afterString[i]);
      i++;
    }
    return Integer.parseInt(digitsString);
  }
  
  @SuppressWarnings("rawtypes")
  private static Vector parseQuantExcelFile(String quantFile, float minusTime, float plusTime, int amountOfIsotopes, int isotopesMustMatch, boolean searchUnknownTime, float basePeakCutoff,
      float rtShift, float lowestRetTime, float highestRetTime) throws IOException,SpectrummillParserException,ExcelInputFileException, ChemicalFormulaException, RulesException, HydroxylationEncodingException{
    return parseQuantExcelFile(quantFile, minusTime, plusTime, amountOfIsotopes, isotopesMustMatch, searchUnknownTime, basePeakCutoff,
        rtShift, lowestRetTime, highestRetTime, true);
  }
  
  
  @SuppressWarnings({ "unchecked", "rawtypes", "resource" })
  public static Vector parseQuantExcelFile(String quantFile, float minusTime, float plusTime, int amountOfIsotopes, int isotopesMustMatch, boolean searchUnknownTime, float basePeakCutoff,
      float rtShift, float lowestRetTime, float highestRetTime, boolean respectMassShift) throws IOException,SpectrummillParserException,ExcelInputFileException, ChemicalFormulaException, RulesException, HydroxylationEncodingException{
    InputStream myxls = new FileInputStream(quantFile);
    Workbook workbook = null;
    if (quantFile.endsWith(".xlsx")) workbook = new XSSFWorkbook(myxls);
    else if (quantFile.endsWith(".xls")) workbook = new HSSFWorkbook(myxls);
    Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();
    LinkedHashMap<String,Integer> classSequence = new LinkedHashMap<String,Integer>();
    Hashtable<String,Boolean> adductInsensitiveRtFilter = new Hashtable<String,Boolean>();
    Hashtable<String,Vector<String>> analyteSequence = new Hashtable<String,Vector<String>>();
    boolean excelOK = false;
    ElementConfigParser aaParser = Settings.getElementParser();
    for (int sheetNumber = 0; sheetNumber!=workbook.getNumberOfSheets(); sheetNumber++){
      Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = new Hashtable<String,Hashtable<String,QuantVO>>();
      Vector<String> analytes = new Vector<String>();
      Sheet sheet = workbook.getSheetAt(sheetNumber);
      boolean rtFilterInsensitive = false;
      int sideChainColumn = -1;
      int doubleBondColumn = -1;

      Hashtable<Integer,String> massOfInterestColumns = new Hashtable<Integer,String>();
      Hashtable<String,Hashtable<String,Integer>> adductComposition = new Hashtable<String,Hashtable<String,Integer>>();
      Hashtable<String,Integer> charges = new Hashtable<String,Integer>();
      Hashtable<String,Integer> multi = new Hashtable<String,Integer>();
      int retTimeColumn = -1;
      boolean foundColumns = false;
      float fixedStartTime = 0;
      float fixedEndTime = Float.MAX_VALUE;
      Hashtable<Integer,String> elementColumns = new  Hashtable<Integer,String>();
      int msLevel = 1;
      int ohNumber = LipidomicsConstants.EXCEL_NO_OH_INFO;
      RangeInteger ohRange = null;
      for (int rowCount=0;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
        Row row = sheet.getRow(rowCount);
        String sideChain = "";
        int doubleBonds = -1;
        
        Hashtable<String,Integer> elementalComposition = new Hashtable<String,Integer>();
        Hashtable <String,Double> massesOfInterest = new Hashtable <String,Double>();
        float retTime = -1;
        Hashtable<Integer,String> possibleElementColumns = new  Hashtable<Integer,String>();
        for (int i=0; row!=null && i!=(row.getLastCellNum()+1);i++){
          Cell cell = row.getCell(i);
          String contents = "";
          Double numeric = null;
          int cellType = -1;
          if (cell!=null) cellType = cell.getCellType();
          if (cellType==Cell.CELL_TYPE_STRING){
            contents = cell.getStringCellValue();
            try{ 
              if (contents!=null)numeric = new Double(contents.replaceAll(",", "."));
            }catch(NumberFormatException nfx){};
          }else if (cellType==Cell.CELL_TYPE_NUMERIC || cellType==Cell.CELL_TYPE_FORMULA){
           numeric = cell.getNumericCellValue();
           contents = String.valueOf(numeric);
          }
          //String contents = sheet.getCell(i,rowCount).getContents();
          if (contents!=null)
            contents = contents.trim();
          if (!foundColumns){
            if (contents.equalsIgnoreCase("Seitenkette")||contents.equalsIgnoreCase("Name")){
              sideChainColumn = i;
            } else if (contents.equalsIgnoreCase("dbs")||contents.equalsIgnoreCase("dbs_TAG")){
              doubleBondColumn = i;
            } 

            else if (contents.startsWith("mass")&&contents.contains("(")&&contents.contains(")")){
              String[] formulaAndName = extractFormulaAndAdductName(contents);
              adductComposition.put(formulaAndName[1],StaticUtils.categorizeFormula(formulaAndName[0]));
              massOfInterestColumns.put(i,formulaAndName[1]);
              charges.put(formulaAndName[1], Integer.parseInt(formulaAndName[2]));
              multi.put(formulaAndName[1], Integer.parseInt(formulaAndName[3]));
            }
            else if (contents.equalsIgnoreCase("tR (min)")){
              retTimeColumn = i;
            }
            else if (contents.startsWith("Start-RT:")){
              try{
                fixedStartTime = Float.parseFloat(contents.substring("Start-RT:".length()).trim().replaceAll(",", "."));
              }catch(NumberFormatException nfx){nfx.printStackTrace();};
            }
            else if (contents.startsWith("Stop-RT:")){
              try{
                fixedEndTime = Float.parseFloat(contents.substring("Stop-RT:".length()).trim().replaceAll(",", "."));
              }catch(NumberFormatException nfx){};
            }
            else if (contents.startsWith("Mass-Trace:")){
              try{
                msLevel = Integer.parseInt(contents.substring("Mass-Trace:".length()).trim().replaceAll(",", "."));
              }catch(NumberFormatException nfx){};  
            }
            else if (contents.trim().length()==1||contents.trim().length()==2){
              boolean ok = false;
              if (Character.isUpperCase(contents.trim().toCharArray()[0])){
                if (contents.trim().length()==2){
                  if (Character.isLowerCase(contents.trim().toCharArray()[1]))
                    ok = true;
                }else{
                  ok = true;
                }
                if (ok){
                  String element = contents.trim();
                  if (aaParser.isElementAvailable(element)){
                    possibleElementColumns.put(i, element);
                  }else{
                    if (sideChainColumn>-1)System.out.println("Warning: The elemental column \""+element+"\" does not seem to be a chemical element!");                
                  }
                }
              }
            } else if (contents.trim().equalsIgnoreCase("adductInsensitiveRtFilter"))
              rtFilterInsensitive = true;
            else if (contents.startsWith("OH-Number:")){
              String ohString = contents.substring("OH-Number:".length()).trim().replaceAll(",", ".");
              try{
                ohNumber = Integer.parseInt(ohString);
              //this value need not necessarily to be a number - it might be the hydroxy encoded string
              }catch(NumberFormatException nfx){
                ohNumber = Settings.getLcbHydroxyEncoding().getHydroxyNumber(ohString);
              }
            } else if (contents.startsWith("OH-Range:")){
              String ohRangeString = contents.substring("OH-Range:".length()).trim().replaceAll(",", ".");
              String[] ohRangeParts = ohRangeString.split("-");
              boolean error = false;
              try{
                if (ohRangeParts.length>2) error = true;
                int start = Integer.parseInt(ohRangeParts[0]);
                int stop = start;
                if (ohRangeParts.length==2)
                  stop = Integer.parseInt(ohRangeParts[1]);
                ohRange = new RangeInteger(start,stop);
              }catch(NumberFormatException nfx){error = true;}
              if (error)
                throw new HydroxylationEncodingException("The value \"OH-Range\" must be a single integer, or a range in the format $lower$-$higher$; the value \""+ohRangeString+"\" in sheet "+sheet.getSheetName()+" does not comply!");
            }
          }else{            
            if (i==sideChainColumn&&contents!=null&contents.length()>0){
//          for Marlene metabolomics implementation - exclude the if - only "sideChain = contents;" must remain 
              if (numeric!=null){
                sideChain =  String.valueOf((int)Math.round(numeric));
              }else
                sideChain = contents;
            }
            if (i==doubleBondColumn&&contents!=null&&contents.length()>0){
              doubleBonds = numeric.intValue();
            }
            if (elementColumns.containsKey(i)&&contents.length()>0){
              int value = 0;
              // this is for columns such as "M" for mass, there the values are float and have
              // to be removed from the chemical formula - however there is no way to figure out
              // directly from Excel if there is an integer entry in the cell or not!
              if (numeric!=null && contents.endsWith(".0")){
                value = (int)Math.round(numeric);
              }else{
                try{
                  value = Integer.parseInt(contents);
                }catch (NumberFormatException nfx3){  
                  System.out.println("Warning: The elemental column \""+elementColumns.get(i)+"\" does not seem to be a chemical element, since float values are there -> I remove it!");
                  elementColumns.remove(i);
                }
              }
              if (value>0)
                elementalComposition.put(elementColumns.get(i), value);
            }

            if (massOfInterestColumns.containsKey(i)&&contents!=null&contents.length()>0){
              double massOfInterest = numeric;
              if (respectMassShift) massOfInterest += LipidomicsConstants.getMassShift();
              massesOfInterest.put(massOfInterestColumns.get(i), massOfInterest);
            }
            if (i==retTimeColumn&&contents!=null&contents.length()>0){
              retTime = numeric.floatValue();
              retTime += rtShift;
            }
          }
        }

        if (sideChainColumn != -1 && possibleElementColumns.size()>0 && 
            massOfInterestColumns.size()>0 && retTimeColumn != -1&&
            !foundColumns){
          foundColumns = true;
          elementColumns = new Hashtable<Integer,String>(possibleElementColumns);
        }
        if (foundColumns&&sideChain!=null&&sideChain.length()>0&&massesOfInterest.size()>0/*&&retTime>0*/){
          float usedMinusTime = new Float(minusTime);
          if (usedMinusTime<0) usedMinusTime = usedMinusTime*-1f;
          float usedPlusTime = new Float(plusTime);
          if (usedPlusTime<0) usedPlusTime = usedPlusTime*-1f;
          if (retTime>0||searchUnknownTime){
            // if a global retention time for a class is set the other retention times have to be recalculated         
            if (fixedStartTime>0 || fixedEndTime<Float.MAX_VALUE && (fixedStartTime<fixedEndTime)){
              if (retTime>0){
                if (fixedStartTime<retTime && retTime<fixedEndTime){
                  if (fixedStartTime>(retTime-usedMinusTime))
                    usedMinusTime = retTime-fixedStartTime;
                  if (fixedEndTime<(retTime+usedPlusTime))
                    usedPlusTime = fixedEndTime-retTime;
                }else{
                // Here I do not know what to do; I guess best is to keep the set RT times, then I can quantify items that are outside
                // the general class range, if necessary!
                }
              }else{
                float startTime = lowestRetTime;
                float stopTime = highestRetTime;
                if (fixedStartTime>startTime)
                  startTime = fixedStartTime;
                if (fixedEndTime<stopTime)
                  stopTime = fixedEndTime;
                retTime = (startTime+stopTime)/2;
                usedMinusTime = retTime-startTime;
                usedPlusTime = stopTime-retTime;
              }
            }
            
            int startOh = 0;
            int stopOh = 0;
            if (ohNumber>LipidomicsConstants.EXCEL_NO_OH_INFO) {
              startOh = ohNumber;
              stopOh = ohNumber;
              if (ohRange!=null) {
                startOh = ohRange.getStart();
                stopOh = ohRange.getStop();
              }
            }
            for (int oh=startOh; oh<(stopOh+1); oh++) {
              Hashtable<String,QuantVO> quantsOfAnalyte = new Hashtable<String,QuantVO>();
              String analEncoded = null;
              Hashtable<String,Integer> correctedElementalComposition = new Hashtable<String,Integer>(elementalComposition);
              double ohDiff = 0d;
              int ohToUse = LipidomicsConstants.EXCEL_NO_OH_INFO;
              if (ohNumber>LipidomicsConstants.EXCEL_NO_OH_INFO) {
                ohToUse = oh;
                if (oh!=ohNumber) {
                  int oxygens = 0;
                  if (correctedElementalComposition.containsKey("O")) oxygens = correctedElementalComposition.get("O");
                  oxygens += (oh-ohNumber);
                  if (oxygens<0) continue;
                  correctedElementalComposition.put("O", oxygens);
                  ohDiff = (oh-ohNumber)*Settings.getElementParser().getElementDetails("O").getMonoMass();
                }
              }
              
              for (String modName : massesOfInterest.keySet()){
                double massOfInterest = massesOfInterest.get(modName)+ohDiff;
                Hashtable<String,Integer> modElements = adductComposition.get(modName);
                Integer charge = charges.get(modName);
                Integer mult = multi.get(modName);
                
                String[] formulas = getFormulasAsString(correctedElementalComposition,modElements,mult);
                String analyteFormula = formulas[0];
                String modificationFormula = formulas[1];
                String chemicalFormula = formulas[2];
                //no negative elements are allowed after an applied modification
                if (chemicalFormula.indexOf("-")!=-1)
                  continue;
                Object[] distris = getTheoreticalIsoDistributions(aaParser,isotopesMustMatch,amountOfIsotopes,chemicalFormula);
                Vector<Double> mustMatchProbabs = (Vector<Double>)distris[0];
                Vector<Double> probabs = (Vector<Double>)distris[1];
                int negativeStartValue = (Integer)distris[2];
              
                QuantVO quantVO = new QuantVO(sheet.getSheetName(), sideChain, doubleBonds,
                    ohToUse,analyteFormula, massOfInterest, charge, modName,
                    modificationFormula, retTime, usedMinusTime, usedPlusTime,
                    mustMatchProbabs, probabs,negativeStartValue);
                analEncoded = quantVO.getAnalyteName();
                quantsOfAnalyte.put(modName, quantVO);
              }
              String analyteName = StaticUtils.generateLipidNameString((analEncoded!=null ? analEncoded : sideChain), doubleBonds);
              analytes.add(analyteName);
              quantsOfClass.put(analyteName, quantsOfAnalyte);
            }
          }
        }
      } 
      quantObjects.put(sheet.getSheetName(), quantsOfClass);
      analyteSequence.put(sheet.getSheetName(), analytes);
      classSequence.put(sheet.getSheetName(),msLevel);
      adductInsensitiveRtFilter.put(sheet.getSheetName(), rtFilterInsensitive);
      if (foundColumns) excelOK = true;
    }
    myxls.close();
    if (!excelOK) throw new ExcelInputFileException("The Excel file is not valid!");
    checkForIsobaricSpecies(classSequence,analyteSequence,quantObjects);
    
    Vector results = new Vector();
    results.add(classSequence);
    results.add(analyteSequence);
    results.add(adductInsensitiveRtFilter);
    results.add(quantObjects);
    return results;
  }
  
  @SuppressWarnings({ "rawtypes" })
  public static Vector getCorrectAnalyteSequence(String filePath, boolean ionMode) throws Exception{
    File quant = new File(filePath);
    Vector quantContent = null;
    if (quant.isFile() && (quant.getName().endsWith(".xls") || quant.getName().endsWith(".xlsx"))){
      quantContent = QuantificationThread.parseQuantExcelFile(filePath, 0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f);
    } else if ((quant.isFile() && quant.getName().endsWith(".txt")) || quant.isDirectory()){
      System.out.println("Ion mode: "+ionMode);
      quantContent = QuantificationThread.parseAlex123TargetList(filePath, 0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f, ionMode);
    }
    if (quantContent!=null)
      return quantContent;
    else
      return null;
  }
  
  private static CellStyle getHeaderStyle(Workbook wb){
    CellStyle arial12style = wb.createCellStyle();
    Font arial12font = wb.createFont();
    arial12font.setBoldweight(Font.BOLDWEIGHT_BOLD);
    arial12font.setFontName("Arial");
    arial12font.setFontHeightInPoints((short)12);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(CellStyle.ALIGN_CENTER);
    return arial12style;
  }
  
  private float[] initThreadMonitors(String[] chromPaths, int numberOfProcessors, float basePeakCutoff) throws CgException{
    availableThreads_ = new Hashtable<Integer,Boolean>();
    analyzers_ = new Hashtable<Integer,LipidomicsAnalyzer>();
    threadToClass_ = new Hashtable<Integer,String> ();
    threadToAnalyte_ = new Hashtable<Integer,String>();
    threadToMod_ = new Hashtable<Integer,String>();
    float[] maxRetTimes = new float[2];
    
    for (int i=0; i!=numberOfProcessors;i++){
      availableThreads_.put(i, true);
      LipidomicsAnalyzer analyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],Settings.useCuda());
      if (i==0){
        float highestRetTime = 0;
        float lowestRetTime = Float.MAX_VALUE;
        for (float retTime : analyzer.getRetentionTimes().values()){
          if (retTime>highestRetTime) highestRetTime = retTime;
          if (retTime<lowestRetTime) lowestRetTime = retTime;
        }
        maxRetTimes[1] = highestRetTime/60;
        maxRetTimes[0] = lowestRetTime/60;
      }
      QuantificationThread.setAnalyzerProperties(analyzer);
      analyzer.setGeneralBasePeakCutoff(basePeakCutoff/3f);
      analyzers_.put(i, analyzer);
    }
    return maxRetTimes;
  }
  
  
  private class ThreadSupervisor extends TimerTask{
    private LinkedHashMap<String,Integer> classSequence_;
    private Hashtable<String,Vector<String>> analyteSequence_;
    private Hashtable<String,Boolean> adductInsensitiveRtFilter_;
    private Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects_;
    private float bpCutoff_;
    private String rsFile_;

    
    @SuppressWarnings({ "rawtypes", "unchecked" })
    public ThreadSupervisor(Vector excelContent, float  basePeakCutoff, String resultFile){
      currentLipidCount_ = 0;
      classSequence_ = (LinkedHashMap<String,Integer>)excelContent.get(0);
      analyteSequence_ = (Hashtable<String,Vector<String>>)excelContent.get(1);
      adductInsensitiveRtFilter_ = (Hashtable<String,Boolean>)excelContent.get(2);
      quantObjects_ = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)excelContent.get(3);
      bpCutoff_ = basePeakCutoff;
      rsFile_ = resultFile;
      initThreadHashes();
    }

    public void run()
    { 
      try{
        if (!finished_)
          handleTimerEvent(classSequence_,analyteSequence_,adductInsensitiveRtFilter_,quantObjects_,bpCutoff_,rsFile_,chromFileName_);
      } catch (Exception ex){
        ex.printStackTrace();
        errorString_ = ex.toString();
        finished_ = true;
        for (Integer analyzer : analyzers_.keySet()){
          if (analyzers_.get(analyzer).getUseCuda()){
            analyzers_.get(analyzer).getSavGolJNI().Frees();
          }
        }
      }
    }
    
    private void initThreadHashes(){
      quantStatus_ = new Hashtable<String,Hashtable<String,Hashtable<String,Integer>>>();
      results_ = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
      ms2Removed_ = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
      unsplittedPeaks_ = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
      threads_ = new  Hashtable<Integer,SingleQuantThread>();
      onlyMS1DataPresent_ = new  Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();
      for (String className : classSequence_.keySet()){
        Hashtable<String,Hashtable<String,QuantVO>> classQuant = quantObjects_.get(className);
        Hashtable<String,Hashtable<String,Integer>> statusClass = new Hashtable<String,Hashtable<String,Integer>>();
        Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> resultsClass = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
        Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> resultsNeg = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
        Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> unsplitted = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
        for (String analyteName : analyteSequence_.get(className)){
          Hashtable<String,QuantVO> analyteQuant = classQuant.get(analyteName);
          Hashtable<String,Integer> status = new Hashtable<String,Integer>(); 
          for (String mod : analyteQuant.keySet()){
            if (analyteQuant.get(mod).isQuantifiedByOtherIsobar()) status.put(mod, STATUS_FINISHED);
            else status.put(mod, STATUS_WAITING);
          }
          statusClass.put(analyteName, status);
          resultsClass.put(analyteName, new Hashtable<String,Hashtable<String,LipidParameterSet>>());
          resultsNeg.put(analyteName, new Hashtable<String,Hashtable<String,LipidParameterSet>>());
          unsplitted.put(analyteName, new Hashtable<String,Hashtable<String,LipidParameterSet>>());
        }
        quantStatus_.put(className, statusClass);
        results_.put(className, resultsClass);
        ms2Removed_.put(className, resultsNeg);
        unsplittedPeaks_.put(className, unsplitted);
      }  
    }

  }
  
  
  private void handleTimerEvent(LinkedHashMap<String,Integer> classSequence,Hashtable<String,Vector<String>> analyteSequence,
      Hashtable<String,Boolean> adductInsensitiveRtFilter, Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects,
      float basePeakCutoff, String resultFile, String chromFile){
    
    boolean allFinished = true;
    boolean error = false;
    boolean stopThread = false;
    // find out if some of the threads are finished, read the results, and make the thread available again
    for (int i=0;i!=availableThreads_.size();i++){
      if (this.availableThreads_.get(i)==false && this.threads_.get(i)!=null && this.threads_.get(i).finished()==true){
        SingleQuantThread singleThread = this.threads_.get(i);
        String className =  threadToClass_.get(i);
        String analyte = threadToAnalyte_.get(i);
        String mod = threadToMod_.get(i);
        if (singleThread.getErrorString()!=null){
          this.errorString_ = singleThread.getErrorString();
          error = true;
          stopThread = true;
        }else{
          Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant = singleThread.getResults();
          Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> ms2RemovedToQuant = singleThread.getMs2Removedhits();
          Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> beforeSplitToQuant = singleThread.getPeaksBeforeSplit();
          
          Hashtable<String,Hashtable<QuantVO,LipidParameterSet>> msnHitsWithSameRt = groupHitsWithSameRt(hitsAccordingToQuant);
          
          for (QuantVO oneQuant : hitsAccordingToQuant.keySet()){
            Hashtable<String,LipidParameterSet> hitsOfOneMod = hitsAccordingToQuant.get(oneQuant);
            Hashtable<String,LipidParameterSet> ms2Removed = new Hashtable<String,LipidParameterSet>();
            if (ms2RemovedToQuant.containsKey(oneQuant)) ms2Removed = ms2RemovedToQuant.get(oneQuant);
            Hashtable<String,LipidParameterSet> beforeSplit = new Hashtable<String,LipidParameterSet>();
            if (beforeSplitToQuant.containsKey(oneQuant)) beforeSplit = beforeSplitToQuant.get(oneQuant);
            if (hitsOfOneMod.size()>0){
// the lines with the 4* are necessary if the prediction is based on consecutive model predictions -
// there is as well a section in the PostQuantificationProcessor that has to be activated
/****            if (msnRoundFinished_ && latestRtPredictions_.containsKey(className) && latestRtPredictions_.get(className).containsKey(mod) && latestRtPredictions_.get(className).get(mod).getCounterModel()!=null){
              RtPredictVO predVO = latestRtPredictions_.get(className).get(mod);
              Hashtable<String,LipidParameterSet> hitsOfOneModCorrected = new Hashtable<String,LipidParameterSet>();
              Hashtable<String,LipidParameterSet> ms2Removed = new Hashtable<String,LipidParameterSet>();
              System.out.println("xxxxxxxxxxxxxxxxxxx");
              for (String rt : hitsOfOneMod.keySet()){*/
// the lines with 1* are necessary if a counter model is necessary           
/*                try {
                  if (PostQuantificationProcessor.fallsOnCounterModelSide(className,mod,analyte,rt,predVO.getModel(),predVO.getCounterModel())){
                    System.out.println("!!!!!!!!!!!!!!!!!!!!!!");
                    ms2Removed.put(rt, hitsOfOneMod.get(rt));
                  }else{*/
////                    hitsOfOneModCorrected.put(rt, hitsOfOneMod.get(rt));
/*                  }
                }
                catch (RulesException | NoRuleException | IOException  | SpectrummillParserException | LMException e) {
                  hitsOfOneModCorrected.put(rt, hitsOfOneMod.get(rt));
                }*/
/****              }
              if (hitsOfOneModCorrected.size()>0) results_.get(className).get(analyte).put(mod,hitsOfOneModCorrected);
              if (ms2Removed.size()>0) ms2Removed_.get(className).get(analyte).put(mod, ms2Removed);
            } else {*/
              Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> resultsClass = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
              if (results_.containsKey(oneQuant.getAnalyteClass())) resultsClass = results_.get(oneQuant.getAnalyteClass());
              Hashtable<String,Hashtable<String,LipidParameterSet>> resultsAnalyte = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
              if (resultsClass.containsKey(oneQuant.getIdString())) resultsAnalyte = resultsClass.get(oneQuant.getIdString());
              //if there are quantitations of isobaric species, the m/z value may shift into not allowed regions -
              //this routine removes highly likely false positive that fall outside the m/z region
              Hashtable<String,LipidParameterSet> newHits = new Hashtable<String,LipidParameterSet>();
              for (String rt : hitsOfOneMod.keySet()){
                LipidParameterSet set = hitsOfOneMod.get(rt);
                boolean insideAllowedMz = false;
                if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE){
                  insideAllowedMz = true;
                } else {
                  for (CgProbe probe : set.getIsotopicProbes().get(0)){
                    if (((float)oneQuant.getAnalyteMass()-LipidomicsConstants.getProfilePeakAcceptanceRange())<=probe.Mz && probe.Mz<=((float)oneQuant.getAnalyteMass()+LipidomicsConstants.getProfilePeakAcceptanceRange())){
                      insideAllowedMz = true;
                      break;
                    }
                  }
                }
                if (insideAllowedMz) newHits.put(rt, set);
              }
              //if there are quantitations of isobaric species, where the the m/z is not completely the same,
              //both species are quantified, but with 3D parameters, it is hard to check if the identifications are the same
              //thus, if the same species is added the first time, all of the analytes are taken as they come from the quantitation processor
              //for the following identifications, only species are allowed that do not overlap with already identifed RT regions
              if (resultsAnalyte.containsKey(oneQuant.getModName())){
                Hashtable<String,LipidParameterSet> hashMods = resultsAnalyte.get(oneQuant.getModName());
                // if there are other hits found, and they are not overlapped by anything -> add them to the hash
                for (LipidParameterSet set : newHits.values()){
                  Vector<CgProbe> probes = set.getIsotopicProbes().get(0);
                  boolean overlap = false;
                  for (LipidParameterSet there : hashMods.values()){
                    if (overlap) break;
                    for (CgProbe probe : probes){
                      if (overlap) break;
                      for (CgProbe other : there.getIsotopicProbes().get(0)){
                        if (overlap) break;
                        if (probe.isCoveredByThisProbe(other)) overlap = true;
                        if (probe.Peak<other.Peak){
                          if (probe.UpperValley>other.LowerValley) overlap = true;
                        }else{
                          if (probe.LowerValley<other.UpperValley) overlap = true;
                        }
                      }
                    }
                  }
                  if (!overlap) hashMods.put(set.getRt(), set);
                }
              } else 
                resultsAnalyte.put(oneQuant.getModName(), newHits);
              resultsClass.put(oneQuant.getIdString(), resultsAnalyte);
              results_.put(oneQuant.getAnalyteClass(), resultsClass); 
              if (beforeSplit.size()>0){
                unsplittedPeaks_.get(oneQuant.getAnalyteClass()).get(oneQuant.getIdString()).put(oneQuant.getModName(), beforeSplit);
              }
/****            }*/
            }else if (singleThread.isMsnFirst() && !singleThread.areMSnSpectraPresent()){
              try{
 //               QuantVO quantVO = singleThread.getQuantSet();
                int msIdentOrder = RulesContainer.getMSIdentificationOrder(StaticUtils.getRuleName(oneQuant.getAnalyteClass(),oneQuant.getModName()));
                if (msIdentOrder==RulesContainer.ORDER_MSN_FIRST){
                  oneQuant.removeOtherIsobaricSpecies();
                  Hashtable<String,Hashtable<String,QuantVO>> quantClass = new Hashtable<String,Hashtable<String,QuantVO>>();
                  if (onlyMS1DataPresent_.containsKey(oneQuant.getAnalyteClass())) quantClass = onlyMS1DataPresent_.get(oneQuant.getAnalyteClass());
                  Hashtable<String,QuantVO> quantAnalyte = new Hashtable<String,QuantVO>();
                  if (quantClass.containsKey(oneQuant.getIdString())) quantAnalyte = quantClass.get(oneQuant.getIdString());
                  quantAnalyte.put(oneQuant.getModName(), oneQuant);
                  quantClass.put(oneQuant.getIdString(), quantAnalyte);
                  onlyMS1DataPresent_.put(oneQuant.getAnalyteClass(), quantClass);
                }
              } catch(Exception ex){
            }
          }
          try{
            if (RulesContainer.isRtPostprocessing(StaticUtils.getRuleName(oneQuant.getAnalyteClass(),oneQuant.getModName())) &&
                RulesContainer.correctRtForParallelModel(StaticUtils.getRuleName(oneQuant.getAnalyteClass(),oneQuant.getModName()))){
              if (ms2Removed.size()>0) ms2Removed_.get(oneQuant.getAnalyteClass()).get(oneQuant.getIdString()).put(oneQuant.getModName(), ms2Removed);
            }
          }catch(Exception ex){}
          }         
          //check if an overlapping MS2 instance has been identified and, after the m/z check, the overlap still exists
          //if not, apply the base peak rules
          for (String rt : msnHitsWithSameRt.keySet()) {
            if (msnHitsWithSameRt.get(rt) .size()<2)
              continue;
            int count = 0;
            LipidParameterSet set = null;
            QuantVO toQuant = null; 
            for (QuantVO quant : msnHitsWithSameRt.get(rt).keySet()) {
              if (results_.containsKey(quant.getAnalyteClass()) && results_.get(quant.getAnalyteClass()).containsKey(quant.getIdString()) &&
                  results_.get(quant.getAnalyteClass()).get(quant.getIdString()).containsKey(quant.getModName()) &&
                  results_.get(quant.getAnalyteClass()).get(quant.getIdString()).get(quant.getModName()).containsKey(rt)) {
                set = results_.get(quant.getAnalyteClass()).get(quant.getIdString()).get(quant.getModName()).get(rt);
                toQuant = quant;
                count++;
              }
            }
            if (count!=1 || !(set instanceof LipidomicsMSnSet))
              continue;
            try {
              MSnAnalyzer msnAnalyzer = new MSnAnalyzer(toQuant.getAnalyteClass(),toQuant.getModName(),set,singleThread.getAnalyzer(),toQuant,true,false);
              if (msnAnalyzer.checkStatus()<=LipidomicsMSnSet.DISCARD_HIT) {
                results_.get(toQuant.getAnalyteClass()).get(toQuant.getIdString()).get(toQuant.getModName()).remove(rt);
                if (results_.get(toQuant.getAnalyteClass()).get(toQuant.getIdString()).get(toQuant.getModName()).size()==0)
                  results_.get(toQuant.getAnalyteClass()).get(toQuant.getIdString()).remove(toQuant.getModName());
                if (results_.get(toQuant.getAnalyteClass()).get(toQuant.getIdString()).size()==0)
                  results_.get(toQuant.getAnalyteClass()).remove(toQuant.getIdString());
              }else {
                results_.get(toQuant.getAnalyteClass()).get(toQuant.getIdString()).get(toQuant.getModName()).put(rt, msnAnalyzer.getResult());
              }
            } catch (RulesException | IOException | SpectrummillParserException | HydroxylationEncodingException | ChemicalFormulaException | CgException | LipidCombinameEncodingException e) {
              e.printStackTrace();
            }  
          }
        }
        quantStatus_.get(className).get(analyte).put(mod, STATUS_FINISHED);
        singleThread = null;
        availableThreads_.put(i,true);
      }
    }
    // find available threads
    Vector<Integer> availableThread = new Vector<Integer>();
    if (!error){
      for (int i=0;i!=availableThreads_.size();i++){
        if (this.availableThreads_.get(i)==true) availableThread.add(i);
      }
    }
    //if there are free threads, assign jobs to them!
    int currentThreadNumber = 0;
    if (!error && availableThread.size()>0){
      for (String className : classSequence.keySet()){
        if (currentThreadNumber>=availableThread.size())
          break;
        Hashtable<String,Hashtable<String,QuantVO>> classQuant = quantObjects.get(className);
        quantObjects.get(className);
        int msLevel = classSequence.get(className);
        for (String analyteName : analyteSequence.get(className)){
          if (currentThreadNumber>=availableThread.size())
            break;
          Hashtable<String,QuantVO> analyteQuant = classQuant.get(analyteName);
          int modCount = 0;
          for (String mod : analyteQuant.keySet()){
            if (currentThreadNumber>=availableThread.size())
              break;
            if (quantStatus_.get(className).get(analyteName).get(mod)==STATUS_WAITING){
              int threadIndex = (availableThread.get(currentThreadNumber));
              availableThreads_.put(threadIndex, false);
              quantStatus_.get(className).get(analyteName).put(mod,STATUS_CALCULATING);
              boolean msnFirst = false;
              if (LipidomicsConstants.isMS2() && !this.msnRoundFinished_){
                Vector<QuantVO> quants = new Vector<QuantVO>();
                quants.add(analyteQuant.get(mod));
                quants.addAll(analyteQuant.get(mod).getOtherIsobaricSpecies());
                for (int i=0; i!=quants.size();i++){
                  QuantVO quant = quants.get(i);
                  try{
                    int msIdentOrder = RulesContainer.getMSIdentificationOrder(StaticUtils.getRuleName(quant.getAnalyteClass(),quant.getModName()));
                    if (i==0 && (msIdentOrder==RulesContainer.ORDER_MSN_FIRST || msIdentOrder==RulesContainer.ORDER_MSN_ONLY)) msnFirst = true;
                    else if (msIdentOrder==RulesContainer.ORDER_MS1_FIRST) msnFirst = false;
                  } catch(Exception ex){
                  }                

                }
              }
              SingleQuantThread thread = new SingleQuantThread(analyzers_.get(threadIndex), analyteQuant.get(mod), msLevel, msnFirst);
              threads_.put(threadIndex, thread);
              threadToClass_.put(threadIndex,className);
              threadToAnalyte_.put(threadIndex,analyteName);
              threadToMod_.put(threadIndex,mod);
              thread.start();
              if (modCount==0){
//                currentLipidCount_++;
                currentLipid_ = className+" "+analyteName;
              }
              currentThreadNumber++;
            }
            modCount++;
          }
        }  
      }
      int currentLipidCount = 0;
      for (String className : quantStatus_.keySet()){
        for (String analyteName : quantStatus_.get(className).keySet()){
          Hashtable<String,Integer> modStati = quantStatus_.get(className).get(analyteName);
          boolean oneDifferent = false;
          for (Integer status : modStati.values()){
            if (status != STATUS_WAITING){
              oneDifferent = true;
              break;
            }
          }
          if (oneDifferent && !(this.onlyMS1DataPresent_.containsKey(className) && this.onlyMS1DataPresent_.get(className).containsKey(analyteName))) currentLipidCount++;
        }
      }
      currentLipidCount_ = currentLipidCount;
    }
    //check if all of the treads are finished
    if (!error && availableThread.size()>0 && currentThreadNumber == 0){
      for (String className : quantStatus_.keySet()){
        for (String analyte : quantStatus_.get(className).keySet()){
          for (String mod : quantStatus_.get(className).get(analyte).keySet()){
            if (quantStatus_.get(className).get(analyte).get(mod)!=STATUS_FINISHED){
              allFinished=false;
              break;
            }
          }
        }
      }
    } else {
      allFinished = false;
    }
    if (allFinished == true) stopThread = true;
    if (stopThread){
      // the following lines are for MSnFirst: here the LM models for the time prediction are calculated
      boolean areThereMS1HitsToProcess = false;
      msnRoundFinished_ = true;
      if (!error){
        if (onlyMS1DataPresent_.size()>0) areThereMS1HitsToProcess = true;
      }
      int countToProcess = 0;
      if (areThereMS1HitsToProcess){
        PostQuantificationProcessor processor = new PostQuantificationProcessor(results_,ms2Removed_,adductInsensitiveRtFilter);
        try {
          latestRtPredictions_ = processor.predictRetentionTimesBasedOnResults(onlyMS1DataPresent_,latestRtPredictions_);
          Hashtable<String,Hashtable<String,Boolean>> predictionFound = new Hashtable<String,Hashtable<String,Boolean>>();
          Hashtable<String,Hashtable<String,Hashtable<String,String>>> toRemove = new Hashtable<String,Hashtable<String,Hashtable<String,String>>>();
          // these reinits the quantification of analytes where an RT prediction was made
          for (String className : onlyMS1DataPresent_.keySet()){
            Hashtable<String,Hashtable<String,QuantVO>> ofClass = onlyMS1DataPresent_.get(className);
            for (String analyteName : ofClass.keySet()){
              Hashtable<String,QuantVO> ofAnalyte = ofClass.get(analyteName);
              for (String modName : ofAnalyte.keySet()){
                QuantVO vo = ofAnalyte.get(modName);
                if (vo.getRetTime()<0) continue;
                quantStatus_.get(className).get(analyteName).put(modName, STATUS_WAITING);
                countToProcess++;
                Hashtable<String,Boolean> modFound = new Hashtable<String,Boolean>();
                if (predictionFound.containsKey(className)) modFound = predictionFound.get(className);
                modFound.put(modName, true);
                predictionFound.put(className, modFound);
                Hashtable<String,Hashtable<String,String>> classToRemove = new Hashtable<String,Hashtable<String,String>>();
                if (toRemove.containsKey(className)) classToRemove = toRemove.get(className);
                Hashtable<String,String> analyteToRemove = new Hashtable<String,String>();
                if (classToRemove.containsKey(analyteName)) analyteToRemove = classToRemove.get(analyteName);
                analyteToRemove.put(modName, modName);
                classToRemove.put(analyteName, analyteToRemove);
                toRemove.put(className, classToRemove);
              }
            }
          }
          // remove hits that are already marked for quantification
          for (String className : toRemove.keySet()){
            for (String analyteName : toRemove.get(className).keySet()){
              for (String modName : toRemove.get(className).get(analyteName).keySet()){
                onlyMS1DataPresent_.get(className).get(analyteName).remove(modName);
              }
              if (onlyMS1DataPresent_.get(className).get(analyteName).size()==0) onlyMS1DataPresent_.get(className).remove(analyteName);
            }
            if (onlyMS1DataPresent_.get(className).size()==0) onlyMS1DataPresent_.remove(className);
          }
          toRemove = new Hashtable<String,Hashtable<String,Hashtable<String,String>>>();
          // these reinits the quantification of analytes where an RT prediction was made
          for (String className : onlyMS1DataPresent_.keySet()){
            Hashtable<String,Hashtable<String,QuantVO>> ofClass = onlyMS1DataPresent_.get(className);
            for (String analyteName : ofClass.keySet()){
              Hashtable<String,QuantVO> ofAnalyte = ofClass.get(analyteName);
              for (String modName : ofAnalyte.keySet()){
                if (predictionFound.containsKey(className) && predictionFound.get(className).containsKey(modName) && predictionFound.get(className).get(modName)) continue;
                quantStatus_.get(className).get(analyteName).put(modName, STATUS_WAITING);
                countToProcess++;
                ofAnalyte.remove(modName);
                Hashtable<String,Hashtable<String,String>> classToRemove = new Hashtable<String,Hashtable<String,String>>();
                if (toRemove.containsKey(className)) classToRemove = toRemove.get(className);
                Hashtable<String,String> analyteToRemove = new Hashtable<String,String>();
                if (classToRemove.containsKey(analyteName)) analyteToRemove = classToRemove.get(analyteName);
                analyteToRemove.put(modName, modName);
                classToRemove.put(analyteName, analyteToRemove);
                toRemove.put(className, classToRemove);
              }
            }
          }
          // remove hits that are already marked for quantification
          for (String className : toRemove.keySet()){
            for (String analyteName : toRemove.get(className).keySet()){
              for (String modName : toRemove.get(className).get(analyteName).keySet()){
                onlyMS1DataPresent_.get(className).get(analyteName).remove(modName);
              }
              if (onlyMS1DataPresent_.get(className).get(analyteName).size()==0) onlyMS1DataPresent_.get(className).remove(analyteName);
            }
            if (onlyMS1DataPresent_.get(className).size()==0) onlyMS1DataPresent_.remove(className);
          }
        }
        catch (Exception e) {
          e.printStackTrace();
          new WarningMessage(new JFrame(), "Warning","Post processing could not be executed because of: "+e.getMessage()+"!");
        }
      }
      //TODO: this is only here to prevent an endless loop
      if (countToProcess>0) stopThread = false;
    }
    if (stopThread){
      if (!error){
        //check here for results that need other adducts to be correct
        if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE) {
          try {
            results_ = OtherAdductChecker.checkTheResultsForOtherAdducts(results_,unsplittedPeaks_,quantObjects,analyzers_.get(0),classSequence);
          }
          catch (CgException | LipidCombinameEncodingException e) {
            e.printStackTrace();
          }
        }
        if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE && LipidomicsConstants.isMS2()){
          PostQuantificationProcessor processor = new PostQuantificationProcessor(results_,ms2Removed_,adductInsensitiveRtFilter);
          try {
            results_ = processor.chooseMoreLikelyOne(quantObjects);
          }
          catch (Exception e) {
            e.printStackTrace();
            new WarningMessage(new JFrame(), "Warning","Choose more likely one could not be executed because of: "+e.getMessage()+"!");
          }
          try {
            results_ = processor.processData();
          }
          catch (Exception e) {
            e.printStackTrace();
            new WarningMessage(new JFrame(), "Warning","Post processing could not be executed because of: "+e.getMessage()+"!");
          }
          results_ = reuniteWronglySeparatedPeaks(results_,quantObjects,unsplittedPeaks_);
        }
        //so that the checkTheResultsForOtherAdducts is really 100% correct, remove all hits that fall below the
        //base peak cutoff, and execute OtherAdductChecker.checkTheResultsForOtherAdducts again
        if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE) {
          try {
            OtherAdductChecker.removePeaksThatFallBelowTheBasepeakCutoff(results_,extractHighestArea()*(basePeakCutoff/1000f));
            results_ = OtherAdductChecker.checkTheResultsForOtherAdducts(results_,unsplittedPeaks_,quantObjects,analyzers_.get(0),classSequence);
          }
          catch (CgException | LipidCombinameEncodingException e) {
            e.printStackTrace();
          }
        }
        executeFinalProcesses(classSequence,analyteSequence,quantObjects,basePeakCutoff,resultFile,chromFile);
      }
      finished_ = true;
      for (Integer analyzer : analyzers_.keySet()){
        if (analyzers_.get(analyzer).getUseCuda()){
          analyzers_.get(analyzer).getSavGolJNI().Frees();
        }
      }
    }
  }
  
  /**
   * this looks if all except one partner were removed from a splitted peak instance -  if so, the original peak is restored and percental splits are removed
   * @param results the results of the quantitation procedure
   * @param quantObjects the objects defining the quantitation parameters
   * @param unsplittedPeaks_ the unsplitted peak instances - if present
   * @return
   */
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> reuniteWronglySeparatedPeaks(Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>results,
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects, Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> unsplittedPeaks_){
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> result = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
    for (String className : results.keySet()){
      Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> classResult = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
      Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = quantObjects.get(className);
//      if (unsplittedPeaks_.containsKey(className)){
      Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> oldClassResult = results.get(className);
//      Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> unsplittedOfClass = unsplittedPeaks_.get(className);
      for (String analyte : oldClassResult.keySet()){
        Hashtable<String,Hashtable<String,LipidParameterSet>> analyteResults = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
        Hashtable<String,QuantVO> quantsOfAnalyte = quantsOfClass.get(analyte);
//        if (unsplittedOfClass.containsKey(analyte)){
        Hashtable<String,Hashtable<String,LipidParameterSet>> oldAnalyteResults = oldClassResult.get(analyte);
//        Hashtable<String,Hashtable<String,LipidParameterSet>> unsplittedOfAnalyte = unsplittedOfClass.get(analyte);
        for (String mod : oldAnalyteResults.keySet()){
          Hashtable<String,LipidParameterSet> modResults = new Hashtable<String,LipidParameterSet>();
          QuantVO quant = quantsOfAnalyte.get(mod);
          Vector<QuantVO> quants = new Vector<QuantVO>();
          quants.add(quant);
          quants.addAll(quant.getOtherIsobaricSpecies());
//          if (unsplittedOfAnalyte.containsKey(mod) && quants.size()>1){
          Hashtable<String,LipidParameterSet> oldModResults = oldAnalyteResults.get(mod);
          Hashtable<String,LipidParameterSet> unsplittedMod = new Hashtable<String,LipidParameterSet>();
          if (unsplittedPeaks_.containsKey(className) && unsplittedPeaks_.get(className).containsKey(analyte) &&
              unsplittedPeaks_.get(className).get(analyte).containsKey(mod))
            unsplittedMod = unsplittedPeaks_.get(className).get(analyte).get(mod);
          for (String rt: oldModResults.keySet()){
            LipidParameterSet set = oldModResults.get(rt);
            if (!splittingPartnerPresent(quants,className,analyte,mod,rt,results)){
//            if (unsplittedMod.containsKey(rt) && !splittingPartnerPresent(quants,className,analyte,mod,rt,results)){
              // this is for MS1 splits
              if (unsplittedMod.containsKey(rt)){
                set = unsplittedMod.get(rt);
              }else{
                set.setPercentalSplit(-1f);
              }
            }
            modResults.put(rt, set);
            
          }
//              }else{
//                modResults = oldAnalyteResults.get(mod);
//              }
            if (modResults.size()>0) analyteResults.put(mod, modResults);
            }
//          }else{
//            analyteResults = oldClassResult.get(analyte);
//          }
          if (analyteResults.size()>0) classResult.put(analyte, analyteResults);
        }
//      }else{
//        classResult = result.get(className);
//      }
      result.put(className, classResult);
    }
    return result;
  }
  
  /**
   * checks for a given rt time of a certain identified peak, if there is a splitting partner present
   * if there is no more partner present, the split will be removed
   * @param quants the potential splitting partners (the objects defes the quantitation parameters)
   * @param className the lipid class
   * @param analyte the analyte name
   * @param mod the adduct of the analyte
   * @param rt the retention time of the peak
   * @param result the results of the quantitation procedure 
   * @return true if there are the splitting partners still present
   */
  private boolean splittingPartnerPresent(Vector<QuantVO> quants, String className, String analyte, String mod,
      String rt, Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> result){
    boolean partnerThere = false;
    for (QuantVO quant: quants){
      String quantAnalName = quant.getIdString();
      if (quant.getAnalyteClass().equalsIgnoreCase(className)&&quantAnalName.equalsIgnoreCase(analyte) &&
          quant.getModName().equalsIgnoreCase(mod)) continue;
//      System.out.println("other");
//      System.out.println(quant.getAnalyteClass()+";"+result.containsKey(quant.getAnalyteClass()));
//      if (result.containsKey(quant.getAnalyteClass())) System.out.println(1);
//      if (result.get(quant.getAnalyteClass()).containsKey(quantAnalName)) System.out.println(2);
//      if (result.get(quant.getAnalyteClass()).get(quantAnalName).containsKey(quant.getModName())) System.out.println(3);
//      if (result.get(quant.getAnalyteClass()).get(quantAnalName).get(quant.getModName()).containsKey(rt)) System.out.println(4);
      if (result.containsKey(quant.getAnalyteClass()) && result.get(quant.getAnalyteClass()).containsKey(quantAnalName)&&
          result.get(quant.getAnalyteClass()).get(quantAnalName).containsKey(quant.getModName())&&
          result.get(quant.getAnalyteClass()).get(quantAnalName).get(quant.getModName()).containsKey(rt)){
        partnerThere = true;
        break;
      }
    }
    return partnerThere;
  }
  
  @SuppressWarnings("unchecked")
  private void executeFinalProcesses(LinkedHashMap<String,Integer> classSequence,Hashtable<String,Vector<String>> analyteSequence,
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects, float basePeakCutoff, String resultFile, String chromFile){
    Hashtable<String,Vector<LipidParameterSet>> sheetParams = new Hashtable<String,Vector<LipidParameterSet>>();
    for (String className : classSequence.keySet()){
      Vector<LipidParameterSet> params = new Vector<LipidParameterSet>();
      Hashtable<String,Hashtable<String,QuantVO>> classQuant = quantObjects.get(className);
      for (String analyteName : analyteSequence.get(className)){
        Hashtable<String,QuantVO> analyteQuant = classQuant.get(analyteName);
        for (String mod : analyteQuant.keySet()){
          if (results_.containsKey(className) && results_.get(className).containsKey(analyteName) && results_.get(className).get(analyteName).containsKey(mod)&&(results_.get(className).get(analyteName).get(mod)!=null)){
            Hashtable<String,LipidParameterSet> hitsOfOneMod = results_.get(className).get(analyteName).get(mod);
            if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE){
              for (String key:hitsOfOneMod.keySet()) params.add(hitsOfOneMod.get(key));
            }else{
              List<DoubleStringVO> keys = new ArrayList<DoubleStringVO>();
              for (String key:hitsOfOneMod.keySet())keys.add(new DoubleStringVO(key,Double.valueOf(key)));
              Collections.sort(keys,new GeneralComparator("at.tugraz.genome.lda.vos.DoubleStringVO", "getValue", "java.lang.Double"));
              for (DoubleStringVO key : keys){
                params.add(hitsOfOneMod.get(key.getKey()));
              }
            }
          }  
        }
      }
      sheetParams.put(className, params);
    }  
    
    
    float highestArea = extractHighestArea();
    float cutOffValue = (basePeakCutoff/1000f)*highestArea;
    Hashtable<String,Vector<LipidParameterSet>> correctedParams = new Hashtable<String,Vector<LipidParameterSet>>();
    for (String sheetName : sheetParams.keySet()){
      Vector<LipidParameterSet> params = sheetParams.get(sheetName);
      Vector<LipidParameterSet> corrected = new Vector<LipidParameterSet>();
      for (LipidParameterSet param : params){
        boolean acceptProbe = false;
        if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE){
          acceptProbe = true;
        }else{
          if (param.getIsotopicProbes().size()>0){
            Vector<CgProbe> newIsoProbes = new Vector<CgProbe>();
            Vector<CgProbe> removedProbes = new Vector<CgProbe>();
            for (CgProbe probe : param.getIsotopicProbes().get(0)){
              if (probe.Area>cutOffValue)
                newIsoProbes.add(probe);
              else
                removedProbes.add(probe);
            }
            if (newIsoProbes.size()==param.getIsotopicProbes().get(0).size()){
              acceptProbe = true;
            }else if (removedProbes.size()==param.getIsotopicProbes().get(0).size()){
            }else{
              Vector<Vector<CgProbe>> isotopicProbes = new Vector<Vector<CgProbe>>();
              isotopicProbes.add(newIsoProbes);
              if (param.getIsotopicProbes().size()>1){
                for (int i=1;i!=param.getIsotopicProbes().size();i++){
                  Vector<CgProbe> newProbes = new Vector<CgProbe>();
                  for (CgProbe probe : param.getIsotopicProbes().get(i)){
                    if (!isWithinRemovedBoundaries(probe,removedProbes))
                      newProbes.add(probe);
                    else {
                      if (isWithinRemovedBoundaries(probe,newIsoProbes))
                        newProbes.add(probe);
                    }
                  }
                  isotopicProbes.add(newProbes);
                }  
              }
              param.setIsotopicProbes(isotopicProbes);
              acceptProbe = true;
            }
          }
        }
        if (acceptProbe)
          corrected.add(param);
      }
      correctedParams.put(sheetName, corrected);
    }
    System.out.println("Required time: "+((System.currentTimeMillis()-startCalcTime_)/(60*1000))+" minutes "+(System.currentTimeMillis()-startCalcTime_)%(60*1000)/1000+" seconds");
    try {
      LipidomicsConstants constants = LipidomicsConstants.getInstance();
      constants.setRawFileName(chromFile);
      boolean isAlexTargetList = false;
      Hashtable<String,Boolean> alexTargetlistUsed = new Hashtable<String,Boolean>();
      if (Settings.useAlex()){
        for (String className : quantObjects.keySet()){
          Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = quantObjects.get(className);
          for (Hashtable<String,QuantVO> quantsOfAnalyte : quantsOfClass.values()){
            for (QuantVO quant : quantsOfAnalyte.values()){
              if (quant instanceof TargetlistEntry){
                isAlexTargetList = true;
                if (((TargetlistEntry)quant).hasAlex123FragmentsForClass() && !alexTargetlistUsed.containsKey(className)){
                  alexTargetlistUsed.put(className, true);
                }
              }
            }
            if (isAlexTargetList) break;
          }
        }
      }
      constants.setAlexTargetlist(isAlexTargetList);
      constants.setAlexTargetlistUsed(alexTargetlistUsed);
      HydroxyEncoding[] encodings = getOnlyUsedHydroxyEncodings(correctedParams,Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding());
      QuantificationResult quantRes = new QuantificationResult(correctedParams,constants,classSequence,encodings[0],encodings[1]);
     
      QuantificationThread.writeResultsToExcel(resultFile,quantRes);
      
      if (isAlexTargetList || cli_){
        String alexResultFile = new String(resultFile);
        if (alexResultFile.endsWith(".xls") || alexResultFile.endsWith(".xlsx"))
          alexResultFile = alexResultFile.substring(0,alexResultFile.lastIndexOf("."));
        alexResultFile += ".tab";
        Vector<QuantificationResult> results = new Vector<QuantificationResult>();
        results.add(quantRes);
        RdbOutputWriter rdbWriter = new RdbOutputWriter();
        rdbWriter.write(alexResultFile, results, classSequence, analyteSequence, quantObjects);
      }
    }
    catch (Exception e) {
      e.printStackTrace();
      this.errorString_ = e.toString();
    }
  }
    
  private float extractHighestArea(){
    float highestArea = 0;
    for (Integer key : analyzers_.keySet()){
      LipidomicsAnalyzer analyzer = analyzers_.get(key);
      if (analyzer.getHighestArea()>highestArea) highestArea = analyzer.getHighestArea();
    }
    return highestArea;
  }
    
  /**
   * checks if the quantitation file contains any isobaric species - if so, this information is stored in the QuantVO
   * @param classSequence the lipid classes in sequence
   * @param analyteSequence the analytes in the lipid classes in sequence
   * @param quantObjects the hash containing the the QuantVO, which pertai information about quantitation
   * @throws RulesException exception thrown if something is wrong with an fragmentation rule
   * @throws IOException if a file is not there
   * @throws SpectrummillParserException if there is something wrong with the elementconfig.xml
   */
  private static void checkForIsobaricSpecies(LinkedHashMap<String,Integer> classSequence, Hashtable<String,Vector<String>> analyteSequence,
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects) throws RulesException, IOException, SpectrummillParserException{
    Vector<String> classes = new Vector<String>(classSequence.keySet());

    // first, the analytes are put in 10Da bins
    Hashtable<Integer,Vector<QuantVO>> tenDaIsobaricClusters = new Hashtable<Integer,Vector<QuantVO>>();
    for (int i=0; i!=classes.size(); i++){
      String lipidClass = classes.get(i);
      Vector<String> analytes = analyteSequence.get(lipidClass);
      Hashtable<String,Hashtable<String,QuantVO>> quantClass = quantObjects.get(lipidClass);
      for (int j=0; j!=analytes.size(); j++){
        String anal = analytes.get(j);
        Hashtable<String,QuantVO> quantAnal = quantClass.get(anal);
        for (String mod : quantAnal.keySet()){
          QuantVO quant = quantAnal.get(mod);
          float mz = (float)quant.getAnalyteMass();
          int clusterId = getTenDaClusterId(mz);
          Vector<QuantVO> vosOfCluster = new Vector<QuantVO>();
          if (tenDaIsobaricClusters.containsKey(clusterId)) vosOfCluster = tenDaIsobaricClusters.get(clusterId);
          vosOfCluster.add(quant);
          tenDaIsobaricClusters.put(clusterId, vosOfCluster);
        }
      }
    }
    Hashtable<String,String> alreadyUsed = new Hashtable<String,String>(); 
    //second the algorithm looks for analytes which are within the tolerance in the 10Da bins
    for (int i=0; i!=classes.size(); i++){
      String lipidClass = classes.get(i);
      Vector<String> analytes = analyteSequence.get(lipidClass);
      Hashtable<String,Hashtable<String,QuantVO>> quantClass = quantObjects.get(lipidClass);
      for (int j=0; j!=analytes.size(); j++){
        String anal = analytes.get(j);
        Hashtable<String,QuantVO> quantAnal = quantClass.get(anal);
        for (String mod : quantAnal.keySet()){
          QuantVO quant1 = quantAnal.get(mod);
          float tol = LipidomicsConstants.getCoarseChromMzTolerance((float)quant1.getAnalyteMass());
          float sameTol = tol/5f;
          boolean noMS21 = false;
          try{ RulesContainer.getAmountOfChains(StaticUtils.getRuleName(quant1.getAnalyteClass(), quant1.getModName())); } catch (NoRuleException nrx){
            noMS21 = true;
          }
          if (noMS21 || quant1.isQuantifiedByOtherIsobar()) continue;
          float lowerMz = (float)quant1.getAnalyteMass()-tol;
          float upperMz = (float)quant1.getAnalyteMass()+tol;
          float lowerSame = (float)quant1.getAnalyteMass()-sameTol;
          float upperSame = (float)quant1.getAnalyteMass()+sameTol;
          int clusterId = getTenDaClusterId(lowerMz);
          Vector<QuantVO> toCompare = new Vector<QuantVO>();
          if (tenDaIsobaricClusters.containsKey(clusterId)) toCompare.addAll(tenDaIsobaricClusters.get(clusterId));
          if (clusterId!=getTenDaClusterId(upperMz)){
            clusterId = getTenDaClusterId(upperMz);
            if (tenDaIsobaricClusters.containsKey(clusterId)) toCompare.addAll(tenDaIsobaricClusters.get(clusterId));
          }
          for (QuantVO quant2 : toCompare){
            if (quant1.equals(quant2) || alreadyUsed.containsKey(getUniqueQuantVOString(quant2))) continue;
            boolean noMS22 = false;
            try{ RulesContainer.getAmountOfChains(StaticUtils.getRuleName(quant2.getAnalyteClass(), quant2.getModName()));}catch (NoRuleException nrx){
              noMS22 = true;
            }
            if (noMS22) continue;
            float mz2 = (float)quant2.getAnalyteMass();
            boolean similarMass = lowerMz<mz2 && mz2<upperMz;
            boolean sameMass = lowerSame<mz2 && mz2<upperSame;
            if (similarMass){
              quant1.addIsobaricSpecies(quant2);
              quant2.addIsobaricSpecies(quant1);
//              if ((quant1.getAnalyteClass().equalsIgnoreCase("P-PE")&&quant1.getIdString().equalsIgnoreCase("36:4"))||
//                  (quant2.getAnalyteClass().equalsIgnoreCase("P-PE")&&quant2.getIdString().equalsIgnoreCase("36:4")))
//                System.out.println(quant1.getAnalyteClass()+quant1.getIdString()+"_"+quant1.getModName()+"/"+quant2.getAnalyteClass()+quant2.getIdString()+"_"+quant2.getModName()+" "+quant1.getAnalyteMass()+"/"+mz2);
              if (sameMass)
                quant2.setQuantifiedByOtherIsobar(true);
            }
          }
          String idString = getUniqueQuantVOString(quant1);
          alreadyUsed.put(idString, idString);
        }
      }
    }      
  }
  
  /**
   * this string is required as unique identifiers, to avoid a repeated comparison
   * @return a unique string separateble from ohter analytes
   */
  private static String getUniqueQuantVOString(QuantVO vo){
    return vo.getAnalyteClass()+"_;%"+vo.getIdString()+"_;%"+vo.getModName();
  }
  
  /**
   * 
   * @param mz mz value
   * @return a unique integer id for a cluster of the size of 10Da
   */
  private static int getTenDaClusterId(float mz){
    return Math.round(Calculator.roundFloat(mz, 0, BigDecimal.ROUND_DOWN));
  }
  
  /**
   * parses an Alex123 target list file or a directory containing Alex123 target list files
   * @param quantFile the Alex123 file or the directory containing the Alex123 files
   * @param minusTime time tolerance in negative direction
   * @param plusTime time tolerance in positive direction
   * @param amountOfIsotopes the number of isotopes to be quantified
   * @param isotopesMustMatch the number of isotopes that have to matcht the theoretical isotopic distribution 
   * @param searchUnknownTime true if analytes have to be searched where no retention time is present
   * @param basePeakCutoff a relative cutoff value
   * @param rtShift shift of retention times in relation to the entered ones
   * @param lowestRetTime lower hard limit for the retention time
   * @param highestRetTime upper hard limit for the retention time
   * @param positiveIonMode should the targets ion positive or in negative ion mode be used for quantitation
   * @return a vector containing the targets for quantitation
   * @throws IOException if a file is not there
   * @throws SpectrummillParserException if there is something wrong with the elementconfig.xml
   * @throws AlexTargetlistParserException if there is something wrong with the target lists
   * @throws ChemicalFormulaException if there is something wrong with the chemical formulae
   * @throws RulesException if there is somehting wrong for the rule parsing
   */
  @SuppressWarnings({ "rawtypes", "unchecked" })
  private static Vector parseAlex123TargetList(String quantFile, float minusTime, float plusTime, int amountOfIsotopes, int isotopesMustMatch, boolean searchUnknownTime, float basePeakCutoff,
      float rtShift, float lowestRetTime, float highestRetTime, boolean positiveIonMode) throws IOException,SpectrummillParserException, AlexTargetlistParserException, ChemicalFormulaException, RulesException{
    File file = new File(quantFile);
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>>> sortedEntries = null;
    if (file.isFile() && file.getName().endsWith(".txt")){
      Vector<Hashtable<Integer,Vector<TargetlistEntry>>> parsedEntries = new Vector<Hashtable<Integer,Vector<TargetlistEntry>>>();
      TargetlistParser parser = new TargetlistParser(quantFile,positiveIonMode);
      parser.parse();
      parsedEntries.add(parser.getResults());
      sortedEntries = TargetlistDirParser.sortEntriesForLDA(parsedEntries);
    }else if (file.isDirectory()){
      TargetlistDirParser dirParser = new TargetlistDirParser(quantFile,positiveIonMode);
      dirParser.parse();
      sortedEntries = dirParser.getResults();
    }
    if (sortedEntries==null || sortedEntries.size()==0)
      throw new AlexTargetlistParserException("There are unusable entries in your target list");
    
    LinkedHashMap<String,Integer> classSequence = new LinkedHashMap<String,Integer>();
    Hashtable<String,Boolean> adductInsensitiveRtFilter = new Hashtable<String,Boolean>();
    Hashtable<String,Vector<String>> analyteSequence = new Hashtable<String,Vector<String>>();
    Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();

    //now generate the corresponding objects
    ElementConfigParser elementParser = Settings.getElementParser();
    for (String className : sortedEntries.keySet()){
      classSequence.put(className, 1);
      LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>> classEntries = sortedEntries.get(className);
      Vector<String> analytes = new Vector<String>();
      Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = new Hashtable<String,Hashtable<String,QuantVO>>();
      for (String analyteOriginalName : classEntries.keySet()){
        LinkedHashMap<String,TargetlistEntry> analyteEntries = classEntries.get(analyteOriginalName);
        String sideChain = "";
        int doubleBonds = -1;
        Hashtable<String,QuantVO> quantsOfAnalyte = new Hashtable<String,QuantVO>();
        for (String mod : analyteEntries.keySet()){
          TargetlistEntry entry = analyteEntries.get(mod);
          sideChain = entry.getAnalyteName();
          doubleBonds = entry.getDbs();
          entry.setTimeConstraints(-1, minusTime, plusTime);
          
          Hashtable<String,Integer> formAnal = StaticUtils.categorizeFormula(entry.getAnalyteFormula());
          Hashtable<String,Integer> formMod = StaticUtils.categorizeFormula(entry.getModFormula());
          String[] formulas = getFormulasAsString(formAnal,formMod,1);
//          String analyteFormula = formulas[0];
//          String modificationFormula = formulas[1];
          String chemicalFormula = formulas[2];
//          System.out.println(className+StaticUtils.generateLipidNameString(sideChain, doubleBonds)+": "+chemicalFormula);
          Object[] distris = getTheoreticalIsoDistributions(elementParser,isotopesMustMatch,amountOfIsotopes,chemicalFormula);
          Vector<Double> mustMatchProbabs = (Vector<Double>)distris[0];
          Vector<Double> probabs = (Vector<Double>)distris[1];
          int negativeStartValue = (Integer)distris[2];
          entry.setDistributionValues(mustMatchProbabs,probabs,negativeStartValue);
          
          
          quantsOfAnalyte.put(entry.getModName(), entry);
        }
        
        String analyteName = StaticUtils.generateLipidNameString(sideChain, doubleBonds);
        analytes.add(analyteName);
        quantsOfClass.put(analyteName, quantsOfAnalyte);
      }
      quantObjects.put(className, quantsOfClass);
      analyteSequence.put(className, analytes);
      adductInsensitiveRtFilter.put(className, false);
    }
    checkForIsobaricSpecies(classSequence,analyteSequence,quantObjects);
    Vector results = new Vector();
    results.add(classSequence);
    results.add(analyteSequence);
    results.add(adductInsensitiveRtFilter);
    results.add(quantObjects);
    //TODO: these few lines are only for testing purposes
//    classSequence = new LinkedHashMap<String,Integer>();
//    classSequence.put("TAG", 1);
//    analyteSequence = new Hashtable<String,Vector<String>>();
//    Vector<String> analytesOfClass = new Vector<String>();
//    analytesOfClass.add("40:0");
//    analyteSequence.put("TAG", analytesOfClass);
//    Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects2 = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();
//    Hashtable<String,Hashtable<String,QuantVO>> ofClass = new Hashtable<String,Hashtable<String,QuantVO>>();
//    ofClass.put("40:0", quantObjects.get("TAG").get("40:0"));   
//    quantObjects2.put("TAG", ofClass);
//    results.add(classSequence);
//    results.add(analyteSequence);
//    results.add(adductInsensitiveRtFilter);
//    results.add(quantObjects2);
    return results;
  }
  
  /**
   * computes chemical formula strings out of the analyte composition and the modification 
   * @param elementalComposition the elemental composition of the analyte
   * @param modElements the elemental composition of the modification
   * @param mult multiplication for dimers, trimers, etc.
   * @return [0] neutral analyte formula; [1] formula of modification; [2] total formula
   */
  private static String[] getFormulasAsString(Hashtable<String,Integer> elementalComposition, Hashtable<String,Integer> modElements,
      int mult){
    Hashtable<String,Integer> chemicalFormula = new Hashtable<String,Integer>(elementalComposition);
    for (String element : modElements.keySet()) {
      int amount = modElements.get(element);
      if (chemicalFormula.containsKey(element))
        amount += chemicalFormula.get(element);
      chemicalFormula.put(element, amount);
    }
    String[] formulas = new String[3];
    formulas[0] = StaticUtils.getFormulaInHillNotation(elementalComposition, true);
    formulas[1] = StaticUtils.getFormulaInHillNotation(modElements, true);
    formulas[2] = StaticUtils.getFormulaInHillNotation(chemicalFormula,true);
    return formulas;
  }
  
  /**
   * calculates the positive and the negative isotopic distribution, and decides which on has to be used
   * @param elementParser parser containint the relative abundances of the elements
   * @param isotopesMustMatch number of isotopes that have to fit the theoretical isotopic distribution
   * @param amountOfIsotopes number of isotopes that shall be quantified by the LDA algorithm
   * @param chemicalFormula the chemical formula of the analyte where the isotopic distribution has to match
   * @return [0] Vector<Double> containing the probabilities of the isotopes that must match; [1] Vector<Double> containing probabilities of all isotopes; [2] if the distribution goes in the negative direction - how negative is the lowest isotope 
   * @throws SpectrummillParserException if there is something wrong with the elementconfig.xml
   */
  private static Object[] getTheoreticalIsoDistributions(ElementConfigParser elementParser, int isotopesMustMatch, int amountOfIsotopes, String chemicalFormula) throws SpectrummillParserException{  
    boolean negativeDistribution = false;
    Vector<Double> probabs = new Vector<Double>();
    if (amountOfIsotopes<isotopesMustMatch)
      amountOfIsotopes = isotopesMustMatch;
    if (amountOfIsotopes>0){
      Vector<Vector<Double>> bothDistris = elementParser.calculateChemicalFormulaIntensityDistribution(chemicalFormula, amountOfIsotopes+1, false);
      probabs = bothDistris.get(0);
      if (bothDistris.size()>1){
        Vector<Double> negDistri = bothDistris.get(1);
        if (StaticUtils.useNegativeDistribution(probabs,negDistri)){
          probabs = negDistri;
          negativeDistribution = true;
        }
      }

    }else{
      probabs.add(1d);
    }
    Vector<Double> mustMatchProbabs = new Vector<Double>();
    if (isotopesMustMatch>0){
      if (amountOfIsotopes == isotopesMustMatch){
        mustMatchProbabs = new Vector<Double>(probabs);
      }else{
        Vector<Vector<Double>> bothDistris = elementParser.calculateChemicalFormulaIntensityDistribution(chemicalFormula, isotopesMustMatch+1, false);
        mustMatchProbabs = bothDistris.get(0);
        if (bothDistris.size()>1){
          Vector<Double> negDistri = bothDistris.get(1);
          if (StaticUtils.useNegativeDistribution(mustMatchProbabs,negDistri)){
            mustMatchProbabs = negDistri;
            negativeDistribution = true;
          }
        }
      }
    }
    int negativeStartValue = 0;
    if (negativeDistribution){
      negativeStartValue = (mustMatchProbabs.size()*-1)+1;
    }
    
    Object[] distris = new Object[3];
    distris[0] = mustMatchProbabs;
    distris[1] = probabs;
    distris[2] = negativeStartValue;
    return distris;
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
    if (columnWidth>255*256)
      columnWidth=255*256;
    sheet.setColumnWidth(column,columnWidth); 
  }
  
  
  /**
   * method that reduced the defined OH encodings to the actually used ones
   * @param correctedParams the found results
   * @param faHydroxyEncoding the actually used OH encodings for FA moiety
   * @param lcbHydroxyEncoding the actually used OH encodings for LCB moiety
   * @return the reduced OH encodings [0] = FA; [1] = LCB
   */
  private HydroxyEncoding[] getOnlyUsedHydroxyEncodings(Hashtable<String,Vector<LipidParameterSet>> correctedParams,
      HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding) {
    Hashtable<String,Short> usedFaEncodings = new Hashtable<String,Short>();
    Hashtable<String,Short> usedLcbEncodings = new Hashtable<String,Short>();
    LipidomicsMSnSet msn;
    Short oh;
    for (Vector<LipidParameterSet> params : correctedParams.values()) {
      for (LipidParameterSet param : params) {
        if (!(param instanceof LipidomicsMSnSet))
          continue;
        msn = (LipidomicsMSnSet)param;
        if (msn.getStatus()<LipidomicsMSnSet.FRAGMENTS_DETECTED)
          continue;
        for (FattyAcidVO fa : msn.getInvolvedFAs().values()) {
          if (fa.getOhNumber()<0)
            continue;
          oh = (short)fa.getOhNumber();
          if (fa.getChainType()<LipidomicsConstants.CHAIN_TYPE_LCB) {
            if (usedFaEncodings.contains(oh))
              continue;
            try {
              usedFaEncodings.put(faHydroxyEncoding.getEncodedPrefix(oh),oh);
            //this exception cannot happen
            }catch (HydroxylationEncodingException e) {}
          }else if (fa.getChainType()==LipidomicsConstants.CHAIN_TYPE_LCB) {
            if (usedLcbEncodings.contains(oh))
              continue;
            try {
              usedLcbEncodings.put(lcbHydroxyEncoding.getEncodedPrefix(oh),oh);
            //this exception cannot happen
            }catch (HydroxylationEncodingException e) {}            
          }
        }
      }
    }
    
    HydroxyEncoding[] encodings = new HydroxyEncoding[2];
    encodings[0] = new HydroxyEncoding(usedFaEncodings);
    encodings[1] = new HydroxyEncoding(usedLcbEncodings);
    return encodings;
  }
  
  /**
   * this reorganizes the hash table returned from the SingleQuantTrhead, so that the retention time is the highest grouping parameter
   * @param hitsAccordingToQuant hash table returned from the SingleQuantTrhead
   * @return the reorganized hash table
   */
  private Hashtable<String,Hashtable<QuantVO,LipidParameterSet>> groupHitsWithSameRt(Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant) {
    Hashtable<String,Hashtable<QuantVO,LipidParameterSet>> sameRt = new Hashtable<String,Hashtable<QuantVO,LipidParameterSet>>();
    for (QuantVO quant : hitsAccordingToQuant.keySet()) {
      for (String rt : hitsAccordingToQuant.get(quant).keySet()) {
        Hashtable<QuantVO,LipidParameterSet> hash = new Hashtable<QuantVO,LipidParameterSet>();
        if (sameRt.containsKey(rt)) hash = sameRt.get(rt);
        hash.put(quant, hitsAccordingToQuant.get(quant).get(rt));
        sameRt.put(rt, hash);
      }
    }
    return sameRt;
  }

}
