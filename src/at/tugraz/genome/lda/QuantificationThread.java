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

package at.tugraz.genome.lda;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.Timer;
import java.util.TimerTask;
import java.util.Vector;

import javax.swing.JFrame;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
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
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.export.QuantificationResultExporter;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.MSnAnalyzer;
import at.tugraz.genome.lda.msn.OtherAdductChecker;
import at.tugraz.genome.lda.msn.PostQuantificationProcessor;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.RtPredictVO;
import at.tugraz.genome.lda.parser.ModificationParser;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.lda.target.export.TargetListExporter;
import at.tugraz.genome.lda.utils.RangeInteger;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleStringVO;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.maspectras.utils.StringUtils;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
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
  
  private long startCalcTime_;
  
  private Timer timer_;
  
  /** in the case of MSnFirst: contains LM-Models and suggestions for the next range for quantitation*/
  private Hashtable<String,Hashtable<String,RtPredictVO>> latestRtPredictions_; 
  
  private final static int STATUS_WAITING = 0;
  private final static int STATUS_CALCULATING = 1;
  private final static int STATUS_FINISHED = 2;
  
    
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
  
  public String getChromFile() {
    return this.chromFile_;
  }
  public String getResultFile() {
    return this.resultFile_;
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
 // LL    Hashtable<String,Vector<String>> analyteSequence = (Hashtable<String,Vector<String>>)quantContent.get(1);
      LinkedHashMap<String,Vector<String>> analyteSequence = (LinkedHashMap<String,Vector<String>>)quantContent.get(1);
    
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
  
  public static void setAnalyzerProperties(LipidomicsAnalyzer analyzer){
    if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE){
      analyzer.setShotgunParameters(LipidomicsConstants.getShogunProcessing(),LipidomicsConstants.getMs2MzTolerance(),LipidomicsConstants.getMs2MzToleranceUnit());
    }else if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_PRM){
      analyzer.setPrmParameters(LipidomicsConstants.getChromSmoothRange(),LipidomicsConstants.getChromSmoothRepeats(),
          LipidomicsConstants.getPeakDiscardingAreaFactor(),LipidomicsConstants.getRelativeAreaCutoff(),
          LipidomicsConstants.getRelativeFarAreaCutoff(),LipidomicsConstants.getRelativeFarAreaTimeSpace(),
          LipidomicsConstants.getMs2MzTolerance(),LipidomicsConstants.getMs2MzToleranceUnit());
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
        LipidomicsConstants.getRelativeIsoInBetweenCutoff(),LipidomicsConstants.getClosePeakTimeTolerance(),
        LipidomicsConstants.getTwinInBetweenCutoff(), LipidomicsConstants.getUnionInBetweenCutoff(),
        LipidomicsConstants.getMs2MzTolerance(),LipidomicsConstants.getMs2MzToleranceUnit());
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
  
  /**
   * @param quantFile File path of the Quant File
   * @param minusTime Retention time before tolerance
   * @param plusTime Retention time after tolerance
   * @param amountOfIsotopes Number of isotopes that shall be quantified
   * @param isotopesMustMatch Number of isotopes that must match the theoretical distribution
   * @param searchUnknownTime Search unknown retention time
   * @param basePeakCutoff The relative cutoff value in per mille
   * @param rtShift Retention time shift
   * @param lowestRetTime Lowest retention time in the chrom files
   * @param highestRetTime Highest retention time in the chrom files
   * @param respectMassShift Take a mass shift range into account
   * @return a vector containing class sequence, analyte sequence, adduct insensitive retention time filter and quantVO objects
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws ExcelInputFileException
   * @throws ChemicalFormulaException
   * @throws RulesException
   * @throws HydroxylationEncodingException
   */
  @SuppressWarnings({ "unchecked", "rawtypes" })
  public static Vector parseQuantExcelFile(String quantFile, float minusTime, float plusTime, int amountOfIsotopes, int isotopesMustMatch, boolean searchUnknownTime, float basePeakCutoff,
      float rtShift, float lowestRetTime, float highestRetTime, boolean respectMassShift) throws IOException,SpectrummillParserException,ExcelInputFileException, ChemicalFormulaException, RulesException, HydroxylationEncodingException{
    InputStream myxls = new FileInputStream(quantFile);
    Workbook workbook = null;
    if (quantFile.endsWith(".xlsx")) workbook = new XSSFWorkbook(myxls);
    else if (quantFile.endsWith(".xls")) workbook = new HSSFWorkbook(myxls);
    Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();
    LinkedHashMap<String,Integer> classSequence = new LinkedHashMap<String,Integer>();
    Hashtable<String,Boolean> adductInsensitiveRtFilter = new Hashtable<String,Boolean>();
    Hashtable<String,Boolean> bestMatchBySpectrumCoverage = new Hashtable<String,Boolean>();
 // LL  Hashtable<String,Vector<String>> analyteSequence = new Hashtable<String,Vector<String>>();
    LinkedHashMap<String,Vector<String>> analyteSequence = new LinkedHashMap<String,Vector<String>>();
    boolean excelOK = false;
    ElementConfigParser aaParser = Settings.getElementParser();
    for (int sheetNumber = 0; sheetNumber!=workbook.getNumberOfSheets(); sheetNumber++){
      Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = new Hashtable<String,Hashtable<String,QuantVO>>();
      Hashtable<String,Hashtable<String,QuantVO>> quantsOfOxClass = new Hashtable<String,Hashtable<String,QuantVO>>();
      
      Vector<String> analytes = new Vector<String>();
      Vector<String> oxAnalytes = new Vector<String>();
      Sheet sheet = workbook.getSheetAt(sheetNumber);
      boolean rtFilterInsensitive = false;
      boolean pickBestMatchBySpectrumCoverage = false;
      int sideChainColumn = -1;
      int doubleBondColumn = -1;
      int molecularSpeciesWithDoubleBondPositionsColumn = -1;
      int oxStateColumn = -1;

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
        String oxState = "";
        Hashtable<String,Integer> elementalComposition = new Hashtable<String,Integer>();
        Hashtable <String,Double> massesOfInterest = new Hashtable <String,Double>();
        float retTime = -1;
        String molecularSpeciesWithDB = null;
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
            } 
            
            else if (contents.equalsIgnoreCase("dbs")||contents.equalsIgnoreCase("dbs_TAG")){
              doubleBondColumn = i;
            } else if (contents.equalsIgnoreCase(LipidomicsConstants.CHAIN_MOD_COLUMN_NAME)){
              oxStateColumn = i;
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
            } else if (contents.trim().equalsIgnoreCase("adductInsensitiveRtFilter")) {
            	rtFilterInsensitive = true;
            } else if (contents.startsWith("OH-Number:")){
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
            } else if (contents.trim().equalsIgnoreCase("pickBestMatchBySpectrumCoverage")) {
            	pickBestMatchBySpectrumCoverage = true;
            } else if (contents.equalsIgnoreCase(TargetListExporter.HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS)){
              molecularSpeciesWithDoubleBondPositionsColumn = i;
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
            if (i==oxStateColumn&&contents!=null&&contents.length()>0){
              oxState = contents;
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
            
            if (i==molecularSpeciesWithDoubleBondPositionsColumn&&contents!=null&contents.length()>0){
              molecularSpeciesWithDB = contents;
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
            
            //create new quantVO objects for each (ox)modified lipid
            String[] oxStates = oxState.split(";");
            if (oxState.length() == 0) {
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
                  Hashtable<String,Integer> modElements = adductComposition.get(modName);
                  Integer charge = charges.get(modName);
                  double massOfInterest = massesOfInterest.get(modName)+ohDiff/((double)charge);
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
                      mustMatchProbabs, probabs,negativeStartValue, "");
                  
                  analEncoded = quantVO.getAnalyteName();
                  quantsOfAnalyte.put(modName, quantVO);
                }
                String analyteName = StaticUtils.generateLipidNameString((analEncoded!=null ? analEncoded : sideChain), doubleBonds, -1, "");
                if (!quantsOfClass.containsKey(analyteName)) {
                  analytes.add(analyteName);
                  quantsOfClass.put(analyteName, quantsOfAnalyte);
                } else if (!hasValidInfoForOmegaAssignment(molecularSpeciesWithDB, retTime)) {
                  System.out.println(String.format("Ignoring duplicate analyte %s in mass list (line %s)!", analyteName, (rowCount+1)));
                }
                if (hasValidInfoForOmegaAssignment(molecularSpeciesWithDB, retTime)) {
                  for (String mod : quantsOfAnalyte.keySet()){
                    try {
                      Vector<FattyAcidVO> chainCombination = StaticUtils.decodeFAsFromHumanReadableName(
                          molecularSpeciesWithDB, Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), false, null);
                      DoubleBondPositionVO doubleBondPositionVO = new DoubleBondPositionVO(
                          chainCombination, retTime);
                      quantsOfClass.get(analyteName).get(mod).addInfoForOmegaAssignment(doubleBondPositionVO);
                    } catch (LipidCombinameEncodingException ex) {
                      System.out.println(ex.getMessage());
                    }
                  }
                }
              }
            	
            } else {
            for (String oxMod : oxStates)
            {
            	oxMod = oxMod.replaceAll("\\s", "");
            	
	            for (int oh=startOh; oh<(stopOh+1); oh++) {
	              Hashtable<String,QuantVO> quantsOfAnalyte = new Hashtable<String,QuantVO>();
	              Hashtable<String,QuantVO> quantsOfOxAnalyte = new Hashtable<String,QuantVO>();
	              
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
	              
	              ModificationParser mp = new ModificationParser(oxMod);
	              mp.parse();
	              correctedElementalComposition = mp.getNewChemicalComposition(correctedElementalComposition);
	              
	              for (String modName : massesOfInterest.keySet()){
	                
	                String analyteClass = sheet.getSheetName();
	                
	                
	                //apply respective (ox)modification to all masses and formulas

                    Hashtable<String,Integer> modElements = adductComposition.get(modName);
                    Integer charge = charges.get(modName);
                    double massOfInterest = mp.getNewMass(massesOfInterest.get(modName)+ohDiff/((double)charge));
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
                	
                	//add ox-prefix to modified lipids
                	if(!oxMod.equals("")) {
                    	analyteClass = "ox" + analyteClass;
                    }
                	
                	QuantVO quantVO = new QuantVO(analyteClass, sideChain, doubleBonds,
                            ohToUse,analyteFormula, massOfInterest, charge, modName,
                            modificationFormula, retTime, usedMinusTime, usedPlusTime,
                            mustMatchProbabs, probabs,negativeStartValue,oxMod);
                            analEncoded = quantVO.getAnalyteName();
                            
                    /*add oxLipids to another sheet (effectively becoming a tab in the GUI)*/
                            //==
                    if(oxMod.equals("")) {
                    	quantsOfAnalyte.put(modName, quantVO);
                    }else {
                    	quantsOfOxAnalyte.put(modName, quantVO);
                    }
	                	
	              }
	              String analyteName = StaticUtils.generateLipidNameString((analEncoded!=null ? analEncoded : sideChain), doubleBonds, -1, oxMod);
	              
	              /*add oxLipids to another sheet (effectively becoming a tab in the GUI)*/
	              if(oxMod.equals("")) {
	            	  analytes.add(analyteName);
	            	  quantsOfClass.put(analyteName, quantsOfAnalyte);
	              }else {
	            	  oxAnalytes.add(analyteName);
	            	  quantsOfOxClass.put(analyteName, quantsOfOxAnalyte);
	              }              
	            }
            }
          	}
          }
        }
      } 
      quantObjects.put(sheet.getSheetName(), quantsOfClass);
      analyteSequence.put(sheet.getSheetName(), analytes);
      classSequence.put(sheet.getSheetName(),msLevel);
      adductInsensitiveRtFilter.put(sheet.getSheetName(), rtFilterInsensitive);
      bestMatchBySpectrumCoverage.put(sheet.getSheetName(), pickBestMatchBySpectrumCoverage);
      
      /*add (if applicable) oxLipids to another sheet (effectively becoming a tab in the GUI)*/
      if(!oxAnalytes.isEmpty())
      {
    	  quantObjects.put("ox"+sheet.getSheetName(), quantsOfOxClass);
          analyteSequence.put("ox"+sheet.getSheetName(), oxAnalytes);
          classSequence.put("ox"+sheet.getSheetName(),msLevel);
          adductInsensitiveRtFilter.put("ox"+sheet.getSheetName(), rtFilterInsensitive);
          bestMatchBySpectrumCoverage.put("ox"+sheet.getSheetName(), pickBestMatchBySpectrumCoverage);
      }
      
      if (foundColumns) excelOK = true;
    }
    myxls.close();
    if (!excelOK) throw new ExcelInputFileException("The Excel file is not valid!");
    checkForIsobaricSpecies(classSequence,analyteSequence,quantObjects);
    
    Vector results = new Vector();
    results.add(classSequence);
    results.add(analyteSequence);
    results.add(adductInsensitiveRtFilter);
    results.add(bestMatchBySpectrumCoverage);
    results.add(quantObjects);
    return results;
  }
  
  /**
   * @param molecularSpeciesWithDB String of the potential molecular species with double bond position assignment
   * @param retTime retention time in minutes
   * @return true if a molecular species with double bond position information and a corresponding retention time was found in the mass list
   */
  private static boolean hasValidInfoForOmegaAssignment(String molecularSpeciesWithDB, float retTime)
  {
    if (molecularSpeciesWithDB != null && retTime>0) return true;
    return false;
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
 // LL  private Hashtable<String,Vector<String>> analyteSequence_;
    private LinkedHashMap<String,Vector<String>> analyteSequence_;
    private Hashtable<String,Boolean> adductInsensitiveRtFilter_;
    private Hashtable<String,Boolean> bestMatchBySpectrumCoverage_;
    private Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects_;
    private float bpCutoff_;
    private String rsFile_;

    
    @SuppressWarnings({ "rawtypes", "unchecked" })
    protected ThreadSupervisor(Vector excelContent, float  basePeakCutoff, String resultFile){
      currentLipidCount_ = 0;
      classSequence_ = (LinkedHashMap<String,Integer>)excelContent.get(0);
   // LL  analyteSequence_ = (Hashtable<String,Vector<String>>)excelContent.get(1);
      analyteSequence_ = (LinkedHashMap<String,Vector<String>>)excelContent.get(1);
      adductInsensitiveRtFilter_ = (Hashtable<String,Boolean>)excelContent.get(2);
      bestMatchBySpectrumCoverage_ = (Hashtable<String,Boolean>)excelContent.get(3);
      quantObjects_ = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)excelContent.get(4);
      bpCutoff_ = basePeakCutoff;
      rsFile_ = resultFile;
      initThreadHashes();
    }

    public void run()
    { 
      try{
        if (!finished_)
          handleTimerEvent(classSequence_,analyteSequence_,adductInsensitiveRtFilter_,bestMatchBySpectrumCoverage_,quantObjects_,bpCutoff_,rsFile_,chromFileName_);
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
   // LL    for (String analyteName : analytes){
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
  
  
  private void handleTimerEvent(LinkedHashMap<String,Integer> classSequence,LinkedHashMap<String,Vector<String>> analyteSequence,
	      Hashtable<String,Boolean> adductInsensitiveRtFilter, Hashtable<String,Boolean> bestMatchBySpectrumCoverage, Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects,
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
                
                set.setOxState(oneQuant.getOxState());           
                
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
                Collections.sort(quants); /// LL
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
/*LL*/        PostQuantificationProcessor processor = new PostQuantificationProcessor(results_,ms2Removed_,adductInsensitiveRtFilter);
        try {
/*LL*/          latestRtPredictions_ = processor.predictRetentionTimesBasedOnResults(onlyMS1DataPresent_,latestRtPredictions_);
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
 /*LL*/         try {
 /*LL*/           results_ = OtherAdductChecker.checkTheResultsForOtherAdducts(results_,unsplittedPeaks_,quantObjects,analyzers_.get(0),classSequence);
 /*LL*/         }
 /*LL*/         catch (CgException | LipidCombinameEncodingException e) {
 /*LL*/           e.printStackTrace();
  /*LL*/        }
        }
        if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE && LipidomicsConstants.isMS2()){
 /*LL*/         PostQuantificationProcessor processor = new PostQuantificationProcessor(results_,ms2Removed_,adductInsensitiveRtFilter);
          try {
 /*LL*/           results_ = processor.chooseMoreLikelyOne(quantObjects);
          }
          catch (Exception e) {
            e.printStackTrace();
            new WarningMessage(new JFrame(), "Warning","Choose more likely one could not be executed because of: "+e.getMessage()+"!");
          }
          try {
/*LL*/            results_ = processor.processData();
          }
          catch (Exception e) {
            e.printStackTrace();
            new WarningMessage(new JFrame(), "Warning","Post processing could not be executed because of: "+e.getMessage()+"!");
          }
/*LL*/          results_ = reuniteWronglySeparatedPeaks(results_,quantObjects,unsplittedPeaks_);
        }
        //so that the checkTheResultsForOtherAdducts is really 100% correct, remove all hits that fall below the
        //base peak cutoff, and execute OtherAdductChecker.checkTheResultsForOtherAdducts again
        if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE) {
  /*LL*/        try {
  /*LL*/          OtherAdductChecker.removePeaksThatFallBelowTheBasepeakCutoff(results_,extractHighestArea()*(basePeakCutoff/1000f));
    /*LL*/        results_ = OtherAdductChecker.checkTheResultsForOtherAdducts(results_,unsplittedPeaks_,quantObjects,analyzers_.get(0),classSequence);
   /*LL*/       }
 /*LL*/         catch (CgException | LipidCombinameEncodingException e) {
  /*LL*/          e.printStackTrace();
  /*LL*/        }
        }
        executeFinalProcesses(classSequence,analyteSequence,quantObjects,basePeakCutoff,resultFile,chromFile,bestMatchBySpectrumCoverage);
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
  private void executeFinalProcesses(LinkedHashMap<String,Integer> classSequence,LinkedHashMap<String,Vector<String>> analyteSequence,
	      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects, float basePeakCutoff, String resultFile, String chromFile, Hashtable<String,Boolean> bestMatchBySpectrumCoverage){
    Hashtable<String,Vector<LipidParameterSet>> sheetParams = new Hashtable<String,Vector<LipidParameterSet>>();
    boolean omegaInfoAvailable = false;
    for (String className : classSequence.keySet()){
      Vector<LipidParameterSet> params = new Vector<LipidParameterSet>();
      Hashtable<String,Hashtable<String,QuantVO>> classQuant = quantObjects.get(className);
      for (String analyteName : analyteSequence.get(className)){
        Hashtable<String,QuantVO> analyteQuant = classQuant.get(analyteName);
        for (String mod : analyteQuant.keySet()){
          if (!omegaInfoAvailable && !analyteQuant.get(mod).getInfoForOmegaAssignment().isEmpty()) omegaInfoAvailable = true;
          if (results_.containsKey(className) && results_.get(className).containsKey(analyteName) && results_.get(className).get(analyteName).containsKey(mod)&&(results_.get(className).get(analyteName).get(mod)!=null)){
            Hashtable<String,LipidParameterSet> hitsOfOneMod = results_.get(className).get(analyteName).get(mod);
            if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE){
              for (String key:hitsOfOneMod.keySet()) params.add(hitsOfOneMod.get(key));
            }else{
              List<DoubleStringVO> keys = new ArrayList<DoubleStringVO>();
              for (String key:hitsOfOneMod.keySet()) {
                keys.add(new DoubleStringVO(key,Double.valueOf(key)));
              }
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
    
    if (omegaInfoAvailable) {
      applyOmegaInfo(correctedParams, classSequence, analyteSequence, quantObjects);
      LipidParameterSet.setOmegaInformationAvailable(true);
    }
    
    long timeMilliSeconds = (System.currentTimeMillis()-startCalcTime_);
    System.out.println(String.format("Required time: %s minutes, %s seconds", timeMilliSeconds/(60*1000), timeMilliSeconds%(60*1000)/1000));
    
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
     
      QuantificationResultExporter.writeResultsToExcel(resultFile,quantRes,bestMatchBySpectrumCoverage);
      
      if (isAlexTargetList || cli_){
        String alexResultFile = new String(resultFile);
        if (alexResultFile.endsWith(".xls") || alexResultFile.endsWith(".xlsx"))
          alexResultFile = alexResultFile.substring(0,alexResultFile.lastIndexOf("."));
        alexResultFile += ".tab";
        Vector<QuantificationResult> results = new Vector<QuantificationResult>();
        results.add(quantRes);
        RdbOutputWriter rdbWriter = new RdbOutputWriter();
        Hashtable<String,Vector<String>> analyteSeq = new Hashtable<String, Vector<String>>();
        analyteSeq.putAll(analyteSequence);
 //       rdbWriter.write(alexResultFile, results, classSequence, analyteSequence, quantObjects);
        rdbWriter.write(alexResultFile, results, classSequence, analyteSeq, quantObjects);
      }
    } catch (ExportException ex) {
      new WarningMessage(new JFrame(), "Error", ex.getMessage());
    }
    catch (Exception e) {
      e.printStackTrace();
      this.errorString_ = e.toString();
    }
  }
  
  //TODO: starting here is the omega stuff
  
  /**
   * Applies information about omega labels where enough evidence is present
   * @param correctedParams Hashtable containing only accepted isotopic probes
   * @param classSequence Class names in correctedParams, analyteSequence and quantObjects
   * @param analyteSequence Analytes in each class
   * @param quantObjects
   * @return correctedParams with added omega information
   */
  private void applyOmegaInfo(
      Hashtable<String,Vector<LipidParameterSet>> correctedParams, 
      LinkedHashMap<String,Integer> classSequence,
      LinkedHashMap<String,Vector<String>> analyteSequence,
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects) 
  {
    //deconvoluting the Hashtables
    for (String className : classSequence.keySet()){
      Hashtable<String,Hashtable<String,QuantVO>> classQuant = quantObjects.get(className);
      Vector<LipidParameterSet> lipidParameterSet = correctedParams.get(className);
      for (String analyteName : analyteSequence.get(className)){
        Hashtable<String,QuantVO> analyteQuant = classQuant.get(analyteName);
        for (String mod : analyteQuant.keySet()){
          Vector<DoubleBondPositionVO> infoForOmegaAssignment = analyteQuant.get(mod).getInfoForOmegaAssignment();
          if (!infoForOmegaAssignment.isEmpty()) {
            try {
              if (Integer.parseInt(RulesContainer.getAmountOfChains(StaticUtils.getRuleName(className, mod))) < 2) {
                addOmegaInformationToParameterSets(lipidParameterSet,analyteName,infoForOmegaAssignment,false);
              } 
            } catch (RulesException | NoRuleException | IOException | SpectrummillParserException ex) {
              System.out.println(ex.getMessage());
            }

            //add omega info to species requiring MSn evidence
            addOmegaInformationToParameterSets(lipidParameterSet,analyteName,infoForOmegaAssignment,true);
          }
        }
      }
    }
  }
  
  /**
   * If the available evidence is sufficient, omega information is added
   * @param lipidParameterSet
   * @param analyteName
   * @param infoForOmegaAssignment
   * @param needsMSnEvidence Lipid species with only one chain do not require MSn evidence for assignment
   */
  private void addOmegaInformationToParameterSets(Vector<LipidParameterSet> lipidParameterSets,
      String analyteName, Vector<DoubleBondPositionVO> infoForOmegaAssignment, boolean needsMSnEvidence) {
    for (int i = 0; i < lipidParameterSets.size(); i++) {
      if ((!needsMSnEvidence || (needsMSnEvidence && (lipidParameterSets.get(i) instanceof LipidomicsMSnSet)) 
          && !(lipidParameterSets.get(i).hasOmegaInformation()))) {
        LipidParameterSet param = lipidParameterSets.get(i);
        if (param.getNameStringWithoutRt().equals(analyteName)) {
          Range[] peakRanges = StaticUtils.determinePeakRanges(param);
          Range peakLimits = peakRanges[0];
          Range mediumAccuracy = peakRanges[1];
          Range highAccuracy = peakRanges[2];
          
          Iterator<DoubleBondPositionVO> it = infoForOmegaAssignment.iterator();
          while(it.hasNext()) {
            DoubleBondPositionVO doubleBondPositionVO = it.next();
            float expectedRetentionTime = doubleBondPositionVO.getExpectedRetentionTime();
            if (peakLimits.insideRange(expectedRetentionTime)) {
              if (mediumAccuracy.insideRange(expectedRetentionTime)) doubleBondPositionVO.setAccuracy(1);
              if (highAccuracy.insideRange(expectedRetentionTime)) doubleBondPositionVO.setAccuracy(2);
              
              if (needsMSnEvidence) {
                String chainCombination = getEquivalentChainCombination((LipidomicsMSnSet)param, doubleBondPositionVO);
                if (chainCombination != null) {
                  orderChainCombination(chainCombination, doubleBondPositionVO);
                  doubleBondPositionVO.setMolecularSpecies(chainCombination);
                  param.addOmegaInformation(doubleBondPositionVO);
                }
              } else if (!needsMSnEvidence) {
                String doubleBondPositionsHumanReadable = doubleBondPositionVO.getDoubleBondPositionsHumanReadable();
                doubleBondPositionVO.setMolecularSpecies(StaticUtils.getHumanReadableWODoubleBondPositions(doubleBondPositionsHumanReadable));
                param.addOmegaInformation(doubleBondPositionVO);
              }
            }
          }
          
          //sort omega information and add assignments
          Vector<DoubleBondPositionVO> paramOmegaInfo = param.getOmegaInformation();
          Collections.sort(paramOmegaInfo);
          Vector<DoubleBondPositionVO> highAccuracyHits = StaticUtils.getHighAccuracyDoubleBondPositions(paramOmegaInfo);
          Vector<DoubleBondPositionVO> assignedHits = StaticUtils.findUnambiguousDoubleBondPositions(highAccuracyHits);
          for (DoubleBondPositionVO assignedHit : assignedHits) {
            assignedHit.setIsAssigned(true);
          }
        }
      }
    }
  }
  
  /**
   * Orders the Vector of FattyAcidVOs in a DoubleBondPositionVO according to a given chain combination, if it contains positional evidence
   * @param chainCombination Human readable String of a molecular species without annotated double bond positions 
   * @param doubleBondPositionVO DoubleBondPositionVO Object describing the same molecular species as chainCombination
   */
  private void orderChainCombination(String chainCombination, DoubleBondPositionVO doubleBondPositionVO) {
    String doubleBondPositionsHumanReadable = doubleBondPositionVO.getDoubleBondPositionsHumanReadable();
    if (chainCombination.contains(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS)) {
      doubleBondPositionsHumanReadable = orderChainsAccordingToTemplate(chainCombination, doubleBondPositionsHumanReadable);
      try {
        Vector<FattyAcidVO> orderedFattyAcids = StaticUtils.decodeFAsFromHumanReadableName(
            doubleBondPositionsHumanReadable, Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), false, null);
        doubleBondPositionVO.setChainCombination(orderedFattyAcids);
      } catch (LipidCombinameEncodingException ex) {
        System.out.println(ex.getMessage());
      }
    }
  }
  
  /**
   * Orders a chain combination according to a template 
   * @param equivalentChainCombination Human readable String of a molecular species without annotated double bond positions
   * @param doubleBondPositionsHumanReadable Human readable String of a molecular species with annotated double bond positions
   * @return reordered String
   */
  private String orderChainsAccordingToTemplate(String equivalentChainCombination, String doubleBondPositionsHumanReadable) {
    String[] templateArray = StaticUtils.splitChainCombinationsAtChainSeparators(equivalentChainCombination);
    String[] doubleBondPositionsArray = StaticUtils.splitChainCombinationsAtChainSeparators(doubleBondPositionsHumanReadable);
    String[] noDoubleBondPositionsArray = new String[templateArray.length];
    String[] nakedPositionsArray = new String[templateArray.length];
    String orderedChains = null;
    
    //split doubleBondPositionsArray into an Array of Strings without double bond position assignment and a second array of just the assignments
    for(int i=0; i<templateArray.length; i++) {
      noDoubleBondPositionsArray[i] = StaticUtils.getHumanReadableWODoubleBondPositions(doubleBondPositionsArray[i]);
      nakedPositionsArray[i] = doubleBondPositionsArray[i].replace(noDoubleBondPositionsArray[i], "");
    }
    
    //create String of ordered chains
    for(int i=0; i<templateArray.length; i++) {
      for (int j=0; j<templateArray.length; j++) {
        if (noDoubleBondPositionsArray[j].equals(templateArray[i])) {
          if (i==0) {
            orderedChains = templateArray[0] + nakedPositionsArray[j];
          } else {
            orderedChains += LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS;
            orderedChains = orderedChains + templateArray[i] + nakedPositionsArray[j];
            break;
          }
        }
      }
    }
    return orderedChains;
  }
  
  /**
   * If available, the human readable String of a molecule with MSn evidence describing the same molecular species as doubleBondPositionVO is returned
   * @param lipidomicsMsnSet
   * @param doubleBondPositionVO
   * @return Human readable String of a molecular species without double bond position information
   */
  private String getEquivalentChainCombination(LipidomicsMSnSet lipidomicsMsnSet, DoubleBondPositionVO doubleBondPositionVO) {
    Set<String> mSnNamesHumanReadable = lipidomicsMsnSet.getHumanReadableNameSet();
    String doubleBondPositionsHumanReadable = doubleBondPositionVO.getDoubleBondPositionsHumanReadable();
    String equivalentChainCombination = null;
    Iterator<String> it = mSnNamesHumanReadable.iterator();
    while(it.hasNext()) {
      String mSnName = it.next();
      if (StaticUtils.isChainCombinationEquivalent(doubleBondPositionsHumanReadable, mSnName)) {
        equivalentChainCombination = mSnName;
      }
    }
    return equivalentChainCombination;
  }
  
  
  //TODO: this concludes the omega stuff
  
    
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
  private static void checkForIsobaricSpecies(LinkedHashMap<String,Integer> classSequence, LinkedHashMap<String,Vector<String>> analyteSequence,
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
          if (quant1.isQuantifiedByOtherIsobar()) continue;
          
          float tol = LipidomicsConstants.getCoarseChromMzTolerance((float)quant1.getAnalyteMass());
          float sameTol = tol/5f;
          try{
            RulesContainer.getAmountOfChains(StaticUtils.getRuleName(quant1.getAnalyteClass(), quant1.getModName()));
          } catch (NoRuleException nrx) {
            if(j==0) { System.out.println(nrx.getMessage()); }
            continue;
          }
          
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
            
            try{
              RulesContainer.getAmountOfChains(StaticUtils.getRuleName(quant2.getAnalyteClass(), quant2.getModName()));
            } catch (NoRuleException nrx) {
              if(j==0) { System.out.println(nrx.getMessage()); }
              continue;
            }
            
            float mz2 = (float)quant2.getAnalyteMass();
            boolean similarMass = lowerMz<mz2 && mz2<upperMz;
            boolean sameMass = lowerSame<mz2 && mz2<upperSame;
            if (similarMass){
              quant1.addIsobaricSpecies(quant2);
              quant2.addIsobaricSpecies(quant1);
//              if ((quant1.getAnalyteClass().equalsIgnoreCase("P-PE")&&quant1.getIdString().equalsIgnoreCase("36:4"))||
//                  (quant2.getAnalyteClass().equalsIgnoreCase("P-PE")&&quant2.getIdString().equalsIgnoreCase("36:4")))
//                System.out.println(quant1.getAnalyteClass()+quant1.getIdString()+"_"+quant1.getModName()+"/"+quant2.getAnalyteClass()+quant2.getIdString()+"_"+quant2.getModName()+" "+quant1.getAnalyteMass()+"/"+mz2);
              if (sameMass) {
                quant2.setQuantifiedByOtherIsobar(true);
              }
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
    Hashtable<String,Boolean> bestMatchBySpectrumCoverage = new Hashtable<String,Boolean>();
// LL   Hashtable<String,Vector<String>> analyteSequence = new Hashtable<String,Vector<String>>();
    LinkedHashMap<String,Vector<String>> analyteSequence = new LinkedHashMap<String,Vector<String>>();
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
        
        String analyteName = StaticUtils.generateLipidNameString(sideChain, doubleBonds, -1, "");
        analytes.add(analyteName);
        quantsOfClass.put(analyteName, quantsOfAnalyte);
      }
      quantObjects.put(className, quantsOfClass);
      analyteSequence.put(className, analytes);
      adductInsensitiveRtFilter.put(className, false);
      bestMatchBySpectrumCoverage.put(className, false);
    }
    checkForIsobaricSpecies(classSequence,analyteSequence,quantObjects);
    Vector results = new Vector();
    results.add(classSequence);
    results.add(analyteSequence);
    results.add(adductInsensitiveRtFilter);
    results.add(bestMatchBySpectrumCoverage);
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
   * @param hitsAccordingToQuant hash table returned from the SingleQuantThread
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

}
