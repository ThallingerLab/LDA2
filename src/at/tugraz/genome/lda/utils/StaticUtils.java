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

package at.tugraz.genome.lda.utils;

import java.awt.Color;
import java.awt.Component;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.export.vos.SummaryVO;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.lda.swing.RangeColor;
import at.tugraz.genome.lda.vos.IntegerStringVO;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;
import at.tugraz.genome.lda.xml.RawToChromThread;
import at.tugraz.genome.maspectras.GlobalConstants;
import at.tugraz.genome.maspectras.chromaviewer.MSMapViewer;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.utils.StringUtils;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
public class StaticUtils
{
  
  /** file extension for a rules file */
  public final static String RULE_FILE_SUFFIX = ".frag.txt";
  
  /** there is no ms2 evidence found for the hit*/
  public final static int NO_MS2 = 0;
  /** ms2 is found, but there is an overlap with an isobar from another species, and the two hits cannot be separated*/
  public final static int PERCENTAL_SPLIT = 1;
  /** ms2 is found, but there is an overlap with an isobar from another species, but the two hits can be separated*/ 
  public final static int SPLIT = 2;
  /** ms2 found and no overlap*/
  public final static int MS2_FULL = 3;
  /** the indicator for an alkylated chain combi name id which includes all information*/ 
  private final static String ALKYL_COMBI_PREFIX = LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_START+LipidomicsConstants.ALKYL_PREFIX+LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_STOP;
  /** the indicator for an alkenylated chain combi name id which includes all information*/ 
  private final static String ALKENYL_COMBI_PREFIX = LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_START+LipidomicsConstants.ALKENYL_PREFIX+LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_STOP;

  public final static String[] physicalMagnitudes_ = {"","m","\u03BC","n","p","f","a"};
  /**char sequence for cleaning empty positions from a chain combination*/
  private final static String EMPTY_AT_BEGINNING = "-"+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS;
  /**char sequence for cleaning empty positions from a chain combination*/
  private final static String EMPTY_IN_MIDDLE = LipidomicsConstants.CHAIN_SEPARATOR_NO_POS+"-"+LipidomicsConstants.CHAIN_SEPARATOR_NO_POS;
  /**char sequence for cleaning empty positions from a chain combination*/
  private final static String EMPTY_AT_END = LipidomicsConstants.CHAIN_SEPARATOR_NO_POS+"-";
  
  public static String[] extractFileNameAndSuffix(String fileNameFull){
    String fileName = "";
    String suffix = "";
    String fileNamePlusSuffix = StaticUtils.extractFileName(fileNameFull);
    int dotIndex = fileNamePlusSuffix.lastIndexOf(".");
    if (dotIndex>-1){
      fileName = fileNamePlusSuffix.substring(0,dotIndex);
      suffix = fileNamePlusSuffix.substring(dotIndex+1);
    }else{
      fileName = fileNamePlusSuffix;
    }
    String[] result = new String[2];
    result[0] = fileName;
    result[1] = suffix;
    return result;
  }
  
  public static String extractDirName(String fileName){
    int idx = 0;
    int slashIdx = fileName.lastIndexOf("/");
    int bkSlashIdx = fileName.lastIndexOf("\\");
    if (slashIdx>bkSlashIdx)
      idx = slashIdx;
    else
      idx = bkSlashIdx;
    String newFileName = fileName.substring(0,idx);
    return newFileName;
  }
  
  /**
   * searches the files for equal parts at the beginning of the file name, and the end of the file name
   * @param fileNames
   * @return [0] amount of equal characters at the beginning; [1] amount of equal characters at the end
   */
  public static int[] detectNameUnequalitiesBeforeAndAfter(Vector<String> fileNames){
    int[] charsToCut = new int[2];
    charsToCut[0] = 0;
    charsToCut[1] = 0;
    if (fileNames.size()>1){
      int smallestFileNameSize = 10000;
      for (String fileName : fileNames){
        if (fileName.length()<smallestFileNameSize)
          smallestFileNameSize = fileName.length();
      }
      for (int i=0;i!=smallestFileNameSize;i++){
        boolean isSameChar = true;
        char[] expName1 = fileNames.get(0).toCharArray();
        for (int j=0;j!=fileNames.size();j++){
          char[] expName2 = fileNames.get(j).toCharArray();
          if (expName1[expName1.length-1-i]!=expName2[expName2.length-1-i])
            isSameChar = false;
        }
        if (isSameChar)
          charsToCut[1]++;
        else
          break;
      }
      for (int i=0;i!=smallestFileNameSize;i++){
        boolean isSameChar = true;
        char[] expName1 = fileNames.get(0).toCharArray();
        for (int j=0;j!=fileNames.size();j++){
          char[] expName2 = fileNames.get(j).toCharArray();
          if (expName1[i]!=expName2[i])
            isSameChar = false;
        }
        if (isSameChar)
          charsToCut[0]++;
        else
          break;
      }
    }
    return charsToCut;
  }
  
  public static String extractChromBaseName(String resultFilePath, String experimentName){
    Vector<String> basePaths = new Vector<String>();
    File resultsFile = new File(resultFilePath);
    if (resultsFile.exists()&&resultsFile.isFile()){
      String directory = StaticUtils.extractDirName(resultFilePath);
      String chromFilePath = directory+File.separator+experimentName+".chrom";
      File chromFile = new File (chromFilePath);
      boolean chromFileExists = false;
      if (chromFile.exists())
        chromFileExists = true;
      else{
        File dirFile = new File(directory);
        if (dirFile.exists()&&dirFile.isDirectory()){
          File[] filesInDir = dirFile.listFiles();
          for (int i=0;i!=filesInDir.length;i++){
            if (filesInDir[i].getName().indexOf(experimentName)==0&&filesInDir[i].getAbsolutePath().endsWith(".chrom")){
//              if (!chromFileExists || filesInDir[i].getName().length()<chromFile.getName().length()){// chromFile = filesInDir[i];
              //  System.out.println("4.");
              basePaths.add(filesInDir[i].getAbsolutePath().substring(0,filesInDir[i].getAbsolutePath().length()-".chrom".length()));
            }
          }
        }
      }
      if (chromFileExists) basePaths.add(chromFile.getAbsolutePath().substring(0,chromFile.getAbsolutePath().length()-".chrom".length()));
    }
    if (basePaths.size()==0)
      return null;
    else if (basePaths.size()==1)
      return basePaths.get(0);
    else{
      String resultFileName = resultsFile.getName();
      String toReturn = basePaths.get(0);
      if (resultFileName.indexOf(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE)!=-1 ||
          resultFileName.indexOf(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE)!=-1){
        for (String basePath : basePaths){
          if ((resultFileName.indexOf(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE)!=-1 &&
              basePath.endsWith(RawToChromThread.FILE_SUFFIX_POLARITY_POSITIVE)) ||
              (resultFileName.indexOf(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE)!=-1 &&
              basePath.endsWith(RawToChromThread.FILE_SUFFIX_POLARITY_NEGATIVE)))
            toReturn = basePath;            
        }
      }else{
        for (String basePath : basePaths){
          if (basePath.length()<resultFileName.length())
            resultFileName = basePath;
        }
      }
      return toReturn;
    }
  }
  
  public static String existChromNecessaryFiles(String chromBaseName){
    String[] chromPaths = StringUtils.getChromFilePaths(chromBaseName+".");
    String oneDoesNotExist = "";
    if (chromPaths[0]==null || !existsFile(chromPaths[0])){
      oneDoesNotExist = "The chrom file is missing: "+chromPaths[0];
    }
    if (chromPaths[1]==null || !existsFile(chromPaths[1])){
      oneDoesNotExist = "The head file is missing: "+chromPaths[1];
    }
    if (chromPaths[2]==null || !existsFile(chromPaths[2])){
      oneDoesNotExist = "The idx file is missing: "+chromPaths[2];
    }
    if (chromPaths[3]==null || !existsFile(chromPaths[3])){
      oneDoesNotExist = "The rtt file is missing: "+chromPaths[3];
    }
    return oneDoesNotExist;    
  }
  
  public static boolean existChromFiles(String chromBaseName){
    String[] chromPaths = StringUtils.getChromFilePaths(chromBaseName+".");
    if (chromPaths[0]==null || !existsFile(chromPaths[0])){
      new WarningMessage(new JFrame(), "Error", "There is no file with the name "+extractFileName(chromBaseName+".chrom")+"!");
      return false;
    }
    if (chromPaths[1]==null || !existsFile(chromPaths[1])){
      new WarningMessage(new JFrame(), "Error", "There is no file with the name "+extractFileName(chromBaseName+".head")+"! A chrom file requires a head file!");
      return false;
    }
    if (chromPaths[2]==null || !existsFile(chromPaths[2])){
      new WarningMessage(new JFrame(), "Error", "There is no file with the name "+extractFileName(chromBaseName+".idx")+"! A chrom file requires an idx file!");
      return false;
    }
    if (chromPaths[3]==null || !existsFile(chromPaths[3])){
      new WarningMessage(new JFrame(), "Error", "There is no file with the name "+extractFileName(chromBaseName+".rtt")+"! A chrom file requires an rtt file!");
      return false;
    }
    return true;
  }
  
  public static String extractFileName(String fileName){
    int idx = 0;
    int slashIdx = fileName.lastIndexOf("/");
    int bkSlashIdx = fileName.lastIndexOf("\\");
    if (slashIdx>bkSlashIdx)
      idx = slashIdx;
    else
      idx = bkSlashIdx;
    String newFileName = fileName.substring(idx+1);
    return newFileName;
  }
  
  public static int getMaxApplicableIsotopeHash(Hashtable<String,Hashtable<String,ResultCompVO>> results, int maxIsotope){
    int isotopes = maxIsotope+1;
    for (Hashtable<String,ResultCompVO> hash : results.values()){
      int isos = StaticUtils.getMaxApplicableIsotope(hash, maxIsotope);
      if (isotopes>isos)
        isotopes = isos;
    }
    return isotopes;
  }
  
  public static int getMaxApplicableIsotope(Hashtable<String,ResultCompVO> results, int maxIsotope){
    int isotopes = maxIsotope+1;
    for (ResultCompVO compVO :results.values()){
      if (compVO.getOriginalArea(0)>0 && isotopes>compVO.getUsedIsotpes())
        isotopes = compVO.getUsedIsotpes();
    }
    return isotopes;
  }
  
  public static void printOutLabelToFile(String label,OutputStream out,Row row, CellStyle cellStyle,int columnCount,boolean excel,boolean tab)throws 
    IOException{
    if (excel){
      StaticUtils.addExcelLabel(cellStyle, row, columnCount, label);
    }else{
      if (tab)
        out.write("\t".getBytes());
      out.write(label.getBytes());
    }
  }
  
  public static void printOutNumberToFile(String number,OutputStream out,Row row, CellStyle cellStyle,int columnCount,boolean excel,boolean tab)throws 
    IOException{
    if (excel){
      StaticUtils.addExcelNumber(cellStyle, row,columnCount, number);
    }else{
      if (tab)
        out.write("\t".getBytes());
      out.write(number.getBytes());
    }
  }

  public static void addExcelLabel(CellStyle  cellStyle, Row row, int columnCount, String toAdd){
    Cell label = row.createCell(columnCount,Cell.CELL_TYPE_STRING);
    label.setCellValue(toAdd);
    if (cellStyle!=null)
      label.setCellStyle(cellStyle);
  }

  private static void addExcelNumber(CellStyle  cellStyle, Row row,int columnCount, String toAdd){
    Cell number = row.createCell(columnCount,Cell.CELL_TYPE_NUMERIC);
    number.setCellValue(Double.valueOf(toAdd));
    if (cellStyle!=null)
      number.setCellStyle(cellStyle);
  }

  public static String getAreaString(double area, ResultDisplaySettingsVO settingVO, String preferredUnit){
    double myArea = getAreaInCorrespondingUnit(area, preferredUnit);
    String unit = getCorrespondingUnit(settingVO, preferredUnit,false);
    String resultString = StaticUtils.getAreaTypeString(settingVO);
    if (resultString!=null && resultString.length()>1)
      resultString+=": ";
    if (!settingVO.getType().equalsIgnoreCase("relative value") /*&& !auStandardized*/){
      resultString += extractDisplayValue(myArea)+unit;
    }
    return resultString;
  }
  
  public static double getAreaInCorrespondingUnit(double area, String preferredUnit){
    double myArea = area;
    if (preferredUnit!=null){
      if (preferredUnit.equalsIgnoreCase("%"))
        myArea = myArea*100d;
      else if (preferredUnit.equalsIgnoreCase("\u2030"))
        myArea = myArea*1000d;
      else if (preferredUnit.equalsIgnoreCase("m"))
        myArea = myArea*1000d;
      else if (preferredUnit.equalsIgnoreCase("\u03BC"))
        myArea = myArea*1000000d;
      else if (preferredUnit.equalsIgnoreCase("n"))
        myArea = myArea*1000000000d;
      else if (preferredUnit.equalsIgnoreCase("p"))
        myArea = myArea*1000000000000d;
      else if (preferredUnit.equalsIgnoreCase("f"))
        myArea = myArea*1000000000000000d;
      else if (preferredUnit.equalsIgnoreCase("a"))
        myArea = myArea*1000000000000000000d;
    }
    return myArea;
  }
  
  public static double getValueDividedByUnit(double inValue, String unit){
    double value = new Double(inValue);
    if (unit!=null&&unit.length()>0){
      if (unit.equalsIgnoreCase("m"))
        value = value/1000d;
      if (unit.equalsIgnoreCase("\u03BC"))
        value = value/1000000d;
      if (unit.equalsIgnoreCase("n"))
        value = value/1000000000d;
      if (unit.equalsIgnoreCase("p"))
        value = value/1000000000000d;
      if (unit.equalsIgnoreCase("f"))
        value = value/1000000000000000d;
      if (unit.equalsIgnoreCase("a"))
        value = value/1000000000000000000d;
    }
    return value;
  }
  
  public static String getCorrespondingUnit(ResultDisplaySettingsVO settingVO, String preferredUnit, boolean addDescription){
    String unit = preferredUnit;
    String divisorUnit = settingVO.getDivisorMagnitude();
    if (addDescription)
      unit = "["+preferredUnit+"]";
    if (settingVO.getType().equalsIgnoreCase("amount end-volume")){
      unit = preferredUnit+"mol";
      if (addDescription)
        unit = "end volume ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("conc. end-volume")){
      unit = preferredUnit+"mol/"+divisorUnit+"L";
      if (addDescription)
        unit = "end conc ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("weight end-volume")){
      unit = preferredUnit+"g";
      if (addDescription)
        unit = "end weight ["+unit+"]";      
    } else if (settingVO.getType().equalsIgnoreCase("amount sample-volume")){
      unit = preferredUnit+"mol";
      if (addDescription)
        unit = "sample volume ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("conc. sample-volume")){
      unit = preferredUnit+"mol/"+divisorUnit+"L";
      if (addDescription)
        unit = "sample conc ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("weight sample-volume")){
      unit = preferredUnit+"g";
      if (addDescription)
        unit = "sample weight ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("relative to sample weight")){
      unit = preferredUnit+"mol/"+divisorUnit+"g";
      if (addDescription)
        unit = "rel. weight ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("relation to protein content")){
      if ((settingVO.getISStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION || settingVO.getESStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION)&&(!settingVO.isAu())){
        unit = preferredUnit+"mol/"+divisorUnit+"g";
        if (addDescription)
          unit = "rel. protein ["+unit+"]";
      } else{
        unit = "AU/"+divisorUnit+"g";
        if (addDescription)
          unit = "rel. protein ["+unit+"]";
//        auStandardized = true;
      }
    } else if (settingVO.getType().equalsIgnoreCase("relation to neutral lipid content")){
      if ((settingVO.getISStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION || settingVO.getESStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION)&&(!settingVO.isAu())){
        unit = preferredUnit+"mol/"+divisorUnit+"g";
        if (addDescription)
          unit = "rel. lipid ["+unit+"]";
      }else {
        unit = "AU/"+divisorUnit+"g";
        if (addDescription)
          unit = "rel. lipid ["+unit+"]";
//        auStandardized = true;
      }
    } else if (settingVO.getType().equalsIgnoreCase("relation to measured neutral lipid")){
      unit = preferredUnit+"mol/"+divisorUnit+"g";
      if (addDescription)
        unit = "rel. to lipid ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("relative to base peak")){
      unit = preferredUnit;
      if (addDescription)
        unit = "base peak ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("relative to measured class amount")){
      unit = preferredUnit;
      if (addDescription)
        unit = "total amount ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("relative to highest total peak")){
      unit = preferredUnit;
      if (addDescription)
        unit = "high tot peak ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("relative to total amount")){
      unit = preferredUnit;
      if (addDescription)
        unit = "rel total int ["+unit+"]";
    } else if (settingVO.getType().equalsIgnoreCase("relative value") && addDescription)
      unit = "area [AU]";
    return unit;
  }
  
  public static String extractDisplayValue(double area){
    return StaticUtils.extractDisplayValue(area, 3);
  }
  
  public static String extractDisplayValue(double area, int commaPositions){
    String value = "";
    if (!(Double.isInfinite(area)||Double.isNaN(area))){
      String doubleString = String.valueOf(StaticUtils.roundDBL(area, commaPositions));
      if (doubleString.indexOf(".")==-1 || doubleString.indexOf("E")!=-1 || doubleString.indexOf(".")+(commaPositions+1)>=doubleString.length())
        value = doubleString;
      else
        value = doubleString.substring(0,doubleString.indexOf(".")+(commaPositions+1));
    }
    return value;
  }
  
  private static double roundDBL(double targetDBL,int decimalPlace){
    BigDecimal bd = new BigDecimal(targetDBL);
    bd = bd.setScale(decimalPlace,BigDecimal.ROUND_HALF_UP);
    return (bd.doubleValue());
  }
  
  public static String getAreaTypeString(ResultDisplaySettingsVO settingVO){
    String resultString = "";
    if (settingVO.getType().equalsIgnoreCase("amount end-volume")){
      resultString +="Amount end-volume";
    } else if (settingVO.getType().equalsIgnoreCase("relative to base peak")){
      resultString +="Base peak";
    } else if (settingVO.getType().equalsIgnoreCase("relative to measured class amount")){
      resultString +="Total class amount";
    } else if (settingVO.getType().equalsIgnoreCase("relative to highest total peak")){
      resultString +="Total base peak";
    } else if (settingVO.getType().equalsIgnoreCase("relative to total amount")){
      resultString +="Total amount";
    } else if (settingVO.getType().equalsIgnoreCase("conc. end-volume")){
      resultString +="Conc. end-volume";
    } else if (settingVO.getType().equalsIgnoreCase("weight end-volume")){
      resultString +="Weight end-volume";
    } else if (settingVO.getType().equalsIgnoreCase("amount sample-volume")){
      resultString +="Amount sample-volume";
    } else if (settingVO.getType().equalsIgnoreCase("conc. sample-volume")){
      resultString +="Conc. sample-volume";
    } else if (settingVO.getType().equalsIgnoreCase("weight sample-volume")){
      resultString +="Weight sample-volume";
    } else if (settingVO.getType().equalsIgnoreCase("relative to sample weight")){
      resultString +="Rel. to sample weight";
    } else if (settingVO.getType().equalsIgnoreCase("relation to protein content")){
      resultString +="Rel. to protein";
    } else if (settingVO.getType().equalsIgnoreCase("relation to neutral lipid content")){
      resultString +="Rel. to neutral lipid";
    } else if (settingVO.getType().equalsIgnoreCase("relation to measured neutral lipid")){
      resultString +="Rel. to neutral lipid (measured)";
    }  
    return resultString;
  }
  
  public static String extractPreferredUnit(double value){
    if (value<=0){
      return null;
    } else if (value>=1)
      return "";
    else if ((value*1000d)>=1)
      return "m";
    else if ((value*1000000d)>=1)
      return "\u03BC";
    else if ((value*1000000000d)>=1)
      return "n";
    else if ((value*1000000000000d)>=1)
      return "p";
    else if ((value*1000000000000000d)>=1)
      return "f";
    else 
      return "a";
  }
  
  /**
   * this splits the chemical formula string in its components
   * the components can be additive '+' or subtractive '-"
   * an empty space is interpreted as '+' -> do not leave any empty spaces between the element and the amount
   * @param formula
   * @return
   * @throws ChemicalFormulaException
   */
  public static Hashtable<String,Integer> categorizeFormula(String formula) throws ChemicalFormulaException {
  	return categorizeFormula(formula, false);
  }
  
  /**
   * this splits the chemical formula string in its components
   * the components can be additive '+' or subtractive '-"
   * an empty space is interpreted as '+' -> do not leave any empty spaces between the element and the amount
   * @param formula									the chemical formula
   * @param allowLowerCase					true if lower case elements should be allowed, e.g. 'h' for a proton
   * @return
   * @throws ChemicalFormulaException
   */
  public static Hashtable<String,Integer> categorizeFormula(String formula, boolean allowLowerCase) throws ChemicalFormulaException {
    Hashtable<String,Integer> categorized = new Hashtable<String,Integer>();
    if (formula==null || formula.length()==0) return categorized;
    boolean add = true;
//    String formulaWithoutSpaces = "";
//    StringTokenizer tokenizer = new StringTokenizer(formula," ");
//    while (tokenizer.hasMoreTokens()){
//      formulaWithoutSpaces+=tokenizer.nextToken();
//    }
    char[] chars = formula.replaceAll(" ", "+").toCharArray();//formulaWithoutSpaces.toCharArray();
    String currentElement = "";
    String currentAmountString = "";
    for (int i=0;i!=chars.length;i++){
      if (chars[i]==' '){
        addFormulaPartToHash(categorized, currentElement, currentAmountString, null, add);
        add = true;
        currentElement = "";
        currentAmountString = "";
      }else if (chars[i]=='+'){
        addFormulaPartToHash(categorized, currentElement, currentAmountString, null, add);
        add = true;
        currentElement = "";
        currentAmountString = "";
      }else if (chars[i]=='-'){
        addFormulaPartToHash(categorized, currentElement, currentAmountString, null, add);
        add = false;
        currentElement = "";
        currentAmountString = "";
      }else if (Character.isLetter(chars[i])){
        if (Character.isUpperCase(chars[i]) || (allowLowerCase && Settings.getElementParser().isElementAvailable(String.valueOf(chars[i])))){
          if (currentElement.length()==0) currentElement+=String.valueOf(chars[i]);
          else{
            addFormulaPartToHash(categorized, currentElement, currentAmountString, null, add);
            currentElement = String.valueOf(chars[i]);
            currentAmountString = "";
          }
        }else{
          if (currentElement.length()==1 && currentAmountString.length()==0){
            currentElement+=String.valueOf(chars[i]);
          }else throw new ChemicalFormulaException("The formula \""+formula+"\" is invalid! An element must not start lower case!");
        }
      } else if (Character.isDigit(chars[i])){
        if (currentElement.length()>0)
          currentAmountString+=String.valueOf(chars[i]);
        else
          throw new ChemicalFormulaException("The formula \""+formula+"\" is invalid! A term cannot start with a digit! Avoid empty spaces between element and amount!");
      }
    }
    addFormulaPartToHash(categorized, currentElement, currentAmountString, null, add);
//    System.out.println("-------------------------------");
//    for (String elem : categorized.keySet()){
//      System.out.println(elem+": "+categorized.get(elem));
//    }
    return categorized;
  }
  
  /**
   * categorizes the formula of an Alex123 adduct in its elemental composition
   * @param formula the adduct formula
   * @return hashtable, key is the element, value is the amount
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   */
  public static Hashtable<String,Integer> categorizeAdduct(String formula) throws ChemicalFormulaException {
    Hashtable<String,Integer> categorized = new Hashtable<String,Integer>();
    if (formula==null || formula.length()==0) return categorized;
    boolean add = true;
    //char[] chars = formula.replaceAll(" ", "+").toCharArray();//formulaWithoutSpaces.toCharArray();
    char[] chars = formula.toCharArray();
    String currentElement = "";
    String currentAmountString = "";
    String currentMultiString = "";
    boolean formulaAlreadyAdded = false;
    for (int i=0;i!=chars.length;i++){
      if (chars[i]==' '){
        addFormulaPartToHash(categorized, currentElement, currentAmountString, currentMultiString, add);
        add = true;
        currentElement = "";
        currentAmountString = "";
        currentMultiString = "";
        formulaAlreadyAdded = true;
      }else if (chars[i]=='+'){
        if (!formulaAlreadyAdded) addFormulaPartToHash(categorized, currentElement, currentAmountString, currentMultiString, add);
        add = true;
        currentElement = "";
        currentAmountString = "";
        currentMultiString = "";
        formulaAlreadyAdded = true;
      }else if (chars[i]=='-'){
        if (!formulaAlreadyAdded) addFormulaPartToHash(categorized, currentElement, currentAmountString, currentMultiString, add);
        add = false;
        currentElement = "";
        currentAmountString = "";
        currentMultiString = "";
        formulaAlreadyAdded = true; 
      }else if (Character.isLetter(chars[i])){
        if (Character.isUpperCase(chars[i])){
          if (currentElement.length()==0){
            currentElement+=String.valueOf(chars[i]);
            formulaAlreadyAdded = false;
          }
          else{
            addFormulaPartToHash(categorized, currentElement, currentAmountString, currentMultiString, add);
            currentElement = String.valueOf(chars[i]);
            currentAmountString = "";
          }
        }else{
          if (currentElement.length()==1 && currentAmountString.length()==0){
            currentElement+=String.valueOf(chars[i]);
          }else throw new ChemicalFormulaException("The formula \""+formula+"\" is invalid! An element must not start lower case!");
        }
      } else if (Character.isDigit(chars[i])){
        if (currentElement.length()>0)
          currentAmountString+=String.valueOf(chars[i]);
        else
          currentMultiString+=String.valueOf(chars[i]); 
      }
    }
    addFormulaPartToHash(categorized, currentElement, currentAmountString, currentMultiString, add);
//    System.out.println("-------------------------------");
//    for (String elem : categorized.keySet()){
//      System.out.println(elem+": "+categorized.get(elem));
//    }
    return categorized;
  }
  
  /**
   * adds a chemical element to the categorized hash
   * @param categorized hashtable, key is the element, value is the amount
   * @param currentElement chemical symbol
   * @param currentAmountString the amount string of this element
   * @param currentMultiString a multiplication factor for elements
   * @param add add or remove these elemental parts
   */
  private static void addFormulaPartToHash(Hashtable<String,Integer> categorized, String currentElement, String currentAmountString, String currentMultiString, boolean add){
    if (currentElement!=null&&currentElement.length()>0){
      int finalAmount = 0;
      if (categorized.containsKey(currentElement))
        finalAmount = categorized.get(currentElement);
      int amountToChange = 1;
      if (currentAmountString!=null && currentAmountString.length()>0) amountToChange = Integer.parseInt(currentAmountString);
      int multi = 1;
      if (currentMultiString!=null && currentMultiString.length()>0) multi = Integer.parseInt(currentMultiString);      
      if (add) finalAmount += amountToChange*multi;
      else finalAmount -= amountToChange*multi;
      categorized.put(currentElement, finalAmount);
    }
  }
  
  /**
   * returns a catogorized formula in its Hill notation
   * @param formAnal the categorized formula
   * @param space make a space between the elements
   * @return formula in hill notation
   */
  public static String getFormulaInHillNotation(Hashtable<String,Integer> formAnal, boolean space){
    return getFormulaInHillNotation(formAnal, space, false);
  }
  
  /**
   * returns a catogorized formula in its Hill notation
   * @param formAnal the categorized formula
   * @param space make a space between the elements
   * @param plusSign puts a plus in front of the added elements 
   * @return formula in hill notation
   */
  private static String getFormulaInHillNotation(Hashtable<String,Integer> formAnal, boolean space, boolean plusSign){
    String hillNotation = "";
    if (formAnal.containsKey("C") && formAnal.get("C")!=0){
      hillNotation += ((plusSign&&formAnal.get("C")>0) ? "+" : "")+(formAnal.get("C")<0 ? "-" : "")+"C"+getHillNotationNumber("C",formAnal)+(space ? " " : "");
    }
    if (formAnal.containsKey("H") && formAnal.get("H")!=0){
      hillNotation += ((plusSign&&formAnal.get("H")>0) ? "+" : "")+(formAnal.get("H")<0 ? "-" : "")+"H"+getHillNotationNumber("H",formAnal)+(space ? " " : "");
    }
    List<String> list = new ArrayList<String>(formAnal.keySet());
    Collections.sort(list);
    for (String element : list){
      if (element.equalsIgnoreCase("C") || element.equalsIgnoreCase("H") || formAnal.get(element)==0) continue;
      hillNotation += ((plusSign&&formAnal.get(element)>0) ? "+" : "")+(formAnal.get(element)<0 ? "-" : "")+element+getHillNotationNumber(element,formAnal)+(space ? " " : "");
    }
    return hillNotation.trim();
  }

  /**
   * returns a catogorized formula in its Hill notation - the positive ones
   * @param formAnal the categorized formula
   * @param space make a space between the elements
   * @return formula in hill notation
   */
  public static String getFormulaInHillNotation_PlusFirst(Hashtable<String,Integer> formAnal, boolean space) {
    Hashtable<String,Integer> added = new Hashtable<String,Integer>();
    Hashtable<String,Integer> subtracted = new Hashtable<String,Integer>();
    String hillNotationPlusFirst = "";
    for (String element : formAnal.keySet()) {
      if (formAnal.get(element)>0)
        added.put(element, formAnal.get(element));
      else if (formAnal.get(element)<0)
        subtracted.put(element, formAnal.get(element));
    }
    if (added.size()>0) {
      hillNotationPlusFirst += getFormulaInHillNotation(added, space, true)+(subtracted.size()>0 ? " " : "");
    }
    if (subtracted.size()>0)
      hillNotationPlusFirst += getFormulaInHillNotation(subtracted, space, true);
    return hillNotationPlusFirst;
  }
  
  /**
   * adds a number to the chemical element in Hill notation; i.e. no number if 0, otherwise add the number
   * @param element the chemical element
   * @param formAnal the categorized formula
   * @return the number in Hill notation
   */
  private static String getHillNotationNumber(String element, Hashtable<String,Integer> formAnal){
    String number = "";
    int amount = Math.abs(formAnal.get(element));
    if (amount>1) number = String.valueOf(amount);
    return number;
  }

  
  public static String generateLipidNameString(String name, Integer doubleBonds, Integer omegaPos, String oxState){
    String nameString = name;
    
    if (doubleBonds!=null && doubleBonds>-1)
    {
      nameString += ":"+String.valueOf(doubleBonds);
      if(oxState!=null && !oxState.equals(""))
        nameString += LipidomicsConstants.CHAIN_MOD_SEPARATOR + oxState;
    }
    if (omegaPos!=null && omegaPos>-1)
      nameString += "(n-"+String.valueOf(omegaPos)+")";
    return nameString;
  }
  
  public static String generateLipidNameString(String name, Integer doubleBonds, String rt, String oxState){
    String nameString = generateLipidNameString(name,doubleBonds,-1,oxState);
    if (rt!=null&&rt.length()>0) nameString+="_"+rt;
    return nameString;
  }
  
  /**
   * encodes a human readable string unique to one chain including all information 
   * @param vo the FattyAcidVO containing all of the required information
   * @param includePrefix include the prefix
   * @param includeOmegaPosition include the omega position
   * @return a human readable string unique to one chain including all information 
   */
  public static String encodeLipidNameForCreatingCombis(FattyAcidVO vo, boolean includePrefix, boolean includeOmegaPosition) {
    return encodeLipidNameForCreatingCombis(vo.getChainType(),includePrefix ? vo.getPrefix() : "",vo.getcAtoms(),vo.getDoubleBonds(),
        vo.getOhNumber(),includeOmegaPosition ? vo.getOmegaPosition() : -1,vo.getOxState());
  }
  
  /**
   * encodes a human readable string unique to one chain including all information 
   * @param chainType LipidomicsConstants.CHAIN_TYPE_FA_ACYL|.CHAIN_TYPE_FA_ALKYL|.CHAIN_TYPE_FA_ALKENYL|.CHAIN_TYPE_LCB
   * @param prefix prefix for any isotopically labeled chains
   * @param cAtoms number of C atoms
   * @param doubleBonds number of double bonds
   * @param ohNumber the number of OH groups
   * @param omega the omega position
   * @param oxState the oxidation state
   * @return a human readable string unique to one chain including all information 
   */
  private static String encodeLipidNameForCreatingCombis(short chainType, String prefix, int cAtoms, int doubleBonds,
      int ohNumber, int omega, String oxState) {
    StringBuilder sb = new StringBuilder();
    if (chainType<LipidomicsConstants.CHAIN_TYPE_LCB)
      sb.append(LipidomicsConstants.CHAIN_TYPE_FA_NAME);
    else if (chainType==LipidomicsConstants.CHAIN_TYPE_LCB)
      sb.append(LipidomicsConstants.CHAIN_TYPE_LCB_NAME);
    //TODO: it would be good to throw an Exception here, if there is an unknown type - since this is not done, I put the ; directly after the else
    else;
    sb.append(LipidomicsConstants.CHAIN_OH_INDEX_SEPARATOR+String.valueOf(ohNumber)+LipidomicsConstants.CHAIN_NAME_TYPE_SEPARATOR);
    if (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ALKYL) 
      sb.append(ALKYL_COMBI_PREFIX);
    else if (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL)
      sb.append(ALKENYL_COMBI_PREFIX);
    sb.append(prefix);
    sb.append(generateLipidNameString(String.valueOf(cAtoms),doubleBonds,omega,oxState));
    return sb.toString();
  }
  
  /**
   * ATTENTION: this method should only be used to bring the lipid combi name in a more accessible format - it does not contain any mass or formula information
   * @param encoded the encoded name of the lipid id used in combinations
   * @return FattyAcidVO to have the information in more accessible format - it does not contain any mass or formula information
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static FattyAcidVO decodeLipidNameForCreatingCombis(String encoded) throws LipidCombinameEncodingException {
    StringTokenizer tokenizer = new StringTokenizer(encoded,LipidomicsConstants.CHAIN_NAME_TYPE_SEPARATOR);
    if (tokenizer.countTokens()!=2)
      throw new LipidCombinameEncodingException("The combiname \""+encoded+"\" cannot be decoded because it does not contain \""+LipidomicsConstants.CHAIN_NAME_TYPE_SEPARATOR+"\" or its usage is incorrect");
    String typeInfo = tokenizer.nextToken();
    String chainInfo = tokenizer.nextToken();
    tokenizer = new StringTokenizer(typeInfo,LipidomicsConstants.CHAIN_OH_INDEX_SEPARATOR);
    if (tokenizer.countTokens()!=2)
      throw new LipidCombinameEncodingException("The combiname \""+encoded+"\" cannot be decoded because it does not contain \""+LipidomicsConstants.CHAIN_OH_INDEX_SEPARATOR+"\" or its usage is incorrect");
    String typeString = tokenizer.nextToken();
    String ohString = tokenizer.nextToken();
    short chainType;
    if (typeString.equalsIgnoreCase(LipidomicsConstants.CHAIN_TYPE_FA_NAME))
      chainType=LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
    else if (typeString.equalsIgnoreCase(LipidomicsConstants.CHAIN_TYPE_LCB_NAME))
      chainType=LipidomicsConstants.CHAIN_TYPE_LCB;
    else
      throw new LipidCombinameEncodingException("The combiname \""+encoded+"\" cannot be decoded because it does not start with \""+LipidomicsConstants.CHAIN_TYPE_FA_NAME+"\" or \""+LipidomicsConstants.CHAIN_TYPE_LCB_NAME+"\"");
    int ohNumber;
    try {
      ohNumber = Integer.parseInt(ohString);
    }catch(NumberFormatException nfx) {
      throw new LipidCombinameEncodingException("The combiname \""+encoded+"\" is irregular, because the value followed by \""+LipidomicsConstants.CHAIN_OH_INDEX_SEPARATOR+"\" must be integer format");
    }
    if (chainInfo.startsWith(LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_START)) {
      if (chainInfo.indexOf(LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_STOP)==-1)
        throw new LipidCombinameEncodingException("The combiname \""+encoded+"\" is irregular, because if it contains a \""+LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_START+"\" must contain \""+LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_STOP+"\"");
      if (chainType>=LipidomicsConstants.CHAIN_TYPE_LCB)
        throw new LipidCombinameEncodingException("The combiname \""+encoded+"\" is irregular, because expressions that contain \""+LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_START+"\" must start with \""+LipidomicsConstants.CHAIN_TYPE_FA_NAME+"\"");
      if (!chainInfo.startsWith(ALKYL_COMBI_PREFIX) && !chainInfo.startsWith(ALKENYL_COMBI_PREFIX))
        throw new LipidCombinameEncodingException("The combiname \""+encoded+"\" is irregular - in the contents of \""+LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_START+"..."+LipidomicsConstants.CHAIN_LINKAGE_INCLUSION_STOP+"\""
            + " is only \""+LipidomicsConstants.ALKYL_PREFIX+"\" or \""+LipidomicsConstants.ALKENYL_PREFIX+"\" permitted");
      if (chainInfo.startsWith(ALKYL_COMBI_PREFIX)) {
        chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
        chainInfo = chainInfo.substring(ALKYL_COMBI_PREFIX.length());
      }else if (chainInfo.startsWith(ALKENYL_COMBI_PREFIX)) {
        chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
        chainInfo = chainInfo.substring(ALKENYL_COMBI_PREFIX.length());        
      }
    }
    StringBuilder prefixBuilder = new StringBuilder();
    char[] chars = chainInfo.toCharArray();
    int i = 0;
    while (i<chars.length && !Character.isDigit(chars[i])) {
      prefixBuilder.append(chars[i]);
      i++;
    }
    String prefix = prefixBuilder.toString();
    if (prefix.length()>0)
      chainInfo = chainInfo.substring(prefix.length());
    String[] cAndDbsAndOx;
    try {
    	cAndDbsAndOx = parseCAndDbsFromChainId(chainInfo);}catch(Exception ex) {
      throw new LipidCombinameEncodingException("The combiname \""+encoded+"\" is irregular. The chain info \""+chainInfo+"\" is not of the format $CAtoms$:$DoubleBonds$");
    }
    
    return new FattyAcidVO(chainType, prefix, Integer.parseInt(cAndDbsAndOx[0]), Integer.parseInt(cAndDbsAndOx[1]), ohNumber, -1d, null,cAndDbsAndOx[2]);
   
  }
  
  /**
   * decodes a combinatorial lipid name string containing chain type and position information to FattyAcidVOs
   * ATTENTION: this method should only be used to bring the lipid combi name in a more accessible format - it does not contain any mass or formula information
   * @param encodedCombi the encoded combinatorial lipid name string
   * @return a vector containing the single fatty acids objects of the combination
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static Vector<FattyAcidVO> decodeLipidNamesFromChainCombi(String encodedCombi) throws LipidCombinameEncodingException {
    Vector<String> encodedParts = splitChainCombiToEncodedStrings(encodedCombi,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
    Vector<FattyAcidVO> decoded = new Vector<FattyAcidVO>();
    for (String encoded : encodedParts) {
      decoded.add(decodeLipidNameForCreatingCombis(encoded));
    }
    return decoded;
  }

  
  /**
   * splits an encoded chain combination to its LDA encoded chain Ids 
   * @param encodedCombi the encoded LDA chain combination
   * @param the combi separator to suse
   * @return a vector containing the single encoded chains
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static Vector<String> splitChainCombiToEncodedStrings(String encodedCombi, String sep) throws LipidCombinameEncodingException {
    String combiPartial = new String(encodedCombi);
    Vector<String> encodedParts = new Vector<String>();
    while (combiPartial.length()>0) {
      String toAdd = "";
      if (combiPartial.indexOf(sep)==-1) {
        toAdd = combiPartial;
        combiPartial = "";
      }else {
        toAdd = combiPartial.substring(0,combiPartial.indexOf(sep));
        combiPartial = combiPartial.substring(combiPartial.indexOf(sep)
            +sep.length());
      }
      encodedParts.add(toAdd);
    }
    return encodedParts;
  }

  
  public static boolean existsFile(String pathName){
    File file = new File(pathName);
    return file.exists();
  }
  
  public static boolean isWithinTolerance(double tolerance, double ref, double value){
    return ((ref-tolerance)<=value && value<=(ref+tolerance));
  }

  public static String[] extractMoleculeRtAndModFromMoleculeName(String moleculeName){
    String[] molRtAndMod = new String[3];
    String mol = null;
    String rt = null;
    String mod = null;
    if (moleculeName.indexOf("_")!=-1){
      try{
        Double.valueOf(moleculeName.substring(moleculeName.lastIndexOf("_")+1));
        rt = moleculeName.substring(moleculeName.lastIndexOf("_")+1);
        mol= moleculeName.substring(0,moleculeName.lastIndexOf("_"));
      }catch(NumberFormatException nfx){
        mod = moleculeName.substring(moleculeName.lastIndexOf("_")+1);
        mol = moleculeName.substring(0,moleculeName.lastIndexOf("_"));
        if (mol.indexOf("_")!=-1){
          try{
            Double.valueOf(mol.substring(mol.lastIndexOf("_")+1));
            rt = mol.substring(mol.lastIndexOf("_")+1);
            mol = moleculeName.substring(0,mol.lastIndexOf("_"));
          }catch(NumberFormatException nfx2){}
        }
      }
    }else{
      mol = moleculeName;
    }
    molRtAndMod[0] = mol;
    molRtAndMod[1] = rt;
    molRtAndMod[2] = mod;
    return molRtAndMod;
  }
  
  public static float calculateDiff (float v1, float v2){
    float diff = v1-v2;
    if (diff<0) diff = diff*-1;
    return diff;
  }
  
  public static Vector<Double> calculateChemicalFormulaIntensityDistribution(ElementConfigParser elementParser, String chemicalFormula, int isosDesired, boolean nh4Corr) throws SpectrummillParserException{
    Vector<Vector<Double>> bothDistris = elementParser.calculateChemicalFormulaIntensityDistribution(chemicalFormula,isosDesired,nh4Corr);
    Vector<Double> distri = bothDistris.get(0);
    if (bothDistris.size()>1){
      Vector<Double> negDistri = bothDistris.get(1);
      if (useNegativeDistribution(distri,negDistri)) distri = negDistri;
    }
    return distri;
  }
  
  public static boolean useNegativeDistribution(Vector<Double> posDistri, Vector<Double> negDistri){
    return (negDistri.get(1)>posDistri.get(1));
  }
  
  
  /**
   * checks if one chain combination key is the permuted version of another combination key
   * @param one combi key one
   * @param two combi key two
   * @param sep the separator to use
   * @return true if they are the permuted combinations of one another
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static boolean isAPermutedVersion(String one, String two, String sep) throws LipidCombinameEncodingException{
    Vector<String> chains = splitChainCombiToEncodedStrings(one,sep);
    for (String combi : getPermutedChainNames(chains, sep)){
      if (combi.equalsIgnoreCase(two)) return true;
    }
    return false;
  }
  
  
  
  /**
   * specifies m/z ranges for coloring from an MSn identification
   * @param paramUncasted the identification not casted in a LipidomicsMSnSet
   * @param selectedMSn shall only a certain MSn identifcation be displayed
   * @param faHydroxyEncoding the OH encodings of the FA moiety
   * @param lcbHydroxyEncoding the OH encodings of the LCB moiety
   * @param isAlex123 do the fragments originate from an Alex123 target list
   * @return m/z ranges for coloring from an MSn identification - the first key is the MS-level
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static Hashtable<Integer,Vector<RangeColor>> createRangeColorVOs(LipidParameterSet paramUncasted, String selectedMSn,
      HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding, boolean isAlex123) throws LipidCombinameEncodingException{
    String selected = null;
    if (selectedMSn!=null && selectedMSn.length()>0){
      selected = new String(selectedMSn);
      selected = (String)cleanEmptyFAPositionsAndEncodeToLDACombiName(selected, faHydroxyEncoding, lcbHydroxyEncoding, isAlex123, null)[0];
    }
    Hashtable<Integer,Vector<RangeColor>> rangeColorLevels = null;
    if (paramUncasted instanceof LipidomicsMSnSet){
      LipidomicsMSnSet param = (LipidomicsMSnSet)paramUncasted;
      rangeColorLevels = new Hashtable<Integer,Vector<RangeColor>>();
      Hashtable<String,CgProbe> headFrags = new Hashtable<String,CgProbe>(param.getHeadGroupFragments());
      while (headFrags.size()>0){
        String highestArea = "";
        float area = 0f;
        for (String name : headFrags.keySet()){
          if (headFrags.get(name).Area>area){
            area = headFrags.get(name).Area;
            highestArea = name;
          }
        }
        CgProbe probe = headFrags.get(highestArea);
        RangeColor vo = new RangeColor(highestArea,MSMapViewer.COLOR_BROWN,probe,probe.Mz-probe.LowerMzBand,probe.Mz+probe.UpperMzBand);
        Vector<RangeColor> rangeColors = new Vector<RangeColor>();
        if (rangeColorLevels.containsKey(probe.getMsLevel()))
          rangeColors = rangeColorLevels.get(probe.getMsLevel());
        rangeColors.add(vo);
        rangeColorLevels.put(probe.getMsLevel(), rangeColors);
        headFrags.remove(highestArea);
      }
      Hashtable<String,Hashtable<String,CgProbe>> chainFrags =  param.getChainFragments();
      Hashtable<String,Double> relAreas = param.getChainCombinationRelativeAreas();
      Vector<String> combiOrdered = new Vector<String>();
      boolean ohInCombi = false;
      if (param.getStatus()<LipidomicsMSnSet.FRAGMENTS_DETECTED) {
        if (relAreas!=null && relAreas.size()>0)
          combiOrdered.add(relAreas.keySet().iterator().next());
      }else {        
        // the bigger areas should come first
        for (String combi : relAreas.keySet()){
          Vector<FattyAcidVO> chains = StaticUtils.decodeLipidNamesFromChainCombi(combi);
          if (areThereOhInCombi(chains))
            ohInCombi = true;
          boolean set = false;
          for (int i=0; i!=combiOrdered.size();i++){
            if (relAreas.get(combi)>relAreas.get(combiOrdered.get(i))){
              combiOrdered.add(i, combi);
              set = true;
              break;
            }
          }
          if (!set) combiOrdered.add(combi);
        }
      }
//      // this is for the sequence of color selection
//      Vector<String> identNames = new Vector<String>();
//      for (Object object : param.getMSnIdentificationNames()){
//        if (object instanceof String) identNames.add((String)object);
//        else identNames.add(((Vector<String>)object).get(0));
//      }
      int colorCount = 1;
      Hashtable<String,String> usedFAs = new Hashtable<String,String>();
      for (String combiOriginal : combiOrdered){
        String combi = combiOriginal;
        if (isAlex123) combi = combi.replaceAll("/",LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
        if (selected!=null && selected.length()>0 && !StaticUtils.isAPermutedVersion(selected,combi,LipidomicsConstants.CHAIN_COMBI_SEPARATOR)) continue;
        String name = combi;
        for (String oneOrder : combiOrdered){
          if (StaticUtils.isAPermutedVersion(oneOrder,combi,LipidomicsConstants.CHAIN_COMBI_SEPARATOR)){
            name = oneOrder;
          }
        }
        Vector<String> fas = splitChainCombiToEncodedStrings(name,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
        for (int i=0; i!= fas.size(); i++){
          String fa = fas.get(i);
          
          Hashtable<String,CgProbe> frags =  new Hashtable<String,CgProbe>();
          String faStored = getStoredFAName(fa,chainFrags);
          if (faStored!=null) frags =  chainFrags.get(faStored);

          Color color = getFragemntColor(colorCount);
          for (String key : frags.keySet()){
            String fragmentName = key;
            String faStoredReadable = StaticUtils.getHumanReadableChainName(decodeLipidNameForCreatingCombis(faStored), faHydroxyEncoding, lcbHydroxyEncoding, ohInCombi,
                isAlex123);
            if (!fragmentName.contains(faStoredReadable)){              
              fragmentName = getChainFragmentDisplayName(fragmentName,
                  StaticUtils.getHumanReadableChainName(decodeLipidNameForCreatingCombis(faStored), faHydroxyEncoding, lcbHydroxyEncoding, ohInCombi));
            }
            if (usedFAs.containsKey(fragmentName))continue;
            usedFAs.put(fragmentName, fragmentName);
            CgProbe probe = frags.get(key);
            RangeColor vo = new RangeColor(fragmentName,color,probe,probe.Mz-probe.LowerMzBand,probe.Mz+probe.UpperMzBand);
            Vector<RangeColor> rangeColors = new Vector<RangeColor>();
            if (rangeColorLevels.containsKey(probe.getMsLevel()))
              rangeColors = rangeColorLevels.get(probe.getMsLevel());
            rangeColors.add(vo);
            rangeColorLevels.put(probe.getMsLevel(), rangeColors);
          }
          colorCount++;
        }
      }
    }
    return rangeColorLevels;
  }
  
  public static String getStoredFAName(String fa, Hashtable<String,Hashtable<String,CgProbe>> chainFrags){
    if (!chainFrags.containsKey(fa) && !chainFrags.containsKey("FA "+fa) && !chainFrags.containsKey("LCB "+fa)) return null;
    String faStored = "";
    if (chainFrags.containsKey(fa)) faStored = fa;
    else if (chainFrags.containsKey("FA "+fa)) faStored = "FA "+fa;
    else if (chainFrags.containsKey("LCB "+fa)) faStored =  "LCB "+fa;
    else if (chainFrags.containsKey("O-"+fa)) faStored =  "O-"+fa;
    return faStored;
  }
  
  /**
   * returns a hash table that can be used by the 3D viewer (input are the VOs returned by createRangeColorVOs)
   * @param rangeColors the RangeColor VOs to be casted in the 3D viewer hash
   * @return hash table that can be used by the 3D viewer
   */
  public static Hashtable<Color,Vector<CgProbe>> get3DColorHash(Vector<RangeColor> rangeColors){
    Hashtable<Color,Vector<CgProbe>> colors = new Hashtable<Color,Vector<CgProbe>>();
    if (rangeColors==null) return colors;
    for (RangeColor color : rangeColors){
      Vector<CgProbe> probes = new Vector<CgProbe>();
      if (colors.containsKey(color.getColor())) probes = colors.get(color.getColor());
      probes.add(color.getProbe());
      colors.put(color.getColor(), probes);
    }
    return colors;
  }

  /**
   * specified sequence in which the colors shall be used
   * @param colorCount current color to use
   * @return the color to use
   */
  private static Color getFragemntColor(int colorCount){
    if (colorCount==1) return MSMapViewer.COLOR_RED;
    else if (colorCount==2) return  MSMapViewer.COLOR_BLUE;
    else if (colorCount==3) return  MSMapViewer.COLOR_GREEN;
    else if (colorCount==4) return  MSMapViewer.COLOR_VIOLET;
    else if (colorCount==5) return  MSMapViewer.COLOR_ORANGE;
    else if (colorCount==6) return  MSMapViewer.COLOR_TURQUOISE;
    else if (colorCount==7) return  MSMapViewer.COLOR_PINK;
    else if (colorCount==8) return  MSMapViewer.COLOR_GREEN_PASTEL;
    else if (colorCount==9) return  MSMapViewer.COLOR_BLUE_PASTEL;
    else if (colorCount==10) return  MSMapViewer.COLOR_GOLD;
    else if (colorCount==11) return  MSMapViewer.COLOR_GREEN_BLUE;
    else return Color.BLACK;
  }
  
  public static String getRuleName(String className, String modName){
    String ruleName = className;
    if (modName!=null && modName.length()>0) ruleName+="_"+modName;
    return ruleName;
  }
  
  /**
   * builds the name for the corresponding rule file name by the lipid class name and the used adduct
   * @param className the lipid class name
   * @param modName the adduct name
   * @return the name of the rule file of this lipid class/adduct combination
   */
  public static String getRuleFileName(String className, String modName){
    return  getRuleName(className,modName)+RULE_FILE_SUFFIX;
  }
  
  
  
  /**
   * cleans the empty positions of the LDA encoded name of the chain combination and returns the amount of positions
   * @param identification the human readable version of the lipid species
   * @param faHydroxyEncoding the OH encodings of the FA moiety
   * @param lcbHydroxyEncoding the OH encodings of the LCB moiety
   * @param isAlexOhEncoding is this a hit of the Alex123 target lists that contains an OH encoding
   * @param lipidomicsConstants constants of the read result file (call trace starting from LDAResultReader), null if called from elsewhere
   * @return [0] the combi name, [1] the number of assignable positions; [2] contains alkyl chains; [3] contains alkenyl chains
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static Object[] cleanEmptyFAPositionsAndEncodeToLDACombiName(String identification, HydroxyEncoding faHydroxyEncoding,
      HydroxyEncoding lcbHydroxyEncoding, boolean isAlexOhEncoding, LipidomicsConstants lipidomicsConstants) throws LipidCombinameEncodingException{
    Object[] result = new Object[4];
    int positions = 0;
    
    String name = identification.replaceAll(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS, LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    while (name.startsWith(EMPTY_AT_BEGINNING) || name.indexOf(EMPTY_IN_MIDDLE)!=-1 || name.endsWith(EMPTY_AT_END)){
      positions++;
      if (name.endsWith(EMPTY_AT_END)){
        name = name.substring(0,name.length()-EMPTY_AT_END.length());
      } else if (name.startsWith(EMPTY_AT_BEGINNING)){
        name = name.substring(EMPTY_AT_BEGINNING.length());
      } else if (name.indexOf(EMPTY_IN_MIDDLE)!=-1){
        name = name.substring(0,name.indexOf(EMPTY_IN_MIDDLE))+name.substring(name.indexOf(EMPTY_IN_MIDDLE)+EMPTY_AT_END.length());
      }
    }
    Vector<FattyAcidVO> chains = decodeFAsFromHumanReadableName(name,faHydroxyEncoding,lcbHydroxyEncoding,isAlexOhEncoding,lipidomicsConstants);
    positions += chains.size();
    boolean containsAlkyl = false;
    boolean containsAlkenyl = false;
    for (FattyAcidVO chain : chains) {
      if (chain.getChainType()==LipidomicsConstants.CHAIN_TYPE_FA_ALKYL)
        containsAlkyl = true;
      if (chain.getChainType()==LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL)
        containsAlkenyl = true;
    }
    result[0] = encodeLipidCombi(chains);
    result[1] = positions;
    result[2] = containsAlkyl;
    result[3] = containsAlkenyl;
    return result;
  }

  
  /**
   * checks if there is an overlap, and if there are overlaps to other species
   * @param set the hit to check
   * @return the status of the hit (NO_MS2, PERCENTAL_SPLIT, SPLIT, MS2_FULL) 
   */
  public static int checkMS2Evidence(LipidParameterSet set){
    int evd = NO_MS2;
    if (set.getPercentalSplit()>=0) evd = PERCENTAL_SPLIT;
    else if (set.getLowerRtHardLimit()>=0 || set.getUpperRtHardLimit()>=0) evd = SPLIT;
    else if (set instanceof LipidomicsMSnSet) evd = MS2_FULL;
    return evd;
  }
  
  @SuppressWarnings({ "unchecked", "rawtypes" })
  public static Vector checkFileStorage (File file, String suffix, Component parentComponent){
    Vector results = new Vector();
    boolean store = true;
    File fileToStore = new File(file.getAbsolutePath());
    if (fileToStore.exists()){
      if (JOptionPane.showConfirmDialog(parentComponent, "The file "+fileToStore.getName()+" exists! Replace existing file?") != JOptionPane.YES_OPTION)
        store = false;
    }else{
      if (fileToStore.getName().indexOf(".")==-1){
        fileToStore = new File(fileToStore.getAbsoluteFile()+"."+suffix);
        if (fileToStore.exists()){
          if (JOptionPane.showConfirmDialog(parentComponent, "The file "+fileToStore.getName()+" exists! Replace existing file?") != JOptionPane.YES_OPTION)
            store = false;
        }  
      }  
    }
    results.add(fileToStore);
    results.add(store);
    return results;
  }
  
  //TODO: this would be faster and more elegant, if FattyAcidVO objects were sorted instead.
  public static String sortFASequenceUnassigned(String key, String separator) throws LipidCombinameEncodingException{
    Vector<String> fas = StaticUtils.splitChainCombiToEncodedStrings(key,separator);
    Hashtable<Integer,Integer> emptyFas = new Hashtable<Integer,Integer>();
    Hashtable<Integer,Integer> unassignedFAs = new Hashtable<Integer,Integer>();
    for (int i=0; i!=fas.size();i++){
      if (fas.get(i).equalsIgnoreCase("-")) emptyFas.put(i, i);
      else unassignedFAs.put(i, i);
    }
    Vector<Integer> assignedFAs = new Vector<Integer>();
    while (unassignedFAs.size()>0){
      Vector<Integer> lowestCarbonNumbers = new Vector<Integer>();
      int lowestCarbonNumber = Integer.MAX_VALUE;
      for (int i : unassignedFAs.keySet()){
        String carbonPart = fas.get(i).substring(0,fas.get(i).indexOf(":"));
        while (carbonPart.length()>0 && !Character.isDigit(carbonPart.toCharArray()[0]))
          carbonPart = carbonPart.substring(1);
        int carbonNumber = Integer.parseInt(carbonPart);
        if (carbonNumber<lowestCarbonNumber ){
          lowestCarbonNumbers = new Vector<Integer>();
          lowestCarbonNumbers.add(i);
          lowestCarbonNumber = carbonNumber;
        }else if (carbonNumber==lowestCarbonNumber){
          lowestCarbonNumbers.add(i);
        }
      }
      if (lowestCarbonNumbers.size()==1){
        assignedFAs.add(lowestCarbonNumbers.get(0));
        unassignedFAs.remove(lowestCarbonNumbers.get(0));
      }else{
      	int containsOmega = 0;
        Vector<Integer> lowestDoubleBonds = new Vector<Integer>();
        int lowestDoubleBondNumber = Integer.MAX_VALUE;
        for (int j=0; j!=lowestCarbonNumbers.size(); j++){
          int indexInFA = lowestCarbonNumbers.get(j);
          String fa = fas.get(indexInFA);
          String dbString = fa.substring(fa.indexOf(":")+1);
          
          if (dbString.indexOf(LipidomicsConstants.OMEGA_POSITION_START)!=-1) {
          	dbString = dbString.substring(0,dbString.indexOf(LipidomicsConstants.OMEGA_POSITION_START));
          	containsOmega++;
          }
          if(dbString.indexOf(LipidomicsConstants.CHAIN_MOD_SEPARATOR)!=-1)
          	dbString = dbString.substring(0,dbString.indexOf(LipidomicsConstants.CHAIN_MOD_SEPARATOR));
          if (dbString.indexOf("(")!=-1)
            dbString = dbString.substring(0,dbString.indexOf("("));
          int doubleBonds = Integer.parseInt(dbString);
          if (doubleBonds<lowestDoubleBondNumber){
            lowestDoubleBonds = new Vector<Integer>();
            lowestDoubleBonds.add(j);
            lowestDoubleBondNumber = doubleBonds;
          } else if (doubleBonds==lowestDoubleBondNumber){
            lowestDoubleBonds.add(j);
          }
        }
        if (lowestDoubleBonds.size()>1 && lowestDoubleBonds.size()==containsOmega){
        	Vector<Integer> lowestOmegas = new Vector<Integer>();
        	int lowestOmegaNumber = Integer.MAX_VALUE;
        	for (int i=0; i<lowestDoubleBonds.size(); i++){
        		int indexInFA = lowestDoubleBonds.get(i);
        		String fa = fas.get(indexInFA);
        		if (fa.indexOf(LipidomicsConstants.OMEGA_POSITION_START)!=-1) {
        			String omegaString = fa.substring(fa.indexOf(LipidomicsConstants.OMEGA_POSITION_START)+3,fa.indexOf(LipidomicsConstants.OMEGA_POSITION_END));
        			int omega = Integer.parseInt(omegaString);
        			if (omega<lowestOmegaNumber){
                lowestOmegas = new Vector<Integer>();
                lowestOmegas.add(i);
                lowestOmegaNumber = omega;
              } else if (omega==lowestOmegaNumber){
              	lowestOmegas.add(i);
              }
        		}
        	}
        	for (Integer nrInLowestDBs : lowestOmegas)
        	{
        		int indexInFA = lowestDoubleBonds.get(nrInLowestDBs);
            assignedFAs.add(indexInFA);
            unassignedFAs.remove(indexInFA);
        	}
        }else{
        	for (Integer nrInLowestCs : lowestDoubleBonds){
            int indexInFA = lowestCarbonNumbers.get(nrInLowestCs);
            assignedFAs.add(indexInFA);
            unassignedFAs.remove(indexInFA);
          }
        }
      }
    }
    for (Integer indexInFA : emptyFas.keySet()){
      assignedFAs.add(indexInFA);
    }
    
    String faString = "";
    for (Integer faNr : assignedFAs){
      if (faString.length()>0) faString += LipidomicsConstants.CHAIN_SEPARATOR_NO_POS;
      faString += fas.get(faNr);
    }
    return faString;
  }
  
  /**
   * creates the corresponding display name - for Alex123, the position of the fatty acid is different
   * @param name name of the fragment
   * @param faName name of the fatty acyl chain
   * @return the display name of the chain fragment
   */
  public static String getChainFragmentDisplayName(String name, String faName){
    String displayName = name+"("+faName+")";
    if (isAnAlex123Fragment(name)){
      String chainPrefix = "FA";
      if (name.equals("LCB") || name.startsWith("LCB ") || name.startsWith("-LCB "))
        chainPrefix = "LCB";
      if (name.indexOf(faName)==-1) {
        if (name.equals(chainPrefix))
          displayName = chainPrefix+" "+faName;
        else if (name.startsWith(chainPrefix+" "))
          displayName = chainPrefix+" "+faName+name.substring((chainPrefix+" ").length());
        else if (name.startsWith("-"+chainPrefix+" "))
          displayName = "-"+chainPrefix+" "+faName+name.substring(("-"+chainPrefix+" ").length());
      }else{      
        displayName = name;
      }
    }
    return displayName;
  }
  
  
  /**
   * checks whether this is an Alex123 fragment encoding
   * @param name the name of the fragment
   * @return true when an Alex123 fragment encoding shall be used
   */
  public static boolean isAnAlex123Fragment(String name) {
    boolean useAlex = false;
    // this is an extension to support the Alex123 nomenclature
    try {
      Class.forName( "at.tugraz.genome.lda.Settings");
      useAlex = Settings.useAlex(); 
    } catch( ClassNotFoundException e ) { }
    if (useAlex  && (name.equals("FA") || name.startsWith("FA ") || name.startsWith("-FA ") ||
        name.equals("LCB") || name.startsWith("LCB ") || name.startsWith("-LCB ") ||
        name.equals("O-") || name.startsWith("O-") || name.startsWith("-O-"))){
      return true;
    }
    return false;
  }
  
  /**
   * extracts the fatty acyl chain name and the fragment name from the readable name
   * @param readableFragmentName the readable name containing fragment name and fatty acyl name
   * @param chainType the type of the chain
   * @param oh the number of OH positions
   * @param useOldEncoding indicates if old (before 2.8.4) or new encoding should be used
   * @return [0] the encoded chain object (FattyAcidVO); [1] fragment name
   * @throws LipidCombinameEncodingException thrown when a human readable chain ID (without type and OH number) cannot be decoded
   */
  public static Object[] parseChainFaAndFragmentNameFromExcel(String readableFragmentName, short chainType, int oh, boolean useOldEncoding) throws LipidCombinameEncodingException{
	  
	  Object[] result = new Object[2];
    String faName = "";
    String fragmentName = "";
    String oxState = "";
    //this is for Alex123 naming
    if (readableFragmentName.startsWith("FA ")||readableFragmentName.startsWith("-FA ")||
        readableFragmentName.startsWith("LCB ")||readableFragmentName.startsWith("-LCB ")||
        readableFragmentName.startsWith("O-")||readableFragmentName.startsWith("-O-")){
      int start = 0;
      String prefix = "FA";
      if (readableFragmentName.startsWith("LCB ")||readableFragmentName.startsWith("-LCB "))
        prefix = "LCB";
      if (readableFragmentName.startsWith("O-")||readableFragmentName.startsWith("-O-"))
        prefix = "O-";
      if (readableFragmentName.startsWith(prefix+" ") || readableFragmentName.startsWith("O-"))
        fragmentName = prefix;
      else
        fragmentName = "-"+prefix;
      if (prefix.equalsIgnoreCase("FA")||prefix.equalsIgnoreCase("LCB"))
        fragmentName += " ";
      start = fragmentName.length();
      boolean isChain = true;
      int stop = start;
      char[] chars = readableFragmentName.toCharArray();
      if (chars.length>(start+1) && chars[start+1]=='-' && (chars[start]=='O' || chars[start]=='P')){
        stop = stop+2;
      }
      while (isChain && stop<readableFragmentName.length()){
        if (Character.isDigit(chars[stop]) || chars[stop]==':' || chars[stop]==LipidomicsConstants.ALEX_OH_SEPARATOR.toCharArray()[0])
          stop++;
        else
          isChain=false;
      }
      faName = readableFragmentName.substring(start,stop);
      if (faName.indexOf(LipidomicsConstants.ALEX_OH_SEPARATOR)!=-1)
        faName = faName.substring(0,faName.indexOf(LipidomicsConstants.ALEX_OH_SEPARATOR));
      fragmentName += readableFragmentName.substring(stop);
    }else{
      faName = readableFragmentName.substring(readableFragmentName.indexOf("(")+1,readableFragmentName.lastIndexOf(")"));
      fragmentName = readableFragmentName.substring(0,readableFragmentName.indexOf("("));
    }
    if(faName.contains(LipidomicsConstants.CHAIN_MOD_SEPARATOR))
    {
    	oxState = faName.substring(faName.indexOf(LipidomicsConstants.CHAIN_MOD_SEPARATOR)+1);
    //	fragmentName = readableFragmentName.substring(0, readableFragmentName.indexOf(LipidomicsConstants.CHAIN_MOD_SEPARATOR1));
    }
    
    Object[] prefixCAndDbs = null;
    try {
    	prefixCAndDbs = parsePrefixCAndDbsFromChainId(faName);
    }catch (Exception e) {
      throw new LipidCombinameEncodingException("The chain-id: \""+readableFragmentName+"\" cannot be decoded");
    }
    result[0] = new FattyAcidVO(chainType, (String)prefixCAndDbs[0], (Integer)prefixCAndDbs[1], (Integer)prefixCAndDbs[2], oh, -1, null,oxState);
    result[1] = fragmentName.trim();
    //for Alex123
    if (readableFragmentName.startsWith("FA ")||readableFragmentName.startsWith("-FA ")||
        readableFragmentName.startsWith("LCB ")||readableFragmentName.startsWith("-LCB ")||
        readableFragmentName.startsWith("O-")||readableFragmentName.startsWith("-O-")){
      result[1] = readableFragmentName;
    }
    return result;
  }

  
  /**
   * returns a list of all permuted name variants of the presented input vector
   * the permuted names are separated by FA_SEPARATOR
   * the method checks and removes duplicate name entries
   * @param singleCombiParts vector containing the names to be permuted
   * @param sep the separator to be used to concatenate the combi - LipidomicsConstants.CHAIN_SEPARATOR_NO_POS is used in the case this is null
   * @return all permuted name variants of the presented input vector - permuted names are separated by FA_SEPARATOR
   */
  public static Vector<String> getPermutedChainNames(Vector<String> singleCombiParts, String sep){
    String currentName = "";
    Vector<String> permutedNames = addPermutedPart(currentName, singleCombiParts,sep);
    return permutedNames;
  }
  
  /**
   * recursive call to assign the individual names to the permuted name string
   * @param currentName value that is assigned to the name until now
   * @param singleCombiParts remaining individual names that have to be assigned
   * @param sep the separator to be used to concatenate the combi - LipidomicsConstants.CHAIN_SEPARATOR_NO_POS is used in the case this is null
   * @return the permuted name
   */
  private static Vector<String> addPermutedPart(String currentName, Vector<String> singleCombiParts, String sep){
    Vector<String> permutedNames = new Vector<String>();
    for (int i=0;i!=singleCombiParts.size();i++){
      String name = new String(currentName);
      name+=singleCombiParts.get(i)+(sep!=null ? sep : LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
      Vector<String> subCombiParts = new Vector<String>(singleCombiParts);
      subCombiParts.remove(i);
      if (subCombiParts.size()>0){
        permutedNames.addAll(addPermutedPart(name,subCombiParts,sep!=null ? sep : LipidomicsConstants.CHAIN_SEPARATOR_NO_POS));
      }else{
        name = name.substring(0,name.length()-(sep!=null ? sep.length() : LipidomicsConstants.CHAIN_SEPARATOR_NO_POS.length()));
        permutedNames.add(name);
      }
    }
    return permutedNames;
  }
  
  
  /**
   * checks whether all of the detected chains have the same number of carbon atoms and double bonds 
   * @param combiName the molecular species name
   * @return true when all of the detected chains have the same number of carbon atoms and double bonds
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static boolean areAllChainsTheSame(String combiName) throws LipidCombinameEncodingException{
    Vector<String> fas = StaticUtils.splitChainCombiToEncodedStrings(combiName.replaceAll("/", "_"),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    boolean theSame = true;
    for (int i=1; i!=fas.size(); i++){
      if (!fas.get(i).equalsIgnoreCase(fas.get(0))){
        theSame = false;
        break;
      }
    }
    return theSame;
  }
  
  public static boolean checkChemicalFormula(String formula){
    if (formula.contains("-")){
      new WarningMessage(new JFrame(), "Error", "The formula "+formula+" must not contain any negative values!");
      return false;
    }
    char[] formulaChars = formula.toCharArray();
    String formulaToCheck = "";
    boolean isPreviousDigit = false;
    for (int i=0;i!=formulaChars.length;i++){
      char currentChar = formulaChars[i];
      if (isPreviousDigit && !Character.isDigit(currentChar)){
        formulaToCheck+=" ";
      }
      formulaToCheck+=String.valueOf(currentChar);
      isPreviousDigit = Character.isDigit(currentChar);
    }
    ElementConfigParser aaParser = Settings.getElementParser();
    try {
      aaParser.calculateTheoreticalMass(formulaToCheck, false);
      return true;
    }
    catch (SpectrummillParserException e) {
      new WarningMessage(new JFrame(), "Error", "The formula "+formula+" is not OK! "+e.getMessage());
      return false;
    } 
  }
  
  /**
   * in a group of peaks, there is one peak with chain information sufficient to split the
   * areas according to this information
   * @param speciesName the name of the lipid species (not the molecular species)
   * @param sets the group of lipid parameter sets
   * @return true when there is chain information present
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static boolean isThereChainInformationAvailable(String speciesName,Vector<LipidParameterSet> sets) throws LipidCombinameEncodingException {
    boolean available = false;
    for (LipidParameterSet set : sets){
      if (isThereChainInformationAvailable(speciesName,set)){
        available = true;
        break;
      }
    }
    return available;
  }

  /**
   * checks whether a peak contains chain information
   * @param set the lipid parameter set
   * @return true when there is chain information present
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static boolean isThereChainInformationAvailable(String speciesName,LipidParameterSet set) throws LipidCombinameEncodingException{
    if (set instanceof LipidomicsMSnSet && ((LipidomicsMSnSet)set).getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED){
      LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
      boolean isSelfIdentification = false;
      Vector<Object> identified = msn.getMSnIdentificationNames();
      if (identified.size()==1 && (identified.get(0) instanceof String) &&
          ((String)identified.get(0)).equalsIgnoreCase(speciesName)){
        isSelfIdentification = true;
      }
      if (!isSelfIdentification){
        return true;
      }
    }
    return false;
  }
  
  /**
   * determines whether this LipidParameterSet contains MS1 evidence only or "clean" MS2 evidence, or there is a split 
   * @param set the LipidParameterSet to vet
   * @return the evidence level according to the definitions in the SummaryVO
   */
  public static short determineEvidenceStateOfHit (LipidParameterSet set) {
    int ev = checkMS2Evidence(set);
    switch (ev) {
    case NO_MS2:
      return SummaryVO.EVIDENCE_MS1_ONLY;
    case PERCENTAL_SPLIT:
      return SummaryVO.EVIDENCE_MS2_NO_SPLIT_POSSIBLE;
    case SPLIT:
      return SummaryVO.EVIDENCE_MS2_SPLIT;
    case MS2_FULL:
      return SummaryVO.EVIDENCE_MS2_UNAMBIGUOUS;
    default:
      return SummaryVO.EVIDENCE_MS1_ONLY;
    }
  }
  
  /**
   * parses an input in fraction, percent and permille format and returns the corresponding double value
   * ATTENTION: the value must be positive
   * @param inValue the String representation of the fraction, percent or permille format
   * @return the double representation of the fraction, percent or permille format
   * @throws NumberFormatException if the value is not a number (percent and permille sign at the end allowed) 
   */
  public static double readPercentPermilleValue(String inValue) throws NumberFormatException {
    boolean percent = false;
    boolean permille = false;
    String value = new String(inValue);
    if (value.endsWith("%")){
      percent = true;
      value = value.substring(0,value.length()-1);
    } else if (value.endsWith("\u2030")){
      permille = true;
      value = value.substring(0,value.length()-1);
    }
    try{
      double doubleValue = Double.parseDouble(value);
      if (percent) doubleValue /= 100;
      else if (permille) doubleValue /= 1000;
      return doubleValue;
    }catch(NumberFormatException nfx){
      throw new NumberFormatException("The value "+inValue+" has not the correct format!");
    }   
  }
  
  /**
   * returns an int[] containing the number of C atoms at int[0] and the number of double bonds at int[1]
   * @param chainId the string representative of a chain
   * @return int[] containing the number of C atoms at int[0] and the number of double bonds at int[1]
   * @throws Exception
   */
  public static String[] parseCAndDbsFromChainId(String chainId) throws Exception {
    String[] cAndDbs = chainId.split(":|\\"+LipidomicsConstants.CHAIN_MOD_SEPARATOR);
    String[] cad = new String[3];
    
    if (cAndDbs.length==2)
    {	
    	cad[2] = "";
    }else if (cAndDbs.length==3) {
    	cad[2] = cAndDbs[2];
    }else {
    	throw new Exception("The chain id \""+chainId+"\" is not valid!");
    }
    
    try {Integer.parseInt(cAndDbs[0]);}catch(NumberFormatException nfx) {throw new Exception("The chain id \""+chainId+"\" is not valid!");}
    try {Integer.parseInt(cAndDbs[1]);}catch(NumberFormatException nfx) {throw new Exception("The chain id \""+chainId+"\" is not valid!");}
    
    cad[0] = cAndDbs[0];
    cad[1] = cAndDbs[1];
    //System.out.println("cAndDbs: " + cad[0] + "-" + cad[1] + "-" + cad[2]);
    return cad;
  }
  
  
  /**
   * returns the display name for an FA/LCB encoded molecular species identification
   * @param combiName the encoded name for the chain combination
   * @param chains the decoded information about the molecular species
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @return the encoded human readable display name for a chain combination
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static String getHumanReadableCombiName(String combiName, HydroxyEncoding faEncoding,
      HydroxyEncoding lcbEncoding) throws LipidCombinameEncodingException {
    return getHumanReadableCombiName(StaticUtils.decodeLipidNamesFromChainCombi(combiName), faEncoding, lcbEncoding);
  }
    
  /**
   * returns the display name for an FA/LCB encoded molecular species identification
   * @param chains the decoded information about the molecular species
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @return the encoded human readable display name for a chain combination
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static String getHumanReadableCombiName(Vector<FattyAcidVO> chains, HydroxyEncoding faEncoding,
        HydroxyEncoding lcbEncoding) throws LipidCombinameEncodingException {
    StringBuilder combi = new StringBuilder();
    boolean ohPresent = areThereOhInCombi(chains);
    for (FattyAcidVO chain : chains) {
      if (combi.length()!=0) combi.append(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
      combi.append(getHumanReadableChainName(chain, faEncoding, lcbEncoding, ohPresent));
    }
    return sortFASequenceUnassigned(combi.toString(),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
  }
  
  
  /**
   * are there any hydroxylation sites in a chain combination
   * @param chains the chains in the combination
   * @return true when there are hydroxylation sites present
   */
  public static boolean areThereOhInCombi(Vector<FattyAcidVO> chains) {
    boolean ohPresent = false;
    for (FattyAcidVO chain : chains) {
      if (chain.getOhNumber()>0) {
        ohPresent = true;
        break;
      }
    }
    return ohPresent;
  }
  

  /**
   * returns the human readable display name for a lipid molecular species, based on the current conventions
   * @param chain the chain object
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @param ohInCombi are there any OH groups in this species, so that the number of hydroxylation sites must be encoded
   * @return the human readable display name for a lipid molecular species
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static String getHumanReadableChainName(FattyAcidVO chain, HydroxyEncoding faEncoding, HydroxyEncoding lcbEncoding, 
      boolean ohInCombi) throws LipidCombinameEncodingException {
    return getHumanReadableChainName(chain,faEncoding,lcbEncoding,ohInCombi,false);
  }

  
  
  /**
   * returns the human readable display name for a lipid molecular species, based on the current conventions
   * @param chain the chain object
   * @param faEncoding the hydroxylation encoding for FA chains
   * @param lcbEncoding the hydroxylation encoding for LCB chains
   * @param ohInCombi are there any OH groups in this species, so that the number of hydroxylation sites must be encoded
   * @return the human readable display name for a lipid molecular species
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static String getHumanReadableChainName(FattyAcidVO chain, HydroxyEncoding faEncoding, HydroxyEncoding lcbEncoding, 
      boolean ohInCombi, boolean useAlexEncoding) throws LipidCombinameEncodingException {
    HydroxyEncoding encoding = null;
    if (chain.getChainType()==LipidomicsConstants.CHAIN_TYPE_FA_ACYL || chain.getChainType()==LipidomicsConstants.CHAIN_TYPE_FA_ALKYL ||
        chain.getChainType()==LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL)
      encoding = faEncoding;
    else if (chain.getChainType()==LipidomicsConstants.CHAIN_TYPE_LCB)
      encoding = lcbEncoding;
    StringBuilder encoded = new StringBuilder();
    if (useAlexEncoding) {
      encoded.append(chain.getCarbonDbsId());
      if (chain.getOhNumber()>0)
        encoded.append(LipidomicsConstants.ALEX_OH_SEPARATOR+chain.getOhNumber());
    }else{
      if (chain.getChainType()==LipidomicsConstants.CHAIN_TYPE_FA_ALKYL)
        encoded.append(LipidomicsConstants.ALKYL_PREFIX);
      else if (chain.getChainType()==LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL)
        encoded.append(LipidomicsConstants.ALKENYL_PREFIX);
      if (ohInCombi) {
        if (encoding==null)
          throw new LipidCombinameEncodingException("For lipid species containing OH groups, the corresponding encoding must be defined. For \""+chain.getChainId()+"\", this is not the case!");
        try {encoded.append(encoding.getEncodedPrefix((short)chain.getOhNumber()));
        }catch (HydroxylationEncodingException e) {
          throw new LipidCombinameEncodingException(e);
        }
      }
      encoded.append(chain.getCarbonDbsId());
    }
    return encoded.toString();
  }
  
  
  /**
   * returns a human readable chain type for the storage in the Excel file
   * @param chainType the encoded LDA chain type
   * @return human readable chain type for the storage in the Excel file
   */
  public static String getHumanReadableChainType(short chainType) {
    if (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ACYL)
      return FragmentRuleVO.CHAIN_NAME.substring(1);
    else if (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ALKYL)
      return FragmentRuleVO.ALKYL_CHAIN_NAME.substring(1);
    else if (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL)
      return FragmentRuleVO.ALKENYL_CHAIN_NAME.substring(1);
    else if (chainType==LipidomicsConstants.CHAIN_TYPE_LCB)
      return FragmentRuleVO.LCB_NAME.substring(1);
    else
      return ("ERROR: the chain type \""+chainType+"\" is not allowed!"); 
  }
  
  
  /**
   * encodes a sequence of chains to the LDA internal presentation of a chain combination
   * @param chains the single chain objects
   * @return LDA internal presentation of a chain combination
   */
  public static String encodeLipidCombi(Vector<FattyAcidVO> chains) {
    Vector<String> chainIds = new Vector<String>();
    for (FattyAcidVO chain : chains) chainIds.add(chain.getChainId());
    return encodeLipidCombiFromIds(chainIds);
  }
  
  /**
   * encodes a sequence of chains to the LDA internal presentation chain ids of a chain combination
   * @param chainIds the single chain objects
   * @return LDA internal presentation of a chain combination
   */
  public static String encodeLipidCombiFromIds(Vector<String> chainIds) {
    StringBuilder encoded = new StringBuilder();
    for (String chainId : chainIds) {
      if (encoded.length()!=0) encoded.append(LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
      encoded.append(chainId);
    }
    return encoded.toString();
  }

  
  /**
   * decodes a human readable chain combination
   * @param humanReadable the human readable lipid molecular species
   * @param faHydroxyEncoding the OH encodings of the FA moiety
   * @param lcbHydroxyEncoding the OH encodings of the LCB moiety
   * @param isAlexOhEncoding is this a hit of the Alex123 target lists that contains an OH encoding
   * @param lipidomicsConstants constants of the read result file (call trace starting from LDAResultReader), null if called from elsewhere
   * @return the decoded chains of the composition
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static Vector<FattyAcidVO> decodeFAsFromHumanReadableName(String humanReadable, HydroxyEncoding faHydroxyEncoding,
      HydroxyEncoding lcbHydroxyEncoding, boolean isAlexOhEncodingname, LipidomicsConstants lipidomicsConstants) throws LipidCombinameEncodingException{
    Vector<FattyAcidVO> chains = new Vector<FattyAcidVO>();
    String toSplit = humanReadable.replaceAll(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS, LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    String[] hrChains = toSplit.split(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    for (String hrChain : hrChains) {
      chains.add(decodeHumanReadableChain(hrChain,faHydroxyEncoding,lcbHydroxyEncoding,isAlexOhEncodingname,lipidomicsConstants));
    }
    return chains;
  }
  
  /**
   * Removes all double bond position annotations from a human readable String denoting a lipid species
   * This method uses a regex pattern including LipidomicsConstants.OMEGA_POSITION_START and LipidomicsConstants.OMEGA_POSITION_END,
   * ensure the pattern still works if these constants change!
   * @param doubleBondPositionsHumanReadable Human readable String
   * @return Human readable String without double bond position annotations
   */
  public static String getHumanReadableWODoubleBondPositions(String doubleBondPositionsHumanReadable) {
    String regex = String.format("\\%s.*?\\%s", LipidomicsConstants.OMEGA_POSITION_START, LipidomicsConstants.OMEGA_POSITION_END);
    return doubleBondPositionsHumanReadable.replaceAll(regex, "");
  }
  
  /**
   * Removes sn position annotations from a human readable String denoting a lipid species
   * sn positions need to be in brackets! ' (sn-' and ')'
   * @param snPositionsHumanReadable Human readable String
   * @return Human readable String without sn position annotations
   */
  public static String removeSNPositions(String snPositionsHumanReadable) {
	  String regex = "\\s\\" + LipidomicsConstants.SN_POSITION_START + ".*?\\" + LipidomicsConstants.SN_POSITION_END;
    return snPositionsHumanReadable.replaceAll(regex, "");
  }
  
  /**
   * Removes modification from a human readable String denoting a lipid species
   * @param snPositionsHumanReadable Human readable String
   * @return Human readable String without sn position annotations
   */
  public static String removeModification(String withModification) {
	  String regex = "\\;.*?[_\\/\\d$]";
	  if (withModification.indexOf(LipidomicsConstants.CHAIN_MOD_SEPARATOR) != -1)
		  return withModification.replaceAll(regex, "");
	  else
		  return withModification;
  }  
  
  /**
   * decodes a human readable chain
   * @param humanReadable the human readable lipid molecular species
   * @param faHydroxyEncoding the OH encodings of the FA moiety
   * @param lcbHydroxyEncoding the OH encodings of the LCB moiety
   * @param isAlexOhEncoding is this a hit of the Alex123 target lists that contains an OH encoding
   * @param lipidomicsConstants constants of the read result file (call trace starting from LDAResultReader), null if called from elsewhere
   * @return the decoded FA chain
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static FattyAcidVO decodeHumanReadableChain(String humanReadable, HydroxyEncoding faHydroxyEncoding,
      HydroxyEncoding lcbHydroxyEncoding, boolean isAlexOhEncodingname, LipidomicsConstants lipidomicsConstants) throws LipidCombinameEncodingException {
	  
  	humanReadable = removeSNPositions(humanReadable);
	  
    short chainType = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
    String prefix = "";
    int cAtoms = -1;
    String dbsPart;
    int dbs = -1;
    int oh = 0;
    int omegaPos = -1;
    String rest = humanReadable;
    String oxState = "";
    if (rest.startsWith(LipidomicsConstants.ALKYL_PREFIX)) {
      rest = rest.substring(LipidomicsConstants.ALKYL_PREFIX.length());
      chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
    }else if (rest.startsWith(LipidomicsConstants.ALKENYL_PREFIX)) {
      rest = rest.substring(LipidomicsConstants.ALKENYL_PREFIX.length());
      chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
    }
    
    char[] charsOfRest = rest.toCharArray();
    int i=0;
    while (!Character.isDigit(charsOfRest[i]) && i < charsOfRest.length)i++;
    String beforeNumberString = rest.substring(0,i);
    rest = rest.substring(i);
    if (beforeNumberString.length()>0) {
      //the beforeNumberString may contain a prefix, an OH encoding or both
      //since the OH encoding comes first, all sizes of the String have to be tested
      //for an existing hydroxy encoding
      int endIndex = beforeNumberString.length();
      String substring;
      while (endIndex!=0) {
        substring = beforeNumberString.substring(0, endIndex);
        try {
          oh = faHydroxyEncoding.getHydroxyNumber(substring);
          //the chain type does not need to be set, since the CHAIN_TYPE_FA_ACYL is anyway the default one
          //and the CHAIN_TYPE_FA_ALKYL and CHAIN_TYPE_FA_ALKENYL are searched for before
          break;
        }
        catch (HydroxylationEncodingException e) {}
        try {
          oh = lcbHydroxyEncoding.getHydroxyNumber(substring);
          if (chainType==LipidomicsConstants.CHAIN_TYPE_FA_ALKYL || chainType==LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL)
            throw new LipidCombinameEncodingException("The chain-id: \""+humanReadable+"\" cannot be decoded");
          chainType = LipidomicsConstants.CHAIN_TYPE_LCB;
          break;
        }
        catch (HydroxylationEncodingException e) {}
        endIndex--;
      }
      prefix = beforeNumberString.substring(endIndex);
    }
    if (isAlexOhEncodingname && rest.indexOf(LipidomicsConstants.ALEX_OH_SEPARATOR)!=-1) {
      try {
        oh = Integer.parseInt(rest.substring(rest.lastIndexOf(LipidomicsConstants.ALEX_OH_SEPARATOR)+1));
        rest = rest.substring(0,rest.lastIndexOf(LipidomicsConstants.ALEX_OH_SEPARATOR));
      }catch(NumberFormatException nfx) {
        throw new LipidCombinameEncodingException("The chain-id: \""+humanReadable+"\" cannot be decoded");
      }
    }
    
    // the chain modification separator changed, so we have to choose depending on version
    String chainModSeparator = LipidomicsConstants.CHAIN_MOD_SEPARATOR;
    if (lipidomicsConstants != null)
    	chainModSeparator = lipidomicsConstants.getModSeparatorDependingOnVersion();
    	
    if (rest.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_DBS)!=-1) {
      try {
    	dbsPart = rest.substring(rest.lastIndexOf(LipidomicsConstants.CHAIN_SEPARATOR_DBS)+1);
		if (rest.indexOf(chainModSeparator)!=-1){
          oxState = dbsPart.substring(dbsPart.indexOf(chainModSeparator)+1);
          dbsPart = dbsPart.substring(0,dbsPart.indexOf(chainModSeparator));
		}
		if (dbsPart.indexOf(LipidomicsConstants.OMEGA_POSITION_START)!=-1 && dbsPart.indexOf(LipidomicsConstants.OMEGA_POSITION_END)!=-1) {
        dbs = Integer.parseInt(dbsPart.substring(0,dbsPart.indexOf(LipidomicsConstants.OMEGA_POSITION_START)));
          omegaPos = Integer.parseInt(dbsPart.substring(
          		dbsPart.indexOf(LipidomicsConstants.OMEGA_POSITION_START)+LipidomicsConstants.OMEGA_POSITION_START.length(),
          		dbsPart.lastIndexOf(LipidomicsConstants.OMEGA_POSITION_END)));
        }else {
          dbs = Integer.parseInt(dbsPart);
        }
        rest = rest.substring(0,rest.lastIndexOf(LipidomicsConstants.CHAIN_SEPARATOR_DBS));
      }catch(NumberFormatException nfx) {
        throw new LipidCombinameEncodingException("The chain-id: \""+humanReadable+"\" cannot be decoded");
      }
    }
    try {
      cAtoms = Integer.parseInt(rest);
    }catch(NumberFormatException nfx) {
      throw new LipidCombinameEncodingException("The chain-id: \""+humanReadable+"\" cannot be decoded");
    }
    FattyAcidVO fa = new FattyAcidVO(chainType,prefix,cAtoms,dbs,oh,-1,null,oxState);
    fa.setOmegaPosition(omegaPos);
    return fa;
  }
  
  
  /**
   * returns an Object[] containing the prefix at Object[0] number of C atoms at Object[1] and the number of double bonds at Object[2]
   * @param chainId the string representative of a chain including the prefix
   * @return Object[] containing the prefix at Object[0] number of C atoms at Object[1] and the number of double bonds at Object[2]
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (not containing type and OH number) cannot be decoded
   */
  public static Object[] parsePrefixCAndDbsFromChainId(String chainId) throws LipidCombinameEncodingException {
    Object[] prefixCAndDbs = new Object[5];
    String rest = chainId;
    String prefix = "";
    //the chain type is extracted for backward compatibility - but it cannot distinguish between acyl and LCB chains
    short chainType = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
    if (chainId.startsWith(LipidomicsConstants.ALKYL_PREFIX)) {
      chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
      rest = rest.substring(LipidomicsConstants.ALKYL_PREFIX.length());
    }else if (chainId.startsWith(LipidomicsConstants.ALKENYL_PREFIX)) {
      chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
      rest = rest.substring(LipidomicsConstants.ALKENYL_PREFIX.length());
    }
    char[] chars = rest.toCharArray();
    int i=0;
    while (!Character.isDigit(chars[i]))i++;
    prefix = rest.substring(0,i);
    rest = rest.substring(i);
    String[] cAndDbAndOx = null;  
    try {cAndDbAndOx = parseCAndDbsFromChainId(rest);
    }catch (Exception e) {
      throw new LipidCombinameEncodingException("The chain-id: \""+chainId+"\" cannot be decoded");
    }
    prefixCAndDbs[0] = prefix;
    prefixCAndDbs[1] = Integer.parseInt(cAndDbAndOx[0]);
    prefixCAndDbs[2] = Integer.parseInt(cAndDbAndOx[1]);
    prefixCAndDbs[3] = chainType;
    prefixCAndDbs[4] = cAndDbAndOx[2];
    return prefixCAndDbs;
  }
  
  /**
   * decodes an Alex123 chain entry to the corresponding LDA representation
   * @param encoded the encoded Alex123 chain entry
   * @param species the the species name - for error handling
   * @return to the corresponding LDA chain representation
   * @throws LipidCombinameEncodingException thrown whenever there is something wrong with the entries
   */
  public static FattyAcidVO decodeAlex123Chain(String encoded, String species) throws LipidCombinameEncodingException {
    short type = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
    int oh = 0;
    int cs = 0;
    int dbs = 0;
    String prefix = "";
    String part = new String(encoded);
    if (part.startsWith(Settings.getInternalStandardDefaultInput())) {
      prefix = Settings.getInternalStandardDefaultInput();
      part = part.substring(Settings.getInternalStandardDefaultInput().length());
    }
    //this is for structure only
    if (part.startsWith("FA ")) {
      type = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
      part = part.substring("FA ".length());
    }
    if (part.startsWith("LCB ")) {
      type = LipidomicsConstants.CHAIN_TYPE_LCB;
      part = part.substring("LCB ".length());
    }
    if (part.startsWith(LipidomicsConstants.ALEX_ALKYL_PREFIX)) {
      type = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
      part = part.substring(LipidomicsConstants.ALEX_ALKYL_PREFIX.length());
    }
    if (part.startsWith(LipidomicsConstants.ALKYL_PREFIX)) {
      type = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
      part = part.substring(LipidomicsConstants.ALKYL_PREFIX.length());
    }
    if (part.startsWith(LipidomicsConstants.ALEX_ALKENYL_PREFIX)) {
      type = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
      part = part.substring(LipidomicsConstants.ALEX_ALKENYL_PREFIX.length());
    }
    if (part.startsWith(LipidomicsConstants.ALKENYL_PREFIX)) {
      type = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
      part = part.substring(LipidomicsConstants.ALKENYL_PREFIX.length());
    }

    if (part.indexOf("+")!=-1) {
      prefix += part.substring(part.indexOf("+"));
      char[] chars = prefix.toCharArray();
      int end = chars.length;
      while (end>0 && Character.isDigit(chars[end-1]))
        end--;
      prefix += prefix.substring(0,end);
      part = part.substring(0,part.indexOf("+"));
    }
    if (part.indexOf(LipidomicsConstants.ALEX_OH_SEPARATOR)!=-1) {
      type = LipidomicsConstants.CHAIN_TYPE_LCB;
      try {
        oh = Integer.parseInt(part.substring(part.indexOf(LipidomicsConstants.ALEX_OH_SEPARATOR)+1)); 
        part = part.substring(0, part.indexOf(LipidomicsConstants.ALEX_OH_SEPARATOR));
      }catch(NumberFormatException nfx) {
        throw new LipidCombinameEncodingException("The molecular species \""+encoded+"\" of the species \""+species+"\" cannot be decoded");
      }
    }
    if (part.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_DBS)!=-1) {
      try {
        dbs = Integer.parseInt(part.substring(part.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_DBS)+1)); 
        part = part.substring(0, part.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_DBS));
      }catch(NumberFormatException nfx) {
        throw new LipidCombinameEncodingException("The molecular species \""+encoded+"\" of the species \""+species+"\" cannot be decoded");
      }
    }
    try {cs = Integer.parseInt(part); 
    }catch(NumberFormatException nfx) {
      throw new LipidCombinameEncodingException("The molecular species \""+encoded+"\" of the species \""+species+"\" cannot be decoded");
    }
    return new FattyAcidVO(type, prefix, cs, dbs, oh, -1, null,"");
  }

  /**
   * encodes the LDA chain objects in an Alex123 compliant name
   * @param className the analyte class name
   * @param chains the chain objects
   * @return the Alex123 readable name
   */
  public static String encodeAlexMolSpeciesName(String className, Vector<FattyAcidVO> chains) {
    StringBuilder sb = new StringBuilder();
    sb.append(className+" ");
    boolean first = true;
    for (FattyAcidVO chain : chains) {
      if (!first) sb.append(LipidomicsConstants.ALEX_CHAIN_SEPARATOR);
      else first = false;
      sb.append(chain.getCarbonDbsId());
      if (chain.getOhNumber()>0)
        sb.append(LipidomicsConstants.ALEX_OH_SEPARATOR+String.valueOf(chain.getOhNumber()));
    }
    
    return sb.toString();
  }
  
  /**
   * generates vectors of potential FattyAcidVO combinations for checking of an intensity rule
   * @param rule the intensity rule
   * @param fas the available chain objects
   * @return vectors of potential FattyAcidVO combinations for checking of an intensity rule
   */
  public static Vector<Vector<FattyAcidVO>> getAllPotentialChainCombinationForThisRule(IntensityRuleVO rule, Vector<FattyAcidVO> fas) {
    Vector<Vector<FattyAcidVO>> combis = new Vector<Vector<FattyAcidVO>>();
    Set<Short> chainTypesOfRule = rule.getAvailableTypes();
    for (Short type : chainTypesOfRule) {
      if (type==LipidomicsConstants.CHAIN_TYPE_NO_CHAIN)
        continue;
      combis = addChainsOfOneType(type,fas,combis);
    }
    return combis;
  }
  
  /**
   * adds chains of one type to the corresponding combinations
   * @param type the chain type to add
   * @param chains the available chains
   * @param oldCombis the combis that were already generated
   * @return the possible combinations
   */
  private static Vector<Vector<FattyAcidVO>> addChainsOfOneType(short type, Vector<FattyAcidVO> chains, Vector<Vector<FattyAcidVO>> oldCombis){
    Vector<Vector<FattyAcidVO>> combis = new Vector<Vector<FattyAcidVO>>();
    Vector<FattyAcidVO> sameType = new Vector<FattyAcidVO>();
    for (FattyAcidVO fa : chains) {
      if (fa.getChainType()==type) {
        boolean add  = true;
        for (FattyAcidVO otherFa : sameType) {
          if (otherFa.getChainId().equalsIgnoreCase(fa.getChainId())) {
            add = false;
            break;
          }
        }
        if (add)
          sameType.add(fa);
      }
    }
    
    //if there is nothing to add for this type, simply return the previous results
    if (sameType.size()==0)
      return oldCombis;
    
    if (oldCombis.size()==0) {
      for (FattyAcidVO fa : sameType) {
        Vector<FattyAcidVO> combiFas = new Vector<FattyAcidVO>();
        combiFas.add(fa);
        combis.add(combiFas);
      }
    } else {
      for (Vector<FattyAcidVO> oldFas : oldCombis) {
        for (FattyAcidVO fa : sameType) {
          Vector<FattyAcidVO> combiFas = new Vector<FattyAcidVO>(oldFas);
          combiFas.add(fa);
          combis.add(combiFas);
        }
      }
    }
    
    return combis;
  }
  
  /**
   * returns all found chain combinations where this rule has an effect on
   * @param validChainCombinations the found chain combinations
   * @param chainsOfRule the chains that are contained in this rule
   * @return all found chain combinations where this rule has an effect on
   * @throws LipidCombinameEncodingException thrown whenever there is something wrong with the entries
   */
  public static Vector<String> getAllAffectedChainCombinations(Vector<String> validChainCombinations, Vector<String> chainsOfRule) throws LipidCombinameEncodingException{
    Hashtable<String,Integer> chainAmounts;
    int amount;
    Vector<String> affectedCombinations = new Vector<String>();
    for (String id : validChainCombinations) {
      //first: build a hash table which and how many chains are present for this combinations
      chainAmounts = new Hashtable<String,Integer>();
      for (String chain : splitChainCombiToEncodedStrings(id,LipidomicsConstants.CHAIN_COMBI_SEPARATOR)) {
        amount = chainAmounts.containsKey(chain) ? chainAmounts.get(chain) : 0;
        amount++;
        chainAmounts.put(chain, amount);
      }
      //second: check whether all of these chains can be found in this combination
      boolean affected = true;
      for (String chain : chainsOfRule) {
        if (chainAmounts.containsKey(chain)) {
          amount = chainAmounts.get(chain);
          amount--;
          if (amount<1)
            chainAmounts.remove(chain);
          else
            chainAmounts.put(chain, amount);
        }else {
          affected = false;
          break;
        }
      }
      if (affected)
        affectedCombinations.add(id);
    }
    return affectedCombinations;  
  }
  
  @SuppressWarnings("unchecked")
  public static Vector<String> sortChainCombinations(Set<String> unsorted) throws LipidCombinameEncodingException{
    if (unsorted.size()<2)
      return new Vector<String>(unsorted);
    Hashtable<String,Vector<FattyAcidVO>> unsortedHash = new Hashtable<String,Vector<FattyAcidVO>>();
    String fa1Encoded;
    Vector<FattyAcidVO> unsortedFAs;
    Vector<FattyAcidVO> sortedFAs;
    int count;
    for (String combi : unsorted) {
      unsortedHash.put(combi, decodeLipidNamesFromChainCombi(combi));
    }
    String encodedFA = null;
    int nrOfChains = unsortedHash.values().iterator().next().size();
    int[][] positions = new int[unsorted.size()][nrOfChains-1];
    Hashtable<Integer,Integer> nrOfOptionsEachLevel = new Hashtable<Integer,Integer>();
    for (int i=0; i!=(nrOfChains-1); i++) {
      Hashtable<String,String> there = new Hashtable<String,String>();
      unsortedFAs = new Vector<FattyAcidVO>();
      for (Vector<FattyAcidVO> fas : unsortedHash.values()) {
        encodedFA = encodeLipidNameForCreatingCombis(fas.get(i),true,true);
        if (there.containsKey(encodedFA))
          continue;
        there.put(encodedFA, encodedFA);
        unsortedFAs.add(fas.get(i));
      }
      sortedFAs = sortChainVOs(unsortedFAs);
      nrOfOptionsEachLevel.put(i, sortedFAs.size());
      count = 0;
      for (Vector<FattyAcidVO> fas : unsortedHash.values()) {
        fa1Encoded = encodeLipidNameForCreatingCombis(fas.get(i),true,true);
        for (int j=0; j!=sortedFAs.size(); j++) {
          if (fa1Encoded.equalsIgnoreCase(encodeLipidNameForCreatingCombis(sortedFAs.get(j),true,true)))
            positions[count][i] = j;
        }
        count++;
      }
    }
    count = 0;
    List<IntegerStringVO> toSort = new ArrayList<IntegerStringVO>();
    for (String combi : unsortedHash.keySet()) {
      int multFactor = 1;
      int totalValue = 0;
      for (int i=(nrOfChains-2); i!=-1; i--) {
        totalValue += positions[count][i]*multFactor;
        multFactor = multFactor*nrOfOptionsEachLevel.get(i);
      }
      toSort.add(new IntegerStringVO(combi,totalValue));
      count++;
    }
    Collections.sort(toSort,new GeneralComparator("at.tugraz.genome.lda.vos.IntegerStringVO", "getValue", "java.lang.Integer"));
//    Hashtable<String,String> keyToCombi = new Hashtable<String,String>();
//    count = 0;
//    for (String combi : unsorted) {
//      StringBuilder key = new StringBuilder();
//      for (int i=0; i!=positions[count].length; i++) {
//        if (key.length()!=0)
//          key.append("-");
//        key.append(String.valueOf(positions[count][i]));
//      }
//      keyToCombi.put(key.toString(), combi);
//      count++;
//    }
//    Vector<String> sorted = new Vector<String>();
//    
//    for (int i=0; i!=(nrOfChains-1); i++) {
//      sdfa
//    }
    Vector<String> sorted = new Vector<String>();
    for (IntegerStringVO vo : toSort)
      sorted.add(vo.getKey());
    return sorted;
  }
  
  
  
  /**
   * takes a list of unsorted chain (or species) names and sorts them in ascending carbon number, followed by ascending double bond, and followed by ascending OH number
   * @param unsorted unsorted chain (or species) names
   * @param ohInCombi are there hydroxylation encodings in the chain names
   * @param excludePrefix if true return sorted chain names without chain prefixes
   * @return sorted list of chain (or species) names
   * @throws LipidCombinameEncodingException thrown whenever there is something wrong with the hydroxylation encodings
   */
  public static Vector<String> sortChainNames(Vector<String> unsorted, boolean ohInCombi, boolean excludePrefix) throws LipidCombinameEncodingException{
    Vector<FattyAcidVO> fas = new Vector<FattyAcidVO>();
    for (String name : unsorted){
    	FattyAcidVO fa = decodeHumanReadableChain(name,Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(),false,null);
    	if (excludePrefix)
    		fa.setPrefix("");
    	fas.add(fa);
    } 
    Vector<FattyAcidVO> sortedChains = sortChainVOs(fas);
    Vector<String> sorted = new Vector<String>();
    for (FattyAcidVO vo : sortedChains)
      sorted.add(StaticUtils.getHumanReadableChainName(vo,Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(),ohInCombi));
    return sorted;
  }
  
  
  /**
   * TODO: this would be more efficient with collections sorting methods
   * takes a list of unsorted chains and sorts them in ascending carbon number, followed by ascending double bond, and followed by ascending OH number
   * @param unsorted vector containing FattyAcidVOs
   * @return sorted list of chains
   */
  public static Vector<FattyAcidVO> sortChainVOs(Vector<FattyAcidVO> unsorted){
    Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<Short,Hashtable<String,FattyAcidVO>>>>>> hash = new Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<Short,Hashtable<String,FattyAcidVO>>>>>>();
    Hashtable<String,Integer> frequencyOfFA = new Hashtable<String,Integer>();
    //for building the hash
    String encoded;
    for (FattyAcidVO fa : unsorted) {
      encoded = encodeLipidNameForCreatingCombis(fa, true, true);
      if (!frequencyOfFA.containsKey(encoded))
        frequencyOfFA.put(encoded, 0);
      frequencyOfFA.put(encoded, (frequencyOfFA.get(encoded)+1));
      Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<Short,Hashtable<String,FattyAcidVO>>>>> sameCarbons = new Hashtable<Integer,Hashtable<Integer,Hashtable<Integer,Hashtable<Short,Hashtable<String,FattyAcidVO>>>>>();
      if (hash.containsKey(fa.getcAtoms())) sameCarbons = hash.get(fa.getcAtoms());
      Hashtable<Integer,Hashtable<Integer,Hashtable<Short,Hashtable<String,FattyAcidVO>>>> sameDbs = new Hashtable<Integer,Hashtable<Integer,Hashtable<Short,Hashtable<String,FattyAcidVO>>>>();
      if (sameCarbons.containsKey(fa.getDoubleBonds())) sameDbs = sameCarbons.get(fa.getDoubleBonds());
      Hashtable<Integer,Hashtable<Short,Hashtable<String,FattyAcidVO>>> sameOh = new Hashtable<Integer,Hashtable<Short,Hashtable<String,FattyAcidVO>>>();
      if (sameDbs.containsKey(fa.getOhNumber())) sameOh = sameDbs.get(fa.getOhNumber());
      Hashtable<Short,Hashtable<String,FattyAcidVO>> sameOmega = new Hashtable<Short,Hashtable<String,FattyAcidVO>>();
      if (sameOh.containsKey(fa.getOmegaPosition())) sameOmega = sameOh.get(fa.getOmegaPosition());
      Hashtable<String,FattyAcidVO> sameType = new Hashtable<String,FattyAcidVO>();
      if (sameOmega.containsKey(fa.getChainType())) sameType = sameOmega.get(fa.getChainType());
      sameType.put(fa.getPrefix(),fa);
      sameOmega.put(fa.getChainType(),sameType);
      sameOh.put(fa.getOmegaPosition(), sameOmega);
      sameDbs.put(fa.getOhNumber(), sameOh);
      sameCarbons.put(fa.getDoubleBonds(), sameDbs);
      hash.put(fa.getcAtoms(), sameCarbons);
    }
    //sorting
    Vector<FattyAcidVO> sorted = new Vector<FattyAcidVO>();
    List<Integer> carbons = new ArrayList<Integer>(hash.keySet());
    Collections.sort(carbons);
    FattyAcidVO fac;
    for (Integer carbon : carbons) {
      List<Integer> dbs = new ArrayList<Integer>(hash.get(carbon).keySet());
      Collections.sort(dbs);
      for (Integer db : dbs) {
        List<Integer> ohs = new ArrayList<Integer>(hash.get(carbon).get(db).keySet());
        Collections.sort(ohs);
        for (Integer oh : ohs) {
          List<Integer> sameOmega = new ArrayList<Integer>(hash.get(carbon).get(db).get(oh).keySet());
          Collections.sort(sameOmega);
          for (Integer omega : sameOmega) {
            List<Short> sameType = new ArrayList<Short>(hash.get(carbon).get(db).get(oh).get(omega).keySet());
            Collections.sort(sameType);
            for (int i=(sameType.size()-1); i!=-1; i--) {
              List<String> prefixes = new ArrayList<String>(hash.get(carbon).get(db).get(oh).get(omega).get(sameType.get(i)).keySet());
              if (prefixes.contains("")) {
                fac = hash.get(carbon).get(db).get(oh).get(omega).get(sameType.get(i)).get("");
                encoded = encodeLipidNameForCreatingCombis(fac, true, true);
                for (int j=0; j!=frequencyOfFA.get(encoded); j++)
                  sorted.add(fac);
                prefixes.remove("");
              }
              Collections.sort(prefixes);
              for (String prefix : prefixes) {
                fac = hash.get(carbon).get(db).get(oh).get(omega).get(sameType.get(i)).get(prefix);
                encoded = encodeLipidNameForCreatingCombis(fac, true, true);
                for (int j=0; j!=frequencyOfFA.get(encoded); j++)
                  sorted.add(fac);
              }
            }
          }
        }
      }
    }
    return sorted;

  }
  
  /**
   * adds the possible isotopic labels stored in the chain libraries to the provided hashes singleLabelLookup, availableLabels and availableSingleLabels;
   * @param labels the labels found in the mass lists or data
   * @param availableSingleLabels the names of the unique single lables
   * @param singleLabelLookup lookup from a (combined) label to the single label
   * @param availableLabels of how many single labels does this label consist of
   */
  public static void extractIsoLabelInformation(Set<String> labels, Set<String> availableSingleLabels, Hashtable<String,String> singleLabelLookup, Hashtable<String,Integer> availableLabels) {
    char[] labelChars;
    StringBuilder sb;
    char ch;
    int current;
    int iter;
    for (String label : labels) {
      if (availableLabels.containsKey(label))
        continue;
      labelChars = label.toCharArray();
      sb = new StringBuilder();
      for (int i=0; i!=label.length(); i++) {
        ch = labelChars[i];
        if (sb.length()==0) {
          sb.append(ch);
        }else {
          if (label.substring(i,label.length()).startsWith(sb.toString())) {
            current = i;
            iter = 1;
            while(label.substring(current,label.length()).startsWith(sb.toString())) {
              current +=  sb.length();
              iter++;
            }
            if (current==label.length()) {
              availableSingleLabels.add(sb.toString());
              singleLabelLookup.put(label,sb.toString());
              availableLabels.put(label, iter);
              break;
            }
          }
          sb.append(ch);
        }
      }
      if (!availableLabels.containsKey(label)) {
        availableSingleLabels.add(sb.toString());
        singleLabelLookup.put(label,sb.toString());
        availableLabels.put(label, 1);
      }
    }
  }

  
  /**
   * checks whether two formulas stored in hash table form are the same
   * TODO: the algorithm asks that both hash table sizes are the same - when one has the amount of zero would be the same as a missing element, but this use case is not covered
   * @param one hash table containing the first formula; key: chemical element; value: the number of elements
   * @param two hash table containing the second formula; key: chemical element; value: the number of elements
   * @return true when both formulas are the same
   */
  public static boolean isChemicalFormulaTheSame(Hashtable<String,Integer> one, Hashtable<String,Integer> two) {
    if (one.size()!=two.size())
      return false;
    for (String element : one.keySet()) {
      if (!two.containsKey(element) || one.get(element)!=two.get(element))
        return false;
    }
    return true;
  }
  

  /**
   * class that translates the Alex123 nomenclature to the one of LDA; return Object array: [0] String: LDA class name; [1] FattyAcidVO: LDA species decoded; [2] String: encoded LDA chain combi name
   * @param alexClass Alex123 class name
   * @param alexSpecies Alex123 species name
   * @param alexMolSpecies Alex123 molecular species name
   * @return object array containing LDA encodings: [0] String: LDA class name; [1] FattyAcidVO: LDA species decoded; [2] String: encoded LDA chain combi name
   * @throws LipidCombinameEncodingException thrown if the ALEX123 names could not be decoded
   */
  public static Object[] translateAlexNomenclatureToLDA(String alexClass, String alexSpecies, String alexMolSpecies) throws LipidCombinameEncodingException {
    Object[] classSpeciesMol = new Object[3];
    String className;
    String species = "";
    String molSpecies = null;
    className = alexClass;
    if (className.startsWith(LipidomicsConstants.ALEX_IS_PREFIX)) {
      className = className.substring(LipidomicsConstants.ALEX_IS_PREFIX.length());
      species = Settings.getInternalStandardDefaultInput();
    }
    if (alexClass.endsWith(" O-")) {
      className = "O-"+className.substring(0,alexClass.length()-" O-".length());
      species += alexSpecies.substring(alexClass.length()-" O-".length()+1);
      if (alexMolSpecies!=null && alexMolSpecies.length()>0)
        molSpecies = alexMolSpecies.substring(alexClass.length()-" O-".length()+1);
    } else if (alexClass.endsWith(" P-")) {
      className = "P-"+className.substring(0,alexClass.length()-" P-".length());
      species += alexSpecies.substring(alexClass.length()-" P-".length()+1);
      if (alexMolSpecies!=null && alexMolSpecies.length()>0)
        molSpecies = alexMolSpecies.substring(alexClass.length()-" P-".length()+1);
    } else {
      species += alexSpecies.substring(alexClass.length()+1);
      if (alexMolSpecies!=null && alexMolSpecies.length()>0)
        molSpecies = alexMolSpecies.substring(alexClass.length()+1);
    }
    FattyAcidVO speciesVO = StaticUtils.decodeAlex123Chain(species, species);
    
//    if (/*alexClass.endsWith(" O-")||alexClass.endsWith(" P-")*/ alexClass.equalsIgnoreCase("Cer"))
//      System.out.println("Species: "+speciesVO);
    classSpeciesMol[0] = className;
    classSpeciesMol[1] = speciesVO;
    if (alexMolSpecies!=null && alexMolSpecies.length()>0) {
      Vector<FattyAcidVO> chains = new Vector<FattyAcidVO>();
      Vector<String> parts = getAlexChainParts(molSpecies);
      for (String part : parts)
        chains.add(StaticUtils.decodeAlex123Chain(part,species));
      chains = StaticUtils.sortChainVOs(chains);
      classSpeciesMol[2] = StaticUtils.encodeLipidCombi(chains);
    }else {
      classSpeciesMol[2] = molSpecies;
    }
    return classSpeciesMol;
  }
  
  /**
   * splits a an Alex123 combination to the individual chains
   * @param combi the Alex123 combination name
   * @return a vector of the single chain identifiers
   */
  private static Vector<String> getAlexChainParts(String combi){
    Vector<String> parts = new Vector<String>();
    String rest = combi.replaceAll("/",LipidomicsConstants.ALEX_CHAIN_SEPARATOR);
    int cut = rest.length();
    while (rest.substring(0,cut).indexOf(LipidomicsConstants.ALEX_CHAIN_SEPARATOR)!=-1) {
      int last = rest.substring(0,cut).lastIndexOf(LipidomicsConstants.ALEX_CHAIN_SEPARATOR);
      if ((last+1)>=LipidomicsConstants.ALEX_ALKYL_PREFIX.length() &&
          LipidomicsConstants.ALEX_ALKYL_PREFIX.equalsIgnoreCase(rest.substring(last-LipidomicsConstants.ALEX_ALKYL_PREFIX.length()+LipidomicsConstants.ALEX_CHAIN_SEPARATOR.length(),last+LipidomicsConstants.ALEX_CHAIN_SEPARATOR.length()))) {
        cut = last;
      }else if ((last+1)>=LipidomicsConstants.ALKYL_PREFIX.length() &&
          LipidomicsConstants.ALKYL_PREFIX.equalsIgnoreCase(rest.substring(last-LipidomicsConstants.ALKYL_PREFIX.length()+LipidomicsConstants.ALEX_CHAIN_SEPARATOR.length(),last+LipidomicsConstants.ALEX_CHAIN_SEPARATOR.length()))) {
        cut = last;
      }else  if ((last+1)>=LipidomicsConstants.ALEX_ALKENYL_PREFIX.length() &&
          LipidomicsConstants.ALEX_ALKYL_PREFIX.equalsIgnoreCase(rest.substring(last-LipidomicsConstants.ALEX_ALKENYL_PREFIX.length()+LipidomicsConstants.ALEX_CHAIN_SEPARATOR.length(),last+LipidomicsConstants.ALEX_CHAIN_SEPARATOR.length()))) {
        cut = last;
      }else  if ((last+1)>=LipidomicsConstants.ALKENYL_PREFIX.length() &&
          LipidomicsConstants.ALKYL_PREFIX.equalsIgnoreCase(rest.substring(last-LipidomicsConstants.ALKENYL_PREFIX.length()+LipidomicsConstants.ALEX_CHAIN_SEPARATOR.length(),last+LipidomicsConstants.ALEX_CHAIN_SEPARATOR.length()))) {
        cut = last;
      }else {
        parts.add(0,rest.substring(last+1));
        rest = rest.substring(0,last);
        cut = rest.length();
      }
    }
    parts.add(0,rest);
    return parts;
  }
  
  public static String[] extractFormulaAndAdductName(String contents) {
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
        String afterChargeString = formNameString.substring(formNameString.indexOf("charge=")+"charge=".length());
        try{
          charge = String.valueOf(extractDigitsAfterEqualSign(afterChargeString));
        }catch (NumberFormatException nfx){
          System.out.println("Warning: The charge entry in the column header \""+contents+"\" is not integer format! Setting it to \"1\"");
        }
      }
      if (formNameString.indexOf("mult=")!=-1){
      	String afterMultiString = formNameString.substring(formNameString.indexOf("mult=")+"mult=".length());
        try{
          multi = String.valueOf(extractDigitsAfterEqualSign(afterMultiString));
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
  
  private static int extractDigitsAfterEqualSign(String afterString){
  	afterString.trim();
    return Math.abs(Integer.parseInt(afterString));
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  //TODO: ALl of this are additions for the omega assignment
  
  
  /**
   * Returns the display name for an FA/LCB encoded molecular species identification.
   * Sorting of the chains is done with the compareTo method in FattyAcidVO; only if @param chainPositionsFixed is false.
   * ATTENTION: the chain sorting will differ in some cases from sortFASequenceUnassigned, as more properties are taken into account.
   * The sorting will be more predictable however.
   * @param combiName 										the encoded name for the chain combination
   * @param faEncoding 										the hydroxylation encoding for FA chains
   * @param lcbEncoding 									the hydroxylation encoding for LCB chains
   * @param chainPositionsFixed 					whether FattyAcidVOs in the Vector are ordered according to known sn positions
   * @return the encoded human readable display name for a chain combination
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static String getHumanReadableCombiName(String combiName, HydroxyEncoding faEncoding,
      HydroxyEncoding lcbEncoding, boolean chainPositionsFixed) throws LipidCombinameEncodingException 
  {
    return getHumanReadableCombiName(StaticUtils.decodeLipidNamesFromChainCombi(combiName), faEncoding, lcbEncoding, chainPositionsFixed);
  }
  
  /**
   * Returns the display name for an FA/LCB encoded molecular species identification.
   * Sorting of the chains is done with the compareTo method in FattyAcidVO; only if @param chainPositionsFixed is false.
   * ATTENTION: the chain sorting will differ in some cases from sortFASequenceUnassigned, as more properties are taken into account.
   * The sorting will be more predictable however.
   * @param chains 									the decoded information about the molecular species
   * @param faEncoding 							the hydroxylation encoding for FA chains
   * @param lcbEncoding 						the hydroxylation encoding for LCB chains
   * @param chainPositionsFixed 		whether FattyAcidVOs in the Vector are ordered according to known sn positions
   * @return the encoded human readable display name for a chain combination
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static String getHumanReadableCombiName(Vector<FattyAcidVO> chains, HydroxyEncoding faEncoding,
      HydroxyEncoding lcbEncoding, boolean chainPositionsFixed) throws LipidCombinameEncodingException
  {
  	if (!chainPositionsFixed)
    {
    	Collections.sort(chains);
    }
  	StringBuilder combi = new StringBuilder();
    boolean ohPresent = areThereOhInCombi(chains);
    for (FattyAcidVO chain : chains) {
      if (combi.length()!=0) {
        if (!chainPositionsFixed) combi.append(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
        else combi.append(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS);
      }
      combi.append(getHumanReadableChainName(chain, faEncoding, lcbEncoding, ohPresent));
    }
    return combi.toString();
  }
  
  /**
   * Determines low, medium and high confidence peak ranges of @param. 
   * First, the lower distance from the peak maximum to the lower and upper border is determined.
   * The medium accuracy cutoff is 2/3 of the shorter distance.
   * The high accuracy cutoff is 1/3 of the shorter distance.
   * Beyond 2/3 of the shorter distance, but within the range of the peak is the low accuracy cutoff.
   * The User provided value for the minimum threshold for high confidence RT matches overrides these predefined cutoffs.
   * @param param
   * @return
   */
  public static Range[] determinePeakRanges(LipidParameterSet param) {
    float lowerValley = param.getIsotopicProbes().get(0).get(0).LowerValley;
    float upperValley = param.getIsotopicProbes().get(0).get(0).UpperValley;
    float peak = param.getIsotopicProbes().get(0).get(0).Peak;
    float shorterDistance = (Math.abs(peak-lowerValley) <= Math.abs(peak-upperValley)) ? Math.abs(peak-lowerValley) : Math.abs(peak-upperValley);
    int n = 3;
    int minimumThreshold = LipidomicsConstants.getMinimumThresholdForHighConfidenceRTMatch();
    float mediumAccuracyCutoff = shorterDistance/n*2;
    float highAccuracyCutoff = shorterDistance/n;
    if (shorterDistance < minimumThreshold*n) { 
    	highAccuracyCutoff = minimumThreshold;
      if (shorterDistance > minimumThreshold) {
        mediumAccuracyCutoff = shorterDistance;
      } else {
      	mediumAccuracyCutoff = minimumThreshold;
      }
    }
    Range[] peakRanges = new Range[3];
    peakRanges[0] = new Range(
        lowerValley/60f,
        upperValley/60f);
    peakRanges[1] = new Range(
        (peak-mediumAccuracyCutoff)/60f, 
        (peak+mediumAccuracyCutoff)/60f);
    peakRanges[2] = new Range(
        (peak-highAccuracyCutoff)/60f, 
        (peak+highAccuracyCutoff)/60f);
    return peakRanges;
  }
  
  /**
   * @param paramOmegaInfo 
   * @return Vector of assigned DoubleBondPositionVOs
   */
  public static Vector<DoubleBondPositionVO> getAssignedDoubleBondPositions(Vector<DoubleBondPositionVO> paramOmegaInfo) {
    Vector<DoubleBondPositionVO> assignedHits = new Vector<DoubleBondPositionVO>();
    for (int i = 0; i < paramOmegaInfo.size(); i++) {
      if (paramOmegaInfo.get(i).getIsAssigned()) {
        assignedHits.add(paramOmegaInfo.get(i));
      }
    }
    return assignedHits;
  }
  
  /**
   * @param paramOmegaInfo 
   * @return Vector of DoubleBondPositionVOs with desired level of accuracy
   */
  public static Vector<DoubleBondPositionVO> getHighAccuracyDoubleBondPositions(Vector<DoubleBondPositionVO> paramOmegaInfo) {
    Vector<DoubleBondPositionVO> highAccuracyHits = new Vector<DoubleBondPositionVO>();
    for (int i = 0; i < paramOmegaInfo.size(); i++) {
      if (paramOmegaInfo.get(i).getAccuracy() > 1) {
        highAccuracyHits.add(paramOmegaInfo.get(i));
      }
    }
    return highAccuracyHits;
  }
  
  /**
   * @param highAccuracyHits sorted DoubleBondPositionVOs with desired level of accuracy
   * @return Vector of DoubleBondPositionVOs which can be unambiguously assigned to the detected peak
   */
  public static Vector<DoubleBondPositionVO> findUnambiguousDoubleBondPositions(Vector<DoubleBondPositionVO> highAccuracyHits) {
    Vector<DoubleBondPositionVO> assignedHits = new Vector<DoubleBondPositionVO>();
    
    int n = highAccuracyHits.size();
    
    if (n == 1) {
      assignedHits.add(highAccuracyHits.get(0));
      return assignedHits;
    } else if (n == 2) {
      if (!highAccuracyHits.get(0).getMolecularSpecies().equals(highAccuracyHits.get(1).getMolecularSpecies())) {
        assignedHits.add(highAccuracyHits.get(0));
        assignedHits.add(highAccuracyHits.get(1));
        return assignedHits;
      }
    } else if (n > 2) {
      if (!highAccuracyHits.get(0).getMolecularSpecies().equals(highAccuracyHits.get(1).getMolecularSpecies())) {
        assignedHits.add(highAccuracyHits.get(0));
      }
      
      for (int i = 1; i < n-1; i++) {
        if (!highAccuracyHits.get(i).getMolecularSpecies().equals(highAccuracyHits.get(i-1).getMolecularSpecies()) &&
            !highAccuracyHits.get(i).getMolecularSpecies().equals(highAccuracyHits.get(i+1).getMolecularSpecies())) {
          assignedHits.add(highAccuracyHits.get(i));
        }
      }
      
      if (!highAccuracyHits.get(n-2).getMolecularSpecies().equals(highAccuracyHits.get(n-1).getMolecularSpecies())) {
        assignedHits.add(highAccuracyHits.get(n-1));
      }
    }
    return assignedHits;
  }
  
  /**
   * Finds C=C to be assigned for each identification
   * First divide into molecular species groups, then subdivide into high, medium and low accuracy hits.
   * If there is only one high accuracy hit: assign. If there are more, do not assign any C=C position.
   * If there are no high accuracy hits and only one other hit, that falls within the appropriate threshold: assign. Else, do not assign any C=C position.
   * 
   * @param param
   * @return
   */
  public static Vector<DoubleBondPositionVO> findUnambiguousDoubleBondPositionsNew(LipidParameterSet param)
  {
  	Vector<DoubleBondPositionVO> assignedHits = new Vector<DoubleBondPositionVO>();
  	LinkedHashMap<String,LinkedHashMap<Integer,ArrayList<DoubleBondPositionVO>>> container = sortDoubleBondPositionVOsOfID(param);
  	
  	for (String molecularSpecies : container.keySet())
  	{
  		LinkedHashMap<Integer,ArrayList<DoubleBondPositionVO>> speciesOfHit = container.get(molecularSpecies);
  		ArrayList<DoubleBondPositionVO> highAccuracy = speciesOfHit.get(DoubleBondPositionVO.ACCURACY_HIGH);
  		if (highAccuracy.size() > 1)
  		{
  			DoubleBondPositionVO combined = combineUnambiguousDoubleBondPositions(highAccuracy, param.getPreciseRT());
  			if (combined != null)
  			{
  				highAccuracy.removeAll(highAccuracy);
  				highAccuracy.add(combined);
  				param.addOmegaInformation(combined);
  			}
  		}
  		if (highAccuracy.size() == 1)
  		{
  			assignedHits.add(highAccuracy.get(0));
  			continue;
  		}
  		
  		if (!highAccuracy.isEmpty()) continue;
  		
  		ArrayList<DoubleBondPositionVO> intermediateConfidence = speciesOfHit.get(DoubleBondPositionVO.ACCURACY_MEDIUM);
  		ArrayList<DoubleBondPositionVO> notHighAccuracy = new ArrayList<DoubleBondPositionVO>(intermediateConfidence);
  		notHighAccuracy.addAll(speciesOfHit.get(DoubleBondPositionVO.ACCURACY_LOW));
  		if (notHighAccuracy.size() > 1 && !intermediateConfidence.isEmpty())
  		{
  			intermediateConfidence.removeAll(intermediateConfidence);
  			DoubleBondPositionVO combined = combineUnambiguousDoubleBondPositions(notHighAccuracy, param.getPreciseRT());
  			if (combined != null)
  			{
  				intermediateConfidence.add(combined);
  				param.addOmegaInformation(combined);
  			}
  		}
  			
  		if (intermediateConfidence.size() == 1)
  		{
  			assignedHits.add(intermediateConfidence.get(0));
  		}
  	}
  	
  	return assignedHits;
  }
  
  private static DoubleBondPositionVO combineUnambiguousDoubleBondPositions(ArrayList<DoubleBondPositionVO> vos, Double rt)
  {
  	ArrayList<Vector<Integer>> patterns = new ArrayList<Vector<Integer>>();
  	for (DoubleBondPositionVO vo : vos)
  	{
  		patterns.add(vo.getPositionAssignmentPattern());
  	}
  	Vector<Integer> combinedPattern = computeConsensusPattern(patterns);
  	boolean anyPositionDefined = false;
  	for (Integer pos : combinedPattern)
  	{
  		if (pos > 0) anyPositionDefined = true;
  	}
  	DoubleBondPositionVO newVO = null;
  	if (anyPositionDefined)
  	{
  		DoubleBondPositionVO closest = vos.get(0);
  		Double closestDiff = Math.abs(closest.getExpectedRetentionTime() - rt);
  		for (DoubleBondPositionVO vo : vos)
  		{
  			if (Math.abs(vo.getExpectedRetentionTime() - rt) < closestDiff) closest = vo;
  		}
  		newVO = new DoubleBondPositionVO(closest);
  		Vector<FattyAcidVO> newChainCombination = new Vector<FattyAcidVO>();
  		for (Integer i=0; i<newVO.getChainCombination().size(); i++)
  		{
  			FattyAcidVO fattyAcid = newVO.getChainCombination().get(i);
  			fattyAcid.setOmegaPosition(combinedPattern.get(i));
  			newChainCombination.add(fattyAcid);
  		}
  		newVO.setChainCombination(newChainCombination);
  	}
  	
  	return newVO;
  }
  
  
  
  private static Vector<Integer> computeConsensusPattern(ArrayList<Vector<Integer>> patterns)
	{
		Vector<Integer> combinedPattern = new Vector<Integer>();
		LinkedHashMap<Integer,HashSet<Integer>> patternOverlap = new LinkedHashMap<Integer,HashSet<Integer>>();
		//compute pattern overlap
		for (Vector<Integer> pattern : patterns)
		{
			for (int j=0; j<pattern.size(); j++)
			{
				if (!patternOverlap.containsKey(j)) patternOverlap.put(j, new HashSet<Integer>());
				patternOverlap.get(j).add(pattern.get(j));
			}
		}
		
		for (int chainPos : patternOverlap.keySet())
		{
			ArrayList<Integer> consensus = new ArrayList<Integer>(patternOverlap.get(chainPos));
			int numberOfPos = 0;
			for (int pos : consensus)
			{
				if (pos > 0) numberOfPos++;
			}
			if (numberOfPos == 1)
			{
				if (consensus.contains(-1)) combinedPattern.add(chainPos, -1);
				else combinedPattern.add(chainPos, consensus.get(0));
			}
			else
			{
				combinedPattern.add(chainPos, -1);
			}
		}
		
		return combinedPattern;
	}
  
  /**
   * Sorts DoubleBondPositionVOs of an identification by their molecular species and their accuracy
   * @param param
   * @return
   */
  private static LinkedHashMap<String,LinkedHashMap<Integer,ArrayList<DoubleBondPositionVO>>> sortDoubleBondPositionVOsOfID(LipidParameterSet param)
  {
  	LinkedHashMap<String,LinkedHashMap<Integer,ArrayList<DoubleBondPositionVO>>> container = new LinkedHashMap<String,LinkedHashMap<Integer,ArrayList<DoubleBondPositionVO>>>();
  	
  	for (DoubleBondPositionVO vo : param.getOmegaInformation())
  	{
  		String molecularSpecies = vo.getMolecularSpecies();
  		if (!container.containsKey(molecularSpecies)) 
  		{
  			container.put(molecularSpecies, new LinkedHashMap<Integer,ArrayList<DoubleBondPositionVO>>());
  			container.get(molecularSpecies).put(DoubleBondPositionVO.ACCURACY_HIGH, new ArrayList<DoubleBondPositionVO>());
  			container.get(molecularSpecies).put(DoubleBondPositionVO.ACCURACY_MEDIUM, new ArrayList<DoubleBondPositionVO>());
  			container.get(molecularSpecies).put(DoubleBondPositionVO.ACCURACY_LOW, new ArrayList<DoubleBondPositionVO>());
  		}
  		Integer accuracy = vo.getAccuracy() == DoubleBondPositionVO.ACCURACY_LOW ? DoubleBondPositionVO.ACCURACY_MEDIUM : vo.getAccuracy();
  		
  		if (accuracy == DoubleBondPositionVO.ACCURACY_HIGH)
  		{
  			container.get(molecularSpecies).get(accuracy).add(vo);
  		}
  		else
  		{
  			if (Math.abs(param.getPreciseRT()-vo.getExpectedRetentionTime())*60 < LipidomicsConstants.getMaximumThresholdForIntermediateConfidenceRTMatch())
    		{
    			container.get(molecularSpecies).get(accuracy).add(vo);
    		}
  			else
  			{
  				container.get(molecularSpecies).get(DoubleBondPositionVO.ACCURACY_LOW).add(vo);
  			}
  		}
  	}
  	return container;
  }
  
  public static Set<String> getMolecularSpeciesSet(Vector<DoubleBondPositionVO> omegaInfo) {
    Set<String> molecularSpeciesSet = new HashSet<String>();
    for (DoubleBondPositionVO doubleBondPositionVO : omegaInfo) {
      molecularSpeciesSet.add(doubleBondPositionVO.getMolecularSpecies());
    }
    return molecularSpeciesSet;
  }
  
  //TODO: use internal position insensitive names, plus this method should only return one species IMO
  public static Vector<DoubleBondPositionVO> getDoubleBondAssignmentsOfMolecularSpecies(Vector<DoubleBondPositionVO> omegaInfo, String molecularSpecies) {
    Vector<DoubleBondPositionVO> doubleBondAssignmentsOfMolecularSpecies = new Vector<DoubleBondPositionVO>();
    for (DoubleBondPositionVO labeledChainCombiVO : omegaInfo) {
      if (labeledChainCombiVO.getMolecularSpecies().equals(molecularSpecies)) {
        doubleBondAssignmentsOfMolecularSpecies.add(labeledChainCombiVO);
      }
    }
    return doubleBondAssignmentsOfMolecularSpecies;
  }
  
  /**
   * Checks whether two human readable Strings, one with annotated double bond positions and the other without, describe the same molecular species 
   * @param doubleBondPositionsHumanReadable Human readable String of a molecular species with annotated double bond positions
   * @param molecularSpecies Human readable String of a molecular species without annotated double bond positions
   * @return boolean informing the caller whether the given chain combinations can be considered equivalent
   */
  public static boolean isChainCombinationEquivalent(String doubleBondPositionsHumanReadable, String molecularSpecies) {
    boolean isChainCombinationEquivalent = false;
    String noDoubleBonds = StaticUtils.getHumanReadableWODoubleBondPositions(doubleBondPositionsHumanReadable);
    String[] arrayNoDoubleBonds = splitChainCombinationsAtChainSeparators(noDoubleBonds);
    String[] arrayIdentifications = splitChainCombinationsAtChainSeparators(molecularSpecies);
    
    if (!(arrayNoDoubleBonds == null || arrayIdentifications == null)) {
      Arrays.sort(arrayNoDoubleBonds);
      Arrays.sort(arrayIdentifications);
      if (Arrays.equals(arrayNoDoubleBonds, arrayIdentifications)) {
        isChainCombinationEquivalent = true;
      }
    }
    return isChainCombinationEquivalent;
  }
  
  /**
   * Splits a String into substrings at CHAIN_SEPARATOR_NO_POS or CHAIN_SEPARATOR_KNOWN_POS as defined in LipidomicsConstants.
   * @param s String to be split
   * @return Array of Strings split at the defined positions
   */
  public static String[] splitChainCombinationsAtChainSeparators(String s) {
    String str = new String(s);
    String[] splitString = {str};
    if (str.contains(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) {
      splitString = str.split(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    } else if (str.contains(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS)) {
      splitString = str.split(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS);
    }
    return splitString;
  }
  
	public static float calculatedMzTolValue(float mz, float tolerance, short tolUnit) {
		float tol = 0f;
		if (tolUnit==LipidomicsConstants.MZUNIT_DA_ID)
			tol = tolerance;
		else if (tolUnit==LipidomicsConstants.MZUNIT_PPM_ID)
			tol = mz*0.000001f*tolerance;
		return tol;
	}
  
}
