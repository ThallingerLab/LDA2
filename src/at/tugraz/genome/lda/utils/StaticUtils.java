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

package at.tugraz.genome.lda.utils;

import java.awt.Color;
import java.awt.Component;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.math.BigDecimal;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.msn.FragmentCalculator;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.swing.RangeColor;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;
import at.tugraz.genome.maspectras.chromaviewer.MSMapViewer;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.utils.StringUtils;

/**
 * 
 * @author Juergen Hartler
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
  

  
  public final static String[] physicalMagnitudes_ = {"","m","\u03BC","n","p","f","a"};
  
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
    String basePath = null;
    File resultsFile = new File(resultFilePath);
    if (resultsFile.exists()&&resultsFile.isFile()){
      String directory = StaticUtils.extractDirName(resultFilePath);
      String chromFilePath = directory+File.separator+experimentName+".chrom";
      File chromFile = new File (chromFilePath);
      boolean chromFileExists = false;
      //This is excluded because the chrom file can now be a directory containing the chrom files
      if (chromFile.exists())//&&chromFile.isFile())
        chromFileExists = true;
      else{
        File dirFile = new File(directory);
        if (dirFile.exists()&&dirFile.isDirectory()){
          File[] filesInDir = dirFile.listFiles();
          for (int i=0;i!=filesInDir.length;i++){
            if (filesInDir[i].getName().indexOf(experimentName)==0&&filesInDir[i].getAbsolutePath().endsWith(".chrom")){
              if (!chromFileExists || filesInDir[i].getName().length()<chromFile.getName().length()) chromFile = filesInDir[i];
              chromFileExists = true;
            }
          }
        }
      }
      if (chromFileExists) basePath = chromFile.getAbsolutePath().substring(0,chromFile.getAbsolutePath().length()-".chrom".length());
    }
    return basePath;
  }
  
  public static String existChromNecessaryFiles(String chromBaseName){
    String[] chromPaths = StringUtils.getChromFilePaths(chromBaseName+".");
    String oneDoesNotExist = "";
    if (chromPaths[0]==null || !existsFile(chromPaths[0])){
      oneDoesNotExist = "The chrom file is missing";
    }
    if (chromPaths[1]==null || !existsFile(chromPaths[1])){
      oneDoesNotExist = "The head file is missing";
    }
    if (chromPaths[2]==null || !existsFile(chromPaths[2])){
      oneDoesNotExist = "The idx file is missing";
    }
    if (chromPaths[3]==null || !existsFile(chromPaths[3])){
      oneDoesNotExist = "The rtt file is missing";
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

  private static void addExcelLabel(CellStyle  cellStyle, Row row, int columnCount, String toAdd){
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
   */
  public static Hashtable<String,Integer> categorizeFormula(String formula) throws ChemicalFormulaException {
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
        addFormulaPartToHash(categorized, currentElement, currentAmountString, add);
        add = true;
        currentElement = "";
        currentAmountString = "";
      }else if (chars[i]=='+'){
        addFormulaPartToHash(categorized, currentElement, currentAmountString, add);
        add = true;
        currentElement = "";
        currentAmountString = "";
      }else if (chars[i]=='-'){
        addFormulaPartToHash(categorized, currentElement, currentAmountString, add);
        add = false;
        currentElement = "";
        currentAmountString = "";
      }else if (Character.isLetter(chars[i])){
        if (Character.isUpperCase(chars[i])){
          if (currentElement.length()==0) currentElement+=String.valueOf(chars[i]);
          else{
            addFormulaPartToHash(categorized, currentElement, currentAmountString, add);
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
    addFormulaPartToHash(categorized, currentElement, currentAmountString, add);
//    System.out.println("-------------------------------");
//    for (String elem : categorized.keySet()){
//      System.out.println(elem+": "+categorized.get(elem));
//    }
    return categorized;
  }
  
  private static void addFormulaPartToHash(Hashtable<String,Integer> categorized, String currentElement, String currentAmountString, boolean add){
    if (currentElement!=null&&currentElement.length()>0){
      int finalAmount = 0;
      if (categorized.containsKey(currentElement))
        finalAmount = categorized.get(currentElement);
      int amountToChange = 1;
      if (currentAmountString!=null && currentAmountString.length()>0) amountToChange = Integer.parseInt(currentAmountString);
      if (add) finalAmount += amountToChange;
      else finalAmount -= amountToChange;
      categorized.put(currentElement, finalAmount);
    }
  }
  
  public static String generateLipidNameString(String name, Integer doubleBonds){
    String nameString = name;
    if (doubleBonds!=null && doubleBonds>-1)
      nameString += ":"+String.valueOf(doubleBonds);
//    nameString += "_";
    return nameString;
  }
  
  public static String generateLipidNameString(String name, Integer doubleBonds, String rt){
    String nameString = generateLipidNameString(name,doubleBonds);
    if (rt!=null&&rt.length()>0) nameString+="_"+rt;
    return nameString;
  }
  
  
  public static boolean existsFile(String pathName){
    File file = new File(pathName);
    return file.exists();
  }
  
  public static boolean isWithinTolerance(double tolerance, double ref, double value){
    return ((ref-tolerance)<value && value<(ref+tolerance));
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
   * @return true if they are the permuted combinations of one another
   */
  public static boolean isAPermutedVersion(String one,String two){
    String[] fas = LipidomicsMSnSet.getFAsFromCombiName(one);
    Vector<String> fasVector = new Vector<String>();
    for (String fa : fas) fasVector.add(fa);
    Vector<String> permuted = FragmentCalculator.getPermutedNames(fasVector);
    for (String combi : permuted){
      if (combi.equalsIgnoreCase(two)) return true;
    }
    return false;  
  }
  
  /**
   * specifies m/z ranges for coloring from an MSn identification
   * @param paramUncasted the identification not casted in a LipidomicsMSnSet
   * @param selectedMSn shall only a certain MSn identifcation be displayed
   * @return m/z ranges for coloring from an MSn identification
   */
  @SuppressWarnings("unchecked")
  public static Vector<RangeColor> createRangeColorVOs(LipidParameterSet paramUncasted, String selectedMSn){
    String selected = null;
    if (selectedMSn!=null && selectedMSn.length()>0){
      selected = new String(selectedMSn);
      selected = (String)cleanEmptyFAPositions(selected)[0];
    }
    Vector<RangeColor> rangeColors = null;
    if (paramUncasted instanceof LipidomicsMSnSet){
      LipidomicsMSnSet param = (LipidomicsMSnSet)paramUncasted;
      rangeColors = new Vector<RangeColor>();
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
        rangeColors.add(vo);
        headFrags.remove(highestArea);
      }
      Hashtable<String,Hashtable<String,CgProbe>> chainFrags =  param.getChainFragments();
      Hashtable<String,Double> relAreas = param.getChainCombinationRelativeAreas();
      Vector<String> combiOrdered = new Vector<String>();
      // the bigger areas should come first
      for (String combi : relAreas.keySet()){
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
      // this is for the sequence of color selection
      Vector<String> identNames = new Vector<String>();
      for (Object object : param.getMSnIdentificationNames()){
        if (object instanceof String) identNames.add((String)object);
        else identNames.add(((Vector<String>)object).get(0));
      }
      int colorCount = 1;
      Hashtable<String,String> usedFAs = new Hashtable<String,String>();
      for (String combi : combiOrdered){
        if (selected!=null && selected.length()>0 && !StaticUtils.isAPermutedVersion(selected.replaceAll("/", "_"),combi)) continue;
        String name = combi;
        for (String oneOrder : identNames){
          if (StaticUtils.isAPermutedVersion(oneOrder.replaceAll("/", "_"),combi)){
            name = oneOrder.replaceAll("/", "_");
          }
        }
        String[] fas = LipidomicsMSnSet.getFAsFromCombiName(name);
        for (int i=0; i!= fas.length; i++){
          String fa = fas[i];
          if (!chainFrags.containsKey(fa)) continue;
          Hashtable<String,CgProbe> frags =  chainFrags.get(fa);
          Color color = getFragemntColor(colorCount);
          for (String key : frags.keySet()){
            String fragmentName = key+" ("+fa+")";
            if (usedFAs.containsKey(fragmentName))continue;
            usedFAs.put(fragmentName, fragmentName);
            CgProbe probe = frags.get(key);
            RangeColor vo = new RangeColor(fragmentName,color,probe,probe.Mz-probe.LowerMzBand,probe.Mz+probe.UpperMzBand);
            rangeColors.add(vo);
            
          }
          colorCount++;
        }
      }
    }
    return rangeColors;
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
   * cleans the empty positions of the final FA name and returns the amount of positions
   * @param faName
   * @return returns the combi name and the number of assignable positions
   */
  public static Object[] cleanEmptyFAPositions(String identification){
    Object[] result = new Object[2];
    int positions = 0;
    String name = identification.replaceAll("/", "_");
    while (name.startsWith("-_") || name.indexOf("_-_")!=-1 || name.endsWith("_-")){
      positions++;
      if (name.endsWith("_-")){
        name = name.substring(0,name.length()-"_-".length());
      } else if (name.startsWith("-_")){
        name = name.substring("-_".length());
      } else if (name.indexOf("_-_")!=-1){
        name = name.substring(0,name.indexOf("_-_"))+name.substring(name.indexOf("_-_")+"_-".length());
      }
    }
    String[] fas = LipidomicsMSnSet.getFAsFromCombiName(name);
    positions += fas.length;
    result[0] = name;
    result[1] = positions;
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

}