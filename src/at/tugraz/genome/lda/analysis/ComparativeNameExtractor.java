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
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.swing.JFrame;

import org.dhatim.fastexcel.reader.Cell;
import org.dhatim.fastexcel.reader.CellType;
import org.dhatim.fastexcel.reader.ReadableWorkbook;
import org.dhatim.fastexcel.reader.Row;
import org.dhatim.fastexcel.reader.Sheet;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.export.QuantificationResultExporter;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleStringVO;
import at.tugraz.genome.lda.vos.ResultFileVO;
import at.tugraz.genome.voutils.GeneralComparator;
import at.tugraz.genome.lda.parser.LDAResultReader;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
public class ComparativeNameExtractor extends ClassNamesExtractor
{
  protected String isSelectionPrefix_;
  protected String esSelectionPrefix_;
  // the first entry is the same from the beginning of the file; just the end is cut away
  // the second entry in the hash is just the difference between the names that should be displayed
  protected Hashtable<String,String> expNames_ = new  Hashtable<String,String>();
  protected Vector<String> expNamesInSequence_ = new  Vector<String>();
  protected Vector<String> groups_;
  protected Hashtable<String,Hashtable<String,String>> allISNames_;
  protected Hashtable<String,Hashtable<String,String>> allESNames_;
  protected Hashtable<String,Vector<File>> filesOfGroup_;
  protected Vector<ResultFileVO> resultFileVO_;
  protected Hashtable<String,File> expNameToFile_ = new Hashtable<String,File>();
  
  private Hashtable<String,Vector<String>> sortedISNames_;
  private Hashtable<String,Vector<String>> sortedESNames_;
  private Hashtable<String,String> lipidClassesHash_;
  
  protected Vector<String> expNameCut1_;
  protected Vector<String> expNameCut2_;
  
  int charsToCut_ = 0;
  int charsToCutPrev_ = 0;
  
  public ComparativeNameExtractor(Vector<File> resultFiles, String isSelectionPrefix, String esSelectionPrefix){
    super(resultFiles);
    if (isSelectionPrefix!=null && isSelectionPrefix.length()>0)
      this.isSelectionPrefix_ = isSelectionPrefix;
    else
      this.isSelectionPrefix_ = null;
    if (esSelectionPrefix!=null &&esSelectionPrefix.length()>0)
      this.esSelectionPrefix_ = esSelectionPrefix;
    else
      this.esSelectionPrefix_ = null;
    groups_ = null;
  }
  
  public ComparativeNameExtractor(Vector<File> resultFiles, String isSelectionPrefix, String esSelectionPrefix,Vector<String> groups, Hashtable<String,Vector<File>> filesOfGroup){
    this(resultFiles, isSelectionPrefix, esSelectionPrefix);
    this.groups_ = groups;
    this.filesOfGroup_ = filesOfGroup;  
  }
    
  protected void extractInformation(int statisticsViewMode, boolean combineOxWithNonOx) throws ExcelInputFileException, LipidCombinameEncodingException{
    expNameCut1_ = new  Vector<String>();
    expNameCut2_ = new  Vector<String>();
    if (resultFiles_.size()>1){
      Vector<String> expNameFull = new Vector<String>();
      for (File resultFile : resultFiles_){
        String fileName = StaticUtils.extractFileName(resultFile.getAbsolutePath());
        expNameFull.add(fileName.substring(0,fileName.lastIndexOf(".")));
      }
      int[] limits = StaticUtils.detectNameUnequalitiesBeforeAndAfter(expNameFull);
      charsToCutPrev_ = limits[0];
      charsToCut_ = limits[1];
      for (String expName : expNameFull){
        expNameCut1_.add(expName.substring(0,expName.length()-charsToCut_));
        expNameCut2_.add(expName.substring(charsToCutPrev_,expName.length()-charsToCut_));
      }
    }  
    allISNames_ = new Hashtable<String,Hashtable<String,String>>();
    allESNames_ = new Hashtable<String,Hashtable<String,String>>();
    sortedISNames_ = new Hashtable<String,Vector<String>>();
    sortedESNames_ = new Hashtable<String,Vector<String>>();
    lipidClassesHash_ = new Hashtable<String,String>();
    lipidClasses_ = new Vector<String>();
    resultFileVO_ = new Vector<ResultFileVO>();
    for (int i=0; i!=resultFiles_.size();i++){
      File resultFile = resultFiles_.get(i);
      String fileName = StaticUtils.extractFileName(resultFile.getAbsolutePath());
      fileName = fileName.substring(0,fileName.lastIndexOf("."));
      if (expNameCut1_.size()>0){
        //expNames_.add(expNameCut_.get(i));
        expNames_.put(expNameCut1_.get(i), expNameCut2_.get(i));
        fileName = expNameCut1_.get(i);
      }else
        expNames_.put(fileName,fileName);
      expNameToFile_.put(fileName, resultFile);
      expNamesInSequence_.add(fileName);
      parseResultFile(resultFile,fileName, statisticsViewMode, combineOxWithNonOx);
    }
    buildResultHashes();
    for (String className : sortedISNames_.keySet()){
      sortedISNames_.put(className,sortLipidNames(sortedISNames_.get(className)));
      sortedESNames_.put(className,sortLipidNames(sortedESNames_.get(className)));
    }
  }
  
  protected void buildResultHashes(){
    
  }
  
  protected void parseResultFile(File resultFile, String fileName, int statisticsViewMode, boolean combineOxWithNonOx) throws ExcelInputFileException, LipidCombinameEncodingException{
  	if (resultFile.getAbsolutePath().endsWith(".xlsx")){
    	parseResultFileFastExcel(resultFile, fileName, statisticsViewMode, combineOxWithNonOx);
    }
  	else 
  	{
    	new WarningMessage(new JFrame(), "ERROR", "The specified file format is not supported!");
      throw new ExcelInputFileException("The specified file format is not supported!");
    }
  }
  
  protected void parseResultFileFastExcel(File resultFile, String fileName, int statisticsViewMode, boolean combineOxWithNonOx) throws ExcelInputFileException, LipidCombinameEncodingException{
    try (InputStream is = new FileInputStream(resultFile);
      ReadableWorkbook wb = new ReadableWorkbook(is);
        Stream<Sheet> sheets = wb.getSheets();)
  	{
      sheets.filter((s) -> (!s.getName().equals(QuantificationResultExporter.SHEET_CONSTANTS)&&
						                !s.getName().endsWith(QuantificationResultExporter.ADDUCT_OMEGA_SHEET)&&
						                !s.getName().endsWith(QuantificationResultExporter.ADDUCT_OVERVIEW_SHEET)&&
						                !s.getName().endsWith(QuantificationResultExporter.ADDUCT_MSN_SHEET)))
            .forEach((s) -> {parseSheet(s);});      
    } catch (IOException ex){
      ex.printStackTrace();
      new WarningMessage(new JFrame(), "ERROR", ex.getMessage());
      throw new ExcelInputFileException(ex);
    }
  }
  //adapted from the original apache poi version
  protected void parseSheet(org.dhatim.fastexcel.reader.Sheet sheet) {
  	Hashtable<String,String> isNames = new Hashtable<String,String>();
    Hashtable<String,String> esNames = new Hashtable<String,String>();
    Vector<String> sortedIS = new Vector<String>();
    Vector<String> sortedES = new Vector<String>();
    String groupName = sheet.getName();
    if (lipidClassesHash_.containsKey(groupName)){
      isNames = allISNames_.get(groupName);
      esNames = allESNames_.get(groupName);
      sortedIS = sortedISNames_.get(groupName);
      sortedES = sortedESNames_.get(groupName);
    }else{
      lipidClasses_.add(groupName);
      lipidClassesHash_.put(groupName, groupName);
    }
    
    List<Row> rows = null;
    try {
      rows = sheet.read();
    } catch (IOException ex) {}
    Row headerRow = rows.get(QuantificationResultExporter.HEADER_ROW);
    List<String> headerTitles = LDAResultReader.readSheetHeaderTitles(headerRow);
    List<Row> contentRows = rows.subList(QuantificationResultExporter.HEADER_ROW+1, rows.size());
    
    for (Row row : contentRows) {
      List<Cell> cells = row.stream().filter((c) -> !(c==null || c.getType().equals(CellType.ERROR))).collect(Collectors.toList());
      int index;
      String rawValue;
      String name = null;
      int dbs = -1;
      String oxState = "";
      
      for (Cell cell : cells) {
        index  = cell.getColumnIndex();
        rawValue = cell.getRawValue();
        
        if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_NAME)) {
          name = rawValue;
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_DBS)) {
          dbs = (int)Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(LipidomicsConstants.CHAIN_MOD_COLUMN_NAME)) {
          oxState = rawValue;
        }
      }
      
      if (name!=null&&name.length()>0)
      {
        Integer doubleBonds = null;
        if (dbs > -1)
        {
        	doubleBonds = dbs;
        }
        else if (!headerTitles.contains(QuantificationResultExporter.HEADER_MODIFICATION)) 
        {
        	Object[] components = ComparativeNameExtractor.splitOldNameStringToComponents(name);
          name = (String)components[0];
          doubleBonds = (Integer)components[1];
        }
        name = StaticUtils.generateLipidNameString(name,doubleBonds,-1, oxState);
        if (isSelectionPrefix_!=null && name.startsWith(isSelectionPrefix_)){
          if (!isNames.containsKey(name)){
            isNames.put(name, name);
            sortedIS.add(name);
          }
        }
        else if (esSelectionPrefix_!=null && name.startsWith(esSelectionPrefix_)){
          if (!esNames.containsKey(name)){
            esNames.put(name, name);
            sortedES.add(name);
          }                
        }
      }
    }
	  allISNames_.put(groupName, isNames);
	  allESNames_.put(groupName, esNames);
	  sortedISNames_.put(groupName, sortedIS);
	  sortedESNames_.put(groupName, sortedES); 
  }
  
//  private void addItemAtSortedPosition(String name, String previousName, Vector<String> sorted){
//    if (previousName!=null&&previousName.length()>0){
//      int position = this.findPositionOfItem(previousName, sorted);
//      if (position>-1&&(position+1)<sorted.size())
//        sorted.add(position,name);
//      else
//        sorted.add(name);
//    }else
//      sorted.add(name);    
//  }
  
  public static String removeChemicalFormula(String name){
    if (name != null && name.contains("_")){
      return name.substring(0,name.lastIndexOf("_"));
      // this is for backward compatibility
    } else if (name!=null && name.contains(":")){
      char[] nameChars = name.toCharArray();
      int idx = name.lastIndexOf(":")+1;
      String returnName = name.substring(0,idx);
      for (int i=idx;i!=nameChars.length;i++){
        if (Character.isDigit(nameChars[i]))
          returnName += String.valueOf(nameChars[i]);
        else
          break;
      }
      return returnName;
    }else
      return name;
  }
  
  protected int findPositionOfItem(String nameToLook,Vector<String> moleculeNames){
    int position = -1;
    for (int i=0;i!=moleculeNames.size();i++){
      if (nameToLook.equalsIgnoreCase(moleculeNames.get(i))){
        position = i;
        break;
      }
    }
    return position;
  }
  
  //TODO: this is very old code. remove as soon as the projects are merged.
  @Deprecated
  public static Object[] splitOldNameStringToComponents(String fullNameString){
    Object[] components = new Object[3];
    String name = "";
    Integer doubleBonds = -1;
    String formula = "";
    String nameAndDoubleBond = removeChemicalFormula(fullNameString);
    formula = fullNameString.substring(nameAndDoubleBond.length());
    if (formula.startsWith("_")) formula = formula.substring(1);
    if (nameAndDoubleBond.lastIndexOf(":")!=-1){
      name = nameAndDoubleBond.substring(0,nameAndDoubleBond.lastIndexOf(":"));
      doubleBonds = new Integer(nameAndDoubleBond.substring(nameAndDoubleBond.lastIndexOf(":")+1));
    }else name = nameAndDoubleBond;
    components[0] = name;
    components[1] = doubleBonds;
    components[2] = formula;
    return components;
  }
  
  @SuppressWarnings("unchecked")
  protected Vector<String> sortLipidNames(Vector<String> moleculeNames){
    boolean containAllDoubleBonds = true;
    boolean containAllRtValues = true;
    if (moleculeNames.size()==0){
      containAllDoubleBonds = false;
      containAllRtValues = false;
    }  
    for (String molName : moleculeNames){
      if (molName.lastIndexOf(":")==-1){
        containAllDoubleBonds = false;
        break;
      }
    }
    for (String molName : moleculeNames){
      boolean isRt = false;
      if (molName.lastIndexOf("_")!=-1){
        try{
          new Double(molName.substring(molName.lastIndexOf("_")+1));
          isRt = true;
        }catch(NumberFormatException nfx){}
      }else if ((isSelectionPrefix_!=null && molName.startsWith(isSelectionPrefix_))||(esSelectionPrefix_!=null && molName.startsWith(esSelectionPrefix_))){
        isRt = true;
      }
      if (!isRt){
        containAllRtValues = false;
        break;        
      }
    }
    if (containAllDoubleBonds){
      Vector<String> moleculeNamesNewOrder = new Vector<String>();
      Hashtable<String,Set<Integer>> mainGroupMembers = new Hashtable<String,Set<Integer>>();
      Hashtable<String,Hashtable<Integer,Vector<DoubleStringVO>>> rtMembers = new  Hashtable<String,Hashtable<Integer,Vector<DoubleStringVO>>>();
      Vector<String> orderMainGroup = new Vector<String>();
      for (String molName : moleculeNames){
        String mainGroupName = molName.substring(0,molName.lastIndexOf(":"));
        String sideChainString = molName.substring(molName.lastIndexOf(":")+1);
        if (containAllRtValues && molName.indexOf("_")!=-1) sideChainString = sideChainString.substring(0,sideChainString.lastIndexOf("_"));
        Integer sideChain = null;
        try{
          sideChain = new Integer(sideChainString);
        } catch (NumberFormatException nfx){
          return moleculeNames;
        }
        DoubleStringVO rt = null;
        if (containAllRtValues && molName.indexOf("_")!=-1){
          String rtAsString = molName.substring(molName.lastIndexOf("_")+1);
          rt = new DoubleStringVO(rtAsString,new Double(rtAsString));
        }
        Set<Integer> groupMembers = new HashSet<Integer>();
        Hashtable<Integer,Vector<DoubleStringVO>> rtHash = new Hashtable<Integer,Vector<DoubleStringVO>>();
        if (!mainGroupMembers.containsKey(mainGroupName)){
          orderMainGroup.add(mainGroupName);          
        }else{
          groupMembers = mainGroupMembers.get(mainGroupName);
          rtHash = rtMembers.get(mainGroupName);
        }
        groupMembers.add(sideChain);
        mainGroupMembers.put(mainGroupName, groupMembers);
        if (containAllRtValues && molName.indexOf("_")!=-1){
          Vector<DoubleStringVO> rts = new Vector<DoubleStringVO>();
          if (rtHash.containsKey(sideChain)) rts = rtHash.get(sideChain);
          rts.add(rt);
          rtHash.put(sideChain, rts);
        }
        rtMembers.put(mainGroupName, rtHash);
      }
      orderMainGroup = sortMainGroup(orderMainGroup);
      for (String mainGroupName: orderMainGroup){
        Vector<Integer> members = new Vector<Integer>(mainGroupMembers.get(mainGroupName));
        Collections.sort(members);
        for (Integer member : members){
          if (containAllRtValues && !(isSelectionPrefix_!=null && mainGroupName.startsWith(isSelectionPrefix_)) && !(esSelectionPrefix_!=null && mainGroupName.startsWith(esSelectionPrefix_))){
            Vector<DoubleStringVO> rts = new Vector<DoubleStringVO>(rtMembers.get(mainGroupName).get(member));
            Collections.sort(rts,new GeneralComparator("at.tugraz.genome.lda.vos.DoubleStringVO", "getValue", "java.lang.Double"));
            for (DoubleStringVO rt : rts) moleculeNamesNewOrder.add(mainGroupName+":"+String.valueOf(member)+"_"+rt.getKey());
          } else moleculeNamesNewOrder.add(mainGroupName+":"+String.valueOf(member));
        }
      }
      moleculeNames = new Vector<String>(moleculeNamesNewOrder);
    }else if (containAllRtValues){
      Vector<String> moleculeNamesNewOrder = new Vector<String>();
      Hashtable<String,Vector<DoubleStringVO>> rtMembers = new  Hashtable<String,Vector<DoubleStringVO>>();
      LinkedHashMap<String,String> mainGroupMembers = new LinkedHashMap<String,String>();
      for (String molName : moleculeNames){
        if (!molName.startsWith(isSelectionPrefix_) && !molName.startsWith(esSelectionPrefix_)){
          String mainGroupName = molName.substring(0,molName.indexOf("_"));
          String rtAsString = molName.substring(molName.lastIndexOf("_")+1);
          DoubleStringVO  rt = new DoubleStringVO(rtAsString,new Double(rtAsString));
          Vector<DoubleStringVO> rts = new Vector<DoubleStringVO>();
          if (!mainGroupMembers.containsKey(mainGroupName)){
            mainGroupMembers.put(mainGroupName, mainGroupName);       
          }else{
            rts = rtMembers.get(mainGroupName);
          }
          rts.add(rt);
          rtMembers.put(mainGroupName, rts);
        } else {
          rtMembers.put(molName, new Vector<DoubleStringVO>());
        }
      }
      for (String mainGroupName: rtMembers.keySet()){
        Vector<DoubleStringVO> rts = new Vector<DoubleStringVO>(rtMembers.get(mainGroupName));
        if (rts.size()>0){
          Collections.sort(rts,new GeneralComparator("at.tugraz.genome.lda.vos.DoubleStringVO", "getValue", "java.lang.Double"));
          for (DoubleStringVO rt : rts) moleculeNamesNewOrder.add(mainGroupName+"_"+rt.getKey());
        } else {
          moleculeNamesNewOrder.add(mainGroupName);
        }
      }
      moleculeNames = new Vector<String>(moleculeNamesNewOrder);
      
    }
    return moleculeNames;
  }
  
  protected Vector<String> sortMainGroup(Vector<String> mainGroup){
    Collections.sort(mainGroup);
    return mainGroup;
  }
  
  public Hashtable<String,String> getExpNames()
  {
    return expNames_;
  }
  
  public Vector<String> getExpNamesInSequence()
  {
    return expNamesInSequence_;
  }
  
  public Vector<String> getGroupNames(){
    return this.groups_;
  }
    
  public Vector<String> getISNames(String className){
    return this.sortedISNames_.get(className);
  }
  
  public Vector<String> getESNames(String className){
    return this.sortedESNames_.get(className);
  }
  
  public File getFullFilePath(String expName){
    return this.expNameToFile_.get(expName);
  }
  
  public Vector<ResultFileVO> getResultFileVOs() {
    return this.resultFileVO_;
  }
  
  /**
   * disable the RT-grouping mechanism (for shotgun data)
   */
  protected void disableRtGrouping(){}
}
