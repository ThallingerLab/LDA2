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

package at.tugraz.genome.lda.msn.parser;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.exception.SheetNotPresentException;
import at.tugraz.genome.lda.msn.FattyAcidsContainer;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.parser.ModificationParser;
import at.tugraz.genome.lda.utils.StaticUtils;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

/**
 * Parser for the Fatty Acid Chain Libraries stored in Excel format
 * @author Juergen Hartler
 *
 */
public class FALibParser
{
  // the Excel file must contain a tab called "FAS"
  private final static String FAS_SHEET_NAME = "FAS";
    
  protected File inputFile_;
  
  /** hash containing the result of the parsing
   * the first String is the encoding of the number of hydroxylation sites
   * the second integer is the amount of carbon atoms
   * the third integer is the amount of double bonds
   * the fourth key is the prefix
   */
  private Hashtable<String,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> result_;
  
  /** the prefixes for a chain - should be the isotopic labels*/
  protected Set<String> availablePrefixes_;
  
  /**
   * Constructor specifying the Excel file to parse
   * @param filePath path to the Excel file
   * @throws IOException
   */
  public FALibParser (String filePath) throws IOException{
	this (new File(filePath));
  }
  
  /**
   * Constructor specifying the Excel file to parse
   * @param file path to the Excel file
   * @throws IOException
   */
  public FALibParser (File file) throws IOException{
    if (!file.exists()) throw new IOException("The file "+file.getAbsolutePath()+" does not exist!");
    inputFile_ = file;
  }
  
  /**
   * Command that starts the parsing of the Excel file specified in the constructor
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException general exception if a file is not there
   * @throws SheetNotPresentException thrown if the tab in the Excel file is not present
   */
  public void parseFile() throws RulesException, IOException, SheetNotPresentException{
    result_ = new Hashtable<String,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>();
    availablePrefixes_ = new HashSet<String>();
    InputStream myxls = null;
    try{
      myxls = new FileInputStream(inputFile_);
      Workbook workbook = null;
      if (inputFile_.getAbsolutePath().endsWith(FattyAcidsContainer.FA_FILE_SUFFIX_NEW))
        workbook = new XSSFWorkbook(myxls);
      else if (inputFile_.getAbsolutePath().endsWith(FattyAcidsContainer.FA_FILE_SUFFIX_OLD))
        workbook = new HSSFWorkbook(myxls);
      Sheet sheet = null;
      for (int sheetNumber = 0; sheetNumber!=workbook.getNumberOfSheets(); sheetNumber++){
        if (workbook.getSheetAt(sheetNumber).getSheetName().equalsIgnoreCase(FAS_SHEET_NAME)) sheet = workbook.getSheetAt(sheetNumber);
      }
      if (sheet==null) throw new SheetNotPresentException("A fatty acid chain library must have a sheet called \""+FAS_SHEET_NAME+"\"! The lib "+inputFile_.getName()+" has not!");
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>> readResult = parseSheet(sheet);
      result_.put(HydroxyEncoding.HYDROXYLATION_ZERO, readResult);
      
      if (Settings.getFaHydroxyEncoding()!=null) {
        //first: find the highest hydroxyNumber
        int highestHydroxyNumber = 0;
        HydroxyEncoding hydroxies = Settings.getFaHydroxyEncoding();
        short hydroxyKey;
        HashSet<String> encodedHydroxyStrings = new HashSet<String>();
        for (Object hydroxyObject : hydroxies.keySet()) {
          hydroxyKey = Short.parseShort((String)hydroxyObject);
          encodedHydroxyStrings.add((String)hydroxies.get(hydroxyObject));
          if (hydroxyKey>highestHydroxyNumber)
            highestHydroxyNumber = hydroxyKey;
        }
        double oxygenMass = Settings.getElementParser().getElementDetails("O").getMonoMass();
        for (int i=1; i<=highestHydroxyNumber; i++) {
          String encoded = null;
          try {
            encoded = hydroxies.getEncodedPrefix((short)i);
          }
          catch (HydroxylationEncodingException e) {
            continue;
          }
          try {
            // now, create a the new hash tables and calculate the mass differences according to the hydroxylation numbers
            Hashtable<Integer,Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> currentHydroxy = new Hashtable<Integer,Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>>();
            for (Integer cAtoms :  result_.get(HydroxyEncoding.HYDROXYLATION_ZERO).keySet()) {
            	Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>> otherCAtoms = result_.get(HydroxyEncoding.HYDROXYLATION_ZERO).get(cAtoms);
            	Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>> cAtomsRes = new Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>();
              for (Integer dbs : otherCAtoms.keySet()) {
                Hashtable<String, Hashtable<String, FattyAcidVO>> otherDbs = otherCAtoms.get(dbs);
                Hashtable<String, Hashtable<String, FattyAcidVO>> dbsRes = new  Hashtable<String, Hashtable<String, FattyAcidVO>>();
                for (String prefix : otherDbs.keySet()) {
                  Hashtable<String, FattyAcidVO> otherOxs = otherDbs.get(prefix);
                  Hashtable<String,FattyAcidVO> OxsRes = new Hashtable<String,FattyAcidVO>();
                  for(String oxState: otherOxs.keySet())
                  {
                	  FattyAcidVO otherFA = otherOxs.get(oxState);
                      Hashtable<String,Integer> formulaCat = StaticUtils.categorizeFormula(otherFA.getFormula());
                      int oxygens = 0;
                      if (formulaCat.containsKey("O")) oxygens = formulaCat.get("O");
                      oxygens += i;
                      formulaCat.put("O", oxygens);
                      FattyAcidVO chainVO = new FattyAcidVO(otherFA.getChainType(), otherFA.getPrefix(), otherFA.getcAtoms(), otherFA.getDoubleBonds(),i, otherFA.getMass()+i*oxygenMass, StaticUtils.getFormulaInHillNotation(formulaCat, true),otherFA.getOxState());
                      OxsRes.put(oxState, chainVO);
                  }
                  dbsRes.put(prefix, OxsRes);	
                }
              cAtomsRes.put(dbs, dbsRes);
              }
              currentHydroxy.put(cAtoms, cAtomsRes);
            }
            result_.put(encoded, currentHydroxy);
          }catch(ChemicalFormulaException cfx) {
            cfx.printStackTrace();
          }
        }
      }
    }finally{
      try{myxls.close();}catch(Exception ex){};
    }
  }


  /**
   * parses the sheet of an Excel file for chain mass lists
   * @param sheet the excel list to parse
   * @return he result hash - first integer is the amount of carbon atoms; second integer is the amount of double bonds; third key is the prefix
   * @throws RulesException specifies in detail which rule has been infringed
   */
  private Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>> parseSheet(Sheet sheet) throws RulesException{
    return this.parseSheet(sheet, 0);
  }

  
  /**
   * parses the sheet of an Excel file for chain mass lists
   * @param sheet the excel list to parse
   * @param ohNumber the number of OH sites
   * @return the result hash - first integer is the amount of carbon atoms; second integer is the amount of double bonds; third key is the prefix
   * @throws RulesException specifies in detail which rule has been infringed
   */
  protected Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>> parseSheet(Sheet sheet, int ohNumber) throws RulesException{
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>> result = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>();
    int cAtomsColumn = -1;
    int dbsColumn = -1;
    int massColumn = -1;
    int oxStateColumn = -1;
    
    boolean foundColumns = false;
    Hashtable<Integer,String> elementColumns = new  Hashtable<Integer,String>();
    Pattern prefixPattern =  Pattern.compile("\\D*(\\d+)");
    for (int rowCount=0;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
      Row row = sheet.getRow(rowCount);
      String prefix = "";
      int cAtoms = -1;
      int dbs = -1;
      double mass = -1d;
      String oxState = "";
    
      
      Hashtable<Integer,String> possibleElementColumns = new  Hashtable<Integer,String>();
      Hashtable<String,Integer> elementalComposition = new Hashtable<String,Integer>();
      for (int i=0; row!=null && i!=(row.getLastCellNum()+1);i++){
        Cell cell = row.getCell(i);
        String contents = "";
        Double numeric = null;
        int cellType = -1;
        if (cell!=null) cellType = cell.getCellType();
        if (cellType==Cell.CELL_TYPE_STRING){
          contents = cell.getStringCellValue();
          if (contents!=null) contents = contents.trim();
          try{ if (contents!=null)numeric = new Double(contents.replaceAll(",", "."));}catch(NumberFormatException nfx){};
        }else if (cellType==Cell.CELL_TYPE_NUMERIC || cellType==Cell.CELL_TYPE_FORMULA){
          numeric = cell.getNumericCellValue();
          contents = String.valueOf(numeric);
        }
        if (contents!=null) contents = contents.trim();
        if (!foundColumns){
          if (contents.equalsIgnoreCase("Name")){
            cAtomsColumn = i;
          } else if (contents.equalsIgnoreCase("dbs")){
            dbsColumn = i;
          } else if (contents.equalsIgnoreCase("mass")){
            massColumn = i;
          } else if (contents.equalsIgnoreCase(LipidomicsConstants.CHAIN_MOD_COLUMN_NAME)){
        	  oxStateColumn = i;
          } else if (contents.trim().length()==1||contents.trim().length()==2){
            boolean ok = false;
            if (Character.isUpperCase(contents.trim().toCharArray()[0])){
              if (contents.trim().length()==2){
                if (Character.isLowerCase(contents.trim().toCharArray()[1]))
                  ok = true;
              }else{
                ok = true;
              }
              if (ok)possibleElementColumns.put(i, contents.trim());
            }
          }

        }else{
          if (i==cAtomsColumn&&contents!=null&&contents.length()>0){
            if (numeric==null){
              Matcher prefixMatcher = prefixPattern.matcher(contents);
              if (prefixMatcher.matches()){
                prefix = prefixMatcher.group(0);
                String cAtomsString = prefixMatcher.group(1);
                prefix = contents.substring(0,contents.indexOf(cAtomsString));
                cAtoms = Integer.parseInt(cAtomsString);
              }
            }else{
              cAtoms = numeric.intValue();
            }
            if (cAtoms<1) throw new RulesException("A fatty acid chain library must have integer values greater or equal 1 in the \"Name\" column! There is the value "+cAtoms+" at line number "+(rowCount+1)+"!");
          }
          if (i==dbsColumn&&contents!=null&&contents.length()>0){
            if (numeric==null) throw new RulesException("A fatty acid chain library must have integer values in the \"dbs\" column! There is the entry \""+contents+"\" at line number "+(rowCount+1)+"!");
            dbs = numeric.intValue();
            if (dbs<0) throw new RulesException("A fatty acid chain library must have integer values greater or equal 0 in the \"dbs\" column! There is the value "+dbs+" at line number "+(rowCount+1)+"!");
          }
          if (i==massColumn&&contents!=null&&contents.length()>0){
            if (numeric==null) throw new RulesException("A fatty acid chain library must have double values in the \"mass\" column! There is the entry \""+contents+"\" at line number "+(rowCount+1)+"!");
            mass = numeric.doubleValue();
            if (mass<=0) throw new RulesException("A fatty acid chain library must have double values greater 0 in the \"mass\" column! There is the value "+mass+" at line number "+(rowCount+1)+"!");
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

        }
      }
      if (!foundColumns && cAtomsColumn != -1 && massColumn != -1){
        foundColumns = true;
        elementColumns = new Hashtable<Integer,String>(possibleElementColumns);
      }
      if (foundColumns && cAtoms>0 && mass>0){
        if (dbsColumn>-1 && dbs<0) continue;
        
        //create new FA objects for each modification
        String[] oxStates = oxState.split(";");
        double newMass = mass;
        for (String oxMod : oxStates)
        {	
        	oxMod = oxMod.replaceAll("\\s", "");
        	//Get new mass and chemical formula for FAS with modifications
        	ModificationParser mp = new ModificationParser(oxMod);
        	mp.parse();
        	newMass = mp.getNewMass(mass);
        	String chemicalFormula = mp.getNewChemicalFormula(elementalComposition);
        	//no negative elements are allowed after an applied modification
            if (chemicalFormula.indexOf("-")!=-1)
              continue;
        	
	       	FattyAcidVO faVO = new FattyAcidVO(LipidomicsConstants.CHAIN_TYPE_FA_ACYL, prefix,cAtoms,dbs,ohNumber,newMass,chemicalFormula,oxMod);
    	    if (prefix!=null && prefix.length()>0)
        	  this.availablePrefixes_.add(prefix);
            Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>> fasWithSameC = new Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>();
            if (result.containsKey(cAtoms)) fasWithSameC = result.get(cAtoms);
            Hashtable<String, Hashtable<String, FattyAcidVO>> fasWithSameDbs = new Hashtable<String,Hashtable<String,FattyAcidVO>>();
            Hashtable<String,FattyAcidVO> fasWithSameOx = new Hashtable<String,FattyAcidVO>();
            
            
            if (fasWithSameC.containsKey(dbs)) fasWithSameDbs = fasWithSameC.get(dbs);
            if (fasWithSameDbs.containsKey(prefix)) fasWithSameOx = fasWithSameDbs.get(prefix);
            if (fasWithSameOx.containsKey(oxMod))
            {
            	throw new RulesException("A fatty acid chain library must not contain the same elements twice! The entry "+faVO.getCarbonDbsId()+" at line number "+(rowCount+1)+" is there for the second time!");
            }
            fasWithSameOx.put(oxMod, faVO);
            fasWithSameDbs.put(prefix, fasWithSameOx);
            fasWithSameC.put(dbs, fasWithSameDbs);
            
            result.put(cAtoms, fasWithSameC);
        }
      }
    }
    if (!foundColumns) throw new RulesException("A fatty acid chain library must have the column headers \"Name\" and \"mass\"! The lib "+inputFile_.getName()+" has not!"); 
    return result;
  }
  
  /**
   * Fetches the results of the parsing - parseFile() has to be called before
   * @return the result hash - first key: encoded hydroxy; second key is the amount of carbon atoms; third key is the amount of double bonds; fourth key is the prefix, fifth key is oxidation state
   */
  public Hashtable<String,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> getFattyAcids(){
    return result_;
  }
  
  /**
   * 
   * @return the isotopic labels
   */
  public Set<String> getAvailableLabels(){
    return this.availablePrefixes_;
  }
}
