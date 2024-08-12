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

package at.tugraz.genome.lda.msn.parser;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Hashtable;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.FattyAcidsContainer;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * Parser for the Long Chain Base Libraries stored in Excel format
 * @author Juergen Hartler
 *
 */
public class LCBLibParser extends FALibParser
{
  
  /** the Excel file must contain a tab ending with "LCB"*/
  private final static String LCB_SHEET_SUFFIX = "LCB";
  /** the highest hydroxylation number available*/
  private short highestHydroxyNumber_ = 0;

  
  /** hash containing the result of the parsing
   * the first String is the encoding of the number of hydroxylation sites
   * the second integer is the amount of carbon atoms
   * the third integer is the amount of double bonds
   */
  private Hashtable<String,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> result_;
  
  
  /**
   * Constructor specifying the Excel file to parse
   * @param file path to the Excel file
   * @throws IOException
   */
  public LCBLibParser (File file) throws IOException{
    super(file);
  }

  /**
   * Command that starts the parsing of the Excel file specified in the constructor
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException general exception if a file is not there
   */
  public void parseFile() throws RulesException, IOException{
    result_ = new Hashtable<String,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>();
    availablePrefixes_ = new HashSet<String>();
    //first, create a hashtable of all available hydroxylation sites
    HydroxyEncoding hydroxies = Settings.getLcbHydroxyEncoding();
    short hydroxyKey;
    Hashtable<String,Short> encodedHydroxyStrings = new Hashtable<String,Short>();
    for (Object hydroxyObject : hydroxies.keySet()) {
      hydroxyKey = Short.parseShort((String)hydroxyObject);
      encodedHydroxyStrings.put((String)hydroxies.get(hydroxyObject),hydroxyKey);
      if (hydroxyKey>highestHydroxyNumber_)
        highestHydroxyNumber_ = hydroxyKey;
    }
    //second, read all sheets whose tabs are named such as $HYDROXY_ENCODING$ followed by LCB_SHEET_SUFFIX, e.g. dLCB
    InputStream myxls = null;
    try{
      myxls = new FileInputStream(inputFile_);
      Workbook workbook = null;
      if (inputFile_.getAbsolutePath().endsWith(FattyAcidsContainer.FA_FILE_SUFFIX_NEW))
        workbook = new XSSFWorkbook(myxls);
      else if (inputFile_.getAbsolutePath().endsWith(FattyAcidsContainer.FA_FILE_SUFFIX_OLD))
        workbook = new HSSFWorkbook(myxls);
      Sheet sheet;
      String sheetName;
      for (int sheetNumber = 0; sheetNumber!=workbook.getNumberOfSheets(); sheetNumber++){
        sheetName = workbook.getSheetAt(sheetNumber).getSheetName();
        sheet = null;
        String hydrKey = null;
        for (String key : encodedHydroxyStrings.keySet()) {
          if (sheetName.equalsIgnoreCase(key+LCB_SHEET_SUFFIX)) {
            sheet = workbook.getSheetAt(sheetNumber);
            hydrKey = key;
            break;
          }
        }
        if (sheet==null) continue;
        result_.put(hydrKey, parseSheet(sheet,encodedHydroxyStrings.get(hydrKey)));
        //correct result for correct chain type
        for (Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>> sameC : result_.get(hydrKey).values()){
          for (Hashtable<String, Hashtable<String, FattyAcidVO>> sameDb : sameC.values()) {
        	for (Hashtable<String, FattyAcidVO> sameOx : sameDb.values()) {  
	            for (FattyAcidVO fa : sameOx.values()) {
	              fa.correctChainType(LipidomicsConstants.CHAIN_TYPE_LCB);
	            }
        	}
          }
        }
//        if (workbook.getSheetAt(sheetNumber).getSheetName().equalsIgnoreCase(FAS_SHEET_NAME)) sheet = workbook.getSheetAt(sheetNumber);
      }
      if (result_.size()==0) throw new RulesException("A long chain base library must contain at least one sheet starting with the hydroxylation code, followed by \""+LCB_SHEET_SUFFIX+"\", e.g. d"+LCB_SHEET_SUFFIX+"! The lib "+inputFile_.getName()+" has not!");
     
      int previousResultNumber = 0;
      Hashtable<Integer,Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> otherHydroxyResults = null;
      // when there are chains entered, fill up the chain values for all possible hydroxylation sites
      double oxygenMass = Settings.getElementParser().getElementDetails("O").getMonoMass();
      String encoded;
      for (int i=0; i<=highestHydroxyNumber_; i++) {
        encoded = null;
        try {encoded=hydroxies.getEncodedPrefix((short)i);}catch(HydroxylationEncodingException hdx) {continue;}
        // when the results already contain this hydroxy number, no further steps are necessary, but we keep track of it as the last found hydroxylation site 
        if (result_.containsKey(encoded)) {
          previousResultNumber = i;
          otherHydroxyResults = result_.get(encoded);
          continue;
        }
        //if there have not been found entries with less hydroxylation sites, look for the next one with more
        if (otherHydroxyResults==null) {
          String otherEncoded;
          for (int j=i+1; j<=highestHydroxyNumber_; j++) {
            try {otherEncoded=hydroxies.getEncodedPrefix((short)j);}catch(HydroxylationEncodingException hdx) {continue;}
            if (result_.containsKey(otherEncoded)) {
              previousResultNumber = j;
              otherHydroxyResults = result_.get(otherEncoded);
              continue;
            }
          }
        }
        // now, create a the new hash tables and calculate the mass differences according to the hydroxylation numbers
        Hashtable<Integer,Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> currentHydroxy = new Hashtable<Integer,Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>>();
        int hydroxyDifference = i-previousResultNumber;
        for (Integer cAtoms : otherHydroxyResults.keySet()) {
          Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>> otherCAtoms = otherHydroxyResults.get(cAtoms);
          Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>> cAtomsRes = new Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>();
          for (Integer dbs : otherCAtoms.keySet()) {
        	Hashtable<String, Hashtable<String, FattyAcidVO>> otherDbs = otherCAtoms.get(dbs);
        	Hashtable<String, Hashtable<String, FattyAcidVO>> dbsRes = new Hashtable<String, Hashtable<String, FattyAcidVO>>();
            for (String prefix : otherDbs.keySet()) {
            	Hashtable<String, FattyAcidVO> otherOxs = otherDbs.get(prefix);
                Hashtable<String,FattyAcidVO> OxsRes = new Hashtable<String,FattyAcidVO>();
                for(String oxState: otherOxs.keySet())
                {
              	    FattyAcidVO otherFA = otherOxs.get(oxState);
                    Hashtable<String,Integer> formulaCat = StaticUtils.categorizeFormula(otherFA.getFormula());
                    int oxygens = 0;
                    if (formulaCat.containsKey("O")) oxygens = formulaCat.get("O");
 // Krettler                   oxygens += i;
                    oxygens += hydroxyDifference;
                    formulaCat.put("O", oxygens);
// Krettler                    FattyAcidVO chainVO = new FattyAcidVO(otherFA.getChainType(), otherFA.getPrefix(), otherFA.getcAtoms(), otherFA.getDoubleBonds(),i, otherFA.getMass()+i*oxygenMass, StaticUtils.getFormulaInHillNotation(formulaCat, true),otherFA.getOxState());
                    FattyAcidVO chainVO = new FattyAcidVO(LipidomicsConstants.CHAIN_TYPE_LCB, otherFA.getPrefix(), otherFA.getcAtoms(), otherFA.getDoubleBonds(),i, otherFA.getMass()+hydroxyDifference*oxygenMass, StaticUtils.getFormulaInHillNotation(formulaCat, true),otherFA.getOxState());
                    OxsRes.put(oxState, chainVO);
                }	
              dbsRes.put(prefix, OxsRes);
            }
            cAtomsRes.put(dbs, dbsRes);
          }
          currentHydroxy.put(cAtoms, cAtomsRes);
        }
        result_.put(encoded, currentHydroxy);
      }
    //the ChemicalFormulaException can never be thrown, since it is generated from an already checked FattyAcidVO
    }catch (ChemicalFormulaException e) { 
      e.printStackTrace();
    }finally{
      try{myxls.close();}catch(Exception ex){};
    }
  }
  
  
  /**
   * @return hash containing the result of the parsing; first key: encoding of the number of hydroxylation sites; second key: number of carbon atoms; third key: number of double bonds
   */
  public Hashtable<String,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>> getResult(){
    return result_;
  }
  
}
