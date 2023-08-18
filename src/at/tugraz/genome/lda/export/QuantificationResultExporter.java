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


package at.tugraz.genome.lda.export;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import org.dhatim.fastexcel.VisibilityState;
import org.dhatim.fastexcel.Workbook;
import org.dhatim.fastexcel.Worksheet;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
//import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.QuantificationThread;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.IntensityChainVO;
import at.tugraz.genome.lda.msn.vos.IntensityPositionVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.StaticUtils;
//import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.Probe3D;
import at.tugraz.genome.lda.utils.Pair;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
//TODO: features for omega assingment have been commented out / not added
public class QuantificationResultExporter
{
  public final static String SHEET_CONSTANTS = "About";
  
  public final static String ADDUCT_MSN_SHEET = " - MSn";
  public final static String ADDUCT_OMEGA_SHEET = " - Omega";
  public final static String ADDUCT_OVERVIEW_SHEET = " - Overview";
  
  //Header constants used in multiple sheets
  public final static int HEADER_ROW = 0;
  public final static String HEADER_SPECIES = "Species";
  public final static String HEADER_AREA = "Area";
  public final static String HEADER_ISOTOPE = "Isotope";
  
  //MS1 Sheet Headers
  public final static String HEADER_INDEX = "Index";
  public final static String HEADER_NAME = "Name";
  public final static String HEADER_DBS = "Dbs";
  public final static String HEADER_OMEGA_POSITION = "MS1 assigned omega-DB position";
  public final static String HEADER_MODIFICATION = "Modification";
  public final static String HEADER_FORMULA = "Formula";
  public final static String HEADER_MOD_FORMULA = "Mod-Formula";
  public final static String HEADER_RT = "RT";
  public final static String HEADER_AREA_ERROR = "AreaError";
  public final static String HEADER_BACKGROUND = "Background";
  public final static String HEADER_CHARGE = "Charge";
  /** different from HEADER_MZ to keep compatibility for older result files */
  public final static String HEADER_MZ_MS1 = "Mz";
  public final static String HEADER_MZ_TOLERANCE = "MzTolerance";
  public final static String HEADER_PEAK = "Peak";
  public final static String HEADER_LOWER_VALLEY = "LowerValley";
  public final static String HEADER_UPPER_VALLEY = "UpperValley";
  public final static String HEADER_LOWER_MZ = "LowMz";
  public final static String HEADER_UPPER_MZ = "UpMz";
  public final static String HEADER_ELL_CENT_TIME = "EllCentTime";
  public final static String HEADER_ELL_CENT_MZ = "EllCentMz";
  public final static String HEADER_ELL_STRETCH_TIME = "EllStretchTime";
  public final static String HEADER_ELL_STRETCH_MZ = "EllStretchMz";
  public final static String HEADER_LOWER_RT_HARD_LIMIT = "LowerRtHardLimit";
  public final static String HEADER_UPPER_RT_HARD_LIMIT = "UpperRtHardLimit";
  public final static String HEADER_PERCENTAL_SPLIT = "PercentalSplit";
  public final static String HEADER_RAW_APEX = "Raw Apex";
  public final static String HEADER_LOWER_VALLEY10PC = "LValley10%";
  public final static String HEADER_LOWER_VALLEY50PC = "LValley50%";
  public final static String HEADER_UPPER_VALLEY50PC = "UValley50%";
  public final static String HEADER_UPPER_VALLEY10PC = "UValley10%";
  public final static String HEADER_LOWER_MZ10PC = "LMz10%";
  public final static String HEADER_LOWER_MZ50PC = "LMz50%";
  public final static String HEADER_UPPER_MZ50PC = "UMz50%";
  public final static String HEADER_UPPER_MZ10PC = "UMz10%";
  public final static String HEADER_MS_LEVEL = "level=";
  
  //MSn Sheet Headers
  public final static String HEADER_ALEX123_MSN_TARGETS_USED = "AlexMSnTargetsUsed";

  //Omega Sheet Headers
  public final static String HEADER_IDENTIFIER = "Identifier";
  public final static String HEADER_MOLECULAR_SPECIES = "Molecular Species";
  public final static String HEADER_DOUBLE_BOND_POSITION_LEVEL = "Double Bond Position Level";
  public final static String HEADER_ASSIGNED = "Assigned";
  public final static String HEADER_EXPECTED_RT = "Expected RT / min";
  public final static String HEADER_ACCURACY = "Accuracy";
  
  //Overview Sheet Headers
  public final static String HEADER_MZ = "m/z";
  
  public final static int MSN_ROW_FRAGMENT_NAME = 0;
  public final static int MSN_ROW_FRAGMENT_FORMULA = 1;
  public final static int MSN_ROW_FRAGMENT_OH = 1;
  public final static int MSN_ROW_FRAGMENT_CHAIN_TYPE = 2;
  
  public final static int MSN_ROW_FRAGMENT_MSLEVEL = 2;
  public final static int MSN_ROW_FRAGMENT_CHARGE = 3;
  public final static int MSN_ROW_FRAGMENT_MZ = 4;
  public final static int MSN_ROW_FRAGMENT_MZ_TOLERANCE = 5;
  public final static int MSN_ROW_FRAGMENT_AREA = 6;
  public final static int MSN_ROW_FRAGMENT_PEAK = 7;
  public final static int MSN_ROW_FRAGMENT_TIME_LOWER = 8;
  public final static int MSN_ROW_FRAGMENT_TIME_UPPER = 9;
  public final static int MSN_ROW_FRAGMENT_MZ_LOWER = 10;
  public final static int MSN_ROW_FRAGMENT_MZ_UPPER = 11;
  public final static int MSN_ROW_FRAGMENT_ELLIPSE_TIME = 12;
  public final static int MSN_ROW_FRAGMENT_ELLIPSE_MZ = 13;
  public final static int MSN_ROW_FRAGMENT_ELLIPSE_TIME_RANGE = 14;
  public final static int MSN_ROW_FRAGMENT_ELLIPSE_MZ_RANGE = 15;

  public final static int MSN_ROW_INTENSITY_RULE = 0;
  public final static int MSN_ROW_INTENSITY_ORIGINAL = 1;
  public final static int MSN_ROW_INTENSITY_VALUES = 2;
  public final static int MSN_ROW_INTENSITY_MISSED = 3; 
  
  private static int amountOfIsotopes_;
  private static boolean hasMSnInformation_;
  private static boolean hasOmegaInformation_;
  private static boolean hasOhInformation_;
  /** total area of all isotopes of one analyte, calculated in writeEvidenceMS1 */
  private static double totalArea_;
  
  /**
   * Write a QuantificationResult to an Excel file
   * @param filePath the path to write the Excel file to
   * @param quantRes the QuantificationResult to write out
   * @throws ExportException
   */
  public static void writeResultsToExcel(String filePath, QuantificationResult quantRes) throws ExportException {
	  writeResultsToExcel(filePath, quantRes, null);
  }
  
  /**
   * Write a QuantificationResult to an Excel file
   * @param filePath the path to write the Excel file to
   * @param quantRes the QuantificationResult to write out
   * @param bestMatchBySpectrumCoverage lipid classes for which the best matches should be picked by spectrum coverage
   * @throws ExportException
   */
  public static void writeResultsToExcel(String filePath, QuantificationResult quantRes, Hashtable<String,Boolean> bestMatchBySpectrumCoverage) throws ExportException {
	  
		
		//remove oxLipids also identified as conventional lipids
		  ArrayList<LipidParameterSet> removeElements = new ArrayList<LipidParameterSet>();
		  for (String key : quantRes.getIdentifications().keySet()){
		  	if (key.contains("ox")) {
		  		for (LipidParameterSet param1 : quantRes.getIdentifications().get(key)) {
		  			boolean remove = false;
		  			if(param1 instanceof LipidomicsMSnSet) {
		  				LipidomicsMSnSet msn_param1 = (LipidomicsMSnSet) param1;
		  				for (int msnp1 : msn_param1.getMsnRetentionTimes().get(2).keySet()) {
		  					for (String key2 : quantRes.getIdentifications().keySet()){
		  						if (!key2.contains("ox")) {
		  							for (LipidParameterSet param2 : quantRes.getIdentifications().get(key2)) {
		  								if(param2 instanceof LipidomicsMSnSet) {
		  									LipidomicsMSnSet msn_param2 = (LipidomicsMSnSet) param2;
		  									for (int msnp2 : msn_param2.getMsnRetentionTimes().get(2).keySet()) {				
		  										//&& msn_param1.getMsnRetentionTimes().get(2).get(msnp1) == msn_param2.getMsnRetentionTimes().get(2).get(msnp2)
		  										if(msnp1==msnp2 ) {
		  											String h1 = msn_param1.getNamePlusModHumanReadable();
		  											String h2 = msn_param2.getNamePlusModHumanReadable();
		  											String mz1 = Float.toString(msn_param1.Mz[0]);
		  											String mz2 = Float.toString(msn_param2.Mz[0]);
		  											remove = true;
		  											//System.out.println(h1 + " " + mz1 + " :: " + h2 + " " + mz2);
		  										}
		  									}
		  								}
		  							}
		  						}
		  						
		  					}
		  				}
		  				
		  				if (remove) {
		  					removeElements.add(param1);
		  				}
		  				
		  			}
		  		}
		  		for(LipidParameterSet element:removeElements) {
		  			quantRes.getIdentifications().get(key).remove(element);
		  		}
		  	}
		  }
		  

		  //version2
		  Hashtable<Integer, Vector<Pair<String, Integer>>> duplicates = new Hashtable<Integer, Vector<Pair<String, Integer>>>();  
		  for (String lipidClass : quantRes.getIdentifications().keySet()){
			  if (bestMatchBySpectrumCoverage != null) {
		  		if (bestMatchBySpectrumCoverage.containsKey(lipidClass) && bestMatchBySpectrumCoverage.get(lipidClass) == true) {
		  	
		  		for (int index = 0; index < quantRes.getIdentifications().get(lipidClass).size()-1; index++) {
		  		//for (LipidParameterSet lms : quantRes.getIdentifications().get(lipidClass)) {
		  			
		  			LipidParameterSet lms = quantRes.getIdentifications().get(lipidClass).get(index);
		  			if (!(lms instanceof LipidomicsMSnSet)){
		  				continue;
		  			}
		  			
		  			LipidomicsMSnSet lmn = (LipidomicsMSnSet)lms;
		  		
		  			for (int key : lmn.getMsnRetentionTimes().keySet()) {
		  				LinkedHashMap<Integer, Float> scan_rt_pair = lmn.getMsnRetentionTimes().get(key);
		  				for (int key2 : scan_rt_pair.keySet()) {
		  					//System.out.println(key2);
		  					
		  					Vector<Pair<String, Integer>> vlps;
		  					if (duplicates.containsKey(key2)) {
		  						vlps = duplicates.get(key2);
		  						vlps.add(new Pair<>(lipidClass, index));
		  					}
		  					else {
		  						vlps = new Vector<Pair<String, Integer>>();
		  						vlps.add(new Pair<>(lipidClass, index));
		  					}
		  					
		  					duplicates.put(key2, vlps);
		  				}
		  				
		  			}
		  		}
		  		}
		  }
		  	}
		//
		//
		//
		  //pick best match from duplicates drop others
		  for (int key : duplicates.keySet()) {
		  	Vector<Pair<String, Integer>> vlps = duplicates.get(key);
		  	if (vlps.size()>1) {
		  		float minDiff = 100;
		  		float maxCoverage = 0;
		  		float minMods = 0;
		  		String bestOut = "";
		  		String out = "";
		  		LipidomicsMSnSet bestMatch = null;
		  		
		  		//find best match
		  		for (Pair<String, Integer> pair: vlps) {
		  			LipidomicsMSnSet lps = (LipidomicsMSnSet)quantRes.getIdentifications().get(pair.getKey()).get(pair.getValue());
		  			float theoMz = lps.Mz[0];
		  			float evMz = lps.getIsotopicProbes().get(0).get(0).Mz;
		  			float diff = Math.abs(theoMz-evMz);
		  			float coverage = lps.getCoverage();
		  			
		  			//figure out total number of modifications
		  			int mods = 0;
		  			String oxState = lps.getOxState();
		  			
		  			if(!oxState.equals("")) {
		  				Pattern p = Pattern.compile("(\\d+)");
			  	        Matcher m = p.matcher(oxState);
			  	        
			  	        if(m.find()) {
			  	        	mods =+ Integer.parseInt(m.group());
			  	        	while(m.find()) {
				  	            mods =+ Integer.parseInt(m.group());
				  	        }
			  	        }
			  	        else {
			  	        	mods =+1 ;
			  	        }
			  	        
		  			}
		  			
		  	        
		  	        out = pair.getKey() + " "+ lps.getNamePlusModHumanReadable() + " theo: " + theoMz + " ev: " + evMz + " diff: " + diff + " coverage: " + coverage;
		  		    
		  	        if (coverage>maxCoverage) {
		  	        	maxCoverage = coverage;
		  	        	bestMatch = lps;
		  	        	minDiff = diff;
		  	        	minMods = mods;
		  	        	bestOut = out;
		  	        }
		  	        else if(coverage==maxCoverage && mods==minMods) {
		  	        	if(diff<minDiff) {
		  	        		bestMatch = lps;
		  	        		minDiff = diff;
		  	        		bestOut = out;
		  	        	}
		  	        }
		  	        else if(coverage==maxCoverage) {
		  	        	if(mods<minMods) {
		  	        		bestMatch = lps;
		  	        		minMods = mods;
		  	        		minDiff = diff;
		  	        		bestOut = out;
		  	        	}
		  	        }
		  	
		  			
		  			System.out.println(out);
		  		}
		  		System.out.println("Best " + bestOut);
		  		
		  		//remove all other matches
		  		for (Pair<String, Integer> pair: vlps) {
		  			LipidomicsMSnSet lps = (LipidomicsMSnSet)quantRes.getIdentifications().get(pair.getKey()).get(pair.getValue());
		  			if (lps!=bestMatch) {
		  				Hashtable<Integer, LinkedHashMap<Integer, Float>> scans = lps.getMsnRetentionTimes();
		  				scans.get(2).remove(key);
		  				lps.setMsnRetentionTimes(scans);
		  				
		  				if(lps.getMsnRetentionTimes().get(2).isEmpty())
		  				{
		  					lps = null;
		  				}
		  				
		  				//new
		  				//lps = null;
		  				
		  				quantRes.getIdentifications().get(pair.getKey()).set(pair.getValue(),lps);
		  			}
		  			
		  		}
		  		
		  		System.out.println("");
		  	}
		  		 
		  } 
		  
	  
    try (BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(filePath));) {
      String s = Settings.VERSION;
      //the constructor can only take a version number in the format xx.yyyy
      Workbook wb = new Workbook(out, "Lipid Data Analyzer", s.substring(0, s.indexOf(".", s.indexOf(".") + 1)));
      
      boolean hasRtInfo = QuantificationThread.hasRtInfo(quantRes.getIdentifications());
      
      LipidomicsConstants constants = quantRes.getConstants();
      if (constants!=null){
        Worksheet ws = wb.newWorksheet(SHEET_CONSTANTS);
        List<Pair<String,String>> propertyRows = constants.getPropertyRowList(quantRes.getFaHydroxyEncoding(),quantRes.getLcbHydroxyEncoding());
        List<String> headerTitles = new ArrayList<String>();
        headerTitles.add(propertyRows.get(0).getKey());
        headerTitles.add(propertyRows.get(0).getValue());
        createHeader(ws,headerTitles);
        writeConstants(ws,propertyRows);
      }
      
      //TODO: simple multithreading by assigning a separate thread to each sheet; only if further speedup is necessary
      for (String sheetName : quantRes.getIdentifications().keySet()) { 
        Vector<LipidParameterSet> params = quantRes.getIdentifications().get(sheetName);
        iterateParamsForClassVariables(params);
        boolean alex123TargetsUsed = constants.getAlexTargetlistUsed().containsKey(sheetName) && constants.getAlexTargetlistUsed().get(sheetName);
        int msLevel = 1;
        //TODO: find out if msLevel always remains 1?
        if (quantRes.getMsLevels()!=null&&quantRes.getMsLevels().containsKey(sheetName)) {
          msLevel = quantRes.getMsLevels().get(sheetName);
        }
        
        Worksheet[] sheetsForClass = createRequiredSheetsForClass(wb, sheetName);
        Worksheet resultSheet = sheetsForClass[0];
        Worksheet resultSheetMSn = sheetsForClass[1];
        Worksheet resultSheetOmega = sheetsForClass[2];
        Worksheet resultSheetOverview = sheetsForClass[3];
        
        List<String> mS1HeaderTitles = createMS1HeaderTitles(hasOhInformation_, hasRtInfo, msLevel);
        List<String> omegaHeaderTitles = null;
        List<String> overviewHeaderTitles = null;
        createHeader(resultSheet, mS1HeaderTitles);
        
        if (resultSheetMSn != null && alex123TargetsUsed) {
          createHeaderCell(resultSheetMSn, HEADER_ROW, 0, String.format("%s=true", HEADER_ALEX123_MSN_TARGETS_USED));
        }
        if (resultSheetOmega != null) {
          omegaHeaderTitles = createOmegaHeaderTitles();
          createHeader(resultSheetOmega, omegaHeaderTitles);
        }
        if (resultSheetOverview != null) {
          overviewHeaderTitles = createOverviewHeaderTitles();
          createHeader(resultSheetOverview, overviewHeaderTitles);
        }
        
        int resultCount = 1; 
        int ms1RowCount = HEADER_ROW;
        int msnRowCount = HEADER_ROW+1; /** leave first row empty for additions like "AlexMSnTargetsUsed" */
//        int omegaRowCount = HEADER_ROW+1; 
        
        for (LipidParameterSet param : params){
        	if (param==null)
        		continue;
       	
        	ms1RowCount = writeEvidenceMS1(ms1RowCount, resultCount, resultSheet, param, mS1HeaderTitles);

        	if (param instanceof LipidomicsMSnSet){
            try {
              msnRowCount = writeEvidenceMSn(msnRowCount,resultSheetMSn,(LipidomicsMSnSet)param,
                  quantRes.getFaHydroxyEncoding(),quantRes.getLcbHydroxyEncoding());
            }catch(Exception ex){
              throw ex;
            }
          }
//          if (param.hasOmegaInformation()) {
//            omegaRowCount = writeEvidenceOmega(omegaRowCount, resultSheetOmega, param, omegaHeaderTitles);
//          }
          if (Settings.isOverviewInExcelDesired()) {
            writeEvidenceOverview(resultCount, resultSheetOverview, param, overviewHeaderTitles);
          }
          
          resultCount++;
        }
      }
      wb.finish();
    } catch (IOException ex) {
      throw new ExportException(ex.getMessage());
    } catch (Exception ex) {
      throw new ExportException(ex.getMessage());
    }
    
  }
  
  
  /**
   * Iterates over a Vector of LipidParameterSets to assign values to the class variables
   * @param params Vector of LipidParameterSets
   */
  private static void iterateParamsForClassVariables(Vector<LipidParameterSet> params) {
    amountOfIsotopes_ = 0;
    hasMSnInformation_ = false;
    hasOmegaInformation_ = false;
    hasOhInformation_ = false;
    for (LipidParameterSet param : params){
    	if (param==null){
      		continue;
      	}
      int maxIso = param.getIsotopicProbes().size();
      if (maxIso>amountOfIsotopes_) amountOfIsotopes_ = maxIso;
      if (param instanceof LipidomicsMSnSet) { hasMSnInformation_ = true;}
//      if (param.hasOmegaInformation()) { hasOmegaInformation_ = true;}
      if (param.getOhNumber()>LipidomicsConstants.EXCEL_NO_OH_INFO) { hasOhInformation_ = true;}
    }
  }
  
  
  /**
   * Initiates required sheets for an analyte class depending on class variable values and settings
   * @param wb the Excel workbook
   * @param className the name of the lipid class
   * @return Array of created worksheets
   */
  private static Worksheet[] createRequiredSheetsForClass(Workbook wb, String className) {
    Worksheet resultSheet = wb.newWorksheet(className);
    Worksheet resultSheetMSn  = null;
    Worksheet resultSheetOmega = null;
    Worksheet resultSheetOverview = null;
    
    if (hasMSnInformation_){
      resultSheetMSn = wb.newWorksheet(className+ADDUCT_MSN_SHEET);
      resultSheetMSn.setVisibilityState(VisibilityState.HIDDEN);
    }
    if (hasOmegaInformation_) {
      resultSheetOmega = wb.newWorksheet(className+ADDUCT_OMEGA_SHEET);
      resultSheetOmega.setVisibilityState(VisibilityState.HIDDEN);
    }
    if (Settings.isOverviewInExcelDesired()) {
      resultSheetOverview = wb.newWorksheet(className+ADDUCT_OVERVIEW_SHEET);
    }
    
    return new Worksheet[] {resultSheet, resultSheetMSn, resultSheetOmega, resultSheetOverview};
  }
  
  
  /**
   * Creates a formatted header row with the given header titles in the row number given by HEADER_ROW
   * @param ws Excel worksheet to write the header to
   * @param headerTitles List of header titles
   */
  private static void createHeader(Worksheet ws, List<String> headerTitles) {
    for (int i=0; i<headerTitles.size(); i++) {
      ws.value(HEADER_ROW, i, headerTitles.get(i));
    }
    ws.range(HEADER_ROW, HEADER_ROW, HEADER_ROW, headerTitles.size()).style()
    .bold().horizontalAlignment("center").fontName("Arial").fontSize(12).set();
  }
  
  
  /**
   * Creates a list of header titles to use for an MS1 sheet
   * @param hasOhInfo if true, a header title for an OH column will be included
   * @param hasRtInfo if true, a header title for an RT info column will be included
   * @param msLevel this integer value will be saved in the last field of the header row
   * @return a list of header titles
   */
  private static List<String> createMS1HeaderTitles(boolean hasOhInfo, boolean hasRtInfo, int msLevel) {
    List<String> headerTitles = new ArrayList<String>();
    headerTitles.add(HEADER_INDEX);
    headerTitles.add(HEADER_NAME);
    headerTitles.add(HEADER_DBS);
    if (hasOhInfo) {
      headerTitles.add(LipidomicsConstants.EXCEL_MS_OH);
    }
    if (hasOmegaInformation_) {
      headerTitles.add(HEADER_OMEGA_POSITION);
    }
    headerTitles.add(HEADER_MODIFICATION);
    headerTitles.add(HEADER_FORMULA);
    headerTitles.add(HEADER_MOD_FORMULA);
    if (hasRtInfo){
      headerTitles.add(HEADER_RT);
    }
    headerTitles.add(HEADER_ISOTOPE);
    headerTitles.add(HEADER_AREA);
    headerTitles.add(HEADER_AREA_ERROR);
    headerTitles.add(HEADER_BACKGROUND);
    headerTitles.add(HEADER_CHARGE);
    headerTitles.add(HEADER_MZ_MS1);
    headerTitles.add(HEADER_MZ_TOLERANCE);
    headerTitles.add(HEADER_PEAK);
    headerTitles.add(HEADER_LOWER_VALLEY);
    headerTitles.add(HEADER_UPPER_VALLEY);
    headerTitles.add(HEADER_LOWER_MZ);
    headerTitles.add(HEADER_UPPER_MZ);
    headerTitles.add(HEADER_ELL_CENT_TIME);
    headerTitles.add(HEADER_ELL_CENT_MZ);
    headerTitles.add(HEADER_ELL_STRETCH_TIME);
    headerTitles.add(HEADER_ELL_STRETCH_MZ);
    headerTitles.add(HEADER_LOWER_RT_HARD_LIMIT);
    headerTitles.add(HEADER_UPPER_RT_HARD_LIMIT);
    headerTitles.add(HEADER_PERCENTAL_SPLIT);
    headerTitles.add(HEADER_RAW_APEX);
    headerTitles.add(HEADER_LOWER_VALLEY10PC);
    headerTitles.add(HEADER_LOWER_VALLEY50PC);
    headerTitles.add(HEADER_UPPER_VALLEY50PC);
    headerTitles.add(HEADER_UPPER_VALLEY10PC);
    headerTitles.add(HEADER_LOWER_MZ10PC);
    headerTitles.add(HEADER_LOWER_MZ50PC);
    headerTitles.add(HEADER_UPPER_MZ50PC);
    headerTitles.add(HEADER_UPPER_MZ10PC);
    headerTitles.add(LipidomicsConstants.CHAIN_MOD_COLUMN_NAME);
    /** keep in header row for backward compatibility with older result files */
    headerTitles.add(String.format("%s%s", HEADER_MS_LEVEL, String.valueOf(msLevel)));
    
    return headerTitles;
  }
  
  
  /**
   * Creates a list of header titles to use for an Omega sheet
   * @return a list of header titles
   */
  private static List<String> createOmegaHeaderTitles() {
    List<String> headerTitles = new ArrayList<String>();
    headerTitles.add(HEADER_IDENTIFIER);
    headerTitles.add(HEADER_MOLECULAR_SPECIES);
    headerTitles.add(HEADER_DOUBLE_BOND_POSITION_LEVEL);
    headerTitles.add(HEADER_EXPECTED_RT);
    headerTitles.add(HEADER_ACCURACY);
    headerTitles.add(HEADER_ASSIGNED);
    
    return headerTitles;
  }
  
  
  /**
   * Creates a list of header titles to use for an Overview sheet
   * @return a list of header titles
   */
  private static List<String> createOverviewHeaderTitles() {
    List<String> headerTitles = new ArrayList<String>();
    headerTitles.add(HEADER_SPECIES);
    headerTitles.add(HEADER_AREA);
    if (Settings.isMassInOverviewExcelDesired()){
      headerTitles.add(HEADER_MZ);
    }
    if (Settings.isIsotopeInOverviewExcelDesired()) {
      for (int i=0;i!=amountOfIsotopes_;i++){
        headerTitles.add(String.format("%s %s", HEADER_ISOTOPE, i));
      }
    }
    
    return headerTitles;
  }
  
  
  /**
   * Writes settings into an Excel sheet
   * @param ws the Excel worksheet for the settings
   * @param propertyRows a list of key/value pairs; keys are setting names and values are settings
   */
  private static void writeConstants(Worksheet ws, List<Pair<String,String>> propertyRows) {
    for (int row = 1; row<propertyRows.size(); row++ ) {
      ws.value(row, 0, propertyRows.get(row).getKey());
      ws.value(row, 1, propertyRows.get(row).getValue());
    }
  }
  
  
  /**
   * Writes MS1 evidence of one analyte identification
   * @param row integer value of the current Excel row where subsequent information can be written
   * @param resultCount the result number for this analyte identification
   * @param ws the Excel worksheet for MS1 evidence for this class
   * @param param the LipidParameterSet object holding the information for this analyte identification
   * @param headerTitles list of header titles
   * @return integer value of the next row where subsequent information can be written
   */
  private static int writeEvidenceMS1(int row, int resultCount, Worksheet ws, LipidParameterSet param, List<String> headerTitles) {
    row++;
    int rowFirst = row;
    totalArea_ = 0;
    
    Vector<Vector<CgProbe>> isotopicProbes = param.getIsotopicProbes();
    for (int j=0;j!=isotopicProbes.size();j++){
      int chargeState = j;
      Vector<CgProbe> probes = isotopicProbes.get(j);
      float totalIsoArea = 0f;
      
      for (int k=0; k<probes.size(); k++) {
        row++;
        
        if (param.getMinIsotope()<0) chargeState*=-1;
        CgProbe probe = probes.get(k);
        boolean areaStatusOK = probe!=null&&probe.AreaStatus==CgAreaStatus.OK;
        Probe3D probe3D = null;
        if (probe instanceof Probe3D){
          probe3D = (Probe3D)probe;
        }
        totalIsoArea+=probe.Area;
        
        if (areaStatusOK) {
          
          ws.value(row, headerTitles.indexOf(HEADER_ISOTOPE), chargeState);
          ws.value(row, headerTitles.indexOf(HEADER_AREA), new Double(probe.Area));
          ws.value(row, headerTitles.indexOf(HEADER_AREA_ERROR), new Double(probe.AreaError));
          ws.value(row, headerTitles.indexOf(HEADER_BACKGROUND), new Double(probe.Background));
          ws.value(row, headerTitles.indexOf(HEADER_CHARGE), probe.Charge);
          ws.value(row, headerTitles.indexOf(HEADER_PEAK), new Double(probe.Peak));
          ws.value(row, headerTitles.indexOf(HEADER_LOWER_VALLEY), new Double(probe.LowerValley));
          ws.value(row, headerTitles.indexOf(HEADER_UPPER_VALLEY), new Double(probe.UpperValley));

          if (probe.getLowerValley10() != null){
            ws.value(row, headerTitles.indexOf(HEADER_RAW_APEX), new Double(probe.getApexIntensity()));
            ws.value(row, headerTitles.indexOf(HEADER_LOWER_VALLEY10PC), new Double(probe.getLowerValley10()));
            ws.value(row, headerTitles.indexOf(HEADER_LOWER_VALLEY50PC), new Double(probe.getLowerValley50()));
            ws.value(row, headerTitles.indexOf(HEADER_UPPER_VALLEY50PC), new Double(probe.getUpperValley50()));
            ws.value(row, headerTitles.indexOf(HEADER_UPPER_VALLEY10PC), new Double(probe.getUpperValley10()));
          }
          
          if (probe3D != null) {
            ws.value(row, headerTitles.indexOf(HEADER_LOWER_MZ), new Double(probe3D.LowerMzBand));
            ws.value(row, headerTitles.indexOf(HEADER_UPPER_MZ), new Double(probe3D.UpperMzBand));
            ws.value(row, headerTitles.indexOf(HEADER_ELL_CENT_TIME), new Double(probe3D.getEllipseTimePosition()));
            ws.value(row, headerTitles.indexOf(HEADER_ELL_CENT_MZ), new Double(probe3D.getEllipseMzPosition()));
            ws.value(row, headerTitles.indexOf(HEADER_ELL_STRETCH_TIME), new Double(probe3D.getEllipseTimeStretch()));
            ws.value(row, headerTitles.indexOf(HEADER_ELL_STRETCH_MZ), new Double(probe3D.getEllipseMzStretch()));
            
            if (probe3D.getLowMz10()>-1) {
              ws.value(row, headerTitles.indexOf(HEADER_LOWER_MZ10PC), new Double(probe3D.getLowMz10()));
              ws.value(row, headerTitles.indexOf(HEADER_LOWER_MZ50PC), new Double(probe3D.getLowMz50()));
              ws.value(row, headerTitles.indexOf(HEADER_UPPER_MZ50PC), new Double(probe3D.getUpMz50()));
              ws.value(row, headerTitles.indexOf(HEADER_UPPER_MZ10PC), new Double(probe3D.getUpMz10()));
            } 
            
          }  
          
        } else {
          ws.value(row, headerTitles.indexOf(HEADER_AREA), new Double(0));
        }
        
        ws.value(row, headerTitles.indexOf(HEADER_MZ_MS1), probe.Mz);
      }
      
      totalArea_+=totalIsoArea;
    }
      
    //the first row follows last, as the value for totalArea_ is not known before
      
    ws.value(rowFirst, headerTitles.indexOf(HEADER_INDEX), new Integer(resultCount));
    ws.value(rowFirst, headerTitles.indexOf(HEADER_NAME), param.Peptide);
    ws.value(rowFirst, headerTitles.indexOf(HEADER_DBS), String.valueOf(param.getDoubleBonds()));
      
    if (headerTitles.indexOf(LipidomicsConstants.EXCEL_MS_OH) >= 0) {
      ws.value(rowFirst, headerTitles.indexOf(LipidomicsConstants.EXCEL_MS_OH), String.valueOf(param.getOhNumber()));
    }
      
      //TODO: I still believe we should not write the double bond positions out here
//      if ((headerTitles.indexOf(HEADER_OMEGA_POSITION) >= 0) && !(param instanceof LipidomicsMSnSet) && param.hasOmegaInformation()) {
//        Vector<DoubleBondPositionVO> assignedDoubleBondPositionVO = StaticUtils.getAssignedDoubleBondPositions(param.getOmegaInformation());
//        
//        /** there can only be one assigned species, which can only have one chain */
//        if (!assignedDoubleBondPositionVO.isEmpty()) {
//          int omegaPosition = assignedDoubleBondPositionVO.get(0).getChainCombination().get(0).getOmegaPosition();
//          ws.value(rowFirst, headerTitles.indexOf(HEADER_OMEGA_POSITION), String.valueOf(omegaPosition));
//        }
//        
//      }
      
    ws.value(rowFirst, headerTitles.indexOf(HEADER_MODIFICATION), param.getModificationName());
    ws.value(rowFirst, headerTitles.indexOf(HEADER_FORMULA), param.getAnalyteFormula());
    ws.value(rowFirst, headerTitles.indexOf(HEADER_MOD_FORMULA), param.getModificationFormula());
    
    if (headerTitles.indexOf(HEADER_RT) >= 0) {
      ws.value(rowFirst, headerTitles.indexOf(HEADER_RT), param.getRt());
    }
      
    ws.value(rowFirst, headerTitles.indexOf(HEADER_AREA), totalArea_);
    ws.value(rowFirst, headerTitles.indexOf(HEADER_CHARGE), param.getCharge());
    ws.value(rowFirst, headerTitles.indexOf(HEADER_MZ_MS1), param.Mz[0]);
    ws.value(rowFirst, headerTitles.indexOf(HEADER_MZ_TOLERANCE), param.LowerMzBand);
    
    if (param.getLowerRtHardLimit()>=0){
      ws.value(rowFirst, headerTitles.indexOf(HEADER_LOWER_RT_HARD_LIMIT), param.getLowerRtHardLimit());
    }
      
    if (param.getUpperRtHardLimit()>=0){
      ws.value(rowFirst, headerTitles.indexOf(HEADER_UPPER_RT_HARD_LIMIT), param.getUpperRtHardLimit());
    }
      
    if (param.getPercentalSplit()>=0){
      ws.value(rowFirst, headerTitles.indexOf(HEADER_PERCENTAL_SPLIT), param.getPercentalSplit());
    }
    
    if (param.getOxState() != null && !param.getOxState().isEmpty())
    	ws.value(rowFirst,  headerTitles.indexOf(LipidomicsConstants.CHAIN_MOD_COLUMN_NAME), param.getOxState());
    
    return row++;
  }
  
  
  /**
   * Writes an overview of the evidence for one analyte identification 
   * @param row integer value of the row for this entry
   * @param ws the Excel worksheet for an overview for this class
   * @param param the LipidParameterSet object holding the information for this analyte identification
   * @param headerTitles list of header titles
   */
  private static void writeEvidenceOverview(int row, Worksheet ws, LipidParameterSet param, List<String> headerTitles) {
    Vector<Vector<CgProbe>> isotopicProbes = param.getIsotopicProbes();
    
    createHeaderCell(ws, row, headerTitles.indexOf(HEADER_SPECIES), param.getNamePlusModHumanReadable());
    ws.value(row, headerTitles.indexOf(HEADER_AREA), totalArea_);
    
    if (Settings.isMassInOverviewExcelDesired()){ 
      ws.value(row, headerTitles.indexOf(HEADER_MZ), param.Mz[0]);
    }
    
    if (Settings.isIsotopeInOverviewExcelDesired()){
      for (int i=0;i<amountOfIsotopes_;i++){
        if (param.getIsotopicProbes().size()>i){
          Vector<CgProbe> probes = isotopicProbes.get(i);
          float totalIsoArea = 0f;
          for (int j=0; j!=probes.size(); j++){
            totalIsoArea+=probes.get(j).Area;
          }
          ws.value(row, headerTitles.indexOf(String.format("%s %s", HEADER_ISOTOPE, i)), totalIsoArea);
        }
      }
    }
  }
  
  
  /**
   * Writes double bond position evidence for an analyte identification
   * @param row integer value of the current Excel row where subsequent information can be written
   * @param ws the Excel worksheet for double bond position evidence for this class
   * @param param the LipidParameterSet object holding the information for this analyte identification
   * @param headerTitles list of header titles
   * @return integer value of the next row where subsequent information can be written
   */
//  public static int writeEvidenceOmega(int row, Worksheet ws, LipidParameterSet param, List<String> headerTitles) {
//    ws.value(row, headerTitles.indexOf(HEADER_IDENTIFIER), param.getNamePlusModHumanReadable());
//    
//    Vector<DoubleBondPositionVO> paramOmegaInfo = param.getOmegaInformation();
//    Iterator<DoubleBondPositionVO> it = paramOmegaInfo.iterator();
//    
//    while(it.hasNext()) {
//      
//      DoubleBondPositionVO doubleBondPositionVO = it.next();
//      ws.value(row, headerTitles.indexOf(HEADER_MOLECULAR_SPECIES), doubleBondPositionVO.getMolecularSpecies());
//      ws.value(row, headerTitles.indexOf(HEADER_DOUBLE_BOND_POSITION_LEVEL), doubleBondPositionVO.getDoubleBondPositionsHumanReadable());
//      ws.value(row, headerTitles.indexOf(HEADER_EXPECTED_RT), doubleBondPositionVO.getExpectedRetentionTime());
//      ws.value(row, headerTitles.indexOf(HEADER_ACCURACY), doubleBondPositionVO.getAccuracy());
//      ws.value(row, headerTitles.indexOf(HEADER_ASSIGNED), doubleBondPositionVO.getIsAssigned());
//      row++;
//      
//    }
//    return row;
//  }

  
  /**
   * Writes the MSn evidence for one identified MSn species
   * @param row the current Excel row where subsequent information can be written
   * @param ws the Excel worksheet for MSn evidence for this class
   * @param param the analyte identification containing MSn evidence
   * @param faHydroxyEncoding the character encoding of the number of hydroxylation sites for the FA
   * @param lcbHydroxyEncoding the character encoding of the number of hydroxylation sites for the LCB
   * @return the next row where subsequent information can be written
   * @throws RulesException
   * @throws LipidCombinameEncodingException thrown when a lipid combi ID (containing type and OH number) cannot be decoded 
   * @throws ChemicalFormulaException thrown when something is wrong with the chemical formula
   */
  @SuppressWarnings("unchecked")
  private static int writeEvidenceMSn(int row, Worksheet ws, LipidomicsMSnSet param, HydroxyEncoding faHydroxyEncoding,
      HydroxyEncoding lcbHydroxyEncoding) throws RulesException, LipidCombinameEncodingException, ChemicalFormulaException {
    int count = row;
    int nameRow = count;
    int areaRow = nameRow+1;
    
    createHeaderCell(ws, nameRow, 0, param.getNamePlusModHumanReadable());
    count++;
    ws.value(areaRow, 0, (double)param.Area);
    count++;
    
    // writing of fragments and rules concerning the head group
    Hashtable<String,CgProbe> headGroupFragments = param.getHeadGroupFragments();
    if (headGroupFragments.size()>0){
      createHeaderCell(ws, count++, 0, LipidomicsConstants.EXCEL_MSN_SECTION_HEAD_FRAGMENTS);
      writeMSnFragmentHeader(ws, count++,LipidomicsConstants.CHAIN_TYPE_NO_CHAIN);
      for (String name : headGroupFragments.keySet()){
        //TODO: ask why number of hydroxylations is -1 at all times?
        writeMSnFragment(ws, count++,name,param.getMSnMzTolerance(),headGroupFragments.get(name),LipidomicsConstants.CHAIN_TYPE_NO_CHAIN,-1);
      }
      Hashtable<String,IntensityRuleVO> headIntRules = param.getHeadIntensityRules();
      if (headIntRules.size()>0){
        createHeaderCell(ws, count++, 0, LipidomicsConstants.EXCEL_MSN_SECTION_HEAD_INTENSITIES);
        writeMSnIntensityHeader(ws, count++);
        Hashtable<String,String> uniqueRules = new Hashtable<String,String>();
        for (IntensityRuleVO ruleVO : headIntRules.values()){
          if (writeMSnIntensity(ws,count,ruleVO,param,uniqueRules,faHydroxyEncoding,lcbHydroxyEncoding)){
            count++;
          }
        }
      }
    }
    
    // writing of fragments and rules concerning the fragments
    Hashtable<String,Hashtable<String,CgProbe>> chainFragments = param.getChainFragments();
    if (chainFragments.size()>0){
      createHeaderCell(ws, count++, 0, LipidomicsConstants.EXCEL_MSN_SECTION_CHAIN_FRAGMENTS);
      writeMSnFragmentHeader(ws, count++, LipidomicsConstants.CHAIN_TYPE_FA_ACYL);
      for (String faName : chainFragments.keySet()){
        Hashtable<String,CgProbe> fragments = chainFragments.get(faName);
        for (String name : fragments.keySet()){
          FattyAcidVO fa = StaticUtils.decodeLipidNameForCreatingCombis(faName);
          String displayName = StaticUtils.getChainFragmentDisplayName(name,fa.getCarbonDbsId());
          writeMSnFragment(ws, count++,displayName,param.getMSnMzTolerance(),fragments.get(name),fa.getChainType(),fa.getOhNumber());
        }
      }
      
      Hashtable<String,Hashtable<String,IntensityChainVO>> chainRules = param.getChainIntensityRules();
      if (chainRules.size()>0){
        createHeaderCell(ws, count++, 0, LipidomicsConstants.EXCEL_MSN_SECTION_CHAIN_INTENSITIES);
        writeMSnIntensityHeader(ws, count++);
        Hashtable<String,String> uniqueRules = new Hashtable<String,String>();
        for (String faName : chainRules.keySet()){
          Hashtable<String,IntensityChainVO> rules = chainRules.get(faName);
          for (IntensityRuleVO rule : rules.values()){
            if (writeMSnIntensity(ws,count,rule,param,uniqueRules,faHydroxyEncoding,lcbHydroxyEncoding)){
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
        createHeaderCell(ws, count++, 0, LipidomicsConstants.EXCEL_MSN_SECTION_POSITION_INTENSITIES+" ("+param.getPositionInsensitiveHumanReadableCombiName(combiName)+")");
        writeMSnIntensityHeader(ws, count++);
        Hashtable<String,String> uniqueRules = new Hashtable<String,String>();
        Vector<IntensityRuleVO> rules = param.getFAsInSequenceAsInRule(combiName); 
        for (IntensityRuleVO rule : rules){
          if (writeMSnIntensity(ws,count,rule,param,uniqueRules,faHydroxyEncoding,lcbHydroxyEncoding)){
            count++;
          }  
        }
      }
    }

    int cellCount = 1;
    
    //TODO: would be nicer to not write out any double bond positions here
//    Vector<DoubleBondPositionVO> assignedDoubleBondPositionVOs = null;
//    if (param.hasOmegaInformation()) {
//      assignedDoubleBondPositionVOs = StaticUtils.getAssignedDoubleBondPositions(param.getOmegaInformation());
//    }
    String identificationString;
    Vector<Object> detected = param.getMSnIdentificationNames();
    Vector<String> newNomenclature = param.getMSnIdentificationNamesWithSNPositions();
    int nameCount = 0;
    for (Object nameObject : detected){
      identificationString = "";
      double area = 0d;
      if (nameObject instanceof Vector){
        area = param.getRelativeIntensity(((Vector<String>)nameObject).get(0))*((double)param.Area);
      // LL         for (String name : (Vector<String>)nameObject){ //
      // LL   identificationString+=name+LipidomicsConstants.CHAIN_COMBI_SEPARATOR_AMBIG_POS_OLD; //
      // LL  }
      // LL  identificationString = identificationString.substring(0,identificationString.length()-1); //
        
        identificationString = newNomenclature.get(nameCount);
        nameCount++;
      }else{
    	    
        String name = (String) nameObject;
      // LL  identificationString = name;
        identificationString = newNomenclature.get(nameCount);
        nameCount++;
        if (param.getStatus()==LipidomicsMSnSet.HEAD_GROUP_DETECTED) area = (double)param.Area;
        else area = param.getRelativeIntensity(name)*((double)param.Area);
     }
      
      /** write double bond position assignment if available */
//      if (assignedDoubleBondPositionVOs != null) {
//        Vector<DoubleBondPositionVO> doubleBondAssignmentsOfMolecularSpecies = StaticUtils.getDoubleBondAssignmentsOfMolecularSpecies(
//            assignedDoubleBondPositionVOs, identificationString);
//        if (doubleBondAssignmentsOfMolecularSpecies.size() == 1) {
//          identificationString = doubleBondAssignmentsOfMolecularSpecies.get(0).getDoubleBondPositionsHumanReadable();
//        } else {
//          /** only one label can be assigned for each species */
//        }
//      }

      createHeaderCell(ws, nameRow, cellCount, identificationString);
      ws.value(areaRow, cellCount, area);
      cellCount++;
    }
    
    //writing the retention times of the used spectra
    List<Integer> msLevels = new ArrayList<Integer>(param.getMsnRetentionTimes().keySet());
    Collections.sort(msLevels);
    for (int msLevel : msLevels){
      LinkedHashMap<Integer,Float> rts = param.getMsnRetentionTimes().get(msLevel);
      String rtString = "";
      for (Integer scanNr : rts.keySet()) rtString += scanNr+"="+rts.get(scanNr)+";";
      if (rtString.length()>0) rtString = rtString.substring(0,rtString.length()-1);
      createHeaderCell(ws, nameRow, cellCount, "MS"+msLevel+" scan RTs");
      ws.value(areaRow, cellCount, rtString);
      cellCount++;
    }
    count++;
    return count;
  }
  
  
  /**
   * Creates a header cell
   * @param ws the Excel worksheet to be written into
   * @param row Excel row that shall be used for writing
   * @param column Excel column that shall be used for writing
   * @param value the value to be written in the header cell
   */
  private static void createHeaderCell(Worksheet ws, int row, int column, String value) {
    ws.value(row, column, value);
    ws.style(row, column).bold().horizontalAlignment("center").fontName("Arial").fontSize(12).set();
  }
  
  
  /**
   * Writes the header row for subsequent fragment evidence
   * @param ws the Excel worksheet to be written into
   * @param row Excel row that shall be used for writing
   * @param chainType allowed chain types: LipidomicsConstants.CHAIN_TYPE_NO_CHAIN|CHAIN_TYPE_FA_ACYL|CHAIN_TYPE_FA_ALKYL|CHAIN_TYPE_FA_ALKENYL|CHAIN_TYPE_LCB
   */
  private static void writeMSnFragmentHeader(Worksheet ws, int row, short chainType){
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_NAME, LipidomicsConstants.EXCEL_MSN_FRAGMENT_NAME);
    int add = 0;
    if (chainType!=LipidomicsConstants.CHAIN_TYPE_NO_CHAIN) {
      createHeaderCell(ws, row, MSN_ROW_FRAGMENT_OH, LipidomicsConstants.EXCEL_MSN_FRAGMENT_OH);
      createHeaderCell(ws, row, MSN_ROW_FRAGMENT_CHAIN_TYPE, LipidomicsConstants.EXCEL_MSN_FRAGMENT_CHAIN_TYPE);      
      add = 2;
    }
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_FORMULA+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_FORMULA);      
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_MSLEVEL+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_MSLEVEL);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_CHARGE+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_CHARGE);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_MZ+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_MZ_TOLERANCE+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_TOLERANCE);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_AREA+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_AREA);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_PEAK+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_PEAK);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_TIME_LOWER+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_TIME_LOWER);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_TIME_UPPER+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_TIME_UPPER);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_MZ_LOWER+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_LOWER);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_MZ_UPPER+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_UPPER);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_ELLIPSE_TIME+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_TIME);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_ELLIPSE_MZ+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_MZ);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_ELLIPSE_TIME_RANGE+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_TIME_RANGE);
    createHeaderCell(ws, row, MSN_ROW_FRAGMENT_ELLIPSE_MZ_RANGE+add, LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_MZ_RANGE);
  }
  
  
  /**
   * Writes a line containing found fragments
   * @param ws the Excel worksheet to be written into
   * @param row Excel row that shall be used for writing
   * @param name display name of the fragment
   * @param mzTolerance m/z tolerance for the identification
   * @param probe VO containing information about the identified fragment
   * @param chainType allowed chain types: LipidomicsConstants.CHAIN_TYPE_NO_CHAIN|CHAIN_TYPE_FA_ACYL|CHAIN_TYPE_FA_ALKYL|CHAIN_TYPE_FA_ALKENYL|CHAIN_TYPE_LCB
   * @param oh the number of hydroxylations on a chain
   * @throws ChemicalFormulaException thrown when something is wrong with the chemical formula
   */
  private static void writeMSnFragment(Worksheet ws, int row, String name, float mzTolerance, CgProbe probe,short chainType, int oh) throws ChemicalFormulaException{
    ws.value(row, MSN_ROW_FRAGMENT_NAME, name);
    int add = 0;
    if (chainType!=LipidomicsConstants.CHAIN_TYPE_NO_CHAIN) {
      ws.value(row, MSN_ROW_FRAGMENT_OH, oh);
      ws.value(row, MSN_ROW_FRAGMENT_CHAIN_TYPE, StaticUtils.getHumanReadableChainType(chainType));
      add = 2;
    }    
    //TODO: this is only here because of a damaged Alex123 file - delete in future version!
    String formula = "";
    if (probe.getFormula()!=null) formula = StaticUtils.getFormulaInHillNotation(StaticUtils.categorizeFormula(probe.getFormula().trim()),true);
    ws.value(row, MSN_ROW_FRAGMENT_FORMULA+add, formula);
    ws.value(row, MSN_ROW_FRAGMENT_MSLEVEL+add, probe.getMsLevel());
    ws.value(row, MSN_ROW_FRAGMENT_CHARGE+add, probe.Charge);
    ws.value(row, MSN_ROW_FRAGMENT_MZ+add, (double)probe.Mz);
    ws.value(row, MSN_ROW_FRAGMENT_MZ_TOLERANCE+add, (double)mzTolerance);
    ws.value(row, MSN_ROW_FRAGMENT_AREA+add, (double)probe.Area);
    ws.value(row, MSN_ROW_FRAGMENT_PEAK+add, (double)probe.Peak);
    ws.value(row, MSN_ROW_FRAGMENT_TIME_LOWER+add, (double)probe.LowerValley);
    ws.value(row, MSN_ROW_FRAGMENT_TIME_UPPER+add, (double)probe.UpperValley);
    
    if (probe instanceof Probe3D){
      Probe3D probe3D = (Probe3D)probe;
      ws.value(row, MSN_ROW_FRAGMENT_MZ_LOWER+add, (double)probe3D.LowerMzBand);
      ws.value(row, MSN_ROW_FRAGMENT_MZ_UPPER+add, (double)probe3D.UpperMzBand);
      ws.value(row, MSN_ROW_FRAGMENT_ELLIPSE_TIME+add, (double)probe3D.getEllipseTimePosition());
      ws.value(row, MSN_ROW_FRAGMENT_ELLIPSE_MZ+add, (double)probe3D.getEllipseMzPosition());
      ws.value(row, MSN_ROW_FRAGMENT_ELLIPSE_TIME_RANGE+add, (double)probe3D.getEllipseTimeStretch());
      ws.value(row, MSN_ROW_FRAGMENT_ELLIPSE_MZ_RANGE+add, (double)probe3D.getEllipseMzStretch());
    }
  }
  
  
  /**
   * Writes the header row for subsequent intensity evidence
   * @param ws the Excel worksheet to be written into
   * @param row Excel row that shall be used for writing
   */
  private static void writeMSnIntensityHeader(Worksheet ws, int row){
    ws.value(row, MSN_ROW_INTENSITY_RULE, LipidomicsConstants.EXCEL_MSN_INTENSITY_RULE);
    ws.value(row, MSN_ROW_INTENSITY_ORIGINAL, LipidomicsConstants.EXCEL_MSN_INTENSITY_ORIGINAL);
    ws.value(row, MSN_ROW_INTENSITY_VALUES, LipidomicsConstants.EXCEL_MSN_INTENSITY_VALUES);
    ws.value(row, MSN_ROW_INTENSITY_MISSED, LipidomicsConstants.EXCEL_MSN_INTENSITY_MISSED);
  }
  
  /**
   * Writes a line containing the intensity evidence
   * @param ws the Excel worksheet to be written into
   * @param row Excel row that shall be used for writing
   * @param ruleVO IntensityRuleVO to be written
   * @param param the lipid identification containing MSn evidence
   * @param uniqueRules Hashtable containing already written rules
   * @param faHydroxyEncoding the character encoding of the number of hydroxylation sites for the
   * @param lcbHydroxyEncoding the character encoding of the number of hydroxylation sites for the LCB
   * @return true if data was written to the Excel sheet
   * @throws RulesException if something is not possible
   * @throws LipidCombinameEncodingException thrown when a lipid combi ID (containing type and OH number) cannot be decoded
   */
  private static boolean writeMSnIntensity(Worksheet ws, int row, IntensityRuleVO ruleVO, LipidomicsMSnSet param, Hashtable<String,String> uniqueRules,
     HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding) throws RulesException, LipidCombinameEncodingException{
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
    
    ws.value(row, MSN_ROW_INTENSITY_RULE, ruleInterpretation);
    ws.value(row, MSN_ROW_INTENSITY_ORIGINAL, rule);
    ws.value(row, MSN_ROW_INTENSITY_VALUES, valueInterpretation);
    ws.value(row, MSN_ROW_INTENSITY_MISSED, missed);
    
    return true;
  }
  
}
