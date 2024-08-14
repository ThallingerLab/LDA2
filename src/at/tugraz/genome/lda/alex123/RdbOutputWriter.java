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

package at.tugraz.genome.lda.alex123;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.alex123.vos.TargetlistEntry;
import at.tugraz.genome.lda.alex123.vos.TargetlistFloatStringVO;
import at.tugraz.genome.lda.analysis.ComparativeAnalysis;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RdbWriterException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.Probe3D;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * 
 * @author Juergen Hartler
 *
 * Class for writing the LDA resuls in ALEX123 format
 */
public class RdbOutputWriter
{
  
  /** these are the column headers */
  public final static String FILENAME_COLUMN = "RAW_ID";
  public final static String MACHINE_COLUMN = "Machine";
  public final static String DETECTOR_COLUMN = "Detector";
  public final static String POLARITY_COLUMN = "Polarity";
  public final static String EVIDENCE_LEVEL_COLUMN = "MS evidence level";
  public final static String SPECIES_COLUMN = "Lipid species";
  public final static String RT_GROUPING_COLUMN = "RT group";
  public final static String MOLECULAR_SPECIES_COLUMN = "Molecular lipid species";
  public final static String MOLECULAR_SPECIES_COLUMN_RANK = "Mol. species rank";
  public final static String LIPID_CLASS_COLUMN = "Lipid class";
  public final static String LEVEL_PREFIX = "MS";
  public final static String FRAGMENT_NAMES_LEVEL_COLUMN_PART = " Fragment names";
  public final static String FRAGMENT_MZ_LEVEL_COLUMN_PART = " Fragment m/z";
  public final static String NUMBER_MS_SCANS_PREFIX = "#";
  public final static String NUMBER_MS_SCANS_SUFFIX= " scans";
  public final static String MSN_RTS_COLUMN_PART = " RTs";
  public final static String PRECURSOR_COLUMN_PART = " precursor m/z";
  public final static String ACTIVATION_COLUMN_PART = " activation";
  public final static String RT_COLUMN = "Retention time";
  public final static String ISOTOPES_COLUMN = "Isotope clusters";
  public final static String LOWER_VALLEY_COLUMN = "Lower valley";
  public final static String UPPER_VALLEY_COLUMN = "Upper valley";
  public final static String LOWER_VALLEY_10_COLUMN = "Lower valley 10%";
  public final static String UPPER_VALLEY_10_COLUMN = "Upper valley 10%";
  public final static String LOWER_VALLEY_50_COLUMN = "Lower valley 50%";
  public final static String UPPER_VALLEY_50_COLUMN = "Upper valley 50%";
  public final static String LOW_MZ_COLUMN = "Low m/z";
  public final static String UP_MZ_COLUMN = "Up m/z";
  public final static String LOW_MZ_10_COLUMN = "Low m/z 10%";
  public final static String UP_MZ_10_COLUMN = "Up m/z 10%";
  public final static String LOW_MZ_50_COLUMN = "Low m/z 50%";
  public final static String UP_MZ_50_COLUMN = "Up m/z 50%";
  public final static String PEAK_AREA_COLUMN = "Peak area (total)";
  public final static String ISO_PEAK_AREA_COLUMN_PART = "Peak area M";
  public final static String INTENSITY_MAX_COLUMN_PART = "Max intensity";
  public final static String TARGET_MZ_COLUMN = "Target m/z";
  public final static String MEASURED_MZ_COLUMN = "Measured m/z";
  public final static String MZ_CORRECTION_COLUMN = "Applied m/z correction";
  public final static String MZ_CORRECTED_COLUMN = "Corrected measured m/z";
  public final static String MZ_OFFSET_COLUMN = "m/z offset";
  public final static String PPM_ERROR_COLUMN = "ppm error";
  public final static String ADDUCT_COLUMN = "Adduct";
  public final static String LIPID_ID_COLUMN = "Lipid ID";
  public final static String LIPID_CATEGORY_COLUMN = "Lipid category";
  public final static String CONFLICTS_COLUMN = "Conflicts";
  public final static String CHARGE_COLUMN = "Charge";
  public final static String C_INDEX_COLUMN = "C index of lipid species";
  public final static String DB_INDEX_COLUMN = "DB index of lipid species";
  public final static String OH_INDEX_COLUMN = "OH index of lipid species";
  public final static String SUM_SPECIES_COLUMN = "Sum composition of lipid species";
  public final static String SUM_FORMULA_COLUMN = "Sum formula of lipid species";
  /**the tab delimiter*/
  public final static String TAB = "\t";

  
  
  /** the prefix for the internal standard*/
//  private String internalStandardPrefix_;
  /** the prefix for the external standard*/
//  private String externalStandardPrefix_;
  
  /**
   * constructor without setting any standard prefixes
   */
  public RdbOutputWriter(){
    this(null,null);
  }
  
  /**
   * constructor setting standard prefixes
   * @param internalStandardPrefix prefix for the internal standard
   * @param externalStandardPrefix prefix for the external standard
   */
  public RdbOutputWriter(String internalStandardPrefix, String externalStandardPrefix){
//    internalStandardPrefix_ = internalStandardPrefix;
//    externalStandardPrefix_ = externalStandardPrefix;
  }
  
  /**
   * starts the writing procedure - this method should be used directly after LDA quantitation
   * @param fileName the name of the RDB file
   * @param results the results to be written - every entry is one LDA results file
   * @param classSequence the sequence in which the classes shall be written
   * @param analyteSequence the sequence in which the analytes of the classes shall be written
   * @param quantObjects the quantification objects providing more information (only for Alex123 target lists)
   * @throws RdbWriterException if something is wrong with the writing
   * @throws ChemicalFormulaException if something is wrong with the chemical formulae
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public void write(String fileName, Vector<QuantificationResult> results, LinkedHashMap<String,Integer> classSequence,
      Hashtable<String,Vector<String>> analyteSequence, Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects) throws RdbWriterException, ChemicalFormulaException, LipidCombinameEncodingException {
    //the third variable is not required - only backward compatibility that is not necessary for this method
    //the sixth and the seventh parameter is not required by this method - all of the molecules and isotopes have to be written
    write(fileName,results,null,classSequence,analyteSequence,null,null,null,quantObjects,false);
  }
  
  /**
   * starts the writing procedure
   * @param fileNameOut name of the RDB file
   * @param results the results to be written - every entry is one LDA results file
   * @param expLookup the name of the excel files - to be used only by export via heat map
   * @param classSequence the sequence in which the classes shall be written
   * @param correctAnalyteSequence the sequence in which the analytes of the classes shall be written
   * @param acceptedMolecules the molecules that are accepted by the GUI - to be used only by export via heat map
   * @param analysisModule the comparative analysis module - containing information about the file comparisons by the heat map (only for exports by LDA files)
   * @param maxIsotopes the maximum isotopes settings - to be used only by export via heat map
   * @param quantObjects the quantification objects providing more information (only for Alex123 target lists)
   * @param alexRtGrouper is this called from the AlexRtGrouper
   * @throws RdbWriterException if something is wrong with the writing
   * @throws ChemicalFormulaException if something is wrong with the chemical formulae
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  @SuppressWarnings("unchecked")
  private void write(String fileNameOut, Vector<QuantificationResult> results, Hashtable<Integer,String> expLookup,
      LinkedHashMap<String,Integer> classSequence, Hashtable<String,Vector<String>> correctAnalyteSequence,
      Hashtable<String,Hashtable<String,String>> acceptedMolecules, ComparativeAnalysis analysisModule,
      Hashtable<String,Integer> maxIsotopes, Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects,
      boolean alexRtGrouper) throws RdbWriterException,
      ChemicalFormulaException, LipidCombinameEncodingException {
    boolean detectorColumn = false;
    boolean polarityColumn = false;
    boolean lipidIdColumn = false;
    boolean lipidCategoryColumn = false;
    boolean conflictsColumn = false;
    boolean ohIndexColumn = false;
    boolean containsMSnInformation = false;
    int highestIsotopeNumber = 0;
    int highestMSLevel = 1;
    boolean isCompModulePresent = analysisModule!=null;
    boolean isRtGrouped = (isCompModulePresent && analysisModule.isRtGrouped());
    
    boolean foundAlexQuantObjects = false;
    //checking if the Alex123 target list specific columns can be filled out
    if (quantObjects!=null){
      for (Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass : quantObjects.values()){
        for (Hashtable<String,QuantVO> quantsOfAnalyte : quantsOfClass.values()){
          for (QuantVO quant : quantsOfAnalyte.values()){
            if (!(quant instanceof TargetlistEntry)) continue;
            if (!isCompModulePresent)
              foundAlexQuantObjects = true;
            TargetlistEntry entry = (TargetlistEntry) quant;
            if (entry.getDetector()!=null) detectorColumn = true;
            if (entry.getPolarity()!=null) polarityColumn = true;
            if (entry.getId()!=null) lipidIdColumn = true;
            if (entry.getCategory()!=null) lipidCategoryColumn = true;
            if (entry.getConflicts()!=null) conflictsColumn = true;
            if (entry.getOhNumber()>=0) ohIndexColumn = true;
            if (!isCompModulePresent) {
              if (entry.getMsnFragments()!=null && entry.getMsnFragments().size()>0) containsMSnInformation = true;
              for (Integer msLevel : entry.getMsnFragments().keySet()){
                if (msLevel>highestMSLevel) highestMSLevel = msLevel;
              }
            } else {
              if (entry.getMsLevel()>highestMSLevel) highestMSLevel = entry.getMsLevel();
            }
//            if (detectorColumn && polarityColumn && lipidIdColumn && lipidCategoryColumn && conflictsColumn && ohIndexColumn && containsMSnInformation)
//              break;
          }
//          if (detectorColumn && polarityColumn && lipidIdColumn && lipidCategoryColumn && conflictsColumn && ohIndexColumn && containsMSnInformation)
//            break;
        }
//        if (detectorColumn && polarityColumn && lipidIdColumn && lipidCategoryColumn && conflictsColumn && ohIndexColumn && containsMSnInformation)
//          break;
      }
    }
    // get the highest MS-level    
    if (!foundAlexQuantObjects){
      for (QuantificationResult result : results){
        for (String className : result.getIdentifications().keySet()){
          for (LipidParameterSet set : result.getIdentifications().get(className)){
            if (!(set instanceof LipidomicsMSnSet) /*|| quantObjects!=null*/) continue;
            if (quantObjects!=null && (detectorColumn||polarityColumn||lipidIdColumn||lipidCategoryColumn||conflictsColumn||ohIndexColumn))
              containsMSnInformation = true;
            LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
            for (CgProbe headFrag : msn.getHeadGroupFragments().values()){
              if (headFrag.getMsLevel()>highestMSLevel) highestMSLevel = headFrag.getMsLevel();
            }
            for (Hashtable<String,CgProbe> chainFrags : msn.getChainFragments().values()){
              for (CgProbe chainFrag : chainFrags.values()){
                if (chainFrag.getMsLevel()>highestMSLevel) highestMSLevel = chainFrag.getMsLevel();
              }
            }
          }

        }
      }
    }
    
    //get the highest isotope number
    for (QuantificationResult result : results){
      for (String className : result.getIdentifications().keySet()){
        if (maxIsotopes!=null){
          int maxIso = maxIsotopes.get(className);
          if (maxIso>highestIsotopeNumber) highestIsotopeNumber = maxIso;
        }else{
          for (LipidParameterSet set : result.getIdentifications().get(className)){
        	  if (set == null)
        		  continue;
        	  if ((set.getIsotopicProbes().size()-1)>highestIsotopeNumber) highestIsotopeNumber = (set.getIsotopicProbes().size()-1);
          }
        }
      }
    }
    
    try {
      BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(fileNameOut));
      StringBuilder headerLine = new StringBuilder();
      headerLine.append(FILENAME_COLUMN+"\t"+MACHINE_COLUMN);
      if (detectorColumn) headerLine.append("\t"+DETECTOR_COLUMN);
      if (polarityColumn) headerLine.append("\t"+POLARITY_COLUMN);
      headerLine.append("\t"+EVIDENCE_LEVEL_COLUMN+"\t"+SPECIES_COLUMN);
      if (isRtGrouped)
        headerLine.append("\t"+RT_GROUPING_COLUMN);
      headerLine.append("\t"+MOLECULAR_SPECIES_COLUMN+"\t"+MOLECULAR_SPECIES_COLUMN_RANK+"\t"+LIPID_CLASS_COLUMN);
      for (int i=2; i<=highestMSLevel; i++){
        String levelName = LEVEL_PREFIX+i;
         headerLine.append("\t"+levelName+FRAGMENT_NAMES_LEVEL_COLUMN_PART+"\t"+levelName+FRAGMENT_MZ_LEVEL_COLUMN_PART+"\t"+NUMBER_MS_SCANS_PREFIX+levelName+NUMBER_MS_SCANS_SUFFIX+"\t"+levelName+MSN_RTS_COLUMN_PART);
         if (containsMSnInformation) headerLine.append("\t"+levelName+PRECURSOR_COLUMN_PART+"\t"+levelName+ACTIVATION_COLUMN_PART);
      }
      headerLine.append("\t"+RT_COLUMN+"\t"+ISOTOPES_COLUMN+"\t"+LOWER_VALLEY_COLUMN+"\t"+UPPER_VALLEY_COLUMN+"\t"+LOWER_VALLEY_10_COLUMN+"\t"+UPPER_VALLEY_10_COLUMN +"\t"+LOWER_VALLEY_50_COLUMN
          + "\t"+UPPER_VALLEY_50_COLUMN+"\t"+LOW_MZ_COLUMN+"\t"+UP_MZ_COLUMN+"\t"+LOW_MZ_10_COLUMN+"\t"+UP_MZ_10_COLUMN+"\t"+LOW_MZ_50_COLUMN+"\t"+UP_MZ_50_COLUMN+"\t"+PEAK_AREA_COLUMN);
      for (int i=0; i!=(highestIsotopeNumber+1); i++) headerLine.append("\t"+ISO_PEAK_AREA_COLUMN_PART+i);
      headerLine.append("\t"+INTENSITY_MAX_COLUMN_PART+"\t"+TARGET_MZ_COLUMN+"\t"+MEASURED_MZ_COLUMN+"\t"+MZ_CORRECTION_COLUMN+"\t"+MZ_CORRECTED_COLUMN+"\t"+MZ_OFFSET_COLUMN+"\t"+PPM_ERROR_COLUMN+"\t"+ADDUCT_COLUMN);
      if (lipidIdColumn)  headerLine.append("\t"+LIPID_ID_COLUMN);
      if (lipidCategoryColumn) headerLine.append("\t"+LIPID_CATEGORY_COLUMN);
      if (conflictsColumn) headerLine.append("\t"+CONFLICTS_COLUMN );
      headerLine.append("\t"+CHARGE_COLUMN+"\t"+C_INDEX_COLUMN+"\t"+DB_INDEX_COLUMN);
      if (ohIndexColumn) headerLine.append("\t"+OH_INDEX_COLUMN);
      headerLine.append("\t"+SUM_SPECIES_COLUMN+"\t"+SUM_FORMULA_COLUMN+"\n");
      stream.write(headerLine.toString().getBytes());      
            
      for (int i=0; i!=results.size(); i++) {
        QuantificationResult result = results.get(i);
        String msMachine = "";
        float mzCorrection = 0f;
        
        String expName = null;
        if (expLookup!=null)
          expName = expLookup.get(result.hashCode());
        String resultsExcelName = null;
        if (result.getConstants()!=null && result.getConstants().getRawFileName()!=null){
          resultsExcelName = result.getConstants().getRawFileName();
          if (resultsExcelName.endsWith(".chrom")) resultsExcelName = resultsExcelName.substring(0,resultsExcelName.length()-".chrom".length());
        //the file lookup is for backward compatibility, and can be called by the heat map only
        }else resultsExcelName = analysisModule.getFullFilePath(expName).getName();
        
        if (result.getConstants()!=null){
          mzCorrection = (float)result.getConstants().getShift();
          msMachine = result.getConstants().getMSMachine();
        }
        //check if all classes are in the class sequence
        LinkedHashMap<String,Integer> classSequenceFull = new LinkedHashMap<String,Integer>();
        if (classSequence!=null)
          classSequenceFull = new LinkedHashMap<String,Integer>(classSequence);
        for (String className : result.getIdentifications().keySet()){
          if (classSequence!=null && !classSequence.containsKey(className)){
            System.out.println("ATTENTION: the class "+className+" is not in your target list!");
            classSequenceFull.put(className, 1);
          }
        }
        for (String className : classSequenceFull.keySet()){
          if (!result.getIdentifications().containsKey(className)) continue;
          Hashtable<String,String> selectedMolHash = null;
          Hashtable<String,Vector<String>> fromSpeciesToSpeciesWithRt = new Hashtable<String,Vector<String>>();
          if (acceptedMolecules!=null){
            if (acceptedMolecules.containsKey(className)) selectedMolHash = acceptedMolecules.get(className);
            else selectedMolHash = new Hashtable<String,String>();
          } else if (isCompModulePresent) {
            selectedMolHash = new Hashtable<String,String>();
            for (String molName : analysisModule.getAllMoleculeNames().get(className))
              selectedMolHash.put(molName, molName);
          }
          if (selectedMolHash!=null) {
            Vector<String> molsWithRt;
            for (String molWithRt : selectedMolHash.keySet()) {
              String molName;
              if (analysisModule.isISorES(className, molWithRt) || LipidomicsConstants.isShotgun() == 1)
                molName = new String(molWithRt);
              else
                molName = molWithRt.substring(0,molWithRt.lastIndexOf("_"));
              molsWithRt = new Vector<String>();
              if (fromSpeciesToSpeciesWithRt.containsKey(molName))
                molsWithRt = fromSpeciesToSpeciesWithRt.get(molName);
              molsWithRt.add(molWithRt);
              fromSpeciesToSpeciesWithRt.put(molName, molsWithRt);
            }            
          }
          
          Vector<LipidParameterSet> classResultsFile = result.getIdentifications().get(className);
          Hashtable<String,Vector<LipidParameterSet>> foundSpecies = new Hashtable<String,Vector<LipidParameterSet>>();
          for (LipidParameterSet set : classResultsFile){
        	if (set == null)
        		continue;
            if (isRtGrouped) {
              if (!fromSpeciesToSpeciesWithRt.containsKey(set.getNameStringWithoutRt()))
                continue;
              Vector<String> possSpeciesDiffRt = fromSpeciesToSpeciesWithRt.get(set.getNameStringWithoutRt());
              for (String molWithRt : possSpeciesDiffRt) {
                ResultAreaVO vo = analysisModule.getResultAreaVO(className, molWithRt, expName);
                if (vo==null)
                  continue;
                if (vo.belongsRtToThisAreaVO(set.getRt(), set.getModificationName())) {
                  Vector<LipidParameterSet> ofOneAnalyte = new Vector<LipidParameterSet>();
                  if (foundSpecies.containsKey(molWithRt)) ofOneAnalyte = foundSpecies.get(molWithRt);
                  ofOneAnalyte.add(set);
                  foundSpecies.put(molWithRt, ofOneAnalyte);                  
                  break;
                }
              }
            }else{
              Vector<LipidParameterSet> ofOneAnalyte = new Vector<LipidParameterSet>();
              if (foundSpecies.containsKey(set.getNameStringWithoutRt())) ofOneAnalyte = foundSpecies.get(set.getNameStringWithoutRt());
              ofOneAnalyte.add(set);
              foundSpecies.put(set.getNameStringWithoutRt(), ofOneAnalyte);
            }
          }
          Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = null;
          if (quantObjects!=null && quantObjects.containsKey(className)) quantsOfClass = quantObjects.get(className);
          

          int maxIsotope = Integer.MAX_VALUE-1;
          if (maxIsotopes!=null) maxIsotope = maxIsotopes.get(className); 
          
          for (String molName: correctAnalyteSequence.get(className)){
            if ((selectedMolHash==null || selectedMolHash.containsKey(molName)) && foundSpecies.containsKey(molName)){
              Hashtable<String,QuantVO> quantsOfAnalyte = null;
              String molNameWoRt = molName;
              if (isRtGrouped && !analysisModule.isISorES(className, molName))
                molNameWoRt = /*className+" "+*/molName.substring(0,molName.lastIndexOf("_"));
              if (quantsOfClass!=null && quantsOfClass.containsKey(molNameWoRt)) quantsOfAnalyte = quantsOfClass.get(molNameWoRt);
              String speciesId = "";
              String classNameToWrite = "";
              int cIndexSpecies = 0;
              int dbIndexSpecies = 0;
              String sumComposition = "";
//              if (internalStandardPrefix_!=null && internalStandardPrefix_.length()>0 && molName.startsWith(internalStandardPrefix_)){
//                speciesId = "IS " + className+" "+ molName.substring(internalStandardPrefix_.length());
//                classNameToWrite = "IS "+className;
//              }else if (externalStandardPrefix_!=null && externalStandardPrefix_.length()>0 && molName.startsWith(externalStandardPrefix_)){
//                speciesId = "ES " + className + " " + molName.substring(externalStandardPrefix_.length());
//                classNameToWrite = "ES "+className;
//              }else{
                speciesId = className+" "+molName;
                classNameToWrite = className;
                int startChar = 0;
                char[] chars = molName.toCharArray();
                for (int j=0; j!=chars.length;j++){
                  startChar = j;
                  if (Character.isDigit(chars[j])) break;
                }
                String indices = molName.substring(startChar);
                String[] cAndDbs = indices.split(":");
                try{
                  cIndexSpecies = Integer.parseInt(cAndDbs[0]);
                  dbIndexSpecies = Integer.parseInt(cAndDbs[1]);
                  sumComposition = indices;
                }catch(NumberFormatException nfx){
                  
                }
//              }
              
              for (LipidParameterSet set : foundSpecies.get(molName)){
                TargetlistEntry quantObject = null;
                if (quantsOfAnalyte!=null && quantsOfAnalyte.containsKey(set.getModificationName()) && (quantsOfAnalyte.get(set.getModificationName()) instanceof TargetlistEntry))
                  quantObject = (TargetlistEntry)quantsOfAnalyte.get(set.getModificationName());
                String msLevel = "ms";
                String molSpeciesId = "";
                float totalArea = 0f;
                Hashtable<Integer,Float> isotopicAreas = new Hashtable<Integer,Float>();
                int maxIterate = maxIsotope+1;
                if (maxIterate>set.getIsotopicProbes().size()) maxIterate = set.getIsotopicProbes().size();
                float measuredMz = 0f;
                float apexIntensity = 0f;
                float lowerValley = Float.MAX_VALUE;
                float lowerValley10Pc = Float.MAX_VALUE;
                float lowerValley50Pc = Float.MAX_VALUE;
                float upperValley = 0f;
                float upperValley10Pc = 0f;
                float upperValley50Pc = 0f;
                float lowerMz = Float.MAX_VALUE;
                float lowerMz10Pc = Float.MAX_VALUE;
                float lowerMz50Pc = Float.MAX_VALUE;
                float upperMz = 0f;
                float upperMz10Pc = 0f;
                float upperMz50Pc = 0f;
                              
                String apexString = "";
                String lowerValleyString = "";
                String lowerValley10PcString = "";
                String lowerValley50PcString = "";
                String upperValleyString = "";
                String upperValley10PcString = "";
                String upperValley50PcString = "";
                String lowerMzString = "";
                String lowerMz10PcString = "";
                String lowerMz50PcString = "";
                String upperMzString = "";
                String upperMz10PcString = "";
                String upperMz50PcString = "";
                  
                String detector = "";
                String polarity = "";
                String speciesToWrite = new String(speciesId);
                String groupingId = "";
                String lipidId = "";
                String lipidCategory = "";
                String conflicts = "";
                int charge = set.getCharge();
                int cIndexSpeciesToWrite = cIndexSpecies;
                int dbIndexSpeciesToWrite = dbIndexSpecies;
                Integer ohIndex = null;
                String sumCompositionToWrite = new String(sumComposition);
                String sumFormula = set.getChemicalFormula();
                if (quantObject!=null){
                  if (quantObject.getDetector()!=null) detector = quantObject.getDetector();
                  if (quantObject.getPolarity()!=null) polarity = quantObject.getPolarity();
                  speciesToWrite = quantObject.getSpecies();
                  classNameToWrite = quantObject.getOriginalClassName();
                  if (quantObject.getId()!=null) lipidId = quantObject.getId();
                  if (quantObject.getCategory()!=null) lipidCategory = quantObject.getCategory();
                  if (quantObject.getConflicts()!=null) conflicts = quantObject.getConflicts();
                  if (quantObject.getPolarity()!=null && quantObject.getPolarity().equalsIgnoreCase("-"))
                    charge = charge*-1;
                  if (quantObject.getOhNumber()>=0) ohIndex = quantObject.getOhNumber();
                  if (quantObject.getCarbonNumber()>=0)
                    cIndexSpeciesToWrite = quantObject.getCarbonNumber();
                  if (quantObject.getDbNumber()>=0)
                    dbIndexSpeciesToWrite = quantObject.getDbNumber();
                  if (quantObject.getSumComposition()!=null)
                    sumCompositionToWrite=quantObject.getSumComposition();
                  else
                    sumCompositionToWrite = "";
                  sumFormula = quantObject.getOriginalSumFormula();
                }
                if (isCompModulePresent && !alexRtGrouper && quantObject==null) {
                  sumFormula = set.getChemicalFormula();
                  Hashtable<String,Integer> categorized = StaticUtils.categorizeFormula(sumFormula);
                  if (Settings.useAlex() && result.getConstants()!=null && result.getConstants().isAlexTargetlist()){
                    Hashtable<String,String> isoLookup = Settings.getAlexIsoLookup();
                    for (String alexIso : isoLookup.keySet()){
                      String ldaIso = isoLookup.get(alexIso);
                      if (categorized.containsKey(ldaIso)){
                        categorized.put(alexIso, categorized.get(ldaIso));
                        categorized.remove(ldaIso);
                      }
                    }
                  }
                  sumFormula = StaticUtils.getFormulaInHillNotation(categorized, false);
                }
                if (isRtGrouped) {
                  groupingId = molName.substring(molName.lastIndexOf("_")+1);
                  if (speciesToWrite.endsWith("_"+groupingId))
                    speciesToWrite = speciesToWrite.substring(0,speciesToWrite.lastIndexOf("_"));
                  if (analysisModule.isISorES(className, molName))
                    groupingId = "";
                }
                  
                
                for (int j=0; j!=maxIterate; j++){
                  for (CgProbe probe : set.getIsotopicProbes().get(j)){
                    totalArea += probe.Area;
                    float isoArea = 0f;
                    if (isotopicAreas.containsKey(j)) isoArea = isotopicAreas.get(j);
                    isoArea += probe.Area;
                    isotopicAreas.put(j, isoArea);
                    if (j==0){
                      measuredMz += probe.Mz;
                      float currentLV = probe.LowerValley; 
                      float currentLV10 = Float.MAX_VALUE; 
                      float currentLV50 = Float.MAX_VALUE; 
                      float currentUV = probe.UpperValley; 
                      float currentUV10 = 0f; 
                      float currentUV50 = 0f; 
                        
                      float currentLMz = probe.Mz-probe.LowerMzBand; 
                      float currentUMz = probe.Mz+probe.UpperMzBand; 

                        
                      if (probe.getLowerValley10()!=null){
                        currentLV10 = probe.getLowerValley10();
                        currentLV50 = probe.getLowerValley50();
                        currentUV10 = probe.getUpperValley10();
                        currentUV50 = probe.getUpperValley50();
                      }
                      if (set.getLowerRtHardLimit()>=0){
                        if (currentLV<set.getLowerRtHardLimit()) currentLV=set.getLowerRtHardLimit();
                        if (currentLV10<set.getLowerRtHardLimit()) currentLV10=set.getLowerRtHardLimit();
                        if (currentLV50<set.getLowerRtHardLimit()) currentLV50=set.getLowerRtHardLimit();
                      }
                      if (set.getUpperRtHardLimit()>=0){
                        if (currentUV>set.getUpperRtHardLimit()) currentUV=set.getUpperRtHardLimit();
                        if (currentUV10>set.getUpperRtHardLimit()) currentUV10=set.getUpperRtHardLimit();
                        if (currentUV50>set.getUpperRtHardLimit()) currentUV50=set.getUpperRtHardLimit();
                      }
                        
                      if (currentLV<lowerValley) lowerValley = currentLV;
                      if (currentUV>upperValley) upperValley = currentUV;
                      if (currentLMz<lowerMz) lowerMz = currentLMz;
                      if (currentUMz>upperMz) upperMz = currentUMz;
                        
                      if (probe.getLowerValley10()!=null){
                        if (probe.getApexIntensity()>apexIntensity) apexIntensity = probe.getApexIntensity();
                        if (currentLV10<lowerValley10Pc) lowerValley10Pc = currentLV10;
                        if (currentLV50<lowerValley50Pc) lowerValley50Pc = currentLV50;
                        if (currentUV10>upperValley10Pc) upperValley10Pc = currentUV10;
                        if (currentUV50>upperValley50Pc) upperValley50Pc = currentUV50;
                      }
                      if (probe instanceof Probe3D){
                        Probe3D probe3D = (Probe3D)probe;
                        float currentLMz10 = Float.MAX_VALUE; 
                        float currentLMz50 = Float.MAX_VALUE; 
                        float currentUMz10 = 0f; 
                        float currentUMz50 = 0f; 
                        if (probe3D.getLowMz10()>0){
                          currentLMz10 = probe3D.getLowMz10();
                          currentLMz50 = probe3D.getLowMz50();
                          currentUMz10 = probe3D.getUpMz10();
                          currentUMz50 = probe3D.getUpMz50();
                        }
                          
                        if (probe3D.getLowMz10()>0){
                          if (currentLMz10<lowerMz10Pc) lowerMz10Pc = currentLMz10;
                          if (currentLMz50<lowerMz50Pc) lowerMz50Pc = currentLMz50;
                          if (currentUMz10>upperMz10Pc) upperMz10Pc = currentUMz10;
                          if (currentUMz50>upperMz50Pc) upperMz50Pc = currentUMz50;
                        }
                          
                      }
                    }
                  }
                  if (j==0){
                    measuredMz = measuredMz/((float)set.getIsotopicProbes().get(j).size());
                    if (apexIntensity>0)  apexString = String.valueOf(Math.round(apexIntensity));
                    lowerValleyString = Calculator.FormatNumberToString(Calculator.roundFloat(lowerValley/60f,2),2d);
                    upperValleyString = Calculator.FormatNumberToString(Calculator.roundFloat(upperValley/60f,2),2d);
                    lowerMzString = String.valueOf(lowerMz);
                    upperMzString = String.valueOf(upperMz);

                    if (lowerValley10Pc<Float.MAX_VALUE){
                      lowerValley10PcString = Calculator.FormatNumberToString(Calculator.roundFloat(lowerValley10Pc/60f,2),2d);
                      upperValley10PcString = Calculator.FormatNumberToString(Calculator.roundFloat(upperValley10Pc/60f,2),2d);
                      lowerValley50PcString = Calculator.FormatNumberToString(Calculator.roundFloat(lowerValley50Pc/60f,2),2d);
                      upperValley50PcString = Calculator.FormatNumberToString(Calculator.roundFloat(upperValley50Pc/60f,2),2d);
                    }
                    if (lowerMz10Pc<Float.MAX_VALUE){
                      lowerMz10PcString = String.valueOf(lowerMz10Pc);
                      upperMz10PcString = String.valueOf(upperMz10Pc);
                      lowerMz50PcString = String.valueOf(lowerMz50Pc);
                      upperMz50PcString = String.valueOf(upperMz50Pc);                       
                    }
                  }
                }
                StringBuilder isotopeClusters = new StringBuilder();
                for (int j=0; j!=isotopicAreas.size(); j++){
                  if (j!=0) isotopeClusters.append("+");
                  isotopeClusters.append("M"+String.valueOf(j));
                }
                float targetMz = set.Mz[0]-mzCorrection;
                float correctedMeasured = measuredMz-mzCorrection;
                float mzOffset = measuredMz-set.Mz[0];
                float ppmError = (mzOffset*1000000f)/targetMz;
                  
                if (set instanceof LipidomicsMSnSet){
                  msLevel = "ms/ms"; 
                  LipidomicsMSnSet msnSet = (LipidomicsMSnSet)set;
                    
                  Hashtable<Integer,LinkedHashMap<Integer,Float>> msnRets = msnSet.getMsnRetentionTimes();
                  Hashtable<Integer,Integer> nrMSnSpectra = new Hashtable<Integer,Integer>();
                  Hashtable<Integer,String> msnRetStrings = new Hashtable<Integer,String>();
                  for (int j=2; j<=highestMSLevel; j++){
                    if (msnRets.containsKey(j)){
                      Vector<Float> rts =  new Vector<Float>(msnRets.get(j).values());
                      String rtsString = "";
                      for (Float rt : rts){
                        if (rtsString.length()>0) rtsString += "|";
                        rtsString += Calculator.FormatNumberToString(Calculator.roundFloat(rt/60f,2),2d);
                      }
                      nrMSnSpectra.put(j, rts.size());
                      msnRetStrings.put(j, rtsString);
                    } else {
                      nrMSnSpectra.put(j, 0);
                      msnRetStrings.put(j, "");
                    }
                  }

                  int rank = 0;
                  Hashtable<Integer,Hashtable<String,Hashtable<String,TargetlistEntry>>> alexFragments = null;
                  Hashtable<String,String> molLookup = null;
                  if (quantObject!=null && containsMSnInformation && quantObject.getMsnFragments()!=null && quantObject.getMsnFragments().size()>0) {
                    alexFragments = quantObject.getMsnFragments();
                    molLookup = quantObject.getMolSpeciesLookup();
                  }

                  for (Object msnNames : msnSet.getMSnIdentificationNames()){
                    String nameString = "";
                    double relativeShare = 1d;
                    String oneCombi = "";
                    rank++;
                    if (msnNames instanceof Vector){
                      Vector<String> names = (Vector<String>)msnNames;
                      for (int j=0; j!= names.size();j++){
                        if (j==0 && msnSet.getMSnIdentificationNames().size()>1) relativeShare = msnSet.getRelativeIntensity(names.get(j));
                        nameString += names.get(j)+"|";
                        oneCombi = names.get(j);
                      }
                      nameString = nameString.substring(0,nameString.length()-1);
                    }else{
                      nameString = (String)msnNames;
                      oneCombi = nameString;
                      if (msnSet.getMSnIdentificationNames().size()>1) relativeShare=msnSet.getRelativeIntensity(nameString);
                    }
                    oneCombi = msnSet.getCombiIdFromHumanReadable(oneCombi);
                      
                    List<TargetlistFloatStringVO> fragments = new ArrayList<TargetlistFloatStringVO>();
                    if (msnSet.getStatus()>=LipidomicsMSnSet.HEAD_GROUP_DETECTED){
                      Hashtable<String,CgProbe> headFragments = msnSet.getHeadGroupFragments();
                      for (String key:headFragments.keySet()){
                        TargetlistEntry ms2Target = getMSnTargetlistEntry(key,null,null,alexFragments);
                        fragments.add(new TargetlistFloatStringVO(key,headFragments.get(key).Mz,headFragments.get(key).getMsLevel(),ms2Target));
                      }
                    }
                    if (msnSet.getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED){
                      //molSpeciesId = classNameToWrite+"("+nameString+")";
                      Hashtable<String,Hashtable<String,CgProbe>> chainFrags =  msnSet.getChainFragments();

//                      Vector<String> fas = StaticUtils.splitChainCombiToEncodedStrings(oneCombi.replaceAll(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS, LipidomicsConstants.CHAIN_SEPARATOR_NO_POS),
//                          LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
                      Vector<FattyAcidVO> fas = StaticUtils.decodeLipidNamesFromChainCombi(oneCombi);
                      molSpeciesId = StaticUtils.encodeAlexMolSpeciesName(classNameToWrite,fas);
                      if (molLookup!=null) {
                        if (molLookup.containsKey(oneCombi))
                          molSpeciesId = molLookup.get(oneCombi);
                        else {
                          for (String combiName : StaticUtils.getPermutedChainNames(StaticUtils.splitChainCombiToEncodedStrings(oneCombi,LipidomicsConstants.CHAIN_COMBI_SEPARATOR),LipidomicsConstants.CHAIN_COMBI_SEPARATOR)) {
                            if (molLookup.containsKey(combiName)) {
                              molSpeciesId = molLookup.get(combiName);
                              break;
                            }
                          }
                        }
                      }
                      Hashtable<String,String> usedFAs = new Hashtable<String,String>();
                      for (int j=0; j!= fas.size(); j++){
                        FattyAcidVO fa = fas.get(j);
                        Hashtable<String,CgProbe> frags =  new Hashtable<String,CgProbe>();
                        //String faStored = StaticUtils.getStoredFAName(fa.getCarbonDbsId(),chainFrags);
                        //if (faStored!=null) frags =  chainFrags.get(faStored);
                        frags =  chainFrags.get(fa.getChainId());
                        if (frags==null)
                          continue;
                        for (String key : frags.keySet()){
                          String fragmentName = key;
                          if (!fragmentName.contains(fa.getChainId())) fragmentName = StaticUtils.getChainFragmentDisplayName(fragmentName,fa.getCarbonDbsId());
                          if (usedFAs.containsKey(fragmentName))continue;
                          usedFAs.put(fragmentName, fragmentName);
                          TargetlistEntry ms2Target = getMSnTargetlistEntry(key,fa.getChainId(),oneCombi,alexFragments);
                          fragments.add(new TargetlistFloatStringVO(fragmentName,frags.get(key).Mz,frags.get(key).getMsLevel(),ms2Target)); 
                        }
                      }
                    } else {
                      molSpeciesId = classNameToWrite+" "+nameString;
                      if (molLookup!=null)
                        molSpeciesId = speciesToWrite;
                    }
                    
                    Collections.sort(fragments,new GeneralComparator("at.tugraz.genome.lda.vos.FloatStringVO", "getValue", "java.lang.Float"));
                    Hashtable<Integer,String> fragmentNameHash = new Hashtable<Integer,String>();
                    Hashtable<Integer,String> fragmentMzHash = new Hashtable<Integer,String>();
                    for (int j=2; j<=highestMSLevel; j++){
                      fragmentNameHash.put(j, "");
                      fragmentMzHash.put(j, "");
                    }
                    String ms2Precursor = "";
                    String ms2Activation = "";
                    String ms3Precursor = "";
                    String ms3Activation = "";
                    for (TargetlistFloatStringVO fragment : fragments){
                      String fragmentNames = fragmentNameHash.get(fragment.getMsLevel());
                      String fragmentMzs = fragmentMzHash.get(fragment.getMsLevel());
                      String mz = Calculator.FormatNumberToString(fragment.getValue(),3);
                      if (fragmentNames.length()>0 ){
                        fragmentNames += "|";
                        fragmentMzs += "|";
                      }
                      fragmentNames += fragment.getKey()+"{"+mz+"}";
                      fragmentMzs+=mz;
                      fragmentNameHash.put(fragment.getMsLevel(), fragmentNames);
                      fragmentMzHash.put(fragment.getMsLevel(), fragmentMzs);
                      if (fragment.getMs2Target()!=null && fragment.getMs2Target().getMs2Precursor()!=null){
                        if (ms2Precursor!=null && ms2Precursor.length()>0){
                          if (!ms2Precursor.equalsIgnoreCase(fragment.getMs2Target().getMs2Precursor()))
                            System.out.println("ATTENTION: The lipid species "+speciesToWrite+" contains fragments that have a different MS-precursors in the ALEX target list");
                          if (!ms2Activation.equalsIgnoreCase(fragment.getMs2Target().getMs2Activation()))
                            System.out.println("ATTENTION: The lipid species "+speciesToWrite+" contains fragments that have a different MS2-activation in the ALEX target list");                                
                        }else{
                          ms2Precursor = fragment.getMs2Target().getMs2Precursor();
                          ms2Activation = fragment.getMs2Target().getMs2Activation();
                        }
                      }
                      if (fragment.getMs2Target()!=null && fragment.getMs2Target().getMs3Precursor()!=null){
                        if (ms3Precursor!=null && ms3Precursor.length()>0){
                          if (!ms3Precursor.contains(fragment.getMs2Target().getMs3Precursor())){
                            ms3Precursor += "|"+fragment.getMs2Target().getMs3Precursor();
                            ms3Activation += "|"+fragment.getMs2Target().getMs3Activation();
                          }
                        }else{
                          ms3Precursor = fragment.getMs2Target().getMs3Precursor();
                          ms3Activation = fragment.getMs2Target().getMs3Activation();
                        }
                      }
                      //TODO: this is a workaround for later RT-grouping
                      if (isCompModulePresent && quantObject!=null){
                        if (quantObject.getMs2Precursor()!=null)
                          ms2Precursor = quantObject.getMs2Precursor();
                        if (quantObject.getMs2Activation()!=null)
                          ms2Activation = quantObject.getMs2Activation();
                        if (quantObject.getMs3Precursor()!=null)
                          ms3Precursor = quantObject.getMs3Precursor();
                        if (quantObject.getMs3Activation()!=null)
                          ms3Activation = quantObject.getMs3Activation();
                      }
                    }
                      
                    //MSn info from quantObject
                      
                    //build the output String for one entry
                    StringBuilder line = new StringBuilder();
                    line.append(resultsExcelName).append("\t").append(msMachine);
                    if (detectorColumn) line.append("\t").append(detector);
                    if (polarityColumn) line.append("\t").append(polarity);
                    line.append("\t").append(msLevel).append("\t").append(speciesToWrite);
                    if (isRtGrouped)
                      line.append("\t").append(groupingId);
                    line.append("\t").append(molSpeciesId);
                    line.append(TAB).append(String.valueOf(rank));
                    line.append("\t").append(classNameToWrite);
                    for (int j=2; j<=highestMSLevel; j++){
                      line.append("\t"+fragmentNameHash.get(j)+"\t"+fragmentMzHash.get(j)+"\t"+nrMSnSpectra.get(j)+"\t"+msnRetStrings.get(j));
                      if (containsMSnInformation){
                        if (j==2)
                          line.append("\t").append(ms2Precursor).append("\t").append(ms2Activation);
                        else if (j==3)
                          line.append("\t").append(ms3Precursor).append("\t").append(ms3Activation);
                      }
                    }
                    line.append("\t").append(set.getRt()).append("\t").append(isotopeClusters).append("\t").append(lowerValleyString);
                    line.append("\t").append(upperValleyString).append("\t").append(lowerValley10PcString).append("\t").append(upperValley10PcString);
                    line.append("\t").append(lowerValley50PcString).append("\t").append(upperValley50PcString).append("\t").append(lowerMzString);
                    line.append("\t").append(upperMzString).append("\t").append(lowerMz10PcString).append("\t").append(upperMz10PcString);
                    line.append("\t").append(lowerMz50PcString).append("\t").append(upperMz50PcString).append("\t").append(String.valueOf(((double)totalArea)*relativeShare));
                    for (int j=0; j!=(highestIsotopeNumber+1); j++){
                      line.append("\t");
                      if (isotopicAreas.containsKey(j)) line.append(String.valueOf(isotopicAreas.get(j)));
                    }
                    line.append("\t").append(apexString).append("\t").append(targetMz).append("\t").append(measuredMz);
                    line.append("\t").append((mzCorrection*-1f)).append("\t").append(correctedMeasured).append("\t").append(mzOffset);
                    line.append("\t").append(ppmError).append("\t").append(set.getModificationName());
                    if (lipidIdColumn) line.append("\t").append(lipidId);
                    if (lipidCategoryColumn) line.append("\t").append(lipidCategory);
                    if (conflictsColumn) line.append("\t").append(conflicts);
                    line.append("\t").append(charge).append("\t").append(cIndexSpeciesToWrite).append("\t").append(dbIndexSpeciesToWrite);
                    if (ohIndexColumn){
                      line.append("\t");
                      if (ohIndex!=null) line.append(ohIndex.intValue());
                    }
                    line.append("\t").append(sumCompositionToWrite).append("\t").append(sumFormula).append("\n");
                    stream.write(line.toString().getBytes());
                  }

                }else{
                  //build the output String for one entry
                  StringBuilder line = new StringBuilder();
                  line.append(resultsExcelName).append("\t").append(msMachine);
                  if (detectorColumn) line.append("\t").append(detector);
                  if (polarityColumn) line.append("\t").append(polarity);
                  line.append("\t").append(msLevel).append("\t").append(speciesToWrite);
                  if (isRtGrouped)
                    line.append("\t").append(groupingId);
                  line.append("\t").append(molSpeciesId+TAB);
                  line.append("\t").append(classNameToWrite);
                  for (int j=2; j<=highestMSLevel; j++){
                    line.append("\t\t\t\t");
                    if (containsMSnInformation){
                      if (j==2 || j==3) line.append("\t\t");
                    }
                  }
                  line.append("\t").append(set.getRt()).append("\t").append(isotopeClusters).append("\t").append(lowerValleyString);
                  line.append("\t").append(upperValleyString).append("\t").append(lowerValley10PcString).append("\t").append(upperValley10PcString);
                  line.append("\t").append(lowerValley50PcString).append("\t").append(upperValley50PcString).append("\t").append(lowerMzString);
                  line.append("\t").append(upperMzString).append("\t").append(lowerMz10PcString).append("\t").append(upperMz10PcString);
                  line.append("\t").append(lowerMz50PcString).append("\t").append(upperMz50PcString).append("\t").append(String.valueOf(totalArea));
                  for (int j=0; j!=(highestIsotopeNumber+1); j++){
                    line.append("\t");
                    if (isotopicAreas.containsKey(j)) line.append(String.valueOf(isotopicAreas.get(j)));
                  }
                  line.append("\t").append(apexString).append("\t").append(targetMz).append("\t").append(measuredMz);
                  line.append("\t").append((mzCorrection*-1f)).append("\t").append(correctedMeasured).append("\t").append(mzOffset);
                  line.append("\t").append(ppmError).append("\t").append(set.getModificationName());
                  if (lipidIdColumn) line.append("\t").append(lipidId);
                  if (lipidCategoryColumn) line.append("\t").append(lipidCategory);
                  if (conflictsColumn) line.append("\t").append(conflicts);
                  line.append("\t").append(charge).append("\t").append(cIndexSpeciesToWrite).append("\t").append(dbIndexSpeciesToWrite);
                  if (ohIndexColumn){
                    line.append("\t");
                    if (ohIndex!=null) line.append(ohIndex.intValue());
                  }
                  line.append("\t").append(sumCompositionToWrite).append("\t").append(sumFormula).append("\n");
                  stream.write(line.toString().getBytes());
                }
              }
            }
          }
        }

      }
      stream.close();
    }
    catch (IOException e) {
      throw new RdbWriterException(e);
    }

  }
  
  /**
   * starts the writing procedure - this method should be used by export by the heat map
   * @param fileName name of the RDB file
   * @param analysisModule the comparative analysis module - containing information about the file comparisons by the heat map
   * @param classSequence the sequence in which the classes shall be written
   * @param analyteSequence the sequence the analytes are stored in the Alex123 target results
   * @param acceptedMolecules the molecules that are accepted by the GUI - to be used only by export via heat map
   * @param maxIsotopes the maximum isotopes settings - to be used only by export via heat map
   * @param quantObjects the quantification objects providing more information (only if there are results from Alex123 target lists)
   * @param allowedExps which experiments (absolute path) are allowed to be written
   * @param alexRtGrouper is this called from the AlexRtGrouper
   * @throws RdbWriterException if something is wrong with the writing
   * @throws ExcelInputFileException if something is wrong with the Excel quantitation files
   * @throws ChemicalFormulaException if something is wrong with the chemical formulae
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public void write(String fileName, ComparativeAnalysis analysisModule, LinkedHashMap<String,Integer> classSequence,
      Hashtable<String,Vector<String>> analyteSequence, Hashtable<String,Hashtable<String,String>> acceptedMolecules, Hashtable<String,Integer> maxIsotopes,
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects, Vector<String> allowedExps, boolean alexRtGrouper)
          throws RdbWriterException, ExcelInputFileException, ChemicalFormulaException, LipidCombinameEncodingException{
    Vector<QuantificationResult> results = new Vector<QuantificationResult>();
    Hashtable<Integer,String> expLookup = new Hashtable<Integer,String>();
    for (int i=0; i!=analysisModule.getExpNamesInSequence().size(); i++) {
      String exp = analysisModule.getExpNamesInSequence().get(i);
      File resultFile = analysisModule.getFullFilePath(exp);
      boolean ok = true;
      if (allowedExps!=null) {
        ok = false;
        for (String excelPath : allowedExps) {
          if (excelPath.equalsIgnoreCase(resultFile.getAbsolutePath())) {
            ok = true;
            break;
          }
        }
      }
      if (ok) {
        QuantificationResult result = LDAResultReader.readResultFile(resultFile.getAbsolutePath(),  new Hashtable<String,Boolean>());
        results.add(result);
        expLookup.put(result.hashCode(), exp);
      }
    }
    Hashtable<String,Vector<String>> molNamesInSequence = analysisModule.getAllMoleculeNames();
    if (analyteSequence!=null) {
      Hashtable<String,Vector<String>> changedOrder = new Hashtable<String,Vector<String>>();
      for (String lClass : molNamesInSequence.keySet()) {
        Vector<String> newOrder = new Vector<String>();
        Vector<String> oldOrder = molNamesInSequence.get(lClass);
        Vector<String> speciesNames = analyteSequence.get(lClass);
        for (String species : speciesNames) {
          for (String nameWRt : oldOrder) {
            if (nameWRt.startsWith(species+"_"))
              newOrder.add(nameWRt);
          }
        }
        changedOrder.put(lClass, newOrder);
      }
      molNamesInSequence = changedOrder;
    }
    this.write(fileName, results, expLookup, analysisModule.getClassSequence()!=null ? analysisModule.getClassSequence() : classSequence,
        molNamesInSequence, acceptedMolecules, analysisModule, maxIsotopes,quantObjects!=null ? quantObjects : analysisModule.getQuantObjects(),
        alexRtGrouper);
  }

  
  private TargetlistEntry getMSnTargetlistEntry(String fragmentName, String faName, String molSpecies, Hashtable<Integer,Hashtable<String,Hashtable<String,TargetlistEntry>>> alexFragments) throws LipidCombinameEncodingException{
    if (alexFragments==null) return null;
    String name = new String(fragmentName);
    if (molSpecies!=null && faName!=null && name.endsWith(" ("+faName+")"))
      name = name.substring(0,name.length()-(" ("+faName+")").length());
    TargetlistEntry entry = null;
    for (Hashtable<String,Hashtable<String,TargetlistEntry>> entriesOfLevel :  alexFragments.values()){
      if (!entriesOfLevel.containsKey(name)){
        String errorMessage = "The fragment \""+name+"\" ";
        if (faName!=null && molSpecies!=null){
          errorMessage += "does not exist in the Alex123 target list for the chain \""+faName+"\" for the molecular species \""+molSpecies+"\"!";
        }else{
          errorMessage += "does not exist in the Alex123 target list!";
        }
        System.out.println(errorMessage);
        return null;
      }
      Hashtable<String,TargetlistEntry> sameFragType = entriesOfLevel.get(name);
      if (molSpecies==null){
        if (sameFragType.size()!=1){
          System.out.println("A head group fragment in an Alex123 target list must not be assigned to more than one lipid molecular species! This is the case for the fragment \""+fragmentName+"\"!");
          return null;
        }
        entry = sameFragType.values().iterator().next();
      }else{
        String molSpeciesKey = molSpecies;
        if (!sameFragType.containsKey(molSpecies)){
          boolean foundPermutedVersion = false;
          for (String key : sameFragType.keySet()) {
            if (StaticUtils.isAPermutedVersion(key,molSpecies,LipidomicsConstants.CHAIN_COMBI_SEPARATOR)) {
              foundPermutedVersion = true;
              molSpeciesKey = key;
            }
          }
          if (!foundPermutedVersion) {
            System.out.println("The molecular species \""+molSpecies+"\" of the chain  fragment cannot be assinged (Alex123 target list)!");
            return null;
          }
        }
        entry = sameFragType.get(molSpeciesKey);
      }
    }
    return entry;
  }
}
