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

import at.tugraz.genome.lda.LDAResultReader;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.alex123.vos.TargetlistEntry;
import at.tugraz.genome.lda.alex123.vos.TargetlistFloatStringVO;
import at.tugraz.genome.lda.analysis.ComparativeAnalysis;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.RdbWriterException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.QuantVO;
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
  
  /** the prefix for the internal standard*/
  private String internalStandardPrefix_;
  /** the prefix for the external standard*/
  private String externalStandardPrefix_;
  
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
    internalStandardPrefix_ = internalStandardPrefix;
    externalStandardPrefix_ = externalStandardPrefix;
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
   */
  public void write(String fileName, Vector<QuantificationResult> results, LinkedHashMap<String,Integer> classSequence,
      Hashtable<String,Vector<String>> analyteSequence, Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects) throws RdbWriterException, ChemicalFormulaException {
    //the third variable is not required - only backward compatibility that is not necessary for this method
    //the sixth and the seventh parameter is not required by this method - all of the molecules and isotopes have to be written
    write(fileName,results,null,classSequence,analyteSequence,null,null,quantObjects);
  }
  
  /**
   * starts the writing procedure
   * @param fileNameOut name of the RDB file
   * @param results the results to be written - every entry is one LDA results file
   * @param fileLookup the name of the excel files - to be used only by export via heat map
   * @param classSequence the sequence in which the classes shall be written
   * @param correctAnalyteSequence the sequence in which the analytes of the classes shall be written
   * @param acceptedMolecules the molecules that are accepted by the GUI - to be used only by export via heat map
   * @param maxIsotopes the maximum isotopes settings - to be used only by export via heat map
   * @param quantObjects the quantification objects providing more information (only for Alex123 target lists)
   * @throws RdbWriterException if something is wrong with the writing
   * @throws ChemicalFormulaException if something is wrong with the chemical formulae
   */
  @SuppressWarnings("unchecked")
  private void write(String fileNameOut, Vector<QuantificationResult> results, Hashtable<Integer,File> fileLookup,
      LinkedHashMap<String,Integer> classSequence, Hashtable<String,Vector<String>> correctAnalyteSequence,
      Hashtable<String,Hashtable<String,String>> acceptedMolecules, Hashtable<String,Integer> maxIsotopes,
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects) throws RdbWriterException, ChemicalFormulaException {
    boolean detectorColumn = false;
    boolean polarityColumn = false;
    boolean lipidIdColumn = false;
    boolean lipidCategoryColumn = false;
    boolean conflictsColumn = false;
    boolean ohIndexColumn = false;
    boolean containsMSnInformation = false;
    int highestIsotopeNumber = 0;
    int highestMSLevel = 1;
    
    //checking if the Alex123 target list specific columns can be filled out
    if (quantObjects!=null){
      for (Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass : quantObjects.values()){
        for (Hashtable<String,QuantVO> quantsOfAnalyte : quantsOfClass.values()){
          for (QuantVO quant : quantsOfAnalyte.values()){
            if (!(quant instanceof TargetlistEntry)) continue;
            TargetlistEntry entry = (TargetlistEntry) quant;
            if (entry.getDetector()!=null) detectorColumn = true;
            if (entry.getPolarity()!=null) polarityColumn = true;
            if (entry.getId()!=null) lipidIdColumn = true;
            if (entry.getCategory()!=null) lipidCategoryColumn = true;
            if (entry.getConflicts()!=null) conflictsColumn = true;
            if (entry.getOhNumber()>=0) ohIndexColumn = true;
            if (entry.getMsnFragments()!=null && entry.getMsnFragments().size()>0) containsMSnInformation = true;
            for (Integer msLevel : entry.getMsnFragments().keySet()){
              if (msLevel>highestMSLevel) highestMSLevel = msLevel;
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
    // get the highest MS-level
    }else{
      for (QuantificationResult result : results){
        for (String className : result.getIdentifications().keySet()){
          for (LipidParameterSet set : result.getIdentifications().get(className)){
            if (!(set instanceof LipidomicsMSnSet) || quantObjects!=null) continue;
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
            if ((set.getIsotopicProbes().size()-1)>highestIsotopeNumber) highestIsotopeNumber = (set.getIsotopicProbes().size()-1);
          }
        }
      }
    }
    
    try {
      BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(fileNameOut));
      StringBuilder headerLine = new StringBuilder();
      headerLine.append("RAW_ID\tMachine");
      if (detectorColumn) headerLine.append("\tDetector");
      if (polarityColumn) headerLine.append("\tPolarity");
      headerLine.append("\tMS evidence level\tLipid species\tMolecular lipid species\tLipid class");
      for (int i=2; i<=highestMSLevel; i++){
        String levelName = "MS"+i;
         headerLine.append("\t"+levelName+" Fragment names\t"+levelName+" Fragment m/z\t#"+levelName+" scans\t"+levelName+" RTs");
         if (containsMSnInformation) headerLine.append("\t"+levelName+" precursor m/z\t"+levelName+" activation");
      }
      headerLine.append("\tRetention time\tIsotope clusters\tLower valley\tUpper valley\tLower valley 10%\tUpper valley 10%\tLower valley 50%"
          + "\tUpper valley 50%\tLow m/z\tUpMz\tLow m/z 10%\tUp m/z 10%\tLow m/z 50%\tUp m/z 50%\tPeak area (total)");
      for (int i=0; i!=(highestIsotopeNumber+1); i++) headerLine.append("\tPeak area M"+i);
      headerLine.append("\tMax intensity\tTarget m/z\tMeasured m/z\tApplied m/z correction\tCorrected measured m/z\tm/z offset\tppm error\tAdduct");
      if (lipidIdColumn)  headerLine.append("\tLipid ID");
      if (lipidCategoryColumn) headerLine.append("\tLipid category");
      if (conflictsColumn) headerLine.append("\tConflicts");
      headerLine.append("\tCharge\tC index of lipid species\tDB index of lipid species");
      if (ohIndexColumn) headerLine.append("\tOH index of lipid species");
      headerLine.append("\tSum composition of lipid species\tSum formula of lipid species\n");
      stream.write(headerLine.toString().getBytes());
      
      boolean useHeatmap = false; 
      
      for (int i=0; i!=results.size(); i++) {
        QuantificationResult result = results.get(i);
        String msMachine = "";
        float mzCorrection = 0f;
        
        String resultsExcelName = null;
        if (result.getConstants()!=null && result.getConstants().getRawFileName()!=null){
          resultsExcelName = result.getConstants().getRawFileName();
          if (resultsExcelName.endsWith(".chrom")) resultsExcelName = resultsExcelName.substring(0,resultsExcelName.length()-".chrom".length());
        //the file lookup is for backward compatibility, and can be called by the heat map only
        }else resultsExcelName = fileLookup.get(result.hashCode()).getName();
        if (result.getConstants()!=null){
          mzCorrection = (float)result.getConstants().getShift();
          msMachine = result.getConstants().getMSMachine();
        }
        //check if all classes are in the class sequence
        LinkedHashMap<String,Integer> classSequenceFull = new LinkedHashMap<String,Integer>(classSequence);
        for (String className : result.getIdentifications().keySet()){
          if (!classSequence.containsKey(className)){
            System.out.println("ATTENTION: the class "+className+" is not in your target list!");
            classSequenceFull.put(className, 1);
          }
        }
        for (String className : classSequenceFull.keySet()){
          if (!result.getIdentifications().containsKey(className)) continue;
          Vector<LipidParameterSet> classResultsFile = result.getIdentifications().get(className);
          Hashtable<String,Vector<LipidParameterSet>> foundSpecies = new Hashtable<String,Vector<LipidParameterSet>>();
          for (LipidParameterSet set : classResultsFile){
            Vector<LipidParameterSet> ofOneAnalyte = new Vector<LipidParameterSet>();
            if (foundSpecies.containsKey(set.getNameStringWithoutRt())) ofOneAnalyte = foundSpecies.get(set.getNameStringWithoutRt());
            ofOneAnalyte.add(set);
            foundSpecies.put(set.getNameStringWithoutRt(), ofOneAnalyte);
          }
          Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = null;
          if (quantObjects!=null && quantObjects.containsKey(className)) quantsOfClass = quantObjects.get(className);
          
          Hashtable<String,String> selectedMolHash = null;
          if (acceptedMolecules!=null){
            if (acceptedMolecules.containsKey(className)) selectedMolHash = acceptedMolecules.get(className);
            else selectedMolHash = new Hashtable<String,String>();
          }

          int maxIsotope = Integer.MAX_VALUE-1;
          if (maxIsotopes!=null) maxIsotope = maxIsotopes.get(className); 
          
          for (String molName: correctAnalyteSequence.get(className)){
            if ((selectedMolHash==null || selectedMolHash.containsKey(molName)) && foundSpecies.containsKey(molName)){
              Hashtable<String,QuantVO> quantsOfAnalyte = null;
              if (quantsOfClass!=null && quantsOfClass.containsKey(molName)) quantsOfAnalyte = quantsOfClass.get(molName);
              String speciesId = "";
              String classNameToWrite = "";
              int cIndexSpecies = 0;
              int dbIndexSpecies = 0;
              String sumComposition = "";
              if (internalStandardPrefix_!=null && internalStandardPrefix_.length()>0 && molName.startsWith(internalStandardPrefix_)){
                speciesId = "IS " + className+" "+ molName.substring(internalStandardPrefix_.length());
                classNameToWrite = "IS "+className;
              }else if (externalStandardPrefix_!=null && externalStandardPrefix_.length()>0 && molName.startsWith(externalStandardPrefix_)){
                speciesId = "ES " + className + " " + molName.substring(externalStandardPrefix_.length());
                classNameToWrite = "ES "+className;
              }else{
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
              }
              
              if  (useHeatmap){
              //TODO: has to be implemented
              } else {
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
                  }else{
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
                    
                    Hashtable<Integer,Vector<Float>> msnRets = msnSet.getMsnRetentionTimes();
                    Hashtable<Integer,Integer> nrMSnSpectra = new Hashtable<Integer,Integer>();
                    Hashtable<Integer,String> msnRetStrings = new Hashtable<Integer,String>();
                    for (int j=2; j<=highestMSLevel; j++){
                      if (msnRets.containsKey(j)){
                        Vector<Float> rts = msnRets.get(j);
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

                    for (Object msnNames : msnSet.getMSnIdentificationNames()){
                      String nameString = "";
                      double relativeShare = 1d;
                      String oneCombi = "";
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
                      
                      Hashtable<Integer,Hashtable<String,Hashtable<String,TargetlistEntry>>> alexFragments = null;
                      if (quantObject!=null && containsMSnInformation && quantObject.getMsnFragments()!=null && quantObject.getMsnFragments().size()>0)
                        alexFragments = quantObject.getMsnFragments();
                      
                      List<TargetlistFloatStringVO> fragments = new ArrayList<TargetlistFloatStringVO>();
                      if (msnSet.getStatus()>=LipidomicsMSnSet.HEAD_GROUP_DETECTED){
                        Hashtable<String,CgProbe> headFragments = msnSet.getHeadGroupFragments();
                        for (String key:headFragments.keySet()){
                          TargetlistEntry ms2Target = getMSnTargetlistEntry(key,null,null,alexFragments);
                          fragments.add(new TargetlistFloatStringVO(key,headFragments.get(key).Mz,headFragments.get(key).getMsLevel(),ms2Target));
                        }
                      }
                      if (msnSet.getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED){
                        molSpeciesId = classNameToWrite+"("+nameString+")";
                        Hashtable<String,Hashtable<String,CgProbe>> chainFrags =  msnSet.getChainFragments();
                        String[] fas = LipidomicsMSnSet.getFAsFromCombiName(oneCombi);
                        Hashtable<String,String> usedFAs = new Hashtable<String,String>();
                        for (int j=0; j!= fas.length; j++){
                          String fa = fas[j];
                          Hashtable<String,CgProbe> frags =  new Hashtable<String,CgProbe>();
                          String faStored = StaticUtils.getStoredFAName(fa,chainFrags);
                          if (faStored!=null) frags =  chainFrags.get(faStored);
                          for (String key : frags.keySet()){
                            String fragmentName = key;
                            if (!fragmentName.contains(faStored)) fragmentName += " ("+faStored+")";
                            if (usedFAs.containsKey(fragmentName))continue;
                            usedFAs.put(fragmentName, fragmentName);
                            TargetlistEntry ms2Target = getMSnTargetlistEntry(key,faStored,oneCombi,alexFragments);
                            fragments.add(new TargetlistFloatStringVO(fragmentName,frags.get(key).Mz,frags.get(key).getMsLevel(),ms2Target)); 
                          }
                        }
                      } else molSpeciesId = classNameToWrite+" "+nameString;
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
                        fragmentNames += fragment.getKey();
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
                      }
                      
                      //MSn info from quantObject
                      
                      //build the output String for one entry
                      StringBuilder line = new StringBuilder();
                      line.append(resultsExcelName).append("\t").append(msMachine);
                      if (detectorColumn) line.append("\t").append(detector);
                      if (polarityColumn) line.append("\t").append(polarity);
                      line.append("\t").append(msLevel).append("\t").append(speciesToWrite).append("\t").append(molSpeciesId);
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
                    line.append("\t").append(msLevel).append("\t").append(speciesToWrite).append("\t").append(molSpeciesId);
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
   * @param acceptedMolecules the molecules that are accepted by the GUI - to be used only by export via heat map
   * @param maxIsotopes the maximum isotopes settings - to be used only by export via heat map
   * @throws RdbWriterException if something is wrong with the writing
   * @throws ExcelInputFileException if something is wrong with the Excel quantitation files
   * @throws ChemicalFormulaException if something is wrong with the chemical formulae
   */
  public void write(String fileName, ComparativeAnalysis analysisModule, LinkedHashMap<String,Integer> classSequence,
      Hashtable<String,Hashtable<String,String>> acceptedMolecules, Hashtable<String,Integer> maxIsotopes)
          throws RdbWriterException, ExcelInputFileException, ChemicalFormulaException{
    Vector<QuantificationResult> results = new Vector<QuantificationResult>();
    Hashtable<Integer,File> fileLookup = new Hashtable<Integer,File>();
    for (int i=0; i!=analysisModule.getExpNamesInSequence().size(); i++) {
      String exp = analysisModule.getExpNamesInSequence().get(i);
      File resultFile = analysisModule.getFullFilePath(exp);
      QuantificationResult result = LDAResultReader.readResultFile(resultFile.getAbsolutePath(),  new Hashtable<String,Boolean>());
      results.add(result);
      fileLookup.put(result.hashCode(), resultFile);
    }
    this.write(fileName, results, fileLookup, analysisModule.getClassSequence()!=null ? analysisModule.getClassSequence() : classSequence,
        analysisModule.getAllMoleculeNames(), acceptedMolecules, maxIsotopes,analysisModule.getQuantObjects());

  }

  
  private TargetlistEntry getMSnTargetlistEntry(String fragmentName, String faName, String molSpecies, Hashtable<Integer,Hashtable<String,Hashtable<String,TargetlistEntry>>> alexFragments){
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
        if (!sameFragType.containsKey(molSpecies)){
          System.out.println("The molecular species \""+molSpecies+"\" of the chain  fragment cannot be assinged (Alex123 target list)!");
          return null;
        }
        entry = sameFragType.get(molSpecies);
      }
    }
    return entry;
  }
}
