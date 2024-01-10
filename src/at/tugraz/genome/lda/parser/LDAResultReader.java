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
package at.tugraz.genome.lda.parser;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Properties;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.swing.JFrame;

import org.dhatim.fastexcel.reader.*;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.analysis.ComparativeNameExtractor;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.exception.SettingsException;
import at.tugraz.genome.lda.export.QuantificationResultExporter;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.msn.vos.IntensityChainVO;
import at.tugraz.genome.lda.msn.vos.IntensityPositionVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.Probe3D;

/**
 * 
 * @author Leonida M. Lamp
 * @author Juergen Hartler
 *
 */
public class LDAResultReader
{  
  
  private static LipidomicsConstants lipidomicsConstants_;
  private static HydroxyEncoding faHydroxyEncoding_;
  private static HydroxyEncoding lcbHydroxyEncoding_;
  private static Hashtable<String,Vector<LipidParameterSet>> resultParameterSets_;
  private static Hashtable<String,Integer> msLevels_;

  /**
   * reads an LDA results file in Excel format
   * @param filePath the absolute path to the Excel file
   * @param showModifications this hash is filled by the method and gives information whether there are more than one modifications present; key: lipid class
   * @return the contents of the Excel file stored in the corresponding value object
   * @throws ExcelInputFileException when there is something wrong with the Excel file
   */
  public static QuantificationResult readResultFile(String filePath, Hashtable<String,Boolean> showModifications)
      throws ExcelInputFileException{
    return LDAResultReader.readResultFile(filePath, showModifications, null);
  }

  /**
   * reads an LDA results file in Excel format
   * @param filePath the absolute path to the Excel file
   * @param showModifications this hash is filled by the method and gives information whether there are more than one modifications present; key: lipid class
   * @param specificClass filter for parsing only the results of one analyte class; enter null when no filter is required
   * @return the contents of the Excel file stored in the corresponding value object
   * @throws ExcelInputFileException when there is something wrong with the Excel file
   */
	public static QuantificationResult readResultFile(String filePath, Hashtable<String,Boolean> showModifications, 
      String specificClass) throws ExcelInputFileException{
    lipidomicsConstants_ = null;
    faHydroxyEncoding_ = null;
    lcbHydroxyEncoding_ = null;
    resultParameterSets_ = new Hashtable<String,Vector<LipidParameterSet>>();
    msLevels_ = new Hashtable<String,Integer>();
    String suffix = "";
    if (filePath!=null && filePath.length()>3)
      suffix = filePath.substring(filePath.lastIndexOf("."));
    if (!(suffix.equalsIgnoreCase(".xlsx"))){
      new WarningMessage(new JFrame(), "ERROR", "The specified file format is not supported!");
      throw new ExcelInputFileException("The specified file format is not supported!");
    } 
    try (InputStream is = new FileInputStream(filePath);
        ReadableWorkbook wb = new ReadableWorkbook(is);
        Stream<Sheet> sheets = wb.getSheets();) {
      //the comparator makes sure MS1 sheets and the lipidomicsConstants are read first (requirement for MSn and double bond position information)
      Comparator<Sheet> sheetComparator = (s1, s2) -> compareBySheetName(s1, s2);
      sheets.filter((s) -> parseSheet(s, specificClass))
            .sorted(sheetComparator)
            .forEach((s) -> {
              try {
                readSheet(s, showModifications);
              } catch (SettingsException | RulesException | LipidCombinameEncodingException | IOException ex) {
                new WarningMessage(new JFrame(), "ERROR", ex.getMessage());
              }
            });
      
    } catch (IOException ex){
      ex.printStackTrace();
      new WarningMessage(new JFrame(), "ERROR", ex.getMessage());
      throw new ExcelInputFileException(ex);
    }
    
    return new QuantificationResult(resultParameterSets_,lipidomicsConstants_,msLevels_,faHydroxyEncoding_,lcbHydroxyEncoding_);
    
  }
  
  /**
   * Calls the corresponding method for the individual Excel sheets
   * @param sheet the Excel sheet
   * @param showModifications this hash is filled by the method and gives information whether there are more than one modifications present; key: lipid class
   * @throws SettingsException
   * @throws RulesException
   * @throws LipidCombinameEncodingException
   */
  private static void readSheet(Sheet sheet, Hashtable<String,Boolean> showModifications) 
      throws SettingsException, RulesException, LipidCombinameEncodingException, IOException {
    String name = sheet.getName();
    if (name.equals(QuantificationResultExporter.SHEET_CONSTANTS)){
      Object[] settings = readSettingsFromExcel(sheet);
      lipidomicsConstants_ = (LipidomicsConstants)settings[0];
      faHydroxyEncoding_ = (HydroxyEncoding)settings[1];
      lcbHydroxyEncoding_ = (HydroxyEncoding)settings[2];
      
    } else if (name.endsWith(QuantificationResultExporter.ADDUCT_MSN_SHEET)) {
      readMSnSheet(sheet);
      
    } else if (name.endsWith(QuantificationResultExporter.ADDUCT_OMEGA_SHEET)) {
      readOmegaSheet(sheet);
      LipidParameterSet.setOmegaInformationAvailable(true);
      
    } else {
      readMS1Sheet(sheet, showModifications);
      
    }
    
  }

  /**
   * Reads the Settings from an Excel sheet
   * @param sheet Constants Excel sheet
   * @throws SettingsException thrown when a settings combination is not possible
   * @return settings: [0] LipidomicsConstants object containing the parameters that were read; [1] FA hydroxylation encoding; [2] LCB hydroxylation encoding
   */
  public static Object[] readSettingsFromExcel(Sheet sheet) throws SettingsException {
    Properties properties = new Properties();
    List<Row> rows = null;
    try {
      rows = sheet.read();
    } catch (IOException ex) { }
    Row headerRow = rows.get(QuantificationResultExporter.HEADER_ROW);
    List<Row> contentRows = rows.subList(QuantificationResultExporter.HEADER_ROW+1, rows.size());
    List<String> headerTitles = LDAResultReader.readSheetHeaderTitles(headerRow);
    Hashtable<String,Short> faOhEncondings = new Hashtable<String,Short>();
    Hashtable<String,Short> lcbOhEncondings = new Hashtable<String,Short>();
    String ohNumberString;

    for (Row row : contentRows) {
      
      String key = row.getCellText(headerTitles.indexOf(LipidomicsConstants.EXCEL_KEY));
      String value = row.getCellText(headerTitles.indexOf(LipidomicsConstants.EXCEL_VALUE));
      
      if (key.startsWith(LipidomicsConstants.EXCEL_HYDROXY_FA_PREFIX) || key.startsWith(LipidomicsConstants.EXCEL_HYDROXY_LCB_PREFIX)) {
        ohNumberString = null;
        if (key.startsWith(LipidomicsConstants.EXCEL_HYDROXY_FA_PREFIX))
          ohNumberString = key.substring(LipidomicsConstants.EXCEL_HYDROXY_FA_PREFIX.length());
        else if (key.startsWith(LipidomicsConstants.EXCEL_HYDROXY_LCB_PREFIX))
          ohNumberString = key.substring(LipidomicsConstants.EXCEL_HYDROXY_LCB_PREFIX.length());
        Short ohNumber = Short.parseShort(ohNumberString);
        if (key.startsWith(LipidomicsConstants.EXCEL_HYDROXY_FA_PREFIX))
          faOhEncondings.put(value, ohNumber);
        else if (key.startsWith(LipidomicsConstants.EXCEL_HYDROXY_LCB_PREFIX))
          lcbOhEncondings.put(value, ohNumber);
      } else {
        properties.put(key, value);
      }
    }
    
    LipidomicsConstants constants = new LipidomicsConstants(false);
    constants.setLDAVersion(properties.getProperty(LipidomicsConstants.LDA_VERSION, Settings.VERSION));
    if (properties.containsKey(LipidomicsConstants.RAW_FILE)) {
      constants.setRawFileName(properties.getProperty(LipidomicsConstants.RAW_FILE));
    }
    if (properties.containsKey(LipidomicsConstants.BASE_PEAK_CUTOFF)) {
      String cutoff = properties.getProperty(LipidomicsConstants.BASE_PEAK_CUTOFF);
      if (cutoff != null) {
        try{
          Double.parseDouble(cutoff);
          constants.setRelativeMS1BasePeakCutoff(cutoff);
        }catch (NumberFormatException nfx){}
      }
    }
    if (properties.containsKey(LipidomicsConstants.MASS_SHIFT)) {
      //TODO: shouldn't we make the 1000d a constant somewhere!?
      properties.put(LipidomicsConstants.MASS_SHIFT, 
          String.valueOf(Double.parseDouble(properties.getProperty(LipidomicsConstants.MASS_SHIFT))*1000d));
    }
    constants.setVariables(properties);
    Object[] returnValues = new Object[3];
    returnValues[0] = constants;
    returnValues[1] = new HydroxyEncoding(faOhEncondings);
    returnValues[2] = new HydroxyEncoding(lcbOhEncondings);
    return returnValues;
  }
  
  
  /**
   * Reads MSn evidence from an Excel sheet. Where applicable, the results are stored in a LipidomicsMSnSet.
   * @param sheet MSn Excel sheet
   * @throws RulesException 
   * @throws LipidCombinameEncodingException thrown when a lipid combi ID (containing type and OH number) cannot be decoded
   */
  private static void readMSnSheet(Sheet sheet) throws RulesException, LipidCombinameEncodingException {
    Hashtable<String,LipidParameterSet> msHash = new Hashtable<String,LipidParameterSet>();
    String lipidClass = sheet.getName().substring(0,sheet.getName().lastIndexOf(QuantificationResultExporter.ADDUCT_MSN_SHEET));
    Vector<LipidParameterSet> resultPrms = resultParameterSets_.get(lipidClass);
    resultPrms.stream().forEach((p) -> msHash.put(p.getNamePlusModHumanReadable(), p));
    
    LipidParameterSet addingMSnEvidence = null;
    Hashtable<Integer,String> columnToIdentification = new Hashtable<Integer,String>();
    Hashtable<String,Double> relativeAreas = new Hashtable<String,Double>();
    String regex = "MS(\\d+) scan RTs";
    Pattern msLevelPattern =  Pattern.compile(regex);
    Hashtable<Integer,LinkedHashMap<Integer,Float>> msnRetentionTimes = new Hashtable<Integer,LinkedHashMap<Integer,Float>>();
    boolean checkMSnAreas = false;
    boolean headGroupFragmentActive = false;
    boolean headGroupRules = false;
    boolean chainFragmentActive = false;
    boolean chainRules = false;
    boolean positionRules = false;
    String combiKey = "";

    // is the header row read and the column indices initialized 
    boolean headerRowRead = false;
    // the header row has to be read
    boolean readFragmentHeaderRow = false;

    // is the intensity header read and the column indices initialized 
    boolean intensityHeaderRead = false;
    // the intensity header row has to be read
    boolean readIntensityHeaderRow = false;
    
    // column indices for the fragment information rows
    int nameColumn = -1;
    int ohColumn = -1;
    int chainTypeColumn = -1;
    int formulaColumn = -1;
    int msLevelColumn = -1;
    int chargeColumn = -1;
    int mzColumn = -1;
    int mzTolColumn = -1;
    int areaColumn = -1;
    int peakColumn = -1;
    int startTimeColumn = -1;
    int stopTimeColumn = -1;
    int startMzColumn = -1;
    int stopMzColumn = -1;
    int ellTimeColumn = -1;
    int ellMzColumn = -1;
    int ellTimeRangeColumn = -1;
    int ellMzRangeColumn = -1;
    
    // column indices for the intensity information
    int ruleColumn = -1;
    int originalRuleColumn = -1;
    int ruleValuesColumn = -1;
    int ruleMissedColumn = -1;

        
    // the information to be written into the returning VO
    int status = LipidomicsMSnSet.NO_MSN_PRESENT;
    float mzTolerance = -1f;
    Hashtable<String,CgProbe> headGroupFragments = new Hashtable<String,CgProbe>();
    Hashtable<String,IntensityRuleVO> headIntensityRules = new Hashtable<String,IntensityRuleVO>();
    Hashtable<String,Hashtable<String,CgProbe>> chainFragments = new Hashtable<String,Hashtable<String,CgProbe>>();
    Hashtable<String,Hashtable<String,IntensityChainVO>> chainIntensityRules = new Hashtable<String,Hashtable<String,IntensityChainVO>>();
    Vector<String> validChainCombinations = new Vector<String>();
    Hashtable<String,Hashtable<Integer,Integer>> positionDefinition = new Hashtable<String,Hashtable<Integer,Integer>>();
    Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> positionEvidence = new Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>();
    Hashtable<Integer,Float> basePeakValues = new Hashtable<Integer,Float>();
    
    Hashtable<String,String> uniqueRules = new Hashtable<String,String>();
    int numberOfPositions = -1;
    boolean containsAlkyl = false;
    boolean containsAlkenyl = false;
    boolean usedAlexMsnTargets = false;
    
    char knownPosSep = LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS.toCharArray()[0];
    char unknownPosSep = LipidomicsConstants.CHAIN_SEPARATOR_NO_POS.toCharArray()[0];
    
    
    List<Row> rows = null;
    try {
      rows = sheet.read();
    } catch (IOException ex) {}
    
    
    int lastRow = rows.get(rows.size()-1).getRowNum();
    //To mimic empty rows as these are delimiters of entries (empty rows are not read in)
    int counter = 0;
    for (int row = 1; row<=lastRow; row++) {
      Hashtable<Integer,String> cellEntries = new Hashtable<Integer,String>();
      
      if (rows.get(counter).getRowNum() == row) {
        rows.get(counter).stream().filter((c) -> !(c==null || c.getType().equals(CellType.ERROR))).forEach((c) -> cellEntries.put(c.getColumnIndex(), c.getRawValue()));
        counter++;
      }
      
      // check if an Alex123 target list was used for the MSn fragments of this class
      if (addingMSnEvidence==null && 
          cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) && 
          ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().startsWith(QuantificationResultExporter.HEADER_ALEX123_MSN_TARGETS_USED)){
        StringTokenizer tokenizer = new StringTokenizer(((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)),"=");
        if (tokenizer.countTokens()!=2) continue;
        tokenizer.nextToken();
        String value = tokenizer.nextToken().trim();
        if (value.equalsIgnoreCase("true") || value.equalsIgnoreCase("yes"))
          usedAlexMsnTargets = true;
      } 
      
      // reading the first row, containing the sum formula, and the individual identifications
      else if (addingMSnEvidence==null && 
          cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) && 
          cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_FORMULA)){
        // for every new lipid MS1 species, the parameters holding the information have to be initialized
        status = LipidomicsMSnSet.NO_MSN_PRESENT;
        mzTolerance = -1f;
        headGroupFragments = new Hashtable<String,CgProbe>();
        headIntensityRules = new Hashtable<String,IntensityRuleVO>();
        chainFragments = new Hashtable<String,Hashtable<String,CgProbe>>();
        chainIntensityRules = new Hashtable<String,Hashtable<String,IntensityChainVO>>();
        validChainCombinations = new Vector<String>();
        positionDefinition = new Hashtable<String,Hashtable<Integer,Integer>>();
        positionEvidence = new Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>>();
        basePeakValues = new Hashtable<Integer,Float>();
        
        columnToIdentification = new Hashtable<Integer,String>();
        String speciesName = (String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME);
        addingMSnEvidence = msHash.get(speciesName);
        int count = QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME+1;
        numberOfPositions = -1;
        containsAlkyl = false;
        containsAlkenyl = false;
        while (cellEntries.containsKey(count)){
          String lipidIdentification = (String)cellEntries.get(count);
          
          if (lipidIdentification!=null && lipidIdentification.length()>0){
            //this put has to be here to store the MSn scan numbers
            columnToIdentification.put(count, lipidIdentification);
            // this if was extended by "lipidIdentification.indexOf(":")!=-1" to support single chains - I hope this has no negative side effects
            if (lipidIdentification.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS)!=-1||lipidIdentification.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)!=-1||
                lipidIdentification.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_DBS)!=-1){
              boolean isAlexOhEncoding = false;
              //if there is a ";", the next character after the following numbers must be a ":"; if there is a "/" or a "_" it is the OH index of an Alex123 notation
              if (usedAlexMsnTargets){
                if (lipidIdentification.indexOf(LipidomicsConstants.CHAIN_COMBI_SEPARATOR_AMBIG_POS_OLD)!=-1){
                  boolean makeSubstring = false;
                  char[] chars = lipidIdentification.toCharArray();
                  for (int i=lipidIdentification.indexOf(lipidomicsConstants_.getChainCombiSeparatorAmbigPosDependingOnVersion())+1; i!=chars.length; i++){
                    //it is the OH index of an Alex123 notation
                    if (chars[i]==knownPosSep || chars[i]==unknownPosSep) {
                      isAlexOhEncoding = true;
                      break;
                    //it is an LDA identification where the position is unknown -> cut
                    }else if (chars[i]==':'){
                      makeSubstring = true;
                      break;
                    }else
                      continue;
                  }
                  if (makeSubstring)
                    lipidIdentification = lipidIdentification.substring(0,lipidIdentification.indexOf(lipidomicsConstants_.getChainCombiSeparatorAmbigPosDependingOnVersion()));
                }
              } else {
            	  if (lipidomicsConstants_.getChainCombiSeparatorAmbigPosDependingOnVersion() != null) {
            		  if (lipidIdentification.indexOf(lipidomicsConstants_.getChainCombiSeparatorAmbigPosDependingOnVersion())!=-1) {
            			  lipidIdentification = lipidIdentification.substring(0,lipidIdentification.indexOf(lipidomicsConstants_.getChainCombiSeparatorAmbigPosDependingOnVersion()));
            		  }
            	  }
              }
              
              // remove sn position annotations if present
              lipidIdentification = StaticUtils.removeSNPositions(lipidIdentification);
              
              // remove omega-DB position annotations if present
              lipidIdentification = StaticUtils.getHumanReadableWODoubleBondPositions(lipidIdentification);
              
              Object[] nameAndNumberOfPos = StaticUtils.cleanEmptyFAPositionsAndEncodeToLDACombiName(lipidIdentification,faHydroxyEncoding_,lcbHydroxyEncoding_,isAlexOhEncoding,lipidomicsConstants_);
              lipidIdentification = (String)nameAndNumberOfPos[0];
              int posNr = (Integer)nameAndNumberOfPos[1];
              if (posNr>numberOfPositions) numberOfPositions = posNr;
              containsAlkyl = (Boolean)nameAndNumberOfPos[2];
              containsAlkenyl = (Boolean)nameAndNumberOfPos[3];
              //when the correct MSn LDA decoded MSn identification name is known, put the correct one to the column
              columnToIdentification.put(count, lipidIdentification);
              validChainCombinations.add(lipidIdentification);
            }
          }
          count++;
        }
        checkMSnAreas = true;
        headGroupFragmentActive = false;
        headGroupRules = false;
        chainFragmentActive = false;
        chainRules = false;
        positionRules = false;

      }
      
      // reading the second row, containing the relative areas of the analytes and the retention times
      else if (checkMSnAreas && addingMSnEvidence!=null){
        relativeAreas = new Hashtable<String,Double>();
        msnRetentionTimes = new Hashtable<Integer,LinkedHashMap<Integer,Float>>();
        double totalArea = Double.parseDouble(cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME));
        for (Integer column : columnToIdentification.keySet()){
          String lipidIdentification = columnToIdentification.get(column);
          //this is for the retention times
          Matcher msLevelMatcher = msLevelPattern.matcher(lipidIdentification);
          if (msLevelMatcher.matches()){
            int msLevel =  Integer.parseInt(msLevelMatcher.group(1));
            LinkedHashMap<Integer,Float> rts = new LinkedHashMap<Integer,Float>();
            StringTokenizer rtTokenizer = null;
            rtTokenizer = new StringTokenizer((String)cellEntries.get(column),";");
            //this is for backward compatibility
            int count = 1;
            while (rtTokenizer.hasMoreTokens()){
              String kvPairString = rtTokenizer.nextToken();
              int scanNr = -1;
              float rt = -1;
              //this is for backward compatibility
              if (kvPairString.indexOf("=")==-1){
                scanNr = count;
                rt = new Float(kvPairString);
                count++;
              }else{
                String[] kvPair = kvPairString.split("=");
                scanNr = Integer.parseInt(kvPair[0]);
                rt = Float.parseFloat(kvPair[1]);
              }
              rts.put(scanNr, rt);
            }
            msnRetentionTimes.put(msLevel, rts);
          //this is for the relative areas
          }else
            relativeAreas.put(lipidIdentification, (Double.parseDouble(cellEntries.get(column)))/totalArea);
        }
        checkMSnAreas = false;
      }
      
      // the following if is for activating the head group fragments section
      else if (!headGroupFragmentActive && addingMSnEvidence!=null && 
          cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) && 
          ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_SECTION_HEAD_FRAGMENTS)){
        headGroupFragmentActive = true;
        headGroupRules = false;
        chainFragmentActive = false;
        chainRules = false;
        positionRules = false;
        readFragmentHeaderRow = true;
        headerRowRead = false;
        intensityHeaderRead = false;
      } 
      
      // the following if is for activating the head group rules section
      else if (!headGroupRules && addingMSnEvidence!=null && 
          cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) && 
          ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_SECTION_HEAD_INTENSITIES)){
        headGroupFragmentActive = false;
        headGroupRules = true;
        chainFragmentActive = false;
        chainRules = false;
        positionRules = false;
        headerRowRead = false;
        readIntensityHeaderRow = true;
        intensityHeaderRead = false;
        uniqueRules = new Hashtable<String,String>();
      }

      // the following if is for activating the chain fragments section
      else if (!chainFragmentActive && addingMSnEvidence!=null && 
          cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) && 
          ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_SECTION_CHAIN_FRAGMENTS)){
        headGroupFragmentActive = false;
        headGroupRules = false;
        chainFragmentActive = true;
        chainRules = false;
        positionRules = false;
        readFragmentHeaderRow = true;
        headerRowRead = false; 
        intensityHeaderRead = false;
      }
      
      // the following if is for activating the chain fragments section
      else if (!chainRules && addingMSnEvidence!=null && 
          cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) && 
          ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_SECTION_CHAIN_INTENSITIES)){
        headGroupFragmentActive = false;
        headGroupRules = false;
        chainFragmentActive = false;
        chainRules = true;
        positionRules = false;
        headerRowRead = false;
        readIntensityHeaderRow = true;
        intensityHeaderRead = false;  
        uniqueRules = new Hashtable<String,String>();
      }
      
      // the following ifs are for activating the position information section
      else if (addingMSnEvidence!=null && 
          cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME) && 
          ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().startsWith(LipidomicsConstants.EXCEL_MSN_SECTION_POSITION_INTENSITIES)){
        headGroupFragmentActive = false;
        headGroupRules = false;
        chainFragmentActive = false;
        chainRules = false;
        positionRules = true;
        headerRowRead = false;
        readIntensityHeaderRow = true;
        intensityHeaderRead = false;
        combiKey = ((String)cellEntries.get(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)).trim().substring(LipidomicsConstants.EXCEL_MSN_SECTION_POSITION_INTENSITIES.length());
        
        combiKey = combiKey.substring(combiKey.indexOf("(")+1,combiKey.indexOf(")"));      
        
        if (lipidomicsConstants_.getChainCombiSeparatorAmbigPosDependingOnVersion() != null) {
        	if (combiKey.indexOf(lipidomicsConstants_.getChainCombiSeparatorAmbigPosDependingOnVersion())!=-1) combiKey = combiKey.substring(0,combiKey.indexOf(lipidomicsConstants_.getChainCombiSeparatorAmbigPosDependingOnVersion()));
        } else {
        	combiKey = StaticUtils.removeSNPositions(combiKey);
        	//combiKey = StaticUtils.removeModification(combiKey);
        }
        
        validChainCombinations = correctUndefinedChainCombinations(combiKey,validChainCombinations,relativeAreas,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
        Hashtable<String,Integer> chainOccurenceInCombi = new Hashtable<String,Integer>();
        combiKey = (String)StaticUtils.cleanEmptyFAPositionsAndEncodeToLDACombiName(combiKey,faHydroxyEncoding_,
            lcbHydroxyEncoding_,  (usedAlexMsnTargets && combiKey.indexOf(lipidomicsConstants_.getChainCombiSeparatorAmbigPosDependingOnVersion())!=-1), lipidomicsConstants_)[0];
        
        combiKey = getPermutationThatIsInValidChainCombinations(combiKey, validChainCombinations);
        
        Vector<String> fas = StaticUtils.splitChainCombiToEncodedStrings(combiKey,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
        for (String fa: fas){
          int count = 0;
          if (chainOccurenceInCombi.containsKey(fa)) count = chainOccurenceInCombi.get(fa);
          count++;
          chainOccurenceInCombi.put(fa, count);
        }
        
        // there is only one type of fatty acid chain
        if (chainOccurenceInCombi.size()==1 && fas.size()==numberOfPositions){
          Hashtable<Integer,Integer> definitions = new Hashtable<Integer,Integer>();
          for (int i=0; i!=fas.size();i++){
            definitions.put(i, i);
          }
          positionDefinition.put(combiKey, definitions);
          positionEvidence.put(combiKey, new Hashtable<Integer,Vector<IntensityPositionVO>>());
          uniqueRules = new Hashtable<String,String>();
          continue;
        }
      }
      
      // the final procedure for creating a LipidomicsMSnSet after the information was read
      else if (addingMSnEvidence!=null && !cellEntries.containsKey(QuantificationResultExporter.MSN_ROW_FRAGMENT_NAME)){
        String speciesName = addingMSnEvidence.getNamePlusModHumanReadable();
//        for (positionDefinition)
        positionDefinition = cleanPositionDefinition(positionDefinition);
        if (isDefinitionPresent(positionDefinition)) status = LipidomicsMSnSet.POSITION_DETECTED;
        addingMSnEvidence = new LipidomicsMSnSet(addingMSnEvidence, status, mzTolerance,headGroupFragments, headIntensityRules,
        chainFragments, chainIntensityRules, validChainCombinations, relativeAreas,  positionDefinition,positionEvidence, numberOfPositions, basePeakValues,
        msnRetentionTimes, faHydroxyEncoding_, lcbHydroxyEncoding_);
        msHash.put(speciesName,addingMSnEvidence);
        addingMSnEvidence=null;
      }
      
      // reading the fragment header row for assigning the columns
      else if (readFragmentHeaderRow){
        for (Integer columnId : cellEntries.keySet()){
          if ((cellEntries.get(columnId)).length()==0) continue;
          String entry = (String)cellEntries.get(columnId);
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_NAME)) nameColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_OH)) ohColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_CHAIN_TYPE)) chainTypeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_FORMULA)) formulaColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_MSLEVEL)) msLevelColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_CHARGE)) chargeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ)) mzColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_TOLERANCE)) mzTolColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_AREA)) areaColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_PEAK)) peakColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_TIME_LOWER)) startTimeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_TIME_UPPER)) stopTimeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_LOWER)) startMzColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_MZ_UPPER)) stopMzColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_TIME)) ellTimeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_MZ)) ellMzColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_TIME_RANGE)) ellTimeRangeColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_FRAGMENT_ELLIPSE_MZ_RANGE)) ellMzRangeColumn = columnId;
        }
        readFragmentHeaderRow = false;
        headerRowRead = true;
      }
      
      // reading the fragment header row for assigning the columns
      else if (readIntensityHeaderRow){
        for (Integer columnId : cellEntries.keySet()){
          if (cellEntries.get(columnId)==null || ((String)cellEntries.get(columnId)).length()==0) continue;
          String entry = (String)cellEntries.get(columnId);
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_INTENSITY_RULE)) ruleColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_INTENSITY_ORIGINAL)) originalRuleColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_INTENSITY_VALUES)) ruleValuesColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_MSN_INTENSITY_MISSED)) ruleMissedColumn = columnId;
        }
        readIntensityHeaderRow = false;
        intensityHeaderRead = true;
      } 
      
      // reading information about the fragment
      else if (addingMSnEvidence!=null && headerRowRead && (headGroupFragmentActive||chainFragmentActive)){ 
        String fragmentName = null;
        if (nameColumn>-1 && cellEntries.containsKey(nameColumn)) fragmentName = (String)cellEntries.get(nameColumn);
        int oh = 0;
        if (ohColumn>-1 && cellEntries.containsKey(ohColumn)) oh = (int)Math.rint(Double.parseDouble(cellEntries.get(ohColumn)));
        // I have to put this to LipidomicsConstants.CHAIN_TYPE_NO_CHAIN for backward compatibility - by this I know that these results were
        // produced by an older LDA version
        short chainType = LipidomicsConstants.CHAIN_TYPE_NO_CHAIN;
        if (chainTypeColumn>-1 && cellEntries.containsKey(chainTypeColumn)) {
          String chainTypeString =  (String)cellEntries.get(chainTypeColumn);
          if (chainTypeString.equalsIgnoreCase(FragmentRuleVO.CHAIN_NAME.substring(1)))
            chainType = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
          else if (chainTypeString.equalsIgnoreCase(FragmentRuleVO.ALKYL_CHAIN_NAME.substring(1)))
            chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
          else if (chainTypeString.equalsIgnoreCase(FragmentRuleVO.ALKENYL_CHAIN_NAME.substring(1)))
            chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
          else if (chainTypeString.equalsIgnoreCase(FragmentRuleVO.LCB_NAME.substring(1)))
            chainType = LipidomicsConstants.CHAIN_TYPE_LCB;
          else
            throw new LipidCombinameEncodingException("The entry \""+chainTypeString+"\" is not allowed in the column \""+LipidomicsConstants.EXCEL_MSN_FRAGMENT_CHAIN_TYPE+"\"!");
        //this is required for backward compatibility
        } else if (chainFragmentActive) {
          chainType = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
          if (containsAlkyl && fragmentName.indexOf("("+LipidomicsConstants.ALKYL_PREFIX)!=-1)
            chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
          if (containsAlkenyl && fragmentName.indexOf("("+LipidomicsConstants.ALKENYL_PREFIX)!=-1)
            chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
        }
        
        String formula = null;
        if (formulaColumn>-1 && cellEntries.containsKey(formulaColumn)) formula = "+"+(cellEntries.get(formulaColumn)).trim().replaceAll(" ", " \\+");
        int msLevel = -1;
        if (msLevelColumn>-1 && cellEntries.containsKey(msLevelColumn)) msLevel = (int)Math.rint(Double.parseDouble(cellEntries.get(msLevelColumn)));
        int charge = -1;
        if (chargeColumn>-1 && cellEntries.containsKey(chargeColumn)) charge = (int)Math.rint(Double.parseDouble(cellEntries.get(chargeColumn)));
        float mz = -1f;
        if (mzColumn>-1 && cellEntries.containsKey(mzColumn)) mz = Float.parseFloat(cellEntries.get(mzColumn));
        if (mzTolerance<0 && mzTolColumn>-1 && cellEntries.containsKey(mzTolColumn)) mzTolerance = Float.parseFloat(cellEntries.get(mzTolColumn));
        float area = -1f;
        if (areaColumn>-1 && cellEntries.containsKey(areaColumn)) area = Float.parseFloat(cellEntries.get(areaColumn));
        float peak = -1f;
        if (peakColumn>-1 && cellEntries.containsKey(peakColumn)) peak = Float.parseFloat(cellEntries.get(peakColumn));
        float startTime = -1f;
        if (startTimeColumn>-1 && cellEntries.containsKey(startTimeColumn)) startTime = Float.parseFloat(cellEntries.get(startTimeColumn));
        float stopTime = -1f;
        if (stopTimeColumn>-1 && cellEntries.containsKey(stopTimeColumn)) stopTime = Float.parseFloat(cellEntries.get(stopTimeColumn));
        float startMz = -1f;
        if (startMzColumn>-1 && cellEntries.containsKey(startMzColumn)) startMz = Float.parseFloat(cellEntries.get(startMzColumn));
        float stopMz = -1f;
        if (stopMzColumn>-1 && cellEntries.containsKey(stopMzColumn)) stopMz = Float.parseFloat(cellEntries.get(stopMzColumn));
        float ellipseTimePosition = -1f;
        if (ellTimeColumn>-1 && cellEntries.containsKey(ellTimeColumn)) ellipseTimePosition = Float.parseFloat(cellEntries.get(ellTimeColumn));
        float ellipseMzPosition = -1f;
        if (ellMzColumn>-1 && cellEntries.containsKey(ellMzColumn)) ellipseMzPosition = Float.parseFloat(cellEntries.get(ellMzColumn));
        float ellipseTimeStretch = -1f;
        if (ellTimeRangeColumn>-1 && cellEntries.containsKey(ellTimeRangeColumn)) ellipseTimeStretch = Float.parseFloat(cellEntries.get(ellTimeRangeColumn));
        float ellipseMzStretch = -1f;
        if (ellMzRangeColumn>-1 && cellEntries.containsKey(ellMzRangeColumn)) ellipseMzStretch = Float.parseFloat(cellEntries.get(ellMzRangeColumn));
        
        //TODO: formula is only excluded for the damaged Alex123 target list
        if (fragmentName!=null && fragmentName.length()>0 && /*formula!=null && formula.length()>0 &&*/ msLevel>0 && charge>0 && mz>0f && mzTolerance>0f
            && area>0f && peak>0f && startTime>-1f && stopTime>0f){
          CgProbe probe = new CgProbe(-1, charge, msLevel, formula);
          probe.AreaStatus = CgAreaStatus.OK;
          probe.Mz = mz;
          probe.Area = area;
          probe.AreaError = 0f;
          probe.Background = 0f;
          probe.Peak = peak;
          probe.LowerValley = startTime;
          probe.UpperValley = stopTime;
          probe.LowerMzBand = mzTolerance;
          probe.UpperMzBand = mzTolerance;
          probe.isotopeNumber = 0;
          if (startMz>0f && stopMz>0f && ellipseTimePosition>0f && ellipseMzPosition>0f && ellipseTimeStretch>0f && ellipseMzStretch>0f){
            probe = new Probe3D(probe,ellipseTimePosition,ellipseMzPosition,
                ellipseTimeStretch,ellipseMzStretch,-1f,-1f,-1f,-1f,-1f,-1f);
            probe.LowerMzBand = startMz;
            probe.UpperMzBand = stopMz;
          }
          if (headGroupFragmentActive){
            headGroupFragments.put(fragmentName, probe);
            status = LipidomicsMSnSet.HEAD_GROUP_DETECTED;
          }else if (chainFragmentActive){
            status = LipidomicsMSnSet.FRAGMENTS_DETECTED;
            Object[] faAndFragment = StaticUtils.parseChainFaAndFragmentNameFromExcel(fragmentName,chainType,oh,lipidomicsConstants_.shouldOldEncodingBeUsed());
            FattyAcidVO faName = (FattyAcidVO)faAndFragment[0];
            fragmentName =  (String)faAndFragment[1];
            Hashtable<String,CgProbe> fragments = new Hashtable<String,CgProbe>();
            if (chainFragments.containsKey(faName.getChainId())) fragments = chainFragments.get(faName.getChainId());
            fragments.put(fragmentName, probe);
            chainFragments.put(faName.getChainId(),fragments);
          }
        }
      }
      
      // reading information about intensity rules
      else if (addingMSnEvidence!=null && intensityHeaderRead && (headGroupRules||chainRules||positionRules)){
        String readableRuleInterpretation = null;
        if (ruleColumn>-1 && cellEntries.containsKey(ruleColumn)) readableRuleInterpretation = (String)cellEntries.get(ruleColumn);
        String rule = null;
        if (originalRuleColumn>-1 && cellEntries.containsKey(originalRuleColumn)) rule = (String)cellEntries.get(originalRuleColumn);
        String ruleValueInterpretation = null;
        if (ruleValuesColumn>-1 && cellEntries.containsKey(ruleValuesColumn)) ruleValueInterpretation = (String)cellEntries.get(ruleValuesColumn);
        if (readableRuleInterpretation!=null && readableRuleInterpretation.length()>0 && rule!=null && rule.length()>0 &&
            ruleValueInterpretation!=null && ruleValueInterpretation.length()>0){
          String uniqueRuleString = readableRuleInterpretation+";"+rule+";"+ruleValueInterpretation;
          String missedString = "";
          if (ruleMissedColumn>-1 && cellEntries.containsKey(ruleMissedColumn)) missedString = (String)cellEntries.get(ruleMissedColumn);
          Hashtable<String,Short> missed = new Hashtable<String,Short>();
          Hashtable<String,Short> missedPosition = new Hashtable<String,Short>();
          StringTokenizer tokenizer = new StringTokenizer(missedString,";");
          while (tokenizer.hasMoreTokens()){
            String token = tokenizer.nextToken();
            //this is necessary for position rules
            if (token.indexOf("[")!=-1) {
              if (positionRules)
                missedPosition.put(token, LipidomicsConstants.CHAIN_TYPE_NO_CHAIN);
              token = token.substring(0, token.indexOf("["));
            }
            short type = LipidomicsConstants.CHAIN_TYPE_NO_CHAIN;
            if (chainRules||positionRules){
              boolean isAChain = false;
              //is it a chain fragment in Alex notation
              if (token.startsWith("FA ") || token.startsWith("-FA ") || token.startsWith("LCB ") || token.startsWith("-LCB ") || token.startsWith("O-") || token.startsWith("-O-")){
                isAChain = true;
              }else {
                String remnant = new String(readableRuleInterpretation);
                while (remnant.indexOf(token)!=-1) {
                  remnant = remnant.substring(remnant.indexOf(token)+token.length());
                  if (remnant.startsWith("(")) {
                    isAChain = true;
                    break;
                  }
                }
              }
              if (isAChain)
                type = LipidomicsConstants.CHAIN_TYPE_MISSED;
            }           
            missed.put(token, type);
          }
          if (uniqueRules.containsKey(uniqueRuleString)) continue;
          uniqueRules.put(uniqueRuleString, uniqueRuleString);
          int currentSection = -1;
          if (headGroupRules) currentSection = FragRuleParser.HEAD_SECTION;
          else if (chainRules) currentSection = FragRuleParser.CHAINS_SECTION;
          else if (positionRules) currentSection = FragRuleParser.POSITION_SECTION;
          Hashtable<String,Short> head = new Hashtable<String,Short>();
          for (String key : headGroupFragments.keySet()) head.put(key,LipidomicsConstants.CHAIN_TYPE_NO_CHAIN);
          Hashtable<String,Short> chain = new Hashtable<String,Short>();
          for (String chainId : chainFragments.keySet()){
            FattyAcidVO fa = StaticUtils.decodeLipidNameForCreatingCombis(chainId);
            Hashtable<String,CgProbe> fragments = chainFragments.get(chainId);
            for (String name : fragments.keySet()) chain.put(name, fa.getChainType());
          }
          IntensityRuleVO ruleVO = FragRuleParser.extractIntensityVOFromEquation(rule, -1, currentSection, head, chain, null, missed);
          if (chainRules || positionRules){
            IntensityRuleVO ruleInst = null;
            if (chainRules) ruleInst = IntensityChainVO.getFattyAcidsFromReadableRule(readableRuleInterpretation,ruleVO,chainFragments,missed,
                faHydroxyEncoding_, lcbHydroxyEncoding_, addingMSnEvidence.getOhNumber()>0,usedAlexMsnTargets);
            else if (positionRules) ruleInst = IntensityPositionVO.getFattyAcidsFromReadableRule(readableRuleInterpretation,ruleVO,chainFragments,missedPosition,
                faHydroxyEncoding_, lcbHydroxyEncoding_, addingMSnEvidence.getOhNumber()>0,usedAlexMsnTargets);
            if (ruleInst!=null) ruleVO = ruleInst;
          }
          if (ruleVO.containsBasePeak()){
            Hashtable<String,Float> values = LipidomicsMSnSet.getFragmentAreas(ruleVO,headGroupFragments,chainFragments);
            for (String name : values.keySet()){
              if (values.get(name)<0) values.put(name, 0f);
            }
            float basePeak = ruleVO.extractBasePeakValue(ruleValueInterpretation,values);
            int msLevel = ruleVO.getMSLevel(headGroupFragments, chainFragments);
            if (msLevel>1) basePeakValues.put(msLevel, basePeak);
          }
          if (ruleVO!=null && headGroupRules) {
            headIntensityRules.put(ruleVO.getRuleIdentifier(), ruleVO);
          }
          else if (ruleVO!=null && chainRules){
            //this line was replaced by the next one since this step was already performed before - I simply casted now the object
            //IntensityChainVO chainVO =  IntensityChainVO.getFattyAcidsFromReadableRule(readableRuleInterpretation,ruleVO,chainFragments,missed);
            IntensityChainVO chainVO = (IntensityChainVO)ruleVO;
            Hashtable<String,IntensityChainVO> rules = new Hashtable<String,IntensityChainVO>();
            Vector<String> ids = new Vector<String>();
            if (chainVO.getChainType()==IntensityChainVO.DIFF_CHAIN_TYPES){
              Vector<String>  chains = new Vector<String>();
              for (FattyAcidVO faVO : chainVO.getParticipatingChains())
                chains.add(faVO.getChainId());
              ids = StaticUtils.getAllAffectedChainCombinations(validChainCombinations,chains);
            } else ids.add(chainVO.getParticipatingChains().get(0).getChainId());
            for (String id : ids) {
              if (chainIntensityRules.containsKey(id)) rules = chainIntensityRules.get(id);
              rules.put(ruleVO.getReadableRuleInterpretation(faHydroxyEncoding_,lcbHydroxyEncoding_), chainVO);
              chainIntensityRules.put(id, rules);
            }
          } else if (ruleVO!=null && positionRules){
            //this line was replaced by the next one since this step was already performed before - I simply casted now the object
            //IntensityPositionVO posVO = IntensityPositionVO.getFattyAcidsFromReadableRule(readableRuleInterpretation,ruleVO,chainFragments,missed);
            IntensityPositionVO posVO = (IntensityPositionVO)ruleVO;
            Hashtable<Integer,Integer> defs = new Hashtable<Integer,Integer>();
            Hashtable<Integer,Vector<IntensityPositionVO>> evidence = new Hashtable<Integer,Vector<IntensityPositionVO>>();
            if (positionDefinition.containsKey(combiKey)) defs = positionDefinition.get(combiKey);
            if (positionEvidence.containsKey(combiKey)) evidence = positionEvidence.get(combiKey);
            Vector<String> fas = StaticUtils.splitChainCombiToEncodedStrings(combiKey,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
            Hashtable<String,String> usedFAs = new Hashtable<String,String>();
            for (int i=0;i!=fas.size();i++){
              String fa = fas.get(i);
              int position = posVO.getPositionByFA(fa);
              if ( position<0 || usedFAs.containsKey(fa)) continue;
              if (defs.containsKey(i)){
                // if the stored assignment is not the same -> no position can be assigned
                if (defs.get(i)!=(position-1)){
                  boolean otherSameFAAtPositionOrUndefined = false;
                  //check if there is the same FA more than once, and this position is assigned or undefined
                  for (int j=0;j!=fas.size();j++){
                    if (i==j || !fa.equalsIgnoreCase(fas.get(j)))continue;
                    if (!defs.containsKey(j)||defs.get(i)==(position-1)||defs.get(j)==(position-1)) otherSameFAAtPositionOrUndefined=true;
                  }
                  if (!otherSameFAAtPositionOrUndefined) defs.remove(i);
                }
              }else{
                //if several rules are fulfilled for the same FA, and this FA exists more than once, the
                //defs would store the same position twice and would not allow another rule with another
                //position for the correct assignment
                boolean sameFaAlreadySamePosition = false;
                for (int j=0; j!=i; j++){
                  if (!fa.equalsIgnoreCase(fas.get(j))) continue;
                  if (defs.get(j)==(position-1)) sameFaAlreadySamePosition = true;
                }
                if (!sameFaAlreadySamePosition){
                  defs.put(i, (position-1));
                  usedFAs.put(fa, fa);
                }
              }
              Vector<IntensityPositionVO> rules = new Vector<IntensityPositionVO>();
              if (evidence.containsKey(position)) rules = evidence.get(position);
              boolean ruleIsThere = false;
              for (IntensityPositionVO other : rules){
                if (other.getReadableRuleInterpretation(faHydroxyEncoding_, lcbHydroxyEncoding_).equalsIgnoreCase(posVO.getReadableRuleInterpretation(faHydroxyEncoding_, lcbHydroxyEncoding_))){
                  ruleIsThere = true;
                  break;
                }
              }
              if (!ruleIsThere)rules.add(posVO);
              
              evidence.put(position, rules);
            }
            positionDefinition.put(combiKey, defs);
            positionEvidence.put(combiKey, evidence);
          }
        }
      }

      // if in the next row the fragment head row has to be read - the column indices have to be initialized
      if (readFragmentHeaderRow){
        nameColumn = -1;
        ohColumn = -1;
        chainTypeColumn = -1;
        formulaColumn = -1;
        msLevelColumn = -1;
        chargeColumn = -1;
        mzColumn = -1;
        mzTolColumn = -1;
        areaColumn = -1;
        peakColumn = -1;
        startTimeColumn = -1;
        stopTimeColumn = -1;
        startMzColumn = -1;
        stopMzColumn = -1;
        ellTimeColumn = -1;
        ellMzColumn = -1;
        ellTimeRangeColumn = -1;
        ellMzRangeColumn = -1;
      }
      
      // if in the next row the intensity head row has to be read - the column indices have to be initialized
      if (readIntensityHeaderRow){
        ruleColumn = -1;
        originalRuleColumn = -1;
        ruleValuesColumn = -1;
      }
    }
    
    if (addingMSnEvidence!=null && (addingMSnEvidence instanceof LipidParameterSet)){
      String speciesName = addingMSnEvidence.getNamePlusModHumanReadable();
      positionDefinition = cleanPositionDefinition(positionDefinition);
      if (isDefinitionPresent(positionDefinition)) status = LipidomicsMSnSet.POSITION_DETECTED;
      addingMSnEvidence = new LipidomicsMSnSet(addingMSnEvidence, status, mzTolerance, headGroupFragments, headIntensityRules,
      chainFragments, chainIntensityRules, validChainCombinations, relativeAreas, positionDefinition,positionEvidence, numberOfPositions, basePeakValues,
      msnRetentionTimes, faHydroxyEncoding_, lcbHydroxyEncoding_);
      msHash.put(speciesName,addingMSnEvidence);
    }

    if (usedAlexMsnTargets) lipidomicsConstants_.getAlexTargetlistUsed().put(lipidClass,true);
    Vector<LipidParameterSet> msnResults = new Vector<LipidParameterSet>();
    for (LipidParameterSet ms1 : resultParameterSets_.get(lipidClass)){
      msnResults.add(msHash.get(ms1.getNamePlusModHumanReadable()));
    }
    resultParameterSets_.put(lipidClass, msnResults);
    
  }

  
  /**
   * Reads double bond position evidence from an Excel sheet. The results are stored in a Vector of DoubleBondPositionVOs for each LipidParameterSet.
   * @param sheet Omega Excel sheet
   */
  private static void readOmegaSheet(Sheet sheet) throws IOException, LipidCombinameEncodingException
  {
    String lipidClass = sheet.getName().replace(QuantificationResultExporter.ADDUCT_OMEGA_SHEET, "");
    Hashtable<String,LipidParameterSet> msHash = new Hashtable<String,LipidParameterSet>();
    for (LipidParameterSet param : resultParameterSets_.get(lipidClass)){
      msHash.put(param.getNamePlusModHumanReadable(), param);
    }
    
    List<Row> rows = null;
    rows = sheet.read();
    Row headerRow = rows.get(QuantificationResultExporter.HEADER_ROW);
    List<String> headerTitles = readSheetHeaderTitles(headerRow);
    List<Row> contentRows = rows.subList(QuantificationResultExporter.HEADER_ROW+1, rows.size());
    
    String identifier = null;
    String molecularSpecies = null;
    String doubleBondPosition = null;
    float expectedRetentionTime = 0f;
    int accuracy = 0;
    boolean isAssigned = false;
    int index;
    String rawValue;
    
    for (Row row : contentRows) 
    {
      List<Cell> cells = row.stream().filter((c) -> !(c==null || c.getType().equals(CellType.ERROR))).collect(Collectors.toList());
      for (Cell cell : cells) {
        index = cell.getColumnIndex();
        rawValue = cell.getRawValue();
        
        if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_IDENTIFIER)) {
          identifier = rawValue;
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_MOLECULAR_SPECIES)) {
          molecularSpecies = rawValue;
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_DOUBLE_BOND_POSITION_LEVEL)) {
          doubleBondPosition = rawValue;
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_EXPECTED_RT)) {
          expectedRetentionTime = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_ACCURACY)) {
          accuracy = Integer.parseInt(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_ASSIGNED)) {
          isAssigned = rawValue.equalsIgnoreCase("1");
        } 
      }
      Vector<FattyAcidVO> chainCombination = StaticUtils.decodeFAsFromHumanReadableName(
          doubleBondPosition, Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), false, lipidomicsConstants_);
      
      DoubleBondPositionVO doubleBondPositionVO = new DoubleBondPositionVO(
          chainCombination, expectedRetentionTime, accuracy, molecularSpecies, isAssigned);
      msHash.get(identifier).addOmegaInformation(doubleBondPositionVO);
    }
  }
  
  
  /**
   * Reads MS1 evidence from an Excel sheet. The results are stored in a Vector of LipidParameterSets.
   * @param sheet MS1 Excel sheet
   * @param showModifications this hash is filled by the method and gives information whether there are more than one modifications present; key: lipid class
   */
  private static void readMS1Sheet(Sheet sheet, Hashtable<String,Boolean> showModifications) throws IOException
  {
    int msLevel=1;
    Vector<LipidParameterSet> resultParams = new Vector<LipidParameterSet>();
    LipidParameterSet params = null;
    boolean showModification = false;
    Hashtable<String,String> analyteNames = new Hashtable<String,String>();
    List<Row> rows = null;
    rows = sheet.read();
    Row headerRow = rows.get(QuantificationResultExporter.HEADER_ROW);
    List<String> headerTitles = readSheetHeaderTitles(headerRow);
    for (String title : headerTitles) {
      if (title.startsWith(QuantificationResultExporter.HEADER_MS_LEVEL)){
        String levelString = title.substring(QuantificationResultExporter.HEADER_MS_LEVEL.length()).trim();
        msLevel = Integer.valueOf(levelString);
      }
    }
    List<Row> contentRows = rows.subList(QuantificationResultExporter.HEADER_ROW+1, rows.size());
    
    for (Row row : contentRows) {
      List<Cell> cells = row.stream()
      		.filter((c) -> !(c==null || c.getType().equals(CellType.ERROR) || c.getType().equals(CellType.EMPTY))).collect(Collectors.toList());
      String name = null;
      int dbs = -1;
      int oh = LipidomicsConstants.EXCEL_NO_OH_INFO;
      int paramCharge = 1;
      int charge = -1;
      String modification = null;
      String formula = null;
      String modFormula = null;
      double preciseRT = 0.0;
      float area = 0f;
      float areaError = 0f;
      float background = 0f; 
      float mz = 0f;
      float mzTolerance = 0f;
      float peak = 0f;
      float lowerValley = 0f;
      float upperValley = 0f;
      float apexIntensity = 0f;
      float lowerValley10Pc = 0f;
      float upperValley10Pc = 0f;
      float lowerValley50Pc = 0f;
      float upperValley50Pc = 0f;
      float lowerMz10Pc = -1f;
      float upperMz10Pc = -1f;
      float lowerMz50Pc = -1f;
      float upperMz50Pc = -1f;
      int isotope = -1;
      float lowMz = -1;
      float upMz = -1;
      float ellipseTimePosition = -1f;
      float ellipseMzPosition = -1f;
      float ellipseTimeStretch = -1f;
      float ellipseMzStretch = -1f;
      float lowerRtHardLimit = -1f;
      float upperRtHardLimit = -1f;
      float percentalSplit = -1f;
      String oxState="";
      
      int index;
      String rawValue;
      
      for (Cell cell : cells) {
        index  = cell.getColumnIndex();
        rawValue = cell.getRawValue();
        
        if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_NAME)) {
          name = rawValue;
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_DBS)) {
          dbs = (int)Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(LipidomicsConstants.EXCEL_MS_OH)) {
          oh = (int)Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_MODIFICATION)) {
          modification = rawValue;
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_FORMULA)) {
          formula = rawValue;
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_MOD_FORMULA)) {
          modFormula = rawValue;
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_RT)) {
          preciseRT = Double.parseDouble(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_ISOTOPE)) {
          isotope = (int)Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_AREA)) {
          area = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_AREA_ERROR)) {
          areaError = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_BACKGROUND)) {
          background = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_CHARGE)) {
          //we left this here for *possible* backward compatibility issues...
          paramCharge = (int)Float.parseFloat(rawValue);
          charge = paramCharge;
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_MZ_MS1)) {
          mz = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_MZ_TOLERANCE)) {
          mzTolerance = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_PEAK)) {
          peak = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_LOWER_VALLEY)) {
          lowerValley = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_UPPER_VALLEY)) {
          upperValley = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_LOWER_MZ)) {
          lowMz = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_UPPER_MZ)) {
          upMz = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_ELL_CENT_TIME)) {
          ellipseTimePosition = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_ELL_CENT_MZ)) {
          ellipseMzPosition = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_ELL_STRETCH_TIME)) {
          ellipseTimeStretch = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_ELL_STRETCH_MZ)) {
          ellipseMzStretch = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_LOWER_RT_HARD_LIMIT)) {
          lowerRtHardLimit = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_UPPER_RT_HARD_LIMIT)) {
          upperRtHardLimit = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_PERCENTAL_SPLIT)) {
          percentalSplit = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_RAW_APEX)) {
          apexIntensity = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_LOWER_VALLEY10PC)) {
          lowerValley10Pc = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_LOWER_VALLEY50PC)) {
          lowerValley50Pc = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_UPPER_VALLEY10PC)) {
          upperValley10Pc = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_UPPER_VALLEY50PC)) {
          upperValley50Pc = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_LOWER_MZ10PC)) {
          lowerMz10Pc = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_LOWER_MZ50PC)) {
          lowerMz50Pc = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_UPPER_MZ10PC)) {
          upperMz10Pc = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(QuantificationResultExporter.HEADER_UPPER_MZ50PC)) {
          upperMz50Pc = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf(LipidomicsConstants.CHAIN_MOD_COLUMN_NAME)) {
        	oxState = rawValue;
        }

      }
      
      if (name!=null&&name.length()>0){
        if (params!=null){
          // this is for backward compatibility
          if (params.ProbeCount()>0)
            params.setCharge(params.Probe(0).Charge);
          resultParams.add(params);
          if (analyteNames.containsKey(params.getNameString())) showModification = true;
          analyteNames.put(params.getNameString(), params.getNameString());
        }
        //this is for backward compatibility TODO: remove after a suitable transition period, written 22.08.2023
        if (headerTitles.indexOf(QuantificationResultExporter.HEADER_MODIFICATION) == -1 ||
            headerTitles.indexOf(QuantificationResultExporter.HEADER_FORMULA) == -1 || 
            headerTitles.indexOf(QuantificationResultExporter.HEADER_MOD_FORMULA) == -1){
          Object[] components = ComparativeNameExtractor.splitOldNameStringToComponents(name);
          name = (String)components[0];
          dbs = (Integer)components[1];
          formula = (String)components[2];
          modification = "";
          modFormula = "";
        }
        params = new LipidParameterSet(mz, name, dbs, modification, preciseRT, formula, modFormula, paramCharge, oh);
        if (lowerRtHardLimit>=0) params.setLowerRtHardLimit(lowerRtHardLimit);
        if (upperRtHardLimit>=0) params.setUpperRtHardLimit(upperRtHardLimit);
        if (percentalSplit>=0) params.setPercentalSplit(percentalSplit);
        params.Area = area;
        params.LowerMzBand = mzTolerance;
        params.UpperMzBand = mzTolerance;
        params.setOxState(oxState); 
      }else{
        if (params!=null){
          //due to this any row without a value in the 'charge' column will be ignored
          if (charge!=-1){
            CgProbe probe = new CgProbe(0,charge);
            if (area>0 || params.getLowerRtHardLimit()>=0 || params.getUpperRtHardLimit()>=0){
              probe.AreaStatus = CgAreaStatus.OK;
              probe.Area = area;
              probe.AreaError = areaError;
              probe.Background = background;
              probe.Peak = peak;
              probe.LowerValley = lowerValley;
              probe.UpperValley = upperValley;
              if (upperValley10Pc>0f){
                probe.setApexIntensity(apexIntensity);
                probe.setLowerValley10(lowerValley10Pc);
                probe.setLowerValley50(lowerValley50Pc);
                probe.setUpperValley10(upperValley10Pc);
                probe.setUpperValley50(upperValley50Pc);
              }
            }else{
              probe.AreaStatus = CgAreaStatus.TooSmall;
            }            
            probe.Mz = mz;
            probe.LowerMzBand = mzTolerance;
            probe.UpperMzBand = mzTolerance;
            probe.isotopeNumber = isotope;
            if (ellipseTimePosition>0&&ellipseMzPosition>0&&ellipseTimeStretch>0&&ellipseMzStretch>0&&
                lowMz>0&&upMz>0){
              Probe3D probe3D = new Probe3D(probe,ellipseTimePosition,ellipseMzPosition,
                  ellipseTimeStretch,ellipseMzStretch,-1f,-1f,lowerMz10Pc,upperMz10Pc,lowerMz50Pc,upperMz50Pc);
              probe3D.LowerMzBand = lowMz;
              probe3D.UpperMzBand = upMz;
              if (params.getLowerRtHardLimit()>=0) probe3D.setLowerHardRtLimit(params.getLowerRtHardLimit());
              if (params.getUpperRtHardLimit()>=0) probe3D.setUpperHardRtLimit(params.getUpperRtHardLimit());
              try {
                params.AddProbe(probe3D);
              } catch (CgException ex) {}
            }else
              try {
                params.AddProbe(probe);
              } catch (CgException ex) {} 
          }
        }
      }
    }
    if (params!=null){
      // this is for backward compatibility
      if (params.ProbeCount()>0)
        params.setCharge(params.Probe(0).Charge);
      resultParams.add(params);
      if (analyteNames.containsKey(params.getNameString())) showModification = true;
      analyteNames.put(params.getNameString(), params.getNameString());
    }
    if (!showModification){
      String modificationString = null;
      for (LipidParameterSet set : resultParams){
        if (modificationString==null) modificationString = set.getModificationName();
        if (!modificationString.equalsIgnoreCase(set.getModificationName())){
          showModification = true;
          break;
        }
      }
    }
    resultParameterSets_.put(sheet.getName(), resultParams);
    showModifications.put(sheet.getName(), showModification);
    msLevels_.put(sheet.getName(), msLevel);
  }
  
  
  /**
   * Parses the header of an Excel sheet
   * @param headerRow row number of the header
   * @return List of the header titles
   */
  public static List<String> readSheetHeaderTitles(Row headerRow) {
    try (Stream<Cell> cells = headerRow.stream();) {
      return cells.map((c) -> (!(c==null || c.getType().equals(CellType.ERROR)) ? c.getText() : "null")).collect(Collectors.toList());
    } 
  }  
  
  
  /**
   * Compares two sheets by their sheet names. The comparison is based on whether the sheet ends with one of the sheet name adducts used by the QuantificationResultExporter.
   * 
   * @param s1 an Excel sheet
   * @param s2 the Excel sheet to be compared to s1
   * @return the value 0 if both strings either end with a sheet adduct or not; a value less than 0 if the first sheet name does not end with a sheet adduct; a value greater than 0 if the first sheet name does end with a sheet adduct.
   */
  private static int compareBySheetName(Sheet s1, Sheet s2) {
    String name1 = s1.getName();
    String name2 = s2.getName();   
    
    if ((endsWithSheetAdduct(name1) && endsWithSheetAdduct(name2)) ||
        (!endsWithSheetAdduct(name1) && !endsWithSheetAdduct(name2))) {
      return 0;
    } else if (!endsWithSheetAdduct(name1)) {
      return -1;
    } else {
      return 1;
    }
  }
  
  
  /**
   * Determines whether a String ends with one of the sheet name adducts used by the QuantificationResultExporter.
   * @param name An Excel sheet name
   * @return true if the name ends with one of the sheet name adducts used by the QuantificationResultExporter.
   */
  private static boolean endsWithSheetAdduct(String name) {
    if (!name.endsWith(QuantificationResultExporter.ADDUCT_OMEGA_SHEET)&&
        !name.endsWith(QuantificationResultExporter.ADDUCT_OVERVIEW_SHEET)&&
        !name.endsWith(QuantificationResultExporter.ADDUCT_MSN_SHEET)) {
      return false;
    }
    return true;
  }
  
  
  /**
   * Determines whether the parameter Excel sheet should be parsed
   * @param sheet an Excel sheet
   * @param specificClass filter for parsing only the results of one analyte class
   * @return true if the sheet should be parsed
   */
  private static boolean parseSheet(Sheet sheet, String specificClass) {
    boolean parseSheet = false;
    String name = sheet.getName();
    if (specificClass != null) {
      if (name.equals(QuantificationResultExporter.SHEET_CONSTANTS) ||
          name.equals(specificClass) || 
          name.equals(specificClass+QuantificationResultExporter.ADDUCT_MSN_SHEET) || 
          name.equals(specificClass+QuantificationResultExporter.ADDUCT_OMEGA_SHEET)) {
        parseSheet = true;
      }
    } else if (!name.endsWith(QuantificationResultExporter.ADDUCT_OVERVIEW_SHEET)) {
      parseSheet = true;
    }
    return parseSheet;
  }
  
  
  /**
   * cleans unidentified position definitions from the position definition hash (a negative integer value corresponds to unidentified)
   * @param posDefs the uncleaned position definition hash
   * @return cleaned position definition hash
   */
  private static Hashtable<String,Hashtable<Integer,Integer>> cleanPositionDefinition (Hashtable<String,Hashtable<Integer,Integer>> posDefs){
    Hashtable<String,Hashtable<Integer,Integer>> cleaned = new Hashtable<String,Hashtable<Integer,Integer>>();
    for (String combi : posDefs.keySet()){
      Hashtable<Integer,Integer> posss = posDefs.get(combi);
      Hashtable<Integer,Integer> possCleaned = new Hashtable<Integer,Integer>();
      Hashtable<Integer,Integer> usedPositions = new Hashtable<Integer,Integer>();
      for (Integer stringPos : posss.keySet()){
        if (posss.get(stringPos)>-1){
          if (usedPositions.containsKey(posss.get(stringPos))){
            //the same position is referenced more than once -> no evidence remove it
            if (possCleaned.containsKey(usedPositions.get(posss.get(stringPos))))possCleaned.remove(usedPositions.get(posss.get(stringPos)));
          }else{
            possCleaned.put(stringPos, posss.get(stringPos));
            usedPositions.put(posss.get(stringPos), stringPos);
          }
        }
      }
      cleaned.put(combi, possCleaned);
    }
    return cleaned;
  }
  
  
  /**
   * the parser adds found chain combinations by replacing the "/" by "_", where the sequence of the fatty acids may
   * be different than in the stored values - this method corrects the sequence to be exactly as in the stored Excel
   * @param combiKey fa combination String as it should be
   * @param validChainCombinations the currently available chain combinations
   * @param relAreas the relative share of each chain combination; first key: chain combination; second key: relative share on MS1 area
   * @param sep the separator for the combination
   * @return the corrected available chain combinations
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private static Vector<String> correctUndefinedChainCombinations(String combiKey, Vector<String> validChainCombinations, Hashtable<String,Double> relAreas, String sep) throws LipidCombinameEncodingException{
    Vector<String> corrected = new Vector<String>();
    for (String combi : validChainCombinations){
      if (StaticUtils.isAPermutedVersion(combiKey,combi,sep)) {
        corrected.add(combiKey);
        if (!combiKey.equalsIgnoreCase(combi)) {
          relAreas.put(combiKey, relAreas.get(combi));
          relAreas.remove(combi);
        }
      }
      else corrected.add(combi);
    }
    return corrected;
  }
  
  
  /**
   * checks if the position definition hash contains any defined position
   * @param posDefs position definition hash
   * @return true if definitions are present
   */
  private static boolean isDefinitionPresent(Hashtable<String,Hashtable<Integer,Integer>> posDefs){
    for (Hashtable<Integer,Integer> posss : posDefs.values()){
      if (posss.size()>0) return true;
    }
    return false;
  }
  
  
  /**
   * selects the permutation that is used in the validChainCombinations
   * @param combiKey an encoded LDA chain combination
   * @param validChainCombinations the available chain combinations in valid validChainCombinations
   * @return the one in validChainCombinations
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private static String getPermutationThatIsInValidChainCombinations(String combiKey, Vector<String> validChainCombinations) throws LipidCombinameEncodingException {
    Vector<String> chains = StaticUtils.splitChainCombiToEncodedStrings(combiKey,LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
    Vector<String> permuts = StaticUtils.getPermutedChainNames(chains, LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
    for (String permut : permuts) {
      if (validChainCombinations.contains(permut))
        return permut;
    }
    return null;
  }
  
}
