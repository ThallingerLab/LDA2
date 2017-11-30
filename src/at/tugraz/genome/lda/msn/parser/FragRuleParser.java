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

package at.tugraz.genome.lda.msn.parser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.FattyAcidsContainer;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.msn.vos.ExpressionForComparisonVO;
import at.tugraz.genome.lda.msn.vos.FragmentMultVO;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.rdi.GeneralSettingsVO;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class FragRuleParser
{
  
  // booleans specifying which attributes are covered by the rules
  boolean foundGeneral_;
  boolean foundHead_;
  boolean foundChains_;
  boolean foundPosition_;
  
  // identifiers for the general sections in the rule file 
  private final static int NO_SECTION = 0;
  private final static int GENERAL_SECTION = 1;
  public final static int HEAD_SECTION = 2;
  public final static int CHAINS_SECTION = 3;
  public final static int POSITION_SECTION = 4;
  
  // identifiers for the general subsections in the rule file 
  private final static int FRAGMENT_SUBSECTION = 1;
  private final static int INTENSITY_SUBSECTION = 2;
  
  // tags for general sections in the rule file
  private final static String GENERAL_SECTION_NAME = "[GENERAL]";
  private final static String HEAD_SECTION_NAME = "[HEAD]";
  private final static String CHAINS_SECTION_NAME = "[CHAINS]";
  private final static String POSITION_SECTION_NAME = "[POSITION]";
  
  //tags for general subsections in the rule file
  private final static String FRAGMENT_SUBSECTION_NAME = "!FRAGMENTS";
  private final static String INTENSITY_SUBSECTION_NAME = "!INTENSITIES";
  
  //general properties of the rule file
  private final static String GENERAL_CHAINS = "AmountOfChains";
  private final static String GENERAL_MS2_LIB = "ChainLibrary";
  public final static String GENERAL_CATOMS_PARSE = "CAtomsFromName";
  public final static String GENERAL_DBOND_PARSE = "DoubleBondsFromName";
  private final static String GENERAL_CUTOFF = "BasePeakCutoff";
  private final static String GENERAL_CHAIN_CUTOFF = "ChainCutoff";
  private final static String GENERAL_SPECTRUM_COVERAGE = "SpectrumCoverage";
  private final static String GENERAL_RT_PROCESSING = "RetentionTimePostprocessing";
  private final static String GENERAL_RT_PARALLEL_SERIES = "RetentionTimeParallelSeries";
  private final static String GENERAL_RT_DEV_MAX = "RetentionTimeMaxDeviation";
  private final static String GENERAL_CHAINS_ALKYL = "AlkylChains";
  private final static String GENERAL_CHAINS_ALKENYL = "AlkenylChains";
  private final static String GENERAL_SINGLE_CHAIN = "SingleChainIdentification";
  private final static String GENERAL_IDENTIFICATION_ORDER = "MSIdentificationOrder";
  private final static String GENERAL_PEAK_UNION_TIME = "EnforcePeakUnionTime";
  private final static String GENERAL_PEAK_UNION_NO_POSITION = "IgnorePositionForUnion";
  private final static String GENERAL_MS1_PEAK_CUTOFF = "ClassSpecificMS1Cutoff";
  private final static String GENERAL_ADD_POSITIONS = "AddChainPositions";
  private final static String GENERAL_ISOBAR_RATIO = "IsobarSCExclusionRatio";
  private final static String GENERAL_ISOBAR_FAR_RATIO = "IsobarSCFarExclusionRatio";
  private final static String GENERAL_ISOBAR_RT = "IsobarRtDiff";

  //properties of the fragment subsection of the rule file
  private final static String FRAGMENT_NAME = "Name";
  private final static String FRAGMENT_FORMULA = "Formula";
  private final static String FRAGMENT_CHARGE = "Charge";
  private final static String FRAGMENT_LEVEL = "MSLevel";
  private final static String FRAGMENT_MANDATORY = "mandatory";
  
  // the various methods for the identification order
  private final static String ORDER_MS1_FIRST = "MS1First";
  private final static String ORDER_MSN_FIRST = "MSnFirst";
  private final static String ORDER_MSN_ONLY = "MSnOnly";
  
  //properties of the intensity subsection of the rule file
  private final static String INTENSITY_EQUATION = "Equation";
  
  //general properties hash - key is are the properties two paragraphs above
  private Hashtable<String,String> generalSettings_;
  // head fragment rules - key is the rule name of the fragment
  private Hashtable<String,FragmentRuleVO> headFragments_;
  // chain fragment rules - key is the name of the fragment
  private Hashtable<String,FragmentRuleVO> chainFragments_;
  // head rules for intensity comparisons
  private Vector<IntensityRuleVO> headIntensities_;
  // chain rules for intensity comparisons
  private Vector<IntensityRuleVO> chainIntensities_;
  // position rules for intensity comparisons
  private Vector<IntensityRuleVO> positionIntensities_;
  
  // requires an ElementConfigParser to evaluate the chemical formulas
  private ElementConfigParser elementParser_;
  
  public final static String NO_HEAD_AND_CHAINS_SECTION = "The rules file must contain a "+HEAD_SECTION_NAME+" or a "+CHAINS_SECTION_NAME+" section!";
  
  /**
   * Constructor requires and ElementConfig parser only - for evaluating the chemical formulas
   * @param elementParser for evaluating the chemical formulas
   */
  public FragRuleParser(ElementConfigParser elementParser){
    elementParser_ = elementParser;
  }
  
  /**
   * Command that starts the parsing of and fragmentation rule file
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException general exception if a file is not there
   */
  public void parseFile(String filePath) throws IOException, RulesException {
    parseFile(new File(filePath));
  }
    
  /**
   * parses the a frag.txt file and stores the values
   * @param file
   * @return
   * @throws IOException
   * @throws RulesException
   */
  public void parseFile(File file) throws IOException, RulesException {
    if (!file.exists()) throw new IOException("The file "+file.getAbsolutePath()+" does not exist!");
    initParseableValues();
    foundGeneral_ = false;
    foundHead_ = false;
    foundChains_ = false;
    foundPosition_ = false;
    LineNumberReader reader = null;
    String line;
    int currentSection = NO_SECTION;
    int subSection = NO_SECTION;
    try{
      reader = new LineNumberReader(new FileReader(file));
      while ((line = reader.readLine()) != null) {
        line = line.trim();
        // this if asks for the current section
        if (line.startsWith("[") && line.endsWith("]")){
          subSection = NO_SECTION;
          // set session specific booleans
          if (line.equalsIgnoreCase(GENERAL_SECTION_NAME)){
            foundGeneral_ = true;
            currentSection = GENERAL_SECTION;
          } else if (line.equalsIgnoreCase(HEAD_SECTION_NAME)){
            foundHead_ = true;
            currentSection = HEAD_SECTION;
          } else if (line.equalsIgnoreCase(CHAINS_SECTION_NAME)){
            foundChains_ = true;
            currentSection = CHAINS_SECTION;
          } else if (line.equalsIgnoreCase(POSITION_SECTION_NAME)){
            foundPosition_ = true;
            currentSection = POSITION_SECTION;
          } else {
            try{reader.close();}catch(Exception ex){}
            throw new RulesException("A section with the name "+line+" is not supported by the fragmentation rules!");
          }
        // what is the subsection
        }else if (line.startsWith("!")){
          if (currentSection==NO_SECTION){
            try{reader.close();}catch(Exception ex){}
            throw new RulesException("A subsection cannot start before a section starts! Error at line number "+reader.getLineNumber()+"!");
          }
          if (currentSection==GENERAL_SECTION){
            try{reader.close();}catch(Exception ex){}
            throw new RulesException("The "+GENERAL_SECTION_NAME+" section must not contain any subsections! Error at line number "+reader.getLineNumber()+"!");
          }
          if (line.equalsIgnoreCase(FRAGMENT_SUBSECTION_NAME)) subSection = FRAGMENT_SUBSECTION;
          else if (line.equalsIgnoreCase(INTENSITY_SUBSECTION_NAME)) subSection = INTENSITY_SUBSECTION;
          else{
            try{reader.close();}catch(Exception ex){}
            throw new RulesException("A section with the name "+line+" is not supported by the fragmentation rules! Error at line number "+reader.getLineNumber()+"!");
          }
          if (currentSection==POSITION_SECTION && subSection==FRAGMENT_SUBSECTION){
            try{reader.close();}catch(Exception ex){}
            throw new RulesException("The section "+POSITION_SECTION_NAME+" does not support the subsection "+FRAGMENT_SUBSECTION_NAME+"! Error at line number "+reader.getLineNumber()+"!");
          }
        // entries are parsed  
        }else{
          if (currentSection==GENERAL_SECTION){
            parseGeneralSectionEntry(line,reader.getLineNumber());
          } else {
            if (subSection==FRAGMENT_SUBSECTION)
              parseFragmentEntry(line,reader.getLineNumber(),currentSection);
            else if (subSection==INTENSITY_SUBSECTION)
              parseIntensityEntry(line,reader.getLineNumber(),currentSection);
          }
        }
      }
      if (!foundGeneral_) throw new RulesException("The rules file does not contain the mandatory "+GENERAL_SECTION_NAME+" section!");
      if (!foundHead_ && !foundChains_) throw new RulesException(NO_HEAD_AND_CHAINS_SECTION);
      if (headFragments_.size()==0 && chainFragments_.size()==0) throw new RulesException("The rules file does not contain any fragments in the "+FRAGMENT_SUBSECTION+"! There must be at least one fragment!");
      checkIfGeneralValuesAreThere();
    }finally{
      try{if (reader!=null)reader.close();}catch(Exception ex){}
    }
  }
  
  /**
   * init variables used for the storage of parameters of the frag.txt file
   */
  private void initParseableValues(){
    generalSettings_ = new Hashtable<String,String>();
    headFragments_ = new Hashtable<String,FragmentRuleVO>();
    chainFragments_ = new Hashtable<String,FragmentRuleVO>();
    headIntensities_ = new Vector<IntensityRuleVO>();
    chainIntensities_ = new Vector<IntensityRuleVO>();
    positionIntensities_ = new Vector<IntensityRuleVO>();

  }
  
  /**
   * parses propertyEntries of the [GENERAL] section and puts the key/value pairs to the generalSettings_ hash
   * @param line the current line to parse
   * @param lineNumber the line number of the current line
   */
  private void parseGeneralSectionEntry(String line, int lineNumber) throws RulesException {
    if (line.length()==0) return;
    if (line.indexOf("=")==-1) throw new RulesException("Properties in the "+GENERAL_SECTION_NAME+" are key/value pairs seperated by \"=\"! There is no \"=\" at line "+lineNumber+"!");
    String key = line.substring(0,line.indexOf("=")).trim();
    String value = line.substring(line.indexOf("=")+1).trim();
    if (key.equalsIgnoreCase(GENERAL_CHAINS)){
      try{
        Integer.parseInt(value);
        generalSettings_.put(GENERAL_CHAINS, value);
      }catch(NumberFormatException nfx){
        throw new RulesException("The value of "+GENERAL_CHAINS+" must be integer, the value \""+value+"\" is not! Error at line number "+lineNumber+"!");
      }
    }else if (key.equalsIgnoreCase(GENERAL_CHAINS_ALKYL)){
      try{
        Integer.parseInt(value);
        generalSettings_.put(GENERAL_CHAINS_ALKYL, value);
      }catch(NumberFormatException nfx){
        throw new RulesException("The value of "+GENERAL_CHAINS_ALKYL+" must be integer, the value \""+value+"\" is not! Error at line number "+lineNumber+"!");
      }
    }else if (key.equalsIgnoreCase(GENERAL_CHAINS_ALKENYL)){
      try{
        Integer.parseInt(value);
        generalSettings_.put(GENERAL_CHAINS_ALKENYL, value);
      }catch(NumberFormatException nfx){
        throw new RulesException("The value of "+GENERAL_CHAINS_ALKENYL+" must be integer, the value \""+value+"\" is not! Error at line number "+lineNumber+"!");
      }
    } else if (key.equalsIgnoreCase(GENERAL_MS2_LIB)){
      if (!value.endsWith(FattyAcidsContainer.FA_FILE_SUFFIX_NEW) && !value.endsWith(FattyAcidsContainer.FA_FILE_SUFFIX_OLD)) throw new RulesException("The value of "+GENERAL_MS2_LIB+" must be an Excel file; the value \""+value+"\" is not! Error at line number "+lineNumber+"!");
      generalSettings_.put(GENERAL_MS2_LIB, value);
    } else if (key.equalsIgnoreCase(GENERAL_CATOMS_PARSE)){
      try{
        Pattern.compile(value);
        if (value.indexOf("(")==-1 || value.indexOf(")")==-1)
          throw new RulesException("The value of "+GENERAL_CATOMS_PARSE+" must be a Java regex containing \"(\" and \")\"; the value \""+value+"\" is not valid! Error at line number "+lineNumber+"!");
        generalSettings_.put(GENERAL_CATOMS_PARSE, value);
      }catch(PatternSyntaxException psx){
        throw new RulesException("The value of "+GENERAL_CATOMS_PARSE+" must be a Java regex; the value \""+value+"\" is not valid! Error at line number "+lineNumber+"!");
      }
    } else if (key.equalsIgnoreCase(GENERAL_DBOND_PARSE)){
      try{
        Pattern.compile(value);
        if (value.indexOf("(")==-1 || value.indexOf(")")==-1)
          throw new RulesException("The value of "+GENERAL_DBOND_PARSE+" must be a Java regex containing \"(\" and \")\"; the value \""+value+"\" is not valid! Error at line number "+lineNumber+"!");
        generalSettings_.put(GENERAL_DBOND_PARSE, value);
      }catch(PatternSyntaxException psx){
        throw new RulesException("The value of "+GENERAL_DBOND_PARSE+" must be a Java regex; the value \""+value+"\" is not valid! Error at line number "+lineNumber+"!");
      }
    } else if (key.equalsIgnoreCase(GENERAL_CUTOFF)){
      readPercentPermilleValue(value,GENERAL_CUTOFF,lineNumber);
      generalSettings_.put(GENERAL_CUTOFF, value);
    } else if (key.equalsIgnoreCase(GENERAL_CHAIN_CUTOFF)){
      readPercentPermilleValue(value,GENERAL_CHAIN_CUTOFF,lineNumber);
      generalSettings_.put(GENERAL_CHAIN_CUTOFF, value);
    } else if (key.equalsIgnoreCase(GENERAL_SPECTRUM_COVERAGE)){  
      readPercentPermilleValue(value,GENERAL_SPECTRUM_COVERAGE,lineNumber);
      generalSettings_.put(GENERAL_SPECTRUM_COVERAGE, value); 
    } else if (key.equalsIgnoreCase(GENERAL_RT_PROCESSING)){
      boolean rtPostProcessing = false;
      if (value!=null && (value.equalsIgnoreCase("true")||value.equalsIgnoreCase("yes"))) rtPostProcessing = true;
      generalSettings_.put(GENERAL_RT_PROCESSING, String.valueOf(rtPostProcessing));
    } else if (key.equalsIgnoreCase(GENERAL_RT_PARALLEL_SERIES)){
      boolean rtParallelSeries = false;
      if (value!=null && (value.equalsIgnoreCase("true")||value.equalsIgnoreCase("yes"))) rtParallelSeries = true;
      generalSettings_.put(GENERAL_RT_PARALLEL_SERIES, String.valueOf(rtParallelSeries));
    } else if (key.equalsIgnoreCase(GENERAL_RT_DEV_MAX)){
      try{
        new Float(value);
        generalSettings_.put(GENERAL_RT_DEV_MAX, value);
      } catch (NumberFormatException nfx){
        throw new RulesException("The value of "+GENERAL_RT_DEV_MAX+" must be float format; the value \""+value+"\" is not valid! Error at line number "+lineNumber+"!");
      }
    } else if (key.equalsIgnoreCase(GENERAL_SINGLE_CHAIN)){
      boolean singleChain = false;
      if (value!=null && (value.equalsIgnoreCase("true")||value.equalsIgnoreCase("yes"))) singleChain = true;
      generalSettings_.put(GENERAL_SINGLE_CHAIN, String.valueOf(singleChain));
    } else if (key.equalsIgnoreCase(GENERAL_IDENTIFICATION_ORDER)){
      int msIdentificationOrder = RulesContainer.ORDER_MS1_FIRST;
      if (!value.equalsIgnoreCase(ORDER_MS1_FIRST)&&!value.equalsIgnoreCase(ORDER_MSN_FIRST)&&!value.equalsIgnoreCase(ORDER_MSN_ONLY))
        throw new RulesException("The value \""+value+"\" is not allowed for \""+GENERAL_IDENTIFICATION_ORDER+"\"! Only "+ORDER_MS1_FIRST+"/"+ORDER_MSN_FIRST+"/"+ORDER_MSN_ONLY+" is allowed! Error at line number "+lineNumber+"!");
      if (value.equalsIgnoreCase(ORDER_MSN_FIRST)) msIdentificationOrder = RulesContainer.ORDER_MSN_FIRST;
      if (value.equalsIgnoreCase(ORDER_MSN_ONLY)) msIdentificationOrder = RulesContainer.ORDER_MSN_ONLY;
      generalSettings_.put(GENERAL_IDENTIFICATION_ORDER, String.valueOf(msIdentificationOrder));
    } else if (key.equalsIgnoreCase(GENERAL_PEAK_UNION_TIME)){
      try{
        new Float(value);
        generalSettings_.put(GENERAL_PEAK_UNION_TIME, value);
      } catch (NumberFormatException nfx){
        throw new RulesException("The value of "+GENERAL_PEAK_UNION_TIME+" must be float format; the value \""+value+"\" is not valid! Error at line number "+lineNumber+"!");
      }
    } else if (key.equalsIgnoreCase(GENERAL_PEAK_UNION_NO_POSITION)){
      boolean noPosition = false;
      if (value!=null && (value.equalsIgnoreCase("true")||value.equalsIgnoreCase("yes"))) noPosition = true;
      generalSettings_.put(GENERAL_PEAK_UNION_NO_POSITION, String.valueOf(noPosition));      
    } else if (key.equalsIgnoreCase(GENERAL_MS1_PEAK_CUTOFF)){
      try{
        new Float(value);
        generalSettings_.put(GENERAL_MS1_PEAK_CUTOFF, value);
      } catch (NumberFormatException nfx){
        throw new RulesException("The value of "+GENERAL_MS1_PEAK_CUTOFF+" must be float format; the value \""+value+"\" is not valid! Error at line number "+lineNumber+"!");
      }
    } else if (key.equalsIgnoreCase(GENERAL_ADD_POSITIONS)){
      try{
        int adds = Integer.parseInt(value);
        if (adds<0) throw new RulesException("The value of "+GENERAL_ADD_POSITIONS+" must be bigger or equal zero, the value \""+value+"\" is not! Error at line number "+lineNumber+"!");
        generalSettings_.put(GENERAL_ADD_POSITIONS, value);
      }catch(NumberFormatException nfx){
        throw new RulesException("The value of "+GENERAL_ADD_POSITIONS+" must be integer, the value \""+value+"\" is not! Error at line number "+lineNumber+"!");
      }      
    } else if (key.equalsIgnoreCase(GENERAL_ISOBAR_RATIO)){
      try{
        new Float(value);
        generalSettings_.put(GENERAL_ISOBAR_RATIO, value);
      } catch (NumberFormatException nfx){
        throw new RulesException("The value of "+GENERAL_ISOBAR_RATIO+" must be float format; the value \""+value+"\" is not valid! Error at line number "+lineNumber+"!");
      }
    } else if (key.equalsIgnoreCase(GENERAL_ISOBAR_FAR_RATIO)){
      try{
        new Float(value);
        generalSettings_.put(GENERAL_ISOBAR_FAR_RATIO, value);
      } catch (NumberFormatException nfx){
        throw new RulesException("The value of "+GENERAL_ISOBAR_FAR_RATIO+" must be float format; the value \""+value+"\" is not valid! Error at line number "+lineNumber+"!");
      }
    } else if (key.equalsIgnoreCase(GENERAL_ISOBAR_RT)){
      try{
        new Float(value);
        generalSettings_.put(GENERAL_ISOBAR_RT, value);
      } catch (NumberFormatException nfx){
        throw new RulesException("The value of "+GENERAL_ISOBAR_RT+" must be float format; the value \""+value+"\" is not valid! Error at line number "+lineNumber+"!");
      }
    } else {
      throw new RulesException("The section "+GENERAL_SECTION_NAME+" does not support the property "+key+"! Error at line number "+lineNumber+"!");
    }
  }
  

  
  /**
   * parses an input in fraction, percent and permille format and returns the corresponding double value
   * ATTENTION: the value must be positive and lower than 100%
   * @param inValue the String representation of the fraction, percent or permille format
   * @param section the type of the value
   * @param lineNumber at which line number of file was this value found.
   * @return the double representaion of the fraction, percent or permille format
   * @throws RulesException if the value is not a number or not in the range of 0-100%, exclusive 100% 
   */
  public static double readPercentPermilleValue(String inValue,String section, int lineNumber) throws RulesException {
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
      if (doubleValue<0){
        throw new RulesException("The "+section+" must not be negative, the value \""+value+"\" is not! Error at line number "+lineNumber+"!");
      } else if (doubleValue>=1){
        if (percent) value+="%";
        else if (permille) value+="\u2030";
        throw new RulesException("The "+section+" must not be bigger than 100%, the value \""+value+"\" is not! Error at line number "+lineNumber+"!");
      }
      return doubleValue;
    }catch(NumberFormatException nfx){
      throw new RulesException("The value of "+section+" must be double format, the value \""+value+"\" is not! Error at line number "+lineNumber+"!");
    }   
  }
  
  /**
   * checks if required values of the [GENERAL] are present
   * @throws RulesException
   */
  private void checkIfGeneralValuesAreThere() throws RulesException {
    if (!generalSettings_.containsKey(GENERAL_CHAINS)) throw new RulesException("The rules file must contain the property \""+GENERAL_CHAINS+"\" in the "+GENERAL_SECTION_NAME+" section!");
    if (!generalSettings_.containsKey(GENERAL_MS2_LIB)) throw new RulesException("The rules file must contain the property \""+GENERAL_MS2_LIB+"\" in the "+GENERAL_SECTION_NAME+" section!");
    if (!generalSettings_.containsKey(GENERAL_CATOMS_PARSE)) throw new RulesException("The rules file must contain the property \""+GENERAL_CATOMS_PARSE+"\" in the "+GENERAL_SECTION_NAME+" section!");
    if (!generalSettings_.containsKey(GENERAL_DBOND_PARSE)) throw new RulesException("The rules file must contain the property \""+GENERAL_DBOND_PARSE+"\" in the "+GENERAL_SECTION_NAME+" section!");
    int amountOfChains =  Integer.parseInt(getAmountOfChains());
    int amountOfAlkylChains = Integer.parseInt(getAmountOfAlkylChains());
    int amountOfAlkenylChains = Integer.parseInt(getAmountOfAlkenylChains());
    if ((amountOfAlkylChains+amountOfAlkenylChains)>amountOfChains)
      throw new RulesException("There must not be more \""+GENERAL_CHAINS_ALKYL+"\" and \""+GENERAL_CHAINS_ALKENYL+"\" than \""+GENERAL_CHAINS+"\"!");
  }
  
  /**
   * parses propertyEntries of the !FRAGMENTS section and puts FragmentRuleVOs to headFragments_ and chainFragments_ hashes
   * @param line one line in the rule file
   * @param lineNumber the line number of the line
   * @param currentSection the section in the file - GENERAL_SECTION, HEAD_SECTION, or CHAINS_SECTION
   * @throws RulesException
   */
  private void parseFragmentEntry(String line, int lineNumber, int currentSection) throws RulesException {
    if (line.length()==0) return;
    int charge = 1;
    int msLevel = 2;
    boolean mandatory = false;
    boolean fromOtherSpecies = false;
    String name = null;
    String formula = null;
    StringTokenizer tokenizer = new StringTokenizer(line,"\t ");
    if (Settings.useAlex()) tokenizer = new StringTokenizer(line,"\t");
    while (tokenizer.hasMoreTokens()){
      String kvPair = tokenizer.nextToken().trim();
      if (kvPair.indexOf("=")==-1) throw new RulesException("The "+FRAGMENT_SUBSECTION_NAME+" are key/value pairs seperated by \"=\"! There is no \"=\" in \""+kvPair+"\"at line "+lineNumber+"!");
      String key = kvPair.substring(0,kvPair.indexOf("="));
      String value = kvPair.substring(kvPair.indexOf("=")+1);
      if (key.equalsIgnoreCase(FRAGMENT_NAME)){
        if (value.indexOf("$")!=-1) throw new RulesException("The key "+FRAGMENT_NAME+" does not support values containing \"$\"! Error at line number "+lineNumber+"!");
        name = value;
      }else if(key.equalsIgnoreCase(FRAGMENT_FORMULA)){
        checkFormula(value, lineNumber);
        formula = value;
      }else if(key.equalsIgnoreCase(FRAGMENT_CHARGE)){
        try{
          charge = Integer.parseInt(value);
        }catch(NumberFormatException nfx){throw new RulesException("The value of "+FRAGMENT_CHARGE+" must be integer, the value \""+value+"\" is not! Error at line number "+lineNumber+"!");}
      }else if(key.equalsIgnoreCase(FRAGMENT_LEVEL)){
        try{
          msLevel = Integer.parseInt(value);
        }catch(NumberFormatException nfx){throw new RulesException("The value of "+FRAGMENT_LEVEL+" must be integer, the value \""+value+"\" is not! Error at line number "+lineNumber+"!");}
      }else if(key.equalsIgnoreCase(FRAGMENT_MANDATORY)){
        if (!(value.equalsIgnoreCase("true")||value.equalsIgnoreCase("false")||value.equalsIgnoreCase("yes")||value.equalsIgnoreCase("no")||value.equalsIgnoreCase("other")))
          throw new RulesException("The value of "+FRAGMENT_MANDATORY+" can contain the values \"true\",\"false\",\"yes\", \"no\" and \"other\" only! Error at line number "+lineNumber+"!");
        if (value.equalsIgnoreCase("true")||value.equalsIgnoreCase("yes"))
          mandatory = true;
        else if (value.equalsIgnoreCase("other"))
          fromOtherSpecies = true;
      }else{
        throw new RulesException("The section "+FRAGMENT_SUBSECTION_NAME+" does not support the key \""+key+"\"! Error at line number "+lineNumber+"!");        
      }
    }
    if (name==null || name.length()==0) throw new RulesException("A "+FRAGMENT_SUBSECTION_NAME+" entry must contain a key called \""+FRAGMENT_NAME+"\"! Error at line number "+lineNumber+"!");
    if (formula==null || formula.length()==0) throw new RulesException("A "+FRAGMENT_SUBSECTION_NAME+" entry must contain a key called \""+FRAGMENT_FORMULA+"\"! Error at line number "+lineNumber+"!");      
    checkNamePresent(name, lineNumber);
    FragmentRuleVO ruleVO = new FragmentRuleVO(name,formula,charge,msLevel,mandatory,fromOtherSpecies,headFragments_, chainFragments_, elementParser_);
    if (currentSection==HEAD_SECTION) headFragments_.put(name, ruleVO);
    if (currentSection==CHAINS_SECTION){
      checkSelectedChainValid(ruleVO,lineNumber);
      chainFragments_.put(name, ruleVO);
    }
  }
  
  /**
   * checks if formula of the !FRAGMENTS section is valid
   * @param formula chemical formula
   * @param lineNumber line that is parsed
   * @throws RulesException
   */
  private void checkFormula(String formula, int lineNumber) throws RulesException {
    try{
      FragmentRuleVO.isFormulaValid(formula, headFragments_, chainFragments_, elementParser_);
    } catch (RulesException rlx) {
      throw new RulesException(rlx.getMessage()+" Error at line number "+lineNumber+"!");
    }
  }
  
  /**
   * checks if the fragment name is unique (or it has been used before)
   * @param name name to check
   * @param lineNumber line that is parsed
   * @throws RulesException
   */
  private void checkNamePresent(String name, int lineNumber) throws RulesException {
    if (headFragments_.containsKey(name)) throw new RulesException("The fragment \""+name+"\" at line "+lineNumber+" was already defined in the "+HEAD_SECTION_NAME+"! Names in rules must be unique all over the file!");
    if (chainFragments_.containsKey(name)) throw new RulesException("The fragment "+name+" at line "+lineNumber+" was already defined in the "+CHAINS_SECTION_NAME+"! Names in rules must be unique all over the file!");
  }
  
  /**
   * parses propertyEntries of the !INTENSITIES section and puts IntensityRuleVOs to headFragments_ and chainFragments_ hashes
   * @param line one line in the rule file
   * @param lineNumber the line number of the line
   * @param currentSection the section in the file - GENERAL_SECTION, HEAD_SECTION, CHAINS_SECTION, or POSITION_SECTION
   * 
   * @throws RulesException
   */
  private void parseIntensityEntry(String line,int lineNumber, int currentSection) throws RulesException{
    if (line.length()==0) return;
    boolean mandatory = false;
    StringTokenizer tokenizer = new StringTokenizer(line,"\t ");
    if (Settings.useAlex()) tokenizer = new StringTokenizer(line,"\t");
    IntensityRuleVO ruleVO = null;
    while (tokenizer.hasMoreTokens()){
      String kvPair = tokenizer.nextToken().trim();
      String key = kvPair.substring(0,kvPair.indexOf("="));
      String value = kvPair.substring(kvPair.indexOf("=")+1);
      if (key.equalsIgnoreCase(INTENSITY_EQUATION)){
        Integer amountOfChains = null;
        if (currentSection == POSITION_SECTION){
          if (!this.generalSettings_.containsKey(GENERAL_CHAINS)) throw new RulesException("If there is a [POSITION] section, "+GENERAL_CHAINS+" has to be declared! Error at line number"+lineNumber+"!");
          amountOfChains = new Integer(generalSettings_.get(GENERAL_CHAINS));
          if (this.generalSettings_.containsKey(GENERAL_ADD_POSITIONS)) amountOfChains += new Integer(generalSettings_.get(GENERAL_ADD_POSITIONS));
        }
        ruleVO = extractIntensityVOFromEquation(value, lineNumber, currentSection, FragmentRuleVO.getStringKeyHash(this.headFragments_),FragmentRuleVO.getStringKeyHash(this.chainFragments_), amountOfChains);
      }else if(key.equalsIgnoreCase(FRAGMENT_MANDATORY)){
        if (!(value.equalsIgnoreCase("true")||value.equalsIgnoreCase("false")||value.equalsIgnoreCase("yes")||value.equalsIgnoreCase("no")))
          throw new RulesException("The value of "+FRAGMENT_MANDATORY+" can contain the values \"true\",\"false\",\"yes\" and \"no\" only! Error at line number "+lineNumber+"!");
        if (value.equalsIgnoreCase("true")||value.equalsIgnoreCase("yes"))
          mandatory = true;
      }else{
        throw new RulesException("The section "+INTENSITY_SUBSECTION_NAME+" does not support the key \""+key+"\"! Error at line number "+lineNumber+"!");
      }
    }
    
    if (ruleVO==null) throw new RulesException("An entry of "+INTENSITY_SUBSECTION+" must contain a an \""+INTENSITY_EQUATION+"=\" attribute! Error at line number "+lineNumber+"!");
    ruleVO.setMandatory(mandatory);
    if (currentSection==HEAD_SECTION) headIntensities_.add(ruleVO);
    if (currentSection==CHAINS_SECTION){
      ruleVO.checkForDiffChainTypes(chainFragments_);
      chainIntensities_.add(ruleVO);
    }
    if (currentSection==POSITION_SECTION) positionIntensities_.add(ruleVO);
  }
  
  /**
   * parses an equation for an intensity rule and returns the corresponding IntensityRuleVO
   * @param equation rule equation
   * @param lineNumber line number of the rule (for error message)
   * @param currentSection the type of intensity rule (HEAD, CHAIN, or position)
   * @param headFragments the possible head fragment names
   * @param chainFragments the possible chain fragment names
   * @param amountOfChains only for POSITION_SECTION the maximum number of chains has to be defined
   * @return the VO representative for the rule
   * @throws RulesException if the rule is not possible (e.g. if a fragment has not been defined before)
   */
  public static IntensityRuleVO extractIntensityVOFromEquation(String equation, int lineNumber, int currentSection, Hashtable<String,String> headFragments, Hashtable<String,String> chainFragments, Integer amountOfChains) throws RulesException{
    return extractIntensityVOFromEquation(equation, lineNumber, currentSection, headFragments, chainFragments, amountOfChains, new Hashtable<String,String>());
  }
  
  /**
   * parses an equation for an intensity rule and returns the corresponding IntensityRuleVO
   * @param equation rule equation
   * @param lineNumber line number of the rule (for error message)
   * @param currentSection the type of intensity rule (HEAD, CHAIN, or position)
   * @param headFragments the possible head fragment names
   * @param chainFragments the possible chain fragment names
   * @param amountOfChains only for POSITION_SECTION the maximum number of chains has to be defined
   * @param missed fragments that were not found in the equation (required for reading of results, since for rules containing "+", not all of the fragments have to be found)
   * @return the VO representative for the rule
   * @throws RulesException if the rule is not possible (e.g. if a fragment has not been defined before)
   */
  public static IntensityRuleVO extractIntensityVOFromEquation(String equation, int lineNumber, int currentSection, Hashtable<String,String> headFragments, Hashtable<String,String> chainFragments, Integer amountOfChains, Hashtable<String,String> missed) throws RulesException{
    Boolean biggerThan = null;
    String originalEquation = "";

    originalEquation = new String (equation);
    int comparatorIdx = 0;
    if (equation.indexOf(">")!=-1) comparatorIdx = equation.indexOf(">");
    if (equation.indexOf("<")!=-1) comparatorIdx = equation.indexOf("<");
    if (comparatorIdx<1 || comparatorIdx>=(equation.length()-1)) throw new RulesException("The value of "+INTENSITY_EQUATION+" must contain a comparator sign (\">\" or \"<\")! Error at line number "+lineNumber+"!");
    String leftSide = equation.substring(0,comparatorIdx);
    String rightSide = equation.substring(comparatorIdx+1);
    String comparator = equation.substring(comparatorIdx,comparatorIdx+1);
    if (comparator.equalsIgnoreCase(">")) biggerThan = true;
    else biggerThan = false;
    if (leftSide.indexOf(">")!=-1 || leftSide.indexOf("<")!=-1 || leftSide.indexOf("=")!=-1)
      throw new RulesException("The value of "+INTENSITY_EQUATION+" must contain only one comparator sign (\">\" or \"<\")! The entry at line number "+lineNumber+" contains more than one!");
    if (rightSide.indexOf(">")!=-1 || rightSide.indexOf("<")!=-1 || rightSide.indexOf("=")!=-1)
      throw new RulesException("The value of "+INTENSITY_EQUATION+" must contain only one comparator sign (\">\" or \"<\")! The entry at line number "+lineNumber+" contains more than one!"); 
    ExpressionForComparisonVO leftExpression = parseIntensityFragment(leftSide,lineNumber,currentSection,comparator,amountOfChains, headFragments, chainFragments, missed);
    ExpressionForComparisonVO rightExpression = parseIntensityFragment(rightSide,lineNumber,currentSection,comparator,amountOfChains, headFragments, chainFragments, missed);
    ExpressionForComparisonVO biggerExpression = rightExpression;
    ExpressionForComparisonVO smallerExpression = leftExpression;
    if (biggerThan){
      biggerExpression = leftExpression;
      smallerExpression = rightExpression;
    } 
    IntensityRuleVO ruleVO = new IntensityRuleVO(currentSection,originalEquation,biggerExpression,smallerExpression);
    return ruleVO;
  }
  
  /**
   * evaluates the intensity expression on one side of the ">" or "<"
   * @param originalValue the intensity expression from the file
   * @param lineNumber the line number of the expression - for error detection
   * @param currentSection the section in the file - GENERAL_SECTION, HEAD_SECTION, CHAINS_SECTION, or POSITION_SECTION
   * @param comp ">" or "<"
   * @param amountOfChains only for POSITION_SECTION the maximum number of chains has to be defined
   * @param headFragments the possible head fragment names
   * @param chainFragments the possible chain fragment names
   * @param missed fragments that were not found in the equation (required for reading of results, since for rules containing "+", not all of the fragments have to be found)
   * @return an object containing all information on one side of the comparator (can be used to calculate the total value)
   * @throws RulesException if the rule is not possible (e.g. if a fragment has not been defined before)
   */
  private static ExpressionForComparisonVO parseIntensityFragment(String originalValue, int lineNumber, int currentSection, String comp, Integer amountOfChains, Hashtable<String,String> headFragments, Hashtable<String,String> chainFragments, Hashtable<String,String> missed) throws RulesException{
    String globalMultiplier = "1";
    String value = new String(originalValue);
    Vector<String> lengthSortedFragmentNames = FragmentRuleVO.getLengthSortedFragmentNames(headFragments.keySet(),chainFragments.keySet(),missed.keySet());
    //if there is a bracket, we can assume that there is a change in the global multiplier
    if (originalValue.indexOf("(")!=-1 || originalValue.indexOf(")")!=-1){
      if (originalValue.indexOf("(")==-1) throw new RulesException("An \""+INTENSITY_EQUATION+"\" containing an opening bracket must contain a closing bracket! The equation \""+value+"\" at line "+lineNumber+" does not!");
      if (originalValue.indexOf(")")==-1) throw new RulesException("An \""+INTENSITY_EQUATION+"\" containing a closing bracket must contain an opening bracket! The equation \""+value+"\" at line "+lineNumber+" does not!");
      if (originalValue.indexOf("(")>originalValue.indexOf(")")) throw new RulesException("An \""+INTENSITY_EQUATION+"\" must not start with a closing bracket before an opening bracket! The equation \""+value+"\" at line "+lineNumber+" does!");
      boolean isMathematicalBracket = true;
      if (Settings.useAlex()){
        for (String fragment : lengthSortedFragmentNames){
          if (fragment.indexOf("(")==-1 || originalValue.indexOf(fragment)==-1) continue;
          if (originalValue.indexOf("(")!=originalValue.indexOf(fragment)+fragment.indexOf("(")) continue;
          isMathematicalBracket=false;
        }
      }
      if (isMathematicalBracket){
        value = originalValue.substring(originalValue.indexOf("(")+1,originalValue.lastIndexOf(")")).trim();
        String prevFragment = originalValue.substring(0,originalValue.indexOf("(")).trim();
        String pastFragment = originalValue.substring(originalValue.lastIndexOf(")")+1).trim();
        String[] results = extractMultiplicationFactor(prevFragment, pastFragment, value, originalValue, lineNumber, currentSection, comp,amountOfChains);
        globalMultiplier = results[0];
      }
    }
    
    Vector<FragmentMultVO> fragments = new Vector<FragmentMultVO>();
    int position = -1;
    
    while (value.length()>0){
      Object[] fragmentAndValue = extractMultiplicativeExpression(value, lengthSortedFragmentNames, originalValue, lineNumber, currentSection, comp, amountOfChains);
      if (fragmentAndValue[0]==null){
        if (fragments.size()==0)
          throw new RulesException("An \""+INTENSITY_EQUATION+"\" must contain one previously declared fragments at each side of the \""+comp+"\"! The value \""+originalValue+"\" at line "+lineNumber+" has not!");
        else
          throw new RulesException("The \""+INTENSITY_EQUATION+"\" contains a part that cannot be interpreted: \""+value+"\"! The problem is the equation \""+originalValue+"\" at line "+lineNumber+"!");
      }
      FragmentMultVO fragMult = (FragmentMultVO)fragmentAndValue[0];
      if (!fragMult.getFragmentName().equalsIgnoreCase(IntensityRuleVO.BASEPEAK_NAME)){
        if (fragments.size()==0 || (fragments.size()==1 && fragments.get(0).getFragmentName().equalsIgnoreCase(IntensityRuleVO.BASEPEAK_NAME))) position = fragMult.getPosition();
        else {
          if (position!=fragMult.getPosition()) throw new RulesException("The \""+INTENSITY_EQUATION+"\" containing position rules must have the same position on one side of the equation! The equation \""+originalValue+"\" at line "+lineNumber+" causes the error!");
        }
      }
      fragments.add((FragmentMultVO)fragmentAndValue[0]);
      value = (String)fragmentAndValue[1];
    }
    ExpressionForComparisonVO exprVO = new ExpressionForComparisonVO(fragments,globalMultiplier);
    return exprVO;
  }
  
  private static Object[] extractMultiplicativeExpression(String inputValue, Vector<String> lengthSortedFragmentNames, String originalValue, int lineNumber, int currentSection, String comp,Integer amountOfChains) throws RulesException{
    Object[] multiplicativeExpressionAndRemainingValue = new Object[2];
    Object[] fragAndStart = getFragmentAndStartPosition(inputValue,lengthSortedFragmentNames);
    if ((Integer)fragAndStart[1]<0){
      multiplicativeExpressionAndRemainingValue[0] = null;
      multiplicativeExpressionAndRemainingValue[1] = null;
      return multiplicativeExpressionAndRemainingValue;
    }
    String name = (String)fragAndStart[0];
    int startPos = (Integer)fragAndStart[1];
    boolean positive = true;
    String beforeString = inputValue.substring(0,startPos).trim();
    char[] charsBefore = beforeString.toCharArray();
    int stopChar = 0;
    for (int i=(charsBefore.length-1);i!=-1;i--){
      char ch = charsBefore[i];
      if (ch=='+' || ch=='-'){
        if (ch=='-') positive=false;
        stopChar = i;
        break;
      }
      if (!(Character.isDigit(ch) || ch=='.' || ch=='*' || ch=='/' || ch==' ' || ch=='[' || ch==']')) break;
      stopChar = i;
    }
    String value = beforeString.substring(0,stopChar).trim();
    String prevFragment = beforeString.substring(stopChar).trim();
    if (prevFragment.startsWith("+")||prevFragment.startsWith("-")) prevFragment = prevFragment.substring(1);  

    stopChar = 0;
    String afterString = inputValue.substring(startPos+name.length()).trim();
    char[] charsAfter = afterString.toCharArray();
    for (int i=0; i!=charsAfter.length; i++){
      char ch = charsAfter[i];
      if (!(Character.isDigit(ch) || ch=='.' || ch=='*' || ch=='/' || ch==' ' || ch=='[' || ch==']')) break;
      stopChar = i+1;
    }
    
    value += afterString.substring(stopChar).trim();
    String pastFragment = afterString.substring(0,stopChar).trim();
    String[] results = extractMultiplicationFactor(prevFragment, pastFragment, value, originalValue, lineNumber, currentSection, comp, amountOfChains);
    String multString = results[0];
    String position = results[1];
    FragmentMultVO multVO = new FragmentMultVO(name, multString, positive, Integer.parseInt(position));
    multiplicativeExpressionAndRemainingValue[0] = multVO;
    multiplicativeExpressionAndRemainingValue[1] = value;
    return multiplicativeExpressionAndRemainingValue;
  }
  
  private static String[] extractMultiplicationFactor(String prevFragment, String pastFragment, String fragmentName, String originalValue, int lineNumber, int currentSection, String comp,Integer amountOfChains) throws RulesException{   
    String multString = "1";
    double mult = 1d;
    int position = 0;
    //now the multiplication factor and the position is extracted
    //this is for multiplication factors before the fragment name
    if (prevFragment.indexOf("*")!=-1 && (pastFragment.indexOf("*")!=-1 || pastFragment.indexOf("/")!=-1))
      throw new RulesException("An \""+INTENSITY_EQUATION+"\" must not contain more than multiplier (\"*\" or \"/\") at each side of the \""+comp+"\"! The value \""+originalValue+"\" at "+lineNumber+" has more than one!");
    if (prevFragment.length()>0){
      if (!prevFragment.endsWith("*")) throw new RulesException("An \""+INTENSITY_EQUATION+"\" allows only a multiplier sign (\"*\") and a number before the declared fragment (\"+name+\")! The expression \""+originalValue+"\" at "+lineNumber+" does not comply this rule!");
      try{
        mult = Double.parseDouble(prevFragment.substring(0,prevFragment.length()-1));
        multString = prevFragment.substring(0,prevFragment.length()-1);
      }catch(NumberFormatException nfx){
        throw new RulesException("An \""+INTENSITY_EQUATION+"\" allows only a multiplier sign (\"*\") and a number before the declared fragment (\"+name+\")! The expression \""+originalValue+"\" at "+lineNumber+" does not comply this rule!");
      }
    }
    //this asks if there is a position identifier
    if (pastFragment.length()>0 && pastFragment.startsWith("[") && pastFragment.indexOf("]")!=-1){
      if (currentSection!=POSITION_SECTION) throw new RulesException("Position specific identifiers in \""+INTENSITY_EQUATION+"\" are allowed in the "+POSITION_SECTION_NAME+" section only! The expression \""+originalValue+"\" at "+lineNumber+" has a position (\"[]\")!");
      try{
        position = Integer.parseInt(pastFragment.substring(1,pastFragment.indexOf("]")));
      }catch(NumberFormatException nfx){
        throw new RulesException("An \""+INTENSITY_EQUATION+"\" must have an integer value for the position! The expression \""+originalValue+"\" at "+lineNumber+" is not integer format!");
      }
      if (position<1 || (amountOfChains!=null && position>amountOfChains)){
        String throwStatement = "The position specific value \""+position+"\" is not allowed!";
        if (amountOfChains!=null && position>amountOfChains)
          throwStatement += " Only values between 1 and "+amountOfChains+" are allowed!";
        throwStatement += " Error at line number "+lineNumber+"!";
        throw new RulesException(throwStatement);
      }
        
      pastFragment = pastFragment.substring(pastFragment.indexOf("]")+1);
    }
    //this is for the multiplication factor after the fragment name
    if (pastFragment.length()>0){
      if (!(pastFragment.startsWith("*")||pastFragment.startsWith("/"))) throw new RulesException("An \""+INTENSITY_EQUATION+"\" allows only a multiplier or divisor sign (\"*\" or \"/\") and a number after the declared fragment ("+fragmentName+")! The expression \""+originalValue+"\" at "+lineNumber+" does not comply this rule!");
      try{
        mult = Double.parseDouble(pastFragment.substring(1,pastFragment.length()));
        multString = pastFragment.substring(1,pastFragment.length());
        if (pastFragment.startsWith("/")){
          mult = 1d/mult;
          multString = String.valueOf(mult);
        }
      }catch(NumberFormatException nfx){
        throw new RulesException("An \""+INTENSITY_EQUATION+"\" allows only a multiplier or divisor sign (\"*\" or \"/\") and a number after the declared fragment ("+fragmentName+")! The expression \""+originalValue+"\" at "+lineNumber+" does not comply this rule!");
      }
    }
    String[] result = new String[2];
    result[0] = multString;
    result[1] = String.valueOf(position);
    return result;
  }
  
  /**
   * returns the fragment name and the index in the expression where it starts
   * @param value the expression
   * @param lengthSortedFragmentNames all available fragment names, sorted by name length in descending order
   * @return Object array of 2 - [0] = name,String ; [1] = position in expression,Integer
   */
  private static Object[] getFragmentAndStartPosition(String value, Vector<String> lengthSortedFragmentNames){
    String fragment = null;
    int start = -1;
    if (value.indexOf(IntensityRuleVO.BASEPEAK_NAME)!=-1){
      fragment = IntensityRuleVO.BASEPEAK_NAME;
      start = value.indexOf(IntensityRuleVO.BASEPEAK_NAME);
    } else {
      for (String frag : lengthSortedFragmentNames){
        if (value.indexOf(frag)!=-1){
          fragment = frag;
          start = (value.indexOf(frag));
          break;
        }
      }
    }
    Object[] result = new Object[2];
    result[0] = fragment;
    result[1] = start;
    return result;
  }
  
  /**
   * checks if the $CHAIN/$ALKYLCHAIN/$ALEKNYLCHAIN is possible in this context
   * @param ruleVO the parsed FragmentRuleVO
   * @param lineNumber the line number of the expression - for error detection
   * @throws RulesException thrown if an invalid chain is used
   */
  private void checkSelectedChainValid(FragmentRuleVO ruleVO, int lineNumber) throws RulesException{
    int amountOfChains = Integer.valueOf(getAmountOfChains());
    int amountOfAlkylChains = Integer.valueOf(getAmountOfAlkylChains());
    int amountOfAlkenylChains = Integer.valueOf(getAmountOfAlkenylChains());
    int amountOfAcylChains = amountOfChains-amountOfAlkylChains-amountOfAlkenylChains;
    String errorMessage = "The chain type ";
    if (ruleVO.getChainType()==FragmentRuleVO.ACYL_CHAIN && amountOfAcylChains<1){
      errorMessage += FragmentRuleVO.CHAIN_NAME+" is not allowed according to the general settings! Only ";
      if (amountOfAlkylChains>0) errorMessage += FragmentRuleVO.ALKYL_CHAIN_NAME;
      if (amountOfAlkylChains>0 && amountOfAlkenylChains>0) errorMessage += "/";
      if (amountOfAlkenylChains>0) errorMessage += FragmentRuleVO.ALKENYL_CHAIN_NAME;
      errorMessage += " is allowed";
      throw new RulesException(errorMessage);
    } else if (ruleVO.getChainType()==FragmentRuleVO.ALKYL_CHAIN && amountOfAlkylChains<1){
      errorMessage += FragmentRuleVO.ALKYL_CHAIN_NAME+" is not allowed according to the general settings! Only ";
      if (amountOfAcylChains>0) errorMessage += FragmentRuleVO.CHAIN_NAME;
      if (amountOfAcylChains>0 && amountOfAlkenylChains>0) errorMessage += "/";
      if (amountOfAlkenylChains>0) errorMessage += FragmentRuleVO.ALKENYL_CHAIN_NAME;
      errorMessage += " is allowed";
      throw new RulesException(errorMessage);
    } else if (ruleVO.getChainType()==FragmentRuleVO.ALKENYL_CHAIN && amountOfAlkenylChains<1){
      errorMessage += FragmentRuleVO.ALKENYL_CHAIN_NAME+" is not allowed according to the general settings! Only ";
      if (amountOfAcylChains>0) errorMessage += FragmentRuleVO.CHAIN_NAME;
      if (amountOfAcylChains>0 && amountOfAlkylChains>0) errorMessage += "/";
      if (amountOfAlkylChains>0) errorMessage += FragmentRuleVO.ALKYL_CHAIN_NAME;
      errorMessage += " is allowed";
      throw new RulesException(errorMessage);
    }
  }
  
  /**
   * retrieves name of chain library- parseFile() has to be called before
   * @return file name of the Excel file containing the fatty acid chains
   */
  public String getChainLibrary(){
    return this.generalSettings_.get(GENERAL_MS2_LIB);
  }

  /**
   * Java regular expression to extract the number C atoms from the analyte name - parseFile() has to be called before 
   * @return regular expression to extract the number C atoms from the analyte name
   */
  public String getCAtomsFromNamePattern(){
    return this.generalSettings_.get(GENERAL_CATOMS_PARSE);
  }

  /**
   * Java regular expression to extract the number of double bonds from the analyte name - parseFile() has to be called before
   * @return regular expression to extract the number of double bonds from the analyte name
   */
  public String getDoubleBondsFromNamePattern(){
    return this.generalSettings_.get(GENERAL_DBOND_PARSE);
  }
  
  /**
   * how many chains does this analyte class have - parseFile() has to be called before
   * @return the amount of chains for this class
   */
  public String getAmountOfChains(){
    return this.generalSettings_.get(GENERAL_CHAINS);
  }

  /**
   * how many alkyl chains does this analyte class have - parseFile() has to be called before
   * @return the amount of chains for this class
   */
  public String getAmountOfAlkylChains(){
    int alkylChains = 0;
    if (generalSettings_.containsKey(GENERAL_CHAINS_ALKYL)) alkylChains = Integer.parseInt(generalSettings_.get(GENERAL_CHAINS_ALKYL));
    return String.valueOf(alkylChains);
  }

  /**
   * how many alkyl chains does this analyte class have - parseFile() has to be called before
   * @return the amount of chains for this class
   */
  public String getAmountOfAlkenylChains(){
    int alkenylChains = 0;
    if (generalSettings_.containsKey(GENERAL_CHAINS_ALKENYL)) alkenylChains = Integer.parseInt(generalSettings_.get(GENERAL_CHAINS_ALKENYL));
    return String.valueOf(alkenylChains);
  }
  
  
  /**
   * cutoff value relative to the base peak - parseFile() has to be called before
   * @return the cutoff value for this class
   * throws RulesException thrown if the base peak cutoff field is invalid - should never be the case in this method, since this is checked in the parsing
   */
  public double getBasePeakCutoff() throws RulesException{
    double cutoff = 0d;
    if (this.generalSettings_.containsKey(GENERAL_CUTOFF)) cutoff = readPercentPermilleValue(getBasePeakCutoffAsString(), GENERAL_CUTOFF, -1);
    return cutoff;
  }

  /**
   * cutoff value relative to the base peak - parseFile() has to be called before
   * @return the cutoff value for this class
   */
  public String getBasePeakCutoffAsString(){
    String cutoff = "0";
    if (this.generalSettings_.containsKey(GENERAL_CUTOFF)) cutoff = this.generalSettings_.get(GENERAL_CUTOFF);
    return cutoff;
  }
  
  /**
   * cutoff value relative to the highest available chain combination - parseFile() has to be called before
   * @return the cutoff value relative to the highest available chain combination
   * @throws RulesException thrown if the chain cutoff value is invalid - should never be the case in this method, since this is checked in the parsing
   */
  public double getChainCutoff() throws RulesException {
    double cutoff = -1d;
    if (this.generalSettings_.containsKey(GENERAL_CHAIN_CUTOFF)) cutoff = readPercentPermilleValue(getChainCutoffAsString(), GENERAL_CHAIN_CUTOFF, -1);
    return cutoff;
  }

  /**
   * cutoff value (in original entered format)cutoff value relative to the highest available chain combination - parseFile() has to be called before
   * @return the cutoff value relative to the highest available chain combination
   */
  public String getChainCutoffAsString(){
    String cutoff = "-1";
    if (this.generalSettings_.containsKey(GENERAL_CHAIN_CUTOFF)) cutoff = this.generalSettings_.get(GENERAL_CHAIN_CUTOFF);
    return cutoff;
  }
  
  /**
   * cutoff value relative to the total intensity of all m/z values of a series of spectra - parseFile() has to be called before
   * @return the spectrum coverage cutoff value for this class
   * @throws RulesException thrown if the spectrum coverage is invalid - should never be the case in this method, since this is checked in the parsing
   */
  public double getSpectrumCoverageMin() throws RulesException{
    double cutoff = 0d;
    if (this.generalSettings_.containsKey(GENERAL_SPECTRUM_COVERAGE)){
      String cutoffString = getSpectrumCoverageMinAsString();
      if (cutoffString!=null) cutoff = readPercentPermilleValue(cutoffString, GENERAL_SPECTRUM_COVERAGE, -1);
    }
    return cutoff;
  }
  
  /**
   * cutoff value (in original entered format) relative to the total intensity of all m/z values of a series of spectra - parseFile() has to be called before
   * @return the originally entered spectrum coverage cutoff value for this class
   */
  public String getSpectrumCoverageMinAsString(){
    String cutoff = null;
    if (this.generalSettings_.containsKey(GENERAL_SPECTRUM_COVERAGE)) cutoff = this.generalSettings_.get(GENERAL_SPECTRUM_COVERAGE);
    return cutoff;
  }

  
  /**
   * shall for this class a post processing by retention time be executed
   * @return true if retention time post processing is desired
   */
  public boolean isRtPostprocessing(){
    boolean rtPostprocessing = false;
    if (this.generalSettings_.containsKey(GENERAL_RT_PROCESSING ) && this.generalSettings_.get(GENERAL_RT_PROCESSING ).equalsIgnoreCase("true")) rtPostprocessing = true;
    return rtPostprocessing;
  }
  
  public boolean correctRtForParallelModel(){
    boolean correctRtParallelModel = false;
    if (this.generalSettings_.containsKey(GENERAL_RT_PARALLEL_SERIES) && this.generalSettings_.get(GENERAL_RT_PARALLEL_SERIES).equalsIgnoreCase("true")) correctRtParallelModel = true;
    return correctRtParallelModel;
  }
  
  public String getRetentionTimeMaxDeviation(){
    String rtMaxDev = null;
    if (this.generalSettings_.containsKey(GENERAL_RT_DEV_MAX)) rtMaxDev = generalSettings_.get(GENERAL_RT_DEV_MAX);
    return rtMaxDev;
  }
  
  public boolean isSingleChainIdentification(){
    boolean singleChain = false;
    if (generalSettings_.containsKey(GENERAL_SINGLE_CHAIN) && generalSettings_.get(GENERAL_SINGLE_CHAIN).equalsIgnoreCase("true")) singleChain = true;
    return singleChain;
  }
  
  public boolean isUnionWithoutPosition(){
    boolean noPosition = false;
    if (generalSettings_.containsKey(GENERAL_PEAK_UNION_NO_POSITION) && generalSettings_.get(GENERAL_PEAK_UNION_NO_POSITION).equalsIgnoreCase("true")) noPosition = true;
    return noPosition;
  }
  
  /**
   * @return the time for uniting peaks sharing the same evidence
   */
  public String getPeakUnionTime(){
    String time = null;
    if (this.generalSettings_.containsKey(GENERAL_PEAK_UNION_TIME)) time = generalSettings_.get(GENERAL_PEAK_UNION_TIME);
    return time;
  }
  
  /**
   * 
   * @return a non-null-value if a class specific cutoff for the MS1 peak detection threshold is set
   */
  public String getMS1PeakCutoff(){
    String cutoff = null;
    if (this.generalSettings_.containsKey(GENERAL_MS1_PEAK_CUTOFF)) cutoff = generalSettings_.get(GENERAL_MS1_PEAK_CUTOFF);
    return cutoff;
  }
   
  /**
   * 
   * @return a non-null-value if a class specific cutoff for isobaric MS1 peaks is set
   */
  public String getIsobarExclusionRatio(){
    String cutoff = null;
    if (this.generalSettings_.containsKey(GENERAL_ISOBAR_RATIO)) cutoff = generalSettings_.get(GENERAL_ISOBAR_RATIO);
    return cutoff;
  }

  /**
   * 
   * @return a non-null-value if a class specific cutoff for isobaric MS1 peaks that are separated a certain time distance from a unique peak
   */
  public String getIsobarFarExclusionRatio(){
    String cutoff = null;
    if (this.generalSettings_.containsKey(GENERAL_ISOBAR_FAR_RATIO)) cutoff = generalSettings_.get(GENERAL_ISOBAR_FAR_RATIO);
    return cutoff;
  }


  /**
   * 
   * @return a non-null-value if a class specific retention time has been defined to use the far exclusion ratio
   */
  public String getIsobarFarRtDifference(){
    String cutoff = null;
    if (this.generalSettings_.containsKey(GENERAL_ISOBAR_RT)) cutoff = generalSettings_.get(GENERAL_ISOBAR_RT);
    return cutoff;
  }
  
  /**
   * 
   * @return in which order the identification by MS should be made (ORDER_MS1_FIRST/ORDER_MSN_FIRST/ORDER_MSN_ONLY)
   */
  public int getMSIdentificationOrder(){
    int msIdentificationOrder = RulesContainer.ORDER_MS1_FIRST;
    if (generalSettings_.containsKey(GENERAL_IDENTIFICATION_ORDER)) msIdentificationOrder = Integer.parseInt(generalSettings_.get(GENERAL_IDENTIFICATION_ORDER));
    return msIdentificationOrder;
  }

  
  /**
   * head fragment rules - key is the rule name of the fragment - parseFile() has to be called before
   * @return head fragment rules
   */
  public Hashtable<String,FragmentRuleVO> getHeadFragmentRules(){
    return this.headFragments_;
  }
  
  /**
   * 
   * @return lowest and highest spectrum level required for MS2 rules
   */
  public int[] getSpectrumLevelRange(){
    int[] range = new int[2];
    int lowestLevel = Integer.MAX_VALUE;
    int highestLevel = 0;
    for (FragmentRuleVO ruleVO : this.headFragments_.values()){
      if (ruleVO.getMsLevel()<lowestLevel) lowestLevel = ruleVO.getMsLevel();
      if (ruleVO.getMsLevel()>highestLevel) highestLevel = ruleVO.getMsLevel();
    }
    for (FragmentRuleVO ruleVO :this.chainFragments_.values()){
      if (ruleVO.getMsLevel()<lowestLevel) lowestLevel = ruleVO.getMsLevel();
      if (ruleVO.getMsLevel()>highestLevel) highestLevel = ruleVO.getMsLevel();      
    }
    range[0] = lowestLevel;
    range[1] = highestLevel;
    return range;
  }
  
  /**
   * how many additional chain positions are possible
   * @return the amount of additional 'empty' chain positions
   */
  public Integer getAddChainPositions(){
    int addChains = 0;
    if (generalSettings_.containsKey(GENERAL_ADD_POSITIONS)) addChains = Integer.parseInt(generalSettings_.get(GENERAL_ADD_POSITIONS));
    return addChains;
  }

  /**
   * the number of possible positions
   * @return the number of possible positions
   */
  public int getAllowedChainPositions(){
    return Integer.parseInt(this.getAmountOfChains())+getAddChainPositions();
  }

  
  /**
   * head rules for intensity comparisons - parseFile() has to be called before
   * @return head rules for intensity comparisons
   */
  public Vector<IntensityRuleVO> getHeadIntensityRules(){
    return this.headIntensities_;
  }
  
  /**
   * chain fragment rules - key is the name of the fragment - parseFile() has to be called before
   * @return chain fragment rules - key is the name of the fragment
   */
  public Hashtable<String,FragmentRuleVO> getChainFragmentRules(){
    return this.chainFragments_;
  }
  
  /**
   * chain rules for intensity comparisons - parseFile() has to be called before
   * @return chain rules for intensity comparisons
   */
  public Vector<IntensityRuleVO> getChainIntensityRules(){
    return this.chainIntensities_;
  }
  
  /**
   * position rules for intensity comparisons - parseFile() has to be called before
   * @return position rules for intensity comparisons
   */
  public Vector<IntensityRuleVO> getPositionIntensityRules(){
    return this.positionIntensities_;
  }
  
  /**
   * writes a complete fragmentation ruleset, by providing the necessary parameters
   * @param dir the directory where the rules shall be written
   * @param lipidClass the lipid class the rules affect
   * @param lipidAdduct the adduct for this class
   * @param generalSettings VO containing the entries for the [GENERAL] section
   * @param headFragments sorted Vector of fragmentation rules for the [HEAD] section
   * @param headIntensityRules sorted Vector of intensity rules for the [HEAD] section
   * @param chainFragments sorted Vector of fragmentation rules for the [CHAINS] section
   * @param chainIntensityRules sorted Vector of intensity rules for the [CHAINS] section
   * @param positionIntensityRules sorted Vector of intensity rules for the [POSITON] section
   * @throws IOException general exception if there is something wrong with the file or the file path
   * @throws RulesException exception that is thrown if there is something wrong with the provided VOs
   */
  public static void writeRules(String dir, String lipidClass, String lipidAdduct, GeneralSettingsVO generalSettings,
      Vector<FragmentRuleVO> headFragments, Vector<IntensityRuleVO> headIntensityRules,
      Vector<FragmentRuleVO> chainFragments, Vector<IntensityRuleVO> chainIntensityRules,
      Vector<IntensityRuleVO> positionIntensityRules) throws IOException, RulesException
  {
    File file;
    String filename = StaticUtils.getRuleFileName(lipidClass, lipidAdduct);
    if (dir == "default")
    {
      file = new File(RulesContainer.currentRulesDir_+"/"+filename);   
    } else {
      file = new File(dir + "/"+filename); 
    }
    if (!file.exists()){
      file.createNewFile();
    } 
    FileWriter fw = new FileWriter(file.getAbsoluteFile());
    BufferedWriter bw = new BufferedWriter(fw);

    bw.write(GENERAL_SECTION_NAME+"\n");        
    bw.write(GENERAL_CHAINS+"=");
    if (generalSettings.getAmountOfChains()!=null) bw.write(generalSettings.getAmountOfChains().toString());             
    bw.write("\n"+GENERAL_MS2_LIB+"=");
    if (generalSettings.getChainLibrary()!=null) bw.write(generalSettings.getChainLibrary());        
    bw.write("\n"+GENERAL_CATOMS_PARSE+"=");
    if (generalSettings.getCarbonAtomsRule()!=null) bw.write(generalSettings.getCarbonAtomsRule());        
    bw.write("\n"+GENERAL_DBOND_PARSE+"=");
    if (generalSettings.getDoubleBondsRule()!=null) bw.write(generalSettings.getDoubleBondsRule());
    if (generalSettings.getAmountOfAlkylChains()!=null && generalSettings.getAmountOfAlkylChains()>0)
      bw.write("\n"+GENERAL_CHAINS_ALKYL+"="+String.valueOf(generalSettings.getAmountOfAlkylChains()));
    if (generalSettings.getAmountOfAlkenylChains()!=null && generalSettings.getAmountOfAlkenylChains()>0)
      bw.write("\n"+GENERAL_CHAINS_ALKENYL+"="+String.valueOf(generalSettings.getAmountOfAlkenylChains()));  
    bw.write("\n"+GENERAL_CUTOFF+"=");
    if (generalSettings.getBasePeakCutoff()!=null) bw.write(generalSettings.getBasePeakCutoff()); 
    if (generalSettings.getChainCutoff()!=null)
      bw.write("\n"+GENERAL_CHAIN_CUTOFF+"="+String.valueOf(generalSettings.getChainCutoff()));    
    if (generalSettings.getSpectrumCoverage()!=null)
      bw.write("\n"+GENERAL_SPECTRUM_COVERAGE+"="+generalSettings.getSpectrumCoverage());
    bw.write("\n");            
    if(generalSettings.isRtPostProcessing())
      bw.write(GENERAL_RT_PROCESSING+"=" + String.valueOf(generalSettings.isRtPostProcessing()) + "\n");
    if(generalSettings.isRtParallelSeries())
      bw.write(GENERAL_RT_PARALLEL_SERIES+"=" + String.valueOf(generalSettings.isRtParallelSeries()) + "\n");
    if(generalSettings.getRtMaxDeviation()!=null)
      bw.write(GENERAL_RT_DEV_MAX+"=" + String.valueOf(generalSettings.getRtMaxDeviation()) + "\n");
    if(generalSettings.isAllowSingleChain())
      bw.write(GENERAL_SINGLE_CHAIN+"=" + String.valueOf(generalSettings.isAllowSingleChain()) + "\n");
    if (generalSettings.getMsIdentificationOrder()!=RulesContainer.ORDER_MS1_FIRST){
      bw.write(GENERAL_IDENTIFICATION_ORDER+"=");
      if (generalSettings.getMsIdentificationOrder()==RulesContainer.ORDER_MSN_FIRST)
        bw.write(ORDER_MSN_FIRST+"\n");
      else if (generalSettings.getMsIdentificationOrder()==RulesContainer.ORDER_MSN_ONLY)
        bw.write(ORDER_MSN_ONLY+"\n");
      else{
        bw.close();
        throw new RulesException("The value "+generalSettings.getMsIdentificationOrder()+" is not allowed for writing "+GENERAL_IDENTIFICATION_ORDER);
      }  
    }
    if (generalSettings.getAddChainPositions()!=null && generalSettings.getAddChainPositions()>0){
      bw.write(GENERAL_ADD_POSITIONS+"=" + String.valueOf(generalSettings.getAddChainPositions()) + "\n");
    }
    
    writeSection(bw, HEAD_SECTION_NAME, headFragments, headIntensityRules);
    writeSection(bw, CHAINS_SECTION_NAME, chainFragments, chainIntensityRules); 
    if(chainFragments.size() != 0 && positionIntensityRules!=null && positionIntensityRules.size()>0){
        writeSection(bw, POSITION_SECTION_NAME, null, positionIntensityRules);
    }
    bw.close();
  }
  
  /**
   * writes a [HEAD], [CHAINS] or [POSITION] section with the provided fragment- and intensity ruels
   * @param bw the buffered writer
   * @param section the section type ([HEAD], [CHAINS] or [POSITION])
   * @param fragments a sorted Vector of fragmentation rules
   * @param intRules a sorted Vector of intensity rules
   * @throws IOException general exception if there is something wrong with the file or the file path
   */
  private static void writeSection(BufferedWriter bw, String section, Vector<FragmentRuleVO> fragments, Vector<IntensityRuleVO> intRules) throws IOException{
    boolean writeIntensities = false;
    boolean fragmentsWritten = false;
    if (section.equalsIgnoreCase(POSITION_SECTION_NAME)){
      bw.write("\n");
      bw.write(section+"\n");
      writeIntensities = true;
    }
    if(fragments!=null && fragments.size() != 0){
      writeIntensities = true;
      bw.write("\n");
      bw.write(section+"\n");
      bw.write(FRAGMENT_SUBSECTION_NAME+"\n");        

      for(int i=0; i<fragments.size(); i++ ){
        writeFragmentRule(bw,fragments.get(i));
      }
      fragmentsWritten = true;
    }
    if(writeIntensities && intRules!=null && intRules.size() != 0)
    {
      if (fragmentsWritten) bw.write("\n");
      bw.write(INTENSITY_SUBSECTION_NAME+"\n");
        
      int headIntensityRulesSize = intRules.size();
      if(headIntensityRulesSize != 0){
        for(int i=0; i < headIntensityRulesSize; i++){
          writeIntensityRule(bw, intRules.get(i));

        }
      }
    }

  }
  
  /**
   * writes the line of a fragment rule
   * @param bw the buffered writer
   * @param rule the rule VO
   * @throws IOException general exception if there is something wrong with the file or the file path
   */
  private static void writeFragmentRule(BufferedWriter bw, FragmentRuleVO rule) throws IOException{
    String mandatory = "false";                
    if(rule.isMandatory()) mandatory = "true";
    bw.write(FRAGMENT_NAME+"=" + rule.getName() + "\t" + FRAGMENT_FORMULA+"=" + 
      rule.getFormula() + "\t" + FRAGMENT_CHARGE+"=" + 
      Integer.toString(rule.getCharge()) + "\t"+ FRAGMENT_LEVEL+"=" + 
      Integer.toString(rule.getMsLevel()) + "\t" + "mandatory=" + mandatory + "\n");     
 
  }
  
  /**
   * writes the line of an intensity rule
   * @param bw the buffered writer
   * @param rule the intensity rule VO
   * @throws IOException general exception if there is something wrong with the file or the file path
   */
  private static void writeIntensityRule(BufferedWriter bw, IntensityRuleVO rule) throws IOException{
    String mandatory = "false";
    if(rule.isMandatory()) mandatory = "true";
    bw.write("Equation=" + rule.getRuleIdentifier() + "\t" + "mandatory=" + mandatory + "\n");            

  }
  
}
