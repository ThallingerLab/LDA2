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

package at.tugraz.genome.lda.msn;

import java.io.File;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.utils.RangeInteger;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;

/**
 * This class caches all fragmentation rules
 * @author Juergen Hartler
 *
 */
public class RulesContainer
{
  /** default directory containing the rules */
  public final static String DEFAULT_RULES_DIR = "fragRules";
  private final static String INTERIM_RULES_DIR = DEFAULT_RULES_DIR+"/interim";
  
  public static String currentRulesDir_ = DEFAULT_RULES_DIR;
  
  /** directory where the rules reside in */
  private String rulesDir_;
  
  /** instances of the classes that cache the rules - the hash is the coresponding rules directory */
  private static Hashtable<String,RulesContainer> instances_;
  
  /** the fragmentation rules - key is the name of the analyte class */
  private Hashtable<String,FragRuleParser> rules_;
  
  public final static int ORDER_MS1_FIRST = 0;
  public final static int ORDER_MSN_FIRST = 1;
  public final static int ORDER_MSN_ONLY = 2;
  
  /**
   * constructor defines directory where the rule files are stored - extraction of rules is started immediately
   * @param rulesDir directory where the rule files are stored
   * @throws RulesException specifies in detail which rules are not valid
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  private RulesContainer(String rulesDir) throws RulesException, IOException, SpectrummillParserException {
    rulesDir_ = rulesDir;
    extractRules();
  }
  
  /**
   * caches an instance of the rules container
   * @param rulesDir directory where the rule files are stored
   * @return instance of the rules container
   * @throws RulesException specifies in detail which rules are not valid
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  private static RulesContainer getInstance(String rulesDir) throws RulesException, IOException, SpectrummillParserException {
    if (instances_==null) clearCache();
    if (!instances_.containsKey(rulesDir)){
      instances_.put(rulesDir, new RulesContainer(rulesDir));
    }
    return instances_.get(rulesDir);
  }
  
  /** tells to RulesContainer to use the permanent rules directory for quantitation*/
  public static void usePermanentRulesDir(){
    currentRulesDir_ = DEFAULT_RULES_DIR;
  }

  /** tells to RulesContainer to use the temporary rules directory for quantitation (if "Apply" is used in the settings)*/
  public static void useInterimRulesDir(){
    currentRulesDir_ = INTERIM_RULES_DIR;
  }

  
  /**
   * parses the rules files and stores them in the "rules_" hash - key is the name of the analyte class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws IOException general exception if there is something wrong about the file
   */
  private void extractRules() throws RulesException, IOException {
    rules_ = new Hashtable<String,FragRuleParser>();
    File rulesDir = new File(rulesDir_);
    if (!rulesDir.exists()) throw new RulesException("The provided fragmentation rules directory does not exist!");
    if (!rulesDir.isDirectory()) throw new RulesException("The provided fragmentation rules directory is a file - not a directory!");
    File[] files = rulesDir.listFiles();
    ElementConfigParser elementParser = Settings.getElementParser();
    for (File file : files){
      if (!file.getAbsolutePath().endsWith(StaticUtils.RULE_FILE_SUFFIX)) continue;
      FragRuleParser parser = new FragRuleParser(elementParser);
      try {
        parser.parseFile(file);
        rules_.put(file.getName().substring(0,file.getName().length()-StaticUtils.RULE_FILE_SUFFIX.length()), parser);
      } catch (RulesException ex){
        throw new RulesException(file.getName()+": "+ex.getMessage());
      }
    }
  }
  
  /**
   * 
   * @param rule the name of the lipid class
   * @return true if rules for this lipid class are present
   */
  private boolean hasRule(String rule){
    return rules_.containsKey(rule);
  }
  
  /**
   * verifies if a rule for this class is in the rules directory
   * @param rule name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  private static RulesContainer checkIfRuleExists(String rule, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    String rDir = rulesDir;
    if (rDir==null || rDir.length()==0) rDir = currentRulesDir_;
    RulesContainer instance = getInstance(rDir);
    if (instance.hasRule(rule)) return instance;
    instance = new RulesContainer(rDir);
    if (!instance.hasRule(rule)) throw new NoRuleException("There is no MS2 rule for the analyte class \""+rule+"\"!");
    instances_.put(rDir, instance);
    return instance;
  }

  /**
   * 
   * @param ruleName name of the lipid class
   * @return the name of the fatty acid chain library
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getChainlibrary(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getChainlibrary(ruleName, currentRulesDir_);
  }
  
  /** 
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the name of the fatty acid chain library
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getChainlibrary(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getChainLibrary();
  }
  
  /**
   * 
   * @param ruleName name of the lipid class
   * @return the name of the long chain base library
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getLcbLibrary(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getLcbLibrary(ruleName, currentRulesDir_);
  }
  
  /** 
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the name of the long chain base library
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getLcbLibrary(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getLcbLibrary();
  }

  
  /**
   * Java regular expression to extract the number C atoms from the analyte name
   * @param ruleName name of the lipid class
   * @return regular expression to extract the number C atoms from the analyte name
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getCAtomsFromNamePattern(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getCAtomsFromNamePattern(ruleName, currentRulesDir_);
  }
  
  /** 
   * Java regular expression to extract the number C atoms from the analyte name
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return regular expression to extract the number C atoms from the analyte name
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getCAtomsFromNamePattern(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getCAtomsFromNamePattern();
  }

  /**
   * Java regular expression to extract the number of double bonds from the analyte name
   * @param ruleName name of the lipid class
   * @return regular expression to extract the number of double bonds from the analyte name
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getDoubleBondsFromNamePattern(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getDoubleBondsFromNamePattern(ruleName, currentRulesDir_);
  }
  
  /** 
   * Java regular expression to extract the number of double bonds from the analyte name
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return regular expression to extract the number of double bonds from the analyte name
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getDoubleBondsFromNamePattern(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getDoubleBondsFromNamePattern();
  }
  
  /**
   * how many chains does this analyte class have
   * @param ruleName name of the lipid class
   * @return the amount of chains for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getAmountOfChains(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getAmountOfChains(ruleName, currentRulesDir_);
  }

  /** 
   * how many chains does this analyte class have
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the amount of chains for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getAmountOfChains(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getAmountOfChains();
  }
  
  /**
   * how many alkyl chains does this analyte class have
   * @param ruleName name of the lipid class
   * @return the amount of alkyl chains for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getAmountOfAlkylChains(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getAmountOfAlkylChains(ruleName, currentRulesDir_);
  }

  /** 
   * how many alkyl chains does this analyte class have
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the amount of alkyl chains for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getAmountOfAlkylChains(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getAmountOfAlkylChains();
  }

  /**
   * how many alkenyl chains does this analyte class have
   * @param ruleName name of the lipid class
   * @return the amount of alkenyl chains for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getAmountOfAlkenylChains(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getAmountOfAlkenylChains(ruleName, currentRulesDir_);
  }

  /** 
   * how many alkenyl chains does this analyte class have
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the amount of alkenyl chains for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getAmountOfAlkenylChains(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getAmountOfAlkenylChains();
  }
  
  /**
   * how many LCBs does this analyte class have
   * @param ruleName name of the lipid class
   * @return the amount of LCB chains for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getAmountOfLCBs(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getAmountOfLCBs(ruleName, currentRulesDir_);
  }

  /** 
   * how many LCBs chains does this analyte class have
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the amount of LCB chains for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getAmountOfLCBs(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getAmountOfLCBs();
  }

  /**
   * the range of possible hydroxylation sites for the FA moiety
   * @param ruleName name of the lipid class
   * @return range of possible hydroxylation sites for the FA moiety
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static RangeInteger getFaHydroxyRange(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getFaHydroxyRange(ruleName, currentRulesDir_);
  }

  /** 
   * the range of possible hydroxylation sites for the FA moiety
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return range of possible hydroxylation sites for the FA moiety
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static RangeInteger getFaHydroxyRange(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getFaHydroxyRange();
  }

  /**
   * the range of possible hydroxylation sites for the LCB moiety
   * @param ruleName name of the lipid class
   * @return range of possible hydroxylation sites for the LCB moiety
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static RangeInteger getLcbHydroxyRange(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getLcbHydroxyRange(ruleName, currentRulesDir_);
  }

  /** 
   * the range of possible hydroxylation sites for the LCB moiety
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return range of possible hydroxylation sites for the LCB moiety
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static RangeInteger getLcbHydroxyRange(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getLcbHydroxyRange();
  }
  

  /**
   * cutoff value relative to the base peak
   * @param ruleName name of the lipid class
   * @return the cutoff value for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static double getBasePeakCutoff(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getBasePeakCutoff(ruleName, currentRulesDir_);
  }

  /** 
   * cutoff value relative to the base peak
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the cutoff value for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */

  public static double getBasePeakCutoff(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getBasePeakCutoff();
  }

  /**
   * cutoff value relative to the base peak (originally entered string)
   * @param ruleName name of the lipid class
   * @return the cutoff value for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getBasePeakCutoffAsString(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getBasePeakCutoffAsString(ruleName, currentRulesDir_);
  }

  /** 
   * cutoff value relative to the base peak (originally entered string)
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the cutoff value for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */

  public static String getBasePeakCutoffAsString(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getBasePeakCutoffAsString();
  }


  /**
   * cutoff value relative to the highest chain combination found
   * @param ruleName name of the lipid class
   * @return the cutoff value for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static double getChainCutoff(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getChainCutoff(ruleName, currentRulesDir_);
  }

  /** 
   * cutoff value relative to the highest chain combination found
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the cutoff value for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static double getChainCutoff(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getChainCutoff();
  }
  
  /**
   * cutoff value (originally entered string) relative to the highest chain combination found
   * @param ruleName name of the lipid class
   * @return the cutoff value (originally entered string) for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getChainCutoffAsString(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getChainCutoffAsString(ruleName, currentRulesDir_);
  }

  /** 
   * cutoff value (originally entered string) relative to the highest chain combination found
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the cutoff value (originally entered string) for this class
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */

  public static String getChainCutoffAsString(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    String cutoffString = cont.rules_.get(ruleName).getChainCutoffAsString();
    if (cutoffString.equalsIgnoreCase("-1")) cutoffString = null;
    return cutoffString;
  }


  
  /** 
   * spectrum coverage that has to be fulfilled by the found fragments
   * @param ruleName name of the lipid class
   * @return minimum spectrum coverage that has to be fulfilled
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static double getSpectrumCoverageMin(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getSpectrumCoverageMin(ruleName, currentRulesDir_);
  }

  
  /** 
   * spectrum coverage that has to be fulfilled by the found fragments
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return minimum spectrum coverage that has to be fulfilled
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static double getSpectrumCoverageMin(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getSpectrumCoverageMin();
  }
  

  /** 
   * spectrum coverage that has to be fulfilled by the found fragments (originally entered string)
   * @param ruleName name of the lipid class
   * @return minimum spectrum coverage that has to be fulfilled (originally entered string)
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getSpectrumCoverageMinAsString(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getSpectrumCoverageMinAsString(ruleName, currentRulesDir_);
  }

  
  /** 
   * spectrum coverage that has to be fulfilled by the found fragments (originally entered string)
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return minimum spectrum coverage that has to be fulfilled (originally entered string)
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getSpectrumCoverageMinAsString(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getSpectrumCoverageMinAsString();
  }

  
  /**
   * shall for this class a post processing by retention time be executed
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean isRtPostprocessing(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return isRtPostprocessing(ruleName, currentRulesDir_);
  }

  
  /**
   * shall for this class a post processing by retention time be executed
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean isRtPostprocessing(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).isRtPostprocessing();
  }

  /**
   * shall the retention time processing take a potential parallel model into account
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean correctRtForParallelModel(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return correctRtForParallelModel(ruleName, currentRulesDir_);
  }

  
  /**
   * shall the retention time processing take a potential parallel model into account
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean correctRtForParallelModel(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).correctRtForParallelModel();
  }

  
  /**
   * returns the maximum retention time which a value may deviate from the hypothetical RT model
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getRetentionTimeMaxDeviation(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getRetentionTimeMaxDeviation(ruleName, currentRulesDir_);
  }

  
  /**
   * returns the maximum retention time which a value may deviate from the hypothetical RT model
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getRetentionTimeMaxDeviation(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getRetentionTimeMaxDeviation();
  }

  
  /**
   * returns if the union of peaks of same and similar MS evidence shall be enforced
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getPeakUnionTime(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getPeakUnionTime(ruleName, currentRulesDir_);
  }

  
  /**
   * returns if the union of peaks of same and similar MS evidence shall be enforced
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getPeakUnionTime(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getPeakUnionTime();
  }
  
  
  /**
   * returns if the union of peaks of same and similar MS should neglect positional evidence
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean isUnionWithoutPosition(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return isUnionWithoutPosition(ruleName, currentRulesDir_);
  }

  
  /**
   * returns if the union of peaks of same and similar MS should neglect positional evidence
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean isUnionWithoutPosition(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).isUnionWithoutPosition();
  }

  
  /**
   * returns a non-null-value if a class specific cutoff for the MS1 peak detection threshold is set
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getMS1PeakCutoff(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getMS1PeakCutoff(ruleName, currentRulesDir_);
  }

  
  /**
   * returns a non-null-value if a class specific cutoff for the MS1 peak detection threshold is set
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getMS1PeakCutoff(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getMS1PeakCutoff();
  }
  
  /**
   * returns a non-null-value if a class specific cutoff for isobaric MS1 peaks is set
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return a non-null-value if a class specific cutoff for isobaric MS1 peaks is set
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getIsobarExclusionRatio(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getIsobarExclusionRatio(ruleName, currentRulesDir_);
  }

  
  /**
   * returns a non-null-value if a class specific cutoff for isobaric MS1 peaks is set
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return a non-null-value if a class specific cutoff for isobaric MS1 peaks is set
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getIsobarExclusionRatio(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getIsobarExclusionRatio();
  }
  
  /**
   * returns a non-null-value if a class specific cutoff for isobaric MS1 peaks that are separated a certain time distance from a unique peak
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return a non-null-value if a class specific cutoff for isobaric MS1 peaks that are separated a certain time distance from a unique peak
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getIsobarFarExclusionRatio(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getIsobarFarExclusionRatio(ruleName, currentRulesDir_);
  }

  
  /**
   * returns a non-null-value if a class specific cutoff for isobaric MS1 peaks that are separated a certain time distance from a unique peak
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return a non-null-value if a class specific cutoff for isobaric MS1 peaks that are separated a certain time distance from a unique peak
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getIsobarFarExclusionRatio(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getIsobarFarExclusionRatio();
  }
  
  /**
   * returns a non-null-value if a class specific retention time has been defined to use the far exclusion ratio
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return a non-null-value if a class specific retention time has been defined to use the far exclusion ratio
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getIsobarFarRtDifference(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getIsobarFarRtDifference(ruleName, currentRulesDir_);
  }

  
  /**
   * returns a non-null-value if a class specific retention time has been defined to use the far exclusion ratio
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return a non-null-value if a class specific retention time has been defined to use the far exclusion ratio
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static String getIsobarFarRtDifference(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getIsobarFarRtDifference();
  }

  
  /**
   * head fragment rules - key is the rule name of the fragment
   * @param ruleName name of the lipid class
   * @return head fragment rules - key is the rule name of the fragment
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Hashtable<String,FragmentRuleVO> getHeadFragmentRules(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getHeadFragmentRules(ruleName, currentRulesDir_);
  }
  
  /** 
   * head fragment rules - key is the rule name of the fragment
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return head fragment rules - key is the rule name of the fragment
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Hashtable<String,FragmentRuleVO> getHeadFragmentRules(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getHeadFragmentRules();
  }

  /**
   * head rules for intensity comparisons
   * @param ruleName name of the lipid class
   * @return head rules for intensity comparisons
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Vector<IntensityRuleVO> getHeadIntensityRules(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getHeadIntensityRules(ruleName, currentRulesDir_);
  }
  
  /** 
   * head rules for intensity comparisons
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return head rules for intensity comparisons
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Vector<IntensityRuleVO> getHeadIntensityRules(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getHeadIntensityRules();
  }

  /**
   * chain fragment rules - key is the rule name of the fragment
   * @param ruleName name of the lipid class
   * @return chain fragment rules - key is the rule name of the fragment
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Hashtable<String,FragmentRuleVO> getChainFragmentRules(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getChainFragmentRules(ruleName, currentRulesDir_);
  }
  
  /** 
   * chain fragment rules - key is the rule name of the fragment
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return chain fragment rules - key is the rule name of the fragment
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Hashtable<String,FragmentRuleVO> getChainFragmentRules(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getChainFragmentRules();
  }

  /**
   * chain rules for intensity comparisons
   * @param ruleName name of the lipid class
   * @return chain rules for intensity comparisons
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Vector<IntensityRuleVO> getChainIntensityRules(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getChainIntensityRules(ruleName, currentRulesDir_);
  }
  
  /** 
   * chain rules for intensity comparisons
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return chain rules for intensity comparisons
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Vector<IntensityRuleVO> getChainIntensityRules(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getChainIntensityRules();
  }

  /**
   * position rules for intensity comparisons
   * @param ruleName name of the lipid class
   * @return position rules for intensity comparisons
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */  
  public static Vector<IntensityRuleVO> getPositionIntensityRules(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getPositionIntensityRules(ruleName, currentRulesDir_);
  }
  
  /** 
   * position rules for intensity comparisons
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return position rules for intensity comparisons
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Vector<IntensityRuleVO> getPositionIntensityRules(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getPositionIntensityRules();
  }
  
  /**
   * is an identification which is based on a single chain already valid?
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean isSingleChainIdentification(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return isSingleChainIdentification(ruleName, currentRulesDir_);
  }

  
  /**
   * is an identification which is based on a single chain already valid?
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean isSingleChainIdentification(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).isSingleChainIdentification();
  }
  
  /**
   * in which order the identification by MS should be made (ORDER_MS1_FIRST/ORDER_MSN_FIRST/ORDER_MSN_ONLY)
   * @param ruleName name of the lipid class
   * @return in which order the identification by MS should be made (ORDER_MS1_FIRST/ORDER_MSN_FIRST/ORDER_MSN_ONLY)
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static int getMSIdentificationOrder(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getMSIdentificationOrder(ruleName, currentRulesDir_);
  }

  
  /**
   * in which order the identification by MS should be made (ORDER_MS1_FIRST/ORDER_MSN_FIRST/ORDER_MSN_ONLY)
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return in which order the identification by MS should be made (ORDER_MS1_FIRST/ORDER_MSN_FIRST/ORDER_MSN_ONLY)
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static int getMSIdentificationOrder(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getMSIdentificationOrder();
  }
  
  /**
   * lowest and highest spectrum level required for MS2 rules
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static int[] getSpectrumLevelRange(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getSpectrumLevelRange(ruleName, currentRulesDir_);
  }

  
  /**
   * lowest and highest spectrum level required for MS2 rules
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return shall for this class a post processing by retention time be executed
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static int[] getSpectrumLevelRange(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getSpectrumLevelRange();
  }
  

  /**
   * how many additional chain positions are possible
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the amount of additional 'empty' chain positions
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Integer getAddChainPositions(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getAddChainPositions(ruleName, currentRulesDir_);
  }

  
  /**
   * how many additional chain positions are possible
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the amount of additional 'empty' chain positions
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Integer getAddChainPositions(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getAddChainPositions();
  }
  
  /**
   * the number of possible positions
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the number of possible positions
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static int getAllowedChainPositions(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getAllowedChainPositions(ruleName, currentRulesDir_);
  }

  
  /**
   * the number of possible positions
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the number of possible positions
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static int getAllowedChainPositions(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getAllowedChainPositions();
  }

  
  /**
   * true when other adducts have to be found to allow for this adduct
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return true when other adducts have to be found to allow for this adduct
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean requiresOtherValidAdduct(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return requiresOtherValidAdduct(ruleName, currentRulesDir_);
  }

  
  /**
   * true when other adducts have to be found to allow for this adduct
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return true when other adducts have to be found to allow for this adduct
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean requiresOtherValidAdduct(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).requiresOtherValidAdduct();
  }


  /**
   * returns the list of required adducts to declare this adduct valid
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return the list of required adducts to declare this adduct valid
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Vector<String> getOtherRequiredAdducts(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getOtherRequiredAdducts(ruleName, currentRulesDir_);
  }

  
  /**
   * the list of required adducts to declare this adduct valid
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return true the list of required adducts to declare this adduct valid
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static Vector<String> getOtherRequiredAdducts(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getOtherRequiredAdducts();
  }


  /**
   * returns true when all of the listed other adducts have to be found
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return true when all of the listed other adducts have to be found
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean areAllOtherAdductsRequired(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return areAllOtherAdductsRequired(ruleName, currentRulesDir_);
  }

  
  /**
   * returns true when all of the listed other adducts have to be found
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return true when all of the listed other adducts have to be found
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean areAllOtherAdductsRequired(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).areAllOtherAdductsRequired();
  }
  

  /**
   * returns the time tolerance to detect another adduct
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return time tolerance to detect another adduct
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static float getOtherTimeTolerance(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return getOtherTimeTolerance(ruleName, currentRulesDir_);
  }

  
  /**
   * returns the time tolerance to detect another adduct
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return time tolerance to detect another adduct
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static float getOtherTimeTolerance(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).getOtherTimeTolerance();
  }
  
  
  /**
   * returns true when other overlapping species have to be removed when other adducts are found
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return true when other overlapping species have to be removed when other adducts are found
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean forceOtherAdductValidity(String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    return forceOtherAdductValidity(ruleName, currentRulesDir_);
  }

  
  /**
   * returns true when other overlapping species have to be removed when other adducts are found
   * @param ruleName name of the lipid class
   * @param rulesDir directory where the rule files are stored
   * @return true when other overlapping species have to be removed when other adducts are found
   * @throws RulesException specifies in detail which rules are not valid
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public static boolean forceOtherAdductValidity(String ruleName, String rulesDir) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    RulesContainer cont = checkIfRuleExists(ruleName,rulesDir);
    return cont.rules_.get(ruleName).forceOtherAdductValidity();
  }


  
  /**
   * removes all stored fragmentation rules
   */
  public static void clearCache(){
    instances_ = new Hashtable<String,RulesContainer>();
  }

  public static void clearCache(String rulesDir){
    if (instances_!=null && instances_.containsKey(rulesDir))
      instances_.remove(rulesDir);
  }

  
}
