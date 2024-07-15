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

package at.tugraz.genome.lda.msn;

import java.io.File;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Set;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.exception.SheetNotPresentException;
import at.tugraz.genome.lda.msn.parser.FALibParser;
import at.tugraz.genome.lda.msn.parser.LCBLibParser;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;

/**
 * Class caches the chain objects
 * @author Juergen Hartler
 *
 */
public class FattyAcidsContainer
{
  /** default directory containing the included fatty acid chain libraries*/
  private final static String DEFAULT_ACIDS_DIR = "fattyAcids";
  /** file suffix new Excel*/
  public final static String FA_FILE_SUFFIX_NEW = ".xlsx";
  /** file suffix old Excel*/
  public final static String FA_FILE_SUFFIX_OLD = ".xls";
  /** prefix that is only present for temporary excel files */
  private final static String EXCEL_TEMPORARY_PREFIX = "~$";
  
  /** the cache instance*/
  private static FattyAcidsContainer instance_;
  
  /** directory containing the fatty acid chain libraries */
  private static String faDir_ = DEFAULT_ACIDS_DIR;
  
  /** hash containing the information
   * first key: Excel library name
   * second key: encoded number of hydroxylation sites
   * third key: amount of carbon atoms of chain
   * fourth key: amount of double bonds of key
   * fifth key: prefix
   */
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>> fattyAcids_;
  
  /**
   * the available isotopic labels
   */
  private Hashtable<String,Set<String>> availableLabels_;
  

  /** hash containing the information
   * first key: Excel library name
   * second key: encoded number of hydroxylation sites
   * third key: amount of carbon atoms of chain
   * fourth key: amount of double bonds of key
   * fifth key: prefix
   */  
  private Hashtable<String,Hashtable<String,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>> lcbs_;

  
  /**
   * constructor defines directory where the Excel libraries are stored
   * @param faDir directory where the Excel libraries are stored
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException general exception if there is something wrong about the file
   */
  private FattyAcidsContainer(String faDir) throws RulesException, IOException {
    faDir_ = faDir;
    fattyAcids_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>>();
    lcbs_ = new Hashtable<String,Hashtable<String,Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>>>>();
    availableLabels_ = new Hashtable<String,Set<String>>();
    extractChains();
  }

  /**
   * caches an instance of the fatty acid container
   * @param faDir directory where the Excel libraries are stored
   * @return instance of the fatty acid container
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException general exception if there is something wrong about the file
   */
  private static FattyAcidsContainer getInstance(String faDir) throws RulesException, IOException {
    if (instance_==null || !faDir_.equalsIgnoreCase(faDir)){
      instance_ = new FattyAcidsContainer(faDir);
    }
    return instance_;
  }

  /**
   * queries if the required FA library exists 
   * @param faLibName name of the FA library
   * @return true if the library is there
   */
  private boolean hasFALib(String faLibName){
    return (fattyAcids_.containsKey(faLibName)||lcbs_.containsKey(faLibName));
  }
  
  /**
   * queries if the required FA library exists and checks if the cache is initialized, before
   * @param faLib name of the FA library
   * @param faLibDir directory where the Excel libraries are stored
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException general exception if there is something wrong about the file
   * @throws NoRuleException thrown if the library is not there
   */
  private static void checkIfFALibExists(String faLib, String faLibDir) throws RulesException, IOException, NoRuleException {
    String faDir = faLibDir;
    if (faDir==null || faDir.length()==0) faDir = DEFAULT_ACIDS_DIR;
    getInstance(faLibDir);
    if (instance_.hasFALib(faLib)) return;
    instance_ = new FattyAcidsContainer(faLibDir);
    if (!instance_.hasFALib(faLib)) throw new NoRuleException("There is no fatty acid lib called \""+faLib+"\"!");
  }

  /**
   * returns the extracted fatty acid chains for the requested library
   * @param faLib name of the FA library
   * @return the extracted fatty acid chains for the requested library; first key: oh encoded; second key: #C-atoms; third key #double bonds; fourth key: prefix 
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<String, Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>>> getAllFattyAcidChains(String faLibName) throws RulesException, NoRuleException, IOException {
    return getAllFattyAcidChains(faLibName, faDir_);
  }

  
  /**
   * returns the extracted fatty acid chains for the requested library
   * @param faLibName name of the FA library
   * @param faLibDir directory where the Excel libraries are stored
   * @return the extracted fatty acid chains for the requested library; first key: oh encoded; second key: #C-atoms; third key #double bonds; fourth key: prefix 
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<String, Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>>> getAllFattyAcidChains(String faLibName, String faLibDir) throws RulesException, NoRuleException, IOException {
    checkIfFALibExists(faLibName,faLibDir);
    return instance_.fattyAcids_.get(faLibName);
  }
  
  /**
   * parses the Excel files and writes the results in the cached hash table
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   */
  private void extractChains() throws RulesException, IOException {
    File faLibDir = new File(faDir_);
    if (!faLibDir.exists()) throw new RulesException("The provided fatty acid lib directory does not exist!");
    if (!faLibDir.isDirectory()) throw new RulesException("The fatty acid lib directory is a file - not a directory!");
    File[] files = faLibDir.listFiles();
    for (File file : files){
    	if (file.getAbsolutePath().contains(EXCEL_TEMPORARY_PREFIX)) continue;
      if (!file.getAbsolutePath().endsWith(FA_FILE_SUFFIX_NEW) && !file.getAbsolutePath().endsWith(FA_FILE_SUFFIX_OLD)) continue;
      try {
        try {
          FALibParser parser = new FALibParser(file);
          parser.parseFile();
          fattyAcids_.put(file.getName().substring(0,file.getName().length()), parser.getResult());
          availableLabels_.put(file.getName().substring(0,file.getName().length()),parser.getAvailableLabels());
        } catch (SheetNotPresentException ex){
          LCBLibParser parser = new LCBLibParser(file);
          parser.parseFile();
          lcbs_.put(file.getName().substring(0,file.getName().length()), parser.getResult());
          availableLabels_.put(file.getName().substring(0,file.getName().length()),parser.getAvailableLabels());
        }
      } catch (RulesException ex){
        throw new RulesException(file.getName()+": "+ex.getMessage());
      }
    }
  }
  
  /**
   * returns the extracted fatty acid chains for the requested library
   * @param faLibName name of the FA library
   * @param encoded the encoded character for this number of hydroxylations
   * @return the result hash - first integer is the amount of carbon atoms; second integer is the amount of double bonds; third key is the prefix
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> getFattyAcidChains(String faLibName, String encoded)
      throws HydroxylationEncodingException, RulesException, NoRuleException, IOException {
    return getFattyAcidChains(faLibName, faDir_, encoded);
  }
  
  
  /**
   * returns the extracted fatty acid chains for the requested library
   * @param faLibName name of the FA library
   * @param faLibDir directory where the Excel libraries are stored
   * @param encoded the encoded character for this number of hydroxylations
   * @return the result hash - first integer is the amount of carbon atoms; second integer is the amount of double bonds; third key is the prefix
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> getFattyAcidChains(String faLibName, String faLibDir, String encoded)
      throws HydroxylationEncodingException, RulesException, NoRuleException, IOException {
    checkIfFALibExists(faLibName, faLibDir);
    if (!instance_.fattyAcids_.get(faLibName).containsKey(encoded))
      throw new HydroxylationEncodingException("The demanded hydroxylation encoding \""+encoded+"\" is not available in your FA-library: "+faLibName+"!");
    return instance_.fattyAcids_.get(faLibName).get(encoded);
  }


  /**
   * returns the available long chain base for the requested library
   * @param faLibName name of the FA library
   * @param hydroxyNr the number of hydroxylation sites
   * @return the result hash - first integer is the amount of carbon atoms; second integer is the amount of double bonds; third key is the prefix
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> getFattyAcidChains(String faLibName, short hydroxyNr)
      throws HydroxylationEncodingException, RulesException, NoRuleException, IOException {
    return getFattyAcidChains(faLibName, faDir_, hydroxyNr);
  }
  
  /**
   * returns the extracted fatty acid chains for the requested library
   * @param faLibName name of the FA library
   * @param faLibDir directory where the Excel libraries are stored
   * @param hydroxyNr the number of hydroxylation sites
   * @return the result hash - first integer is the amount of carbon atoms; second integer is the amount of double bonds; third key is the prefix
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> getFattyAcidChains(String faLibName, String faLibDir, short hydroxyNr)
      throws HydroxylationEncodingException, RulesException, NoRuleException, IOException {
    String encoding = Settings.getFaHydroxyEncoding().getEncodedPrefix(hydroxyNr);
    try {
      return getFattyAcidChains(faLibName, faLibDir, encoding);
    }catch(HydroxylationEncodingException hdx) {
      throw new HydroxylationEncodingException("The demanded hydroxylation number \""+hydroxyNr+"\" is not available in your FA-library!");
    }
  }

  
  /**
   * returns the available long chain base for the requested library
   * @param lcbLibName name of the LCB library
   * @return the extracted LCB chains for the requested library; first key: oh encoded; second key: #C-atoms; third key #double bonds; fourth key: prefix 
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<String, Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>>> getAllLCBs(String lcbLibName) throws RulesException, NoRuleException, IOException {
    return getAllLCBs(lcbLibName, faDir_);
  }

  
  /**
   * returns the available long chain base for the requested library
   * @param lcbLibName name of the LCB library
   * @param lcbLibDir directory where the Excel libraries are stored
   * @return the extracted fatty acid chains for the requested library; first key: oh encoded; second key: #C-atoms; third key #double bonds; fourth key: prefix 
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<String, Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>>> getAllLCBs(String lcbLibName, String lcbLibDir) throws RulesException, NoRuleException, IOException {
    checkIfFALibExists(lcbLibName,lcbLibDir);
    return instance_.lcbs_.get(lcbLibName);
  }


  
  /**
   * returns the available long chain base for the requested library
   * @param lcbLibName name of the LCB library
   * @param encoded the encoded character for this number of hydroxylations
   * @return the result hash - first integer is the amount of carbon atoms; second integer is the amount of double bonds; third key is the prefix
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> getLCBs(String lcbLibName, String encoded)
      throws HydroxylationEncodingException, RulesException, NoRuleException, IOException {
    return getLCBs(lcbLibName, faDir_, encoded);
  }
  
  
  /**
   * returns the available long chain base for the requested library
   * @param lcbLibName name of the LCB library
   * @param lcbLibDir directory where the Excel libraries are stored
   * @param encoded the encoded character for this number of hydroxylations
   * @return the result hash - first integer is the amount of carbon atoms; second integer is the amount of double bonds; third key is the prefix
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> getLCBs(String lcbLibName, String lcbLibDir, String encoded)
      throws HydroxylationEncodingException, RulesException, NoRuleException, IOException {
    checkIfFALibExists(lcbLibName, lcbLibDir);
    if (!instance_.lcbs_.get(lcbLibName).containsKey(encoded))
      throw new HydroxylationEncodingException("The demanded hydroxylation encoding \""+encoded+"\" is not available in your LCB-library: "+lcbLibName+"!");
    return instance_.lcbs_.get(lcbLibName).get(encoded);
  }


  /**
   * returns the available long chain base for the requested library
   * @param lcbLibName name of the FA library
   * @param lcbLibDir directory where the Excel libraries are stored
   * @param hydroxyNr the number of hydroxylation sites
   * @return the result hash - first integer is the amount of carbon atoms; second integer is the amount of double bonds; third key is the prefix
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> getLCBs(String lcbLibName, short hydroxyNr)
      throws HydroxylationEncodingException, RulesException, NoRuleException, IOException {
    return getLCBs(lcbLibName, faDir_, hydroxyNr);
  }
  
  /**
   * returns the available long chain base for the requested library
   * @param lcbLibName name of the FA library
   * @param lcbLibDir directory where the Excel libraries are stored
   * @param hydroxyNr the number of hydroxylation sites
   * @return the result hash - first integer is the amount of carbon atoms; second integer is the amount of double bonds; third key is the prefix
   * @throws HydroxylationEncodingException thrown if the encoding does not exist
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Hashtable<Integer, Hashtable<Integer, Hashtable<String, Hashtable<String, FattyAcidVO>>>> getLCBs(String lcbLibName, String lcbLibDir, short hydroxyNr)
      throws HydroxylationEncodingException, RulesException, NoRuleException, IOException {
    String encoding = Settings.getLcbHydroxyEncoding().getEncodedPrefix(hydroxyNr);
    try {
      return getLCBs(lcbLibName, lcbLibDir, encoding);
    }catch(HydroxylationEncodingException hdx) {
      throw new HydroxylationEncodingException("The demanded hydroxylation number \""+hydroxyNr+"\" is not available in your LCB-library!");
    }
  }
  
  /**
   * returns the available isotopic labels for the requested library
   * @param libName name of the LCB library
   * @return the extracted isotopic labels for the requested library 
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Set<String> getAvailableLabels(String libName) throws RulesException, NoRuleException, IOException {
    return getAvailableLabels(libName, faDir_);
  }

  
  /**
   * returns the available isotopic labels for the requested library
   * @param libName name of the LCB library
   * @param libDir directory where the Excel libraries are stored
   * @return the extracted isotopic labels for the requested library
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException exception if there is something wrong about the file
   */
  public static Set<String>  getAvailableLabels(String libName, String libDir) throws RulesException, NoRuleException, IOException {
    checkIfFALibExists(libName,libDir);
    return instance_.availableLabels_.get(libName);
  }

}
