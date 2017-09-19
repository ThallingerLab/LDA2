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

package at.tugraz.genome.lda;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.LineNumberReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.Vector;

import javax.swing.UIManager;

import at.tugraz.genome.lda.exception.SettingsException;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.swing.RuleDefinitionInterface;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class Settings
{
  private static Settings instance_ = null;

  public final static String VERSION = "2.5.2_6";
  
  public final static String SETTINGS_FILE = ".settings";
  
  private static String operatingSystem_;
  private static String ldaUserHomePath_;
  private static String readWPath_;
  private static String msConvertPath_;
  private static String massWolfPath_;
  private static String massPlusPlusPath_;
  private static String elementConfigPath_;
  private static boolean overviewExcelWorkbook_;
  private static boolean overviewExcelMass_;
  private static boolean overviewExcelIsotope_;
  /** for Waters QTOF, several traces for MS/MS are generated, each one leads to a separate mzXML
   *  if true the mzXML files are merged to one file - this setting should be false for MSE*/
  private static boolean mergeMultipleMSMSFiles_;
  /** only applicable if "MergeMultipleMSMSFiles" is true;
   * the last trace is a quality trace for detection stability - this should not be merged into the data*/
  private static boolean skipLastMzXML_;

  private static String isDefaultInput_;
  private static String esDefaultInput_;
  
  /** the latest stored fragmentation selection 1*/
  private static String fragmentationSelection1_;
  /** the latest stored fragmentation selection 2*/
  private static String fragmentationSelection2_;
  
  /** +/- time tolerance [min] for 3D display of MSn spectra*/
  private static float msn3DDisplayTolerance_;
  
  public final static String PROPERTIES_DIRECTORY = "properties";
  public final static String PROPERTYFILE_PREFIX = "LipidDataAnalyzer_";
  public final static String PROPERTYFILE_SUFFIX = "properties";
  
  public final static String FRAG_SELECTION_NONE = "none";
  public final static String FRAG_SELECTION_NO_INTENSITY = "noIntensity";
  
  public final static String FRAG_SETTINGS_FILE = ".selected";
  
  private LinkedHashMap<String,File> propertiesFiles_;
  
  private static Settings getInstance() {
    if (instance_ == null) {
      instance_ = new Settings();
      operatingSystem_ = System.getProperty("os.name");
      ldaUserHomePath_=System.getProperty("user.home")+File.separator+".lda";
      readSettingsFile();
      readFragSelectedFile();
      instance_.propertiesFiles_ = instance_.getPropertiesFiles();
    }
    return instance_;
  }
  
  private static void readSettingsFile(){
    try{
      File file = new File(SETTINGS_FILE);
      FileInputStream inNew = new FileInputStream(file);
      Properties properties = new Properties();
      properties.load(inNew);
      inNew.close();

      readWPath_ = properties.getProperty("ReadWPath", null);
      msConvertPath_ = properties.getProperty("MsconvertPath", null);
      massWolfPath_ = properties.getProperty("MassWolfPath", null);
      massPlusPlusPath_ = properties.getProperty("MassPPPath", null);
      elementConfigPath_ = properties.getProperty("ElementConfig", null);
      String overviewExcelString = properties.getProperty("OverviewExcelWorkbook", null);
      if (overviewExcelString!=null&&(overviewExcelString.equalsIgnoreCase("true")||overviewExcelString.equalsIgnoreCase("yes"))){
        overviewExcelWorkbook_ = true;
      }
      String overviewExcelMassString = properties.getProperty("OverviewExcelMass", null);
      if (overviewExcelMassString!=null&&(overviewExcelMassString.equalsIgnoreCase("true")||overviewExcelMassString.equalsIgnoreCase("yes"))){
        overviewExcelMass_ = true;
      }
      String overviewExcelIsotopeString = properties.getProperty("OverviewIsotopicValues", null);
      if (overviewExcelIsotopeString!=null&&(overviewExcelIsotopeString.equalsIgnoreCase("true")||overviewExcelIsotopeString.equalsIgnoreCase("yes"))){
        overviewExcelIsotope_ = true;
      }
      mergeMultipleMSMSFiles_ = false;
      String mergeMultipleMSMSFilesString = properties.getProperty("MergeMultipleMSMSFiles", null);
      if (mergeMultipleMSMSFilesString!=null&&(mergeMultipleMSMSFilesString.equalsIgnoreCase("true")||mergeMultipleMSMSFilesString.equalsIgnoreCase("yes"))){
        mergeMultipleMSMSFiles_ = true;
      }
      skipLastMzXML_ = false;
      String skipLastMzXMLString = properties.getProperty("OverviewIsotopicValues", null);
      if (skipLastMzXMLString!=null&&(skipLastMzXMLString.equalsIgnoreCase("true")||skipLastMzXMLString.equalsIgnoreCase("yes"))){
        skipLastMzXML_ = true;
      }
      
      isDefaultInput_ = properties.getProperty("ISDefaultInput", "");
      esDefaultInput_ = properties.getProperty("ESDefaultInput", "Ex-IS");
      msn3DDisplayTolerance_ = 60f;
      String msn3DDisplayToleranceString = properties.getProperty("MSn3DDisplayTimeWindow","1");
      if (msn3DDisplayToleranceString!=null){
        try{ msn3DDisplayTolerance_ = Float.parseFloat(msn3DDisplayToleranceString)*60f;}catch(NumberFormatException ex){}
      }
    }catch(Exception e){
      e.printStackTrace();
    }
  }
  
  /**
   * 
   * @return the path to the file containing the selected fragmentation settings
   */
  private static String getSettingsFilePath(){
    return RulesContainer.currentRulesDir_+"/"+FRAG_SETTINGS_FILE;
  }
  
  /**
   * reads the selected fragmentation settings and keeps it in the cache
   */
  private static void readFragSelectedFile(){
    fragmentationSelection1_ = FRAG_SELECTION_NONE;
    fragmentationSelection2_ = FRAG_SELECTION_NONE;
    File file = new File(getSettingsFilePath());
    if (!file.exists() || !file.isFile()) return;
    try {
      LineNumberReader reader = new LineNumberReader(new FileReader(file));
      String line;
      boolean firstSet = false;
      while ((line = reader.readLine()) != null) {
        if (line==null || line.trim().length()==0) continue;
        if (firstSet){
          fragmentationSelection2_ = line.trim();
          break;
        } else {
          fragmentationSelection1_ = line.trim();
          firstSet = true;
        }
      }
      reader.close();
    }
    catch (IOException e) {
      e.printStackTrace();
    }
  }
  
  public static void rereadSettings(){
    Settings.getInstance();
    readSettingsFile();
    readFragSelectedFile();
  }
  
  public static String getOperatingSystem(){
    Settings.getInstance();
    return operatingSystem_;
  }

  public static String getReadWPath(){
    Settings.getInstance();
    return Settings.readWPath_;
  }
  
  public static String getMsConvertPath(){
    Settings.getInstance();
    return Settings.msConvertPath_;
  }
  
  public static String getMassWolfPath(){
    Settings.getInstance();
    return Settings.massWolfPath_;
  }
  
  public static String getMassPlusPlusPath(){
    Settings.getInstance();
    return Settings.massPlusPlusPath_;
  }
  
  public static String getElementConfigPath(){
    Settings.getInstance();
    return Settings.elementConfigPath_;
  }
  
  public static boolean isOverviewInExcelDesired(){
    Settings.getInstance();
    return Settings.overviewExcelWorkbook_;
  }
  
  public static boolean isMassInOverviewExcelDesired(){
    Settings.getInstance();
    return Settings.overviewExcelMass_;
  }

  public static boolean isIsotopeInOverviewExcelDesired(){
    Settings.getInstance();
    return Settings.overviewExcelIsotope_;
  }


  /**
   * 
   * @return true the mzXML files are merged to one file - this setting should be false for MSE
   */
  public static boolean mergeMultipleMSMSFiles(){
    Settings.getInstance();
    return Settings.mergeMultipleMSMSFiles_;
  }

  /**
   * 
   * @return true if the last trace of the Waters files should be neglected
   */
  public static boolean skipLastMzXML(){
    Settings.getInstance();
    return Settings.skipLastMzXML_;
  }

  public static String getInternalStandardDefaultInput()
  {
    Settings.getInstance();
    return isDefaultInput_;
  }
  
  public static String getExternalStandardDefaultInput()
  {
    Settings.getInstance();
    return esDefaultInput_;
  }
  
  public static String getLdaUserHomePath(){
    Settings.getInstance();
    return ldaUserHomePath_;
  }
  
  public static String getLicensePath(){
    Settings.getInstance();
    return Settings.getLdaUserHomePath()+File.separator+"License";
  }
  
  public static boolean isWindows(){
    if (Settings.getOperatingSystem().startsWith("Windows"))
      return true;
    else return false;
  }
  
  public static boolean isOSMac(){
    if (Settings.getOperatingSystem().startsWith("Mac OS"))
      return true;
    else return false;
  }
  
  public static boolean isOSMacAndJavaLookAndFeel(){
    Settings.getInstance();
    if (Settings.getOperatingSystem().startsWith("Mac OS")&&UIManager.getLookAndFeel().getName().equalsIgnoreCase("Metal"))
      return true;
    else return false;
  }
  
  /**
   * 
   * @return the stored selected fragmentation settings of the first selection
   */
  public static String getFragmentationSelection1(){
    Settings.getInstance();
    return fragmentationSelection1_;
  }

  /**
   * 
   * @return the stored selected fragmentation settings of the second selection
   */
  public static String getFragmentationSelection2(){
    Settings.getInstance();
    return fragmentationSelection2_;
  }

  /**
   * 
   * @return a String showing the selected fragmentation rules - for display in the title
   */
  public static String getFragmentSettingsString(){
    Settings.getInstance();
    String settingsString = "";
    if (fragmentationSelection1_!=null && fragmentationSelection1_.length()>0 && !fragmentationSelection1_.equalsIgnoreCase(FRAG_SELECTION_NONE)){
      settingsString += fragmentationSelection1_;
      if (!fragmentationSelection1_.equalsIgnoreCase(FRAG_SELECTION_NO_INTENSITY) && fragmentationSelection2_!=null && fragmentationSelection2_.length()>0 && !fragmentationSelection2_.equalsIgnoreCase(FRAG_SELECTION_NONE))
        settingsString += " "+fragmentationSelection2_;
    }
    return settingsString;
  }
  
  /**
   * 
   * @return +/- time tolerance [min] for 3D display of MSn spectra
   */
  public static float getMsn3DDisplayTolerance()
  {
    Settings.getInstance();
    return msn3DDisplayTolerance_;
  }

  /**
   * applies the fragmentation rule settings without overwriting the default ones
   * @param machine the selected MS machine settings
   * @throws SettingsException if there is something wrong
   */
  public static void applySettings(String machine) throws SettingsException{
    if (instance_.propertiesFiles_.containsKey(machine)){
      LipidomicsConstants.switchToOtherConfFile(instance_.propertiesFiles_.get(machine));
    } else {
      throw new SettingsException("There are no settings for the machine \""+machine+"\"!");
    }
  }

  public static void saveMachineSettings(String machine) throws SettingsException{
    if (instance_.propertiesFiles_.containsKey(machine)){
      try{
        LipidomicsConstants.switchToOtherDefaultConfFile(instance_.propertiesFiles_.get(machine));
      } catch (IOException iox){throw new SettingsException(iox);}
    } else {
      throw new SettingsException("There are no settings for the machine \""+machine+"\"!");
    }
  }

  
  public static Set<String> getPropertyFileNames(){
    if (instance_ == null) Settings.getInstance();
    return instance_.propertiesFiles_.keySet();
  }
  
  /**
   * 
   * @param currentMSMachine
   * @return the possible fragmentation options for a specified MS machine
   */
  public static Vector<String> getFragmentationSettings(String currentMSMachine){
    if (instance_ == null) Settings.getInstance();
    Vector<String> fragSettings = new Vector<String>();
    fragSettings.add(FRAG_SELECTION_NONE);
    fragSettings.add(FRAG_SELECTION_NO_INTENSITY);
    List<String> fragOptions = getFragmentationOptions(currentMSMachine);
    for (String option : fragOptions) fragSettings.add(option);
    return fragSettings;
  }
  
  /**
   * 
   * @param currentMSMachine
   * @return the possible fragmentation options for a specified MS machine
   */
  private static List<String> getFragmentationOptions(String currentMSMachine){
    List<String> fragOptions = new ArrayList<String>();
    File fragDir = new File(RulesContainer.DEFAULT_RULES_DIR);
    File[] files = fragDir.listFiles();
    for (File file : files){
      if (file.isDirectory()&&file.getName().equalsIgnoreCase(currentMSMachine)){
        for (File optionDir : file.listFiles()){
          if (optionDir.isDirectory()) fragOptions.add(optionDir.getName());
        }
      }
    }
    Collections.sort(fragOptions);
    return fragOptions;
  }
  
  private LinkedHashMap<String,File> getPropertiesFiles(){
    File propDir = new File(Settings.PROPERTIES_DIRECTORY);
    File[] files = propDir.listFiles();
    LinkedHashMap<String,File> propFiles = new  LinkedHashMap<String,File>();
    for (int i=0; i!=files.length;i++){
      String fileName = StaticUtils.extractFileName(files[i].getAbsolutePath());
      if (fileName.length()>(Settings.PROPERTYFILE_PREFIX.length()+Settings.PROPERTYFILE_SUFFIX.length()+1) &&
          fileName.startsWith(Settings.PROPERTYFILE_PREFIX) && fileName.endsWith("."+Settings.PROPERTYFILE_SUFFIX)){
        String machineName = fileName.substring(Settings.PROPERTYFILE_PREFIX.length(),fileName.indexOf("."+Settings.PROPERTYFILE_SUFFIX));
        propFiles.put(machineName, files[i]);
      }  
    }
    return propFiles;
  }
  
  /**
   * applies the selected fragmentation settings for quantitation without overwriting the default settings
   * @param settings1 the first selected settings
   * @param settings2 the second selected settings
   * @throws SettingsException if there is something wrong
   */
  public static void applyFragmentationSettings(String settings1, String settings2) throws SettingsException {
    RulesContainer.useInterimRulesDir();
    copyFragmentationSettings(settings1,settings2);
    readFragSelectedFile();
  }

  /**
   * stores the selected fragmentation settings as default
   * @param settings1 the first selected settings
   * @param settings2 the second selected settings
   * @throws SettingsException if there is something wrong
   */
  public static void saveFragmentationSettings(String settings1, String settings2) throws SettingsException {
    RulesContainer.usePermanentRulesDir();
    copyFragmentationSettings(settings1,settings2);
    readFragSelectedFile();
  }
  
  /**
   * copies permanently stored fragmentation rules to an active directory for quantification
   * @param settings1 the first selected settings
   * @param settings2 the second selected settings
   * @throws SettingsException if there is something wrong
   */
  private static void copyFragmentationSettings(String settings1, String settings2) throws SettingsException {
    String selectedFileContents = "";
    if (settings1==null || settings1.length()==0) return;
    removeOldRules();
    if (!settings1.equalsIgnoreCase(FRAG_SELECTION_NONE)){
      selectedFileContents += settings1;
    }
    boolean settings2Selected = false;
    if (settings2!=null && settings2.length()>0 && !settings2.equalsIgnoreCase(FRAG_SELECTION_NONE))
      settings2Selected = true;
    try{
      if (!settings1.equalsIgnoreCase(FRAG_SELECTION_NONE) && !settings1.equalsIgnoreCase(FRAG_SELECTION_NO_INTENSITY) && settings2Selected){
        selectedFileContents += "\n"+settings2;
        copyFragSettingsFiles(getFragSettingsDirectory(settings2), RulesContainer.currentRulesDir_);
      }
      if (!settings1.equalsIgnoreCase(FRAG_SELECTION_NONE)){
        copyFragSettingsFiles(getFragSettingsDirectory(settings1), RulesContainer.currentRulesDir_);
        OutputStream out = new BufferedOutputStream(new FileOutputStream(getSettingsFilePath()));
        out.write(selectedFileContents.getBytes());
        out.close();
      }
    } catch (IOException iox){
      throw new SettingsException(iox);
    }
  }
  
  /**
   * 
   * @param fragSettings
   * @return the directory name of a selected fragmentation setting value
   */
  private static String getFragSettingsDirectory(String fragSettings){
    String dir = RulesContainer.DEFAULT_RULES_DIR+"/";
    if (!fragSettings.equalsIgnoreCase(FRAG_SELECTION_NO_INTENSITY)) dir += LipidomicsConstants.getCurrentMSMachine()+"/";
    dir += fragSettings;
    return dir;
  }
  
  /**
   * copies the fragmentation rules from the permanent option directory to the active directory
   * @param fromPath path containing the fragmentation rules
   * @param toPath active directory for the fragmentation rules
   * @throws SettingsException thrown if there is something wrong
   * @throws IOException thrown if an I/O error occurs
   */
  private static void copyFragSettingsFiles(String fromPath, String toPath) throws SettingsException, IOException{
    File fromDir = new File(fromPath);
    File toDir = new File (toPath);
    if (!fromDir.exists())  throw new SettingsException("The path \""+fromPath+"\" does not exist!");
    if (!toDir.exists())  throw new SettingsException("The path \""+toPath+"\" does not exist!");
    if (!fromDir.isDirectory()) throw new SettingsException("The path \""+fromPath+"\" does not point to a directory!");
    if (!toDir.isDirectory()) throw new SettingsException("The path \""+toPath+"\" does not point to a directory!");
    
    for (File file: fromDir.listFiles()){
      if (!file.isFile() || !file.getName().endsWith(StaticUtils.RULE_FILE_SUFFIX)) continue;
      String outPath = getActiveRulesFilename(file.getName());
      int chunkSize = 1024;
      InputStream in = new BufferedInputStream(new FileInputStream(file));
      OutputStream out = new BufferedOutputStream(new FileOutputStream(outPath));
      byte[] buffer = new byte[chunkSize];
      int len;
      while ((len = in.read(buffer)) > 0) {
        out.write(buffer, 0, len);
      }
      in.close();
      out.close();
    }
    RulesContainer.clearCache();
    RuleDefinitionInterface.clearCacheDir();
    RulesContainer.clearCache(RuleDefinitionInterface.CACHE_DIR);
  }
  
  /**
   * 
   * @param fileName the name of the file to be stored
   * @return the path name of the file to be stored
   */
  private static String getActiveRulesFilename(String fileName){
    return RulesContainer.currentRulesDir_+"/"+fileName;
  }
  
  /**
   * removes files in the active rules directory
   */
  private static void removeOldRules(){
    File rulesDir = new File(RulesContainer.currentRulesDir_);
    for (File file : rulesDir.listFiles()){
      if (file.isFile() && (file.getName().equalsIgnoreCase(FRAG_SETTINGS_FILE) || file.getName().endsWith(StaticUtils.RULE_FILE_SUFFIX)))
        file.delete();
    }
  }
}