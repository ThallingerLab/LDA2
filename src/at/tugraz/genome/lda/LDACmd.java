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

package at.tugraz.genome.lda;
import java.io.File;
import java.util.Hashtable;
import java.util.Timer;
import java.util.TimerTask;
import java.util.Vector;
import java.util.logging.Logger;

import javax.swing.JLabel;
import javax.swing.JProgressBar;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import at.tugraz.genome.lda.swing.BatchQuantificationTable;
import at.tugraz.genome.lda.swing.BatchQuantificationTableModel;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.RawQuantificationPairVO;


/**
 * 
 * @author Juergen Hartler
 *
 */
public class LDACmd
{
  
  /** for logging error messages*/
  private static Logger log_ = Logger.getLogger(LDACmd.class.getName());
  
  /** the quantification thread*/
  private BatchQuantThread batchQuantThread_ = null;
  
  /** timer for checking whether the thread is finished*/
  private Timer timer_;

  /**
   * class for command line interface
   */
  public LDACmd() {
  }

  /**
   * command line interfaces using the arguments as shown in printUsage - method checks input parameters for validity
   * @param args the input arguments from the command line interface
   */
  public static void main(String[] args)
  {
    Options parameters = new Options();
    Option option=new Option("i","input", true, "directory containing raw data");
    option.setRequired(true);
    parameters.addOption(option);
    option=new Option("q","quant", true, "directory containing quant files");
    option.setRequired(true);
    parameters.addOption(option);
    option=new Option("pc","processorsChrom", true, "number of processors for translation");
    option.setRequired(true);
    parameters.addOption(option);
    option=new Option("p","processors", true, "number of processors for quantification");
    option.setRequired(true);
    parameters.addOption(option);
    option=new Option("c","cutoff", true, "intensity cutoff in per mille");
    option.setRequired(false);
    parameters.addOption(option);
    option=new Option("iso1","isotopes1", true, "isotopes that must match");
    option.setRequired(false);
    parameters.addOption(option);
    option=new Option("iso2","isotopes2", true, "isotopes that shall be quantified");
    option.setRequired(false);
    parameters.addOption(option);
    option=new Option("rtTolBef","rtToleranceBefore", true, "RT before tolerance");
    option.setRequired(false);
    parameters.addOption(option);
    option=new Option("rtTolAft","rtToleranceAfter", true, "RT after tolerance");
    option.setRequired(false);
    parameters.addOption(option);
    option=new Option("rtSh","rtShift", true, "retention time shift");
    option.setRequired(false);
    parameters.addOption(option);
    option=new Option("u","unknown", true, "search unknown retention time");
    option.setRequired(false);
    parameters.addOption(option);
    option=new Option("v","version", false, "version of Lipid Data Analyzer");
    option.setRequired(false);
    parameters.addOption(option);    
    
    CommandLineParser cmdParser = new DefaultParser();
    CommandLine command;
    try{
      command = cmdParser.parse(parameters, args);
    }catch (ParseException pe)
    { 
        pe.printStackTrace();
        printUsage(parameters,pe.getMessage()); 
        return; 
    }
    //TODO: I do not think that version will work - since it will throw an error because of the required values
    if (command.hasOption("v")){
      System.out.println("Version: "+Settings.VERSION);
      System.exit(0);
    }
    
    String rawDirString = command.getOptionValue("i");
    String quantDirString = command.getOptionValue("q");
    int nrChromProcessors  = 1;
    try{
      if (command.getOptionValue("pc")!=null) nrChromProcessors = Integer.parseInt(command.getOptionValue("pc"));
      if (nrChromProcessors<1) {
        log_.severe("Number of Processors match must be greater 0");
        System.exit(1);
      }
    }catch (NumberFormatException nfx){
      log_.severe("Number of Processors is not integer format");
      System.exit(1);
    }
    
    int nrProcessors  = 1;
    try{
      if (command.getOptionValue("p")!=null) nrProcessors = Integer.parseInt(command.getOptionValue("p"));
      if (nrProcessors<1) {
        log_.severe("Number of Processors match must be greater 0");
        System.exit(1);
      }
    }catch (NumberFormatException nfx){
      log_.severe("Number of Processors is not integer format");
      System.exit(1);
    }

    float cutoff = Float.parseFloat(LipidomicsConstants.getBasePeakDefaultCutoff());
    LipidomicsConstants.getInstance().setRelativeMS1BasePeakCutoff(LipidomicsConstants.getBasePeakDefaultCutoff());
    try{
      if (command.getOptionValue("c")!=null){
        cutoff = Float.parseFloat(command.getOptionValue("c"));
        if (cutoff<0f) {
          log_.severe("Base peak cutoff value must not be negative");
          System.exit(1);
        }
        if (cutoff>1000f) {
          log_.severe("Base peak cutoff value must not be greater than 1000 per mille");
          System.exit(1);
        }
        LipidomicsConstants.getInstance().setRelativeMS1BasePeakCutoff(command.getOptionValue("c"));
      }
    }catch (NumberFormatException nfx){
      log_.severe("Base peak cutoff value is not float format");
      System.exit(1);
    }
    int isotopesMustMatch = 1;
    int amountOfIsotopes = 2;
    try{
      if (command.getOptionValue("iso1")!=null) isotopesMustMatch = Integer.parseInt(command.getOptionValue("iso1"));
      if (isotopesMustMatch<0) {
        log_.severe("Isotopes that must match must not be negative");
        System.exit(1);
      }
      if (isotopesMustMatch>10) {
        log_.severe("Isotopes that must match must not be greater than 10");
        System.exit(1);
      }
      if (isotopesMustMatch>amountOfIsotopes)
        amountOfIsotopes = isotopesMustMatch;
    }catch (NumberFormatException nfx){
      log_.severe("Isotopes that must match is not integer format");
      System.exit(1);
    }
    try{
      if (command.getOptionValue("iso2")!=null) amountOfIsotopes = Integer.parseInt(command.getOptionValue("iso2"));
      if (amountOfIsotopes<0) {
        log_.severe("Isotopes that are quantified must not be negative");
        System.exit(1);
      }
      if (amountOfIsotopes>10) {
        log_.severe("Isotopes that are quantified must not be greater than 10");
        System.exit(1);
      }
    }catch (NumberFormatException nfx){
      log_.severe("Isotopes that are quantified is not integer format");
      System.exit(1);
    }
    if (amountOfIsotopes<isotopesMustMatch) {
      log_.severe("Isotopes that are quantified must not be smaller than isotopes that must match");
      System.exit(1);
    }  
    float minusTimeTol = 0f;
    try{
      if (command.getOptionValue("rtTolBef")!=null) {
        minusTimeTol = Float.parseFloat(command.getOptionValue("rtTolBef"));
        if (minusTimeTol<0f) {
          log_.severe("RT before tolerance must not be negative");
          System.exit(1);
        }
      }
    }catch (NumberFormatException nfx){
      log_.severe("RT before tolerance is not float format");
      System.exit(1);
    }
    float plusTimeTol = 0f;
    try{
      if (command.getOptionValue("rtTolAft")!=null) plusTimeTol = Float.parseFloat(command.getOptionValue("rtTolAft"));
      if (plusTimeTol<0f) {
        log_.severe("RT after tolerance must not be negative");
        System.exit(1);
      }
    }catch (NumberFormatException nfx){
      log_.severe("RT after tolerance is not float format");
      System.exit(1);
    }
    float rtShift = 0f;
    try{
      if (command.getOptionValue("rtSh")!=null) rtShift = Float.parseFloat(command.getOptionValue("rtSh"));
    }catch (NumberFormatException nfx){
      log_.severe("retention time shift is not float format");
      System.exit(1);
    }
    boolean searchUnknownBatchTime = true;
    if (command.getOptionValue("u")!=null) {
      if (!(command.getOptionValue("u").equalsIgnoreCase("true") || command.getOptionValue("u").equalsIgnoreCase("false") ||
            command.getOptionValue("u").equalsIgnoreCase("yes") || command.getOptionValue("u").equalsIgnoreCase("no"))) {
        log_.severe("For search unknown retention time is only true/yes/false/no allowed and not \""+command.getOptionValue("u")+"\"");
        System.exit(1);   
      }
      if (command.getOptionValue("u").equalsIgnoreCase("false") || command.getOptionValue("u").equalsIgnoreCase("no"))
        searchUnknownBatchTime = false;
    }
    
    System.out.println("rawDirString: "+rawDirString);
    System.out.println("quantDirString: "+quantDirString);
    System.out.println("nrProcessors: "+nrProcessors);
    System.out.println("nrChromProcessors: "+nrChromProcessors);
    System.out.println("cutoff: "+cutoff);
    System.out.println("isotopesMustMatch: "+isotopesMustMatch);
    System.out.println("amountOfIsotopes: "+amountOfIsotopes);
    System.out.println("minusTimeTol: "+minusTimeTol);
    System.out.println("plusTimeTol: "+plusTimeTol);
    System.out.println("rtShift: "+rtShift);
    System.out.println("searchUnknownBatchTime: "+searchUnknownBatchTime);
    System.out.println("------------------------------------------");
    LDACmd cmd = new LDACmd();
    cmd.prepareAndStartQuatification(rawDirString, quantDirString, nrChromProcessors, nrProcessors, cutoff,
        isotopesMustMatch, amountOfIsotopes, minusTimeTol, plusTimeTol, rtShift,
        searchUnknownBatchTime);
  }
  
  /**
   * checks for the existence of directories, detects quantifiable files and quantification files, and starts the batch quantification
   * @param rawDirString directory containing the MS data
   * @param quantDirString directory containing quantification files
   * @param nrProcessors amount of processors/threads to be used for file translation
   * @param nrProcessors amount of processors/threads to be used for quantification
   * @param cutoff the relative cutoff value in per mille
   * @param isotopesMustMatch number of isotopes that must match the theoretical distribution
   * @param amountOfIsotopes number of isotopes that shall be quantified
   * @param minusTimeTol retention time before tolerance
   * @param plusTimeTol retention time after tolerance
   * @param rtShift retention time shift
   * @param searchUnknownBatchTime search unknown retention time
   */
  private void prepareAndStartQuatification(String rawDirString, String quantDirString,
      int nrChromProcessors, int nrProcessors, float cutoff, int isotopesMustMatch, int amountOfIsotopes,
      float minusTimeTol, float plusTimeTol, float rtShift, boolean searchUnknownBatchTime) {
    BatchQuantificationTableModel batchQuantTableModel = new BatchQuantificationTableModel();
    BatchQuantificationTable batchQuantTable = new BatchQuantificationTable(batchQuantTableModel);
    
    //TODO: these are visual components that are not required for command line
    JLabel quantifyingBatchLabel = new JLabel("Quantifying");
    JProgressBar progressBatchBar = new JProgressBar();
    progressBatchBar.setMaximum(100);

    
    Vector<File> rawFiles = new Vector<File>();
    Vector<File> quantFiles = new Vector<File>();
    if (rawDirString!=null && rawDirString.length()>0 && quantDirString!=null&&quantDirString.length()>0){
      File rawDir = new File(rawDirString );
      File quantDir = new File(quantDirString);
      if (rawDir.exists()&&rawDir.isDirectory()&&quantDir.exists()&&quantDir.isDirectory()){
        File[] rawFileCandidates = rawDir.listFiles();
        Hashtable<String,Vector<File>> avoidDuplication = new Hashtable<String,Vector<File>>();
        boolean mzXMLOrChromPresent = false;
        for (int i=0; i!=rawFileCandidates.length;i++){
          if (rawFileCandidates[i].isFile()){
            String[] fileNameAndSuffix = StaticUtils.extractFileNameAndSuffix(rawFileCandidates[i].getAbsolutePath()); 
            String suffix = fileNameAndSuffix[1];
            String fileName = fileNameAndSuffix[0];
            if (suffix.equalsIgnoreCase("mzxml")||suffix.equalsIgnoreCase("raw")||suffix.equalsIgnoreCase("chrom")||suffix.equalsIgnoreCase("wiff")){
              if (suffix.equalsIgnoreCase("mzxml")||suffix.equalsIgnoreCase("chrom")) mzXMLOrChromPresent = true;
              Vector<File> theFiles = new Vector<File>();
              if (avoidDuplication.containsKey(fileName)){
                
                theFiles = avoidDuplication.get(fileName);
                }
                theFiles.add(rawFileCandidates[i]);
                avoidDuplication.put(fileName, theFiles);
              }
            }
            if (rawFileCandidates[i].isDirectory()){
              String[] fileNameAndSuffix = StaticUtils.extractFileNameAndSuffix(rawFileCandidates[i].getAbsolutePath()); 
              String suffix = fileNameAndSuffix[1];
              String fileName = fileNameAndSuffix[0];
              if (suffix.equalsIgnoreCase("raw")|| suffix.equalsIgnoreCase("d") ||suffix.equalsIgnoreCase("chrom")){
                if (suffix.equalsIgnoreCase("chrom")) mzXMLOrChromPresent = true;
                Vector<File> theFiles = new Vector<File>();
                if (avoidDuplication.containsKey(fileName)){
                  theFiles = avoidDuplication.get(fileName);
                }
                theFiles.add(rawFileCandidates[i]);
                avoidDuplication.put(fileName, theFiles);
              }
            }
          }
          for (String key : avoidDuplication.keySet()){
            Vector<File> theFiles = avoidDuplication.get(key);
            if (theFiles.size()==1){
              String suffix = StaticUtils.extractFileNameAndSuffix(theFiles.get(0).getAbsolutePath())[1];
              if (!mzXMLOrChromPresent || !suffix.equalsIgnoreCase("wiff"))
                rawFiles.add(theFiles.get(0));
            }else{
              int selectedIndex = -1;
              for (int i=0; i!=theFiles.size();i++){
                File file = theFiles.get(i);
                String suffix = file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf(".")+1);
                if (mzXMLOrChromPresent && suffix.equalsIgnoreCase("wiff")) continue;
                if (suffix.equalsIgnoreCase("chrom")){
                  selectedIndex = i;
                }
              }
              if (selectedIndex>-1){
                rawFiles.add(theFiles.get(selectedIndex));
              }else{
                for (int i=0; i!=theFiles.size();i++){
                  File file = theFiles.get(i);
                  String suffix = file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf(".")+1);
                  if (mzXMLOrChromPresent && suffix.equalsIgnoreCase("wiff")) continue;
                  if (suffix.equalsIgnoreCase("mzXML")){
                    rawFiles.add(theFiles.get(i));
                  }
                }  
              }
            }
          }
          File[] quantificationFileCandidates = quantDir.listFiles();
          boolean containsTxtFiles = false;
          for (int i=0; i!=quantificationFileCandidates.length;i++){
            String suffix = quantificationFileCandidates[i].getAbsolutePath().substring(quantificationFileCandidates[i].getAbsolutePath().lastIndexOf(".")+1);
            if (suffix.equalsIgnoreCase("xls")||suffix.equalsIgnoreCase("xlsx")){
              quantFiles.add(quantificationFileCandidates[i]);
            } else if (suffix.equalsIgnoreCase("txt"))
              containsTxtFiles = true;
          }
          if (Settings.useAlex() && containsTxtFiles)
            quantFiles.add(quantDir);
          Vector<RawQuantificationPairVO> pairs = new Vector<RawQuantificationPairVO>();
          
          if (rawFiles.size()>0 && quantFiles.size()>0 && (pairs = LipidDataAnalyzer.generateQuantificationPairVOs(rawFiles,quantFiles)).size()>0){
            boolean ionMode = false;
//            if (this.ionModeBatch_!=null && ((String)ionModeBatch_.getSelectedItem()).equalsIgnoreCase("+"))
//              ionMode = true;
            batchQuantTableModel.clearFiles();
            batchQuantTableModel.addFiles(pairs);
            progressBatchBar.setValue(0);
            try {
              batchQuantThread_ = new BatchQuantThread(batchQuantTable, batchQuantTableModel,progressBatchBar, 
              quantifyingBatchLabel, minusTimeTol,plusTimeTol,amountOfIsotopes,isotopesMustMatch,searchUnknownBatchTime, cutoff, 
                rtShift, nrChromProcessors,nrProcessors,ionMode,true);
              batchQuantThread_.start();
              this.initTimer();
            }catch(Exception ex) {
              ex.printStackTrace();
              System.exit(1);
            }
          }else{
            if (rawFiles.size()==0){
              log_.severe("In the specified raw directory are no quantifyable files");
              System.exit(1);
            }
            if (quantFiles.size()==0){
              log_.severe("In the specified quant directory are no quantifyable files");
              System.exit(1);
            }
            if (rawFiles.size()>0 && quantFiles.size()>0) {
              log_.severe("In the specified directories are no quantifyable raw/quant pairs");
              System.exit(1);
            }
          }
        }else{
          if (!rawDir.exists()||!rawDir.isDirectory()){
            log_.severe("The raw directory does not exist");
            System.exit(1);
          }
          if (!quantDir.exists()||!quantDir.isDirectory()){
            log_.severe("The quantification directory does not exist");
            System.exit(1);
          }
        }
    }
  }
  
  
  /**
   * prints the usage of the command line interface
   * @param options the options specified
   * @param message the message of the exception
   */
  private static void printUsage(Options options, String message)
  {
    HelpFormatter formatter = new HelpFormatter();
    formatter.printHelp( "LDA command line interface... automated quantitation and identification of lipids\n\nERROR: "+message+"\n\n", options );
  }
  
  private void handleTimerEvent(){
    if (this.batchQuantThread_!=null && this.batchQuantThread_.finished()){
      this.batchQuantThread_ = null;
      System.exit(0);
    }
  }

  /**
   * checks the batch quantification thread whether it is finished
   * @author Juergen Hartler
   *
   */
  private class ThreadSupervisor extends TimerTask{

    public void run()
    {
      handleTimerEvent();
    }
  }
  
  /**
   * timer for activating the thread supervisor
   */
  private void initTimer(){
    timer_ = new java.util.Timer();
    timer_.schedule(new ThreadSupervisor(), 10, 1000);
  }

}
