/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2019 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.tugraz.genome;

import java.io.File;
import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Vector;

import at.tugraz.genome.lda.alex123.RdbOutputWriter;
import at.tugraz.genome.lda.alex123.RdbParser;
import at.tugraz.genome.lda.analysis.ComparativeAnalysis;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.QuantVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class AlexRtGrouper
{
  
  /** the directory where the files are located in*/
  private String dir_;
  /** the selected grouping time in minutes*/
  private double groupingTime_;

  
  /**
   * constructor for the Alex123 RT-grouper
   * @param dir directory where the files are located in
   * @param groupingTime selected grouping time in minutes
   */
  public AlexRtGrouper(String dir, double groupingTime) {
    this.dir_ = dir;
    this.groupingTime_ = groupingTime;
  }
  
  
  /**
   * starts the RT grouping process, which includes the rewriting of the Alex123 RDB file by adding a column called
   * RT group, which includes the grouping RT time
   * @throws Exception
   */
  public void groupTheEntries() throws Exception {
    Hashtable<String,File> avoidDuplicates = new Hashtable<String,File>();
    Vector<File> resultFiles = new Vector<File>();
    File resultsDir = new File(dir_);
    if (resultsDir.exists() && resultsDir.isDirectory()){
      File[] resultFileCandidates = resultsDir.listFiles();
      for (int i=0; i!=resultFileCandidates.length;i++){
        if (resultFileCandidates[i].isFile() && !avoidDuplicates.containsKey(resultFileCandidates[i].getAbsolutePath())){
          String fileName = StaticUtils.extractFileName(resultFileCandidates[i].getAbsolutePath()); 
          String suffix = fileName.substring(fileName.lastIndexOf(".")+1);
          if (suffix.equalsIgnoreCase("xls")||suffix.equalsIgnoreCase("xlsx")){
            avoidDuplicates.put(resultFileCandidates[i].getAbsolutePath(),resultFileCandidates[i]);
            resultFiles.add(resultFileCandidates[i]);
          }
        }
      }
      resultFiles = sortFilesByName(resultFiles);
      
      
      ComparativeAnalysis analysisModule = new ComparativeAnalysis(resultFiles,"IS","Ex-IS", null, null, -1, 
        null,null,null,groupingTime_);
      //TODO: using default values for statistics view mode: 0 and not combining ox with non ox.
      analysisModule.parseInput(0,false);
      analysisModule.calculateStatistics();
        
      //TODO: this is a fixed value - might be taken as input parameter in future
      RdbOutputWriter rdbWriter = new RdbOutputWriter("IS","Ex-IS");
      for (String exp : analysisModule.getExpNamesInSequence()) {
        File excelFile = analysisModule.getFullFilePath(exp);
        String alexFile = getAlexFileNameFromExcelFile(excelFile);
        File file = new File(alexFile);
        System.out.println("Rt-grouping: "+file.getName());
        LinkedHashMap<String,Integer> classSequence = null;
        Hashtable<String,Vector<String>> analyteSequence = null;
        Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> targetLists = null;
        if (file.exists() && file.isFile()) {
          RdbParser parser = new RdbParser(alexFile);
          parser.parse();
          classSequence = parser.getClassSequence();
          targetLists = parser.getTargetlistInfo();
          analyteSequence = parser.getAnalyteSequence();
        }
        Vector<String> excelFiles = new Vector<String>();
        excelFiles.add(excelFile.getAbsolutePath());
        rdbWriter.write(alexFile, analysisModule, classSequence,analyteSequence,
            null, null,targetLists,excelFiles,true);   
      }
    }
  }
  
  private String getAlexFileNameFromExcelFile(File excelFile) {
    return excelFile.getAbsolutePath().substring(0,excelFile.getAbsolutePath().lastIndexOf("."))+".tab";
  }
  
  private static Vector<File> sortFilesByName(Vector<File> filesToSort){
    Vector<File> newOrder = new Vector<File>(filesToSort.size());
    Hashtable<String,File> hash = new Hashtable<String,File>();
    for (File file : filesToSort) hash.put(file.getAbsolutePath(), file);
    Vector<String> pathVect = new Vector<String>(hash.keySet());
    Collections.sort(pathVect);
    for (String path : pathVect) newOrder.add(hash.get(path));
    return newOrder;
  }
  
  /**
   * 
   * @param args
   */
  public static void main(String[] args)
  {
    if (args.length!=2){
      if (args.length!=0){
        System.out.println("ERROR: you must provide 2 arguements! You provided just: "+args.length);
        System.out.println();
      }
      printUsage();
      return;
    }
    
    
    File dir = new File(args[0]);
    if (!dir.exists()){
      System.out.println("ERROR: The directory \""+dir.getAbsolutePath()+"\" does not exist!");
      System.out.println();
      return;
    }
    if (dir.isFile()){
      System.out.println("The directory \""+dir.getAbsolutePath()+"\" is a file! Only directories are allowed!");
      System.out.println();
      return;
    }
    double groupingRt = -1d;
    try{
      groupingRt=Double.parseDouble(args[1].replace(",", "."));
    }catch (NumberFormatException nfx){
      System.out.println("grouping retention time is not double format");
      System.out.println();
      return;
    }
    if (groupingRt<=0d) {
      System.out.println("ERROR: Only rt grouping times > 0 are allowed!");
      System.out.println();
      return;
    }
    AlexRtGrouper grouper = new AlexRtGrouper(dir.getAbsolutePath(),groupingRt);
    try {
      grouper.groupTheEntries();
      System.out.println("Finished");
    }
    catch (Exception e) {
      System.out.println("ERROR:");
      e.printStackTrace();
    }
  }
  
  /**
   * prints the usage of the command line tool to the console
   */
  private static void printUsage()
  {
    System.out.println("\nAlex123 Rt-grouper for LDA files:");
    System.out.println();
    System.out.println("AlexRtGrouper $dirToGroupingFiles $groupingTimeInMinutes");
  }

  
  

}
