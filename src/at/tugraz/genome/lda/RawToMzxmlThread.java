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
import java.io.File;

import at.tugraz.genome.lda.parser.MzXMLMergerForWaters;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class RawToMzxmlThread extends Thread
{
  
  String[] params_;
  boolean isMassPlusPlus_;
  boolean finished_ = false;
  String errorString_;
  
  public RawToMzxmlThread(String[] params, boolean isMassPlusPlus){
    params_ = params;
    isMassPlusPlus_ = isMassPlusPlus;
  }
  
  private static void startRawToMzxmlTranslation(String[] params) throws Exception{
    Process process = Runtime.getRuntime().exec(params);
    if (params.length>2&&params[2].endsWith(".d")){
      BufferedInputStream in = new BufferedInputStream(process.getInputStream());
      in.close();
    }
    process.waitFor();
////    while (isAlive(process)) Thread.sleep( 100 );
  }
  
  public static boolean isAlive( Process p ) {
    try
    {
        p.exitValue();
        return false;
    } catch (IllegalThreadStateException e) {
        return true;
    }
  }
  
  public void run(){
    try{
      RawToMzxmlThread.startRawToMzxmlTranslation(params_);
      if (isMassPlusPlus_){
        File[] files = (new File(params_[2])).listFiles();
        int msLevels = 0;
        for (int i=0;i!=files.length;i++){
          if (files[i].getName().startsWith("_FUN") && files[i].getAbsolutePath().endsWith(".DAT"))
            msLevels++;
        }
        String outputBaseFile = new String(params_[5]);
        for (int i=2;i<=msLevels;i++){
          System.out.println("MS-Level: "+i);
          params_[5] = outputBaseFile+String.valueOf(i);
          params_[7] = String.valueOf((i-1));
          RawToMzxmlThread.startRawToMzxmlTranslation(params_);
        }
        if (msLevels>1 && Settings.mergeMultipleMSMSFiles()){
          int mergingLevels = msLevels;
          if (Settings.skipLastMzXML()) mergingLevels--;
          MzXMLMergerForWaters merger = new MzXMLMergerForWaters(outputBaseFile,mergingLevels);
          merger.merge();
          // here I have to implement the deletion of the other mzXML files and the renaming of the merged file
          for (int i=1; i<=msLevels;i++){
            String fileName = new String(outputBaseFile);
            if (i>1) fileName += String.valueOf(i);
            (new File(fileName)).delete();
          }
          File mergedFile = new File(merger.getMergedFileName());
          mergedFile.renameTo(new File(outputBaseFile));
        }
        //here I have to add the merging of mzXML Files
      }
    } catch (Exception ex){
      ex.printStackTrace();
      errorString_ = ex.toString();     
    }
    finished_ = true;
  }
  
  public boolean finished(){
    return this.finished_;
  }
  
  public String getErrorString(){
    return this.errorString_;
  }
  
  public static void deleteMzXMLFiles(String baseMzxmlPath){
    File mzXMLFile = new File(baseMzxmlPath);
    mzXMLFile.delete();
    int count=2;
    mzXMLFile = new File(baseMzxmlPath+String.valueOf(count));
    while (mzXMLFile.exists()){
      mzXMLFile.delete();
      count++;
      mzXMLFile = new File(baseMzxmlPath+String.valueOf(count));
    }
  }
}
