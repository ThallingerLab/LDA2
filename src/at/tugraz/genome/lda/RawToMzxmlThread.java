/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2021 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp 
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
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
public class RawToMzxmlThread extends Thread
{
  
  String[] params_;
  boolean isMassPlusPlus_;
  boolean watersMsConvert_;
  boolean finished_ = false;
  String errorString_;
  
  public RawToMzxmlThread(String[] params, boolean isMassPlusPlus, boolean watersMsConvert){
    params_ = params;
    isMassPlusPlus_ = isMassPlusPlus;
    watersMsConvert_ = watersMsConvert;
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
      if (isMassPlusPlus_ || watersMsConvert_){
        File[] files = (new File(params_[2])).listFiles();
        int msLevels = 0;
        for (int i=0;i!=files.length;i++){
          if (files[i].getName().startsWith("_FUN") && files[i].getAbsolutePath().endsWith(".DAT"))
            msLevels++;
        }
        String outputBaseFile = "";
        if (isMassPlusPlus_) {
          outputBaseFile = new String(params_[5]);
          for (int i=2;i<=msLevels;i++){
//            System.out.println("MS-Level: "+i);
            params_[5] = outputBaseFile+String.valueOf(i);
            params_[7] = String.valueOf((i-1));
            RawToMzxmlThread.startRawToMzxmlTranslation(params_);
          }
        }
        if (watersMsConvert_) {
          outputBaseFile = params_[2].substring(0,params_[2].lastIndexOf("."))+"."+LipidomicsConstants.getIntermediateFileFormat();
          for (int i=2;i<=msLevels;i++){
            params_[params_.length-1] = "msLevel "+i;              
            params_[4] = StaticUtils.extractDirName(params_[2])+"/"+i;
            RawToMzxmlThread.startRawToMzxmlTranslation(params_);
            String fileName = StaticUtils.extractFileName(params_[2]);
            fileName = fileName.substring(0,fileName.lastIndexOf("."))+"."+LipidomicsConstants.getIntermediateFileFormat();
            File outputFile = new File(params_[4]+"/"+fileName);
            outputFile.renameTo(new File(outputBaseFile+String.valueOf(i)));
            File oldDir = new File(params_[4]);
            oldDir.delete();
          }
        }
        if (msLevels>1 && Settings.mergeMultipleMSMSFiles()){
          int mergingLevels = msLevels;
          if (Settings.skipLastMzXML()) mergingLevels--;
          MzXMLMergerForWaters merger = new MzXMLMergerForWaters(outputBaseFile,mergingLevels);
          if (mergingLevels>1) {
            merger.merge();
          }else {
            File mergedFile = new File(outputBaseFile);
            mergedFile.renameTo(new File(merger.getMergedFileName()));            
          }
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
