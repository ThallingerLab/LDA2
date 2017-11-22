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

package at.tugraz.genome;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.StringTokenizer;

import at.tugraz.genome.voutils.GeneralComparator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MgfSorter
{

  public static void main(String[] args)
  {
    
    if (args.length!=1){
      System.out.println("ERROR: you must provide 1 arguements! The directory that contains the mgf files");
      System.out.println();
      return;
    }
    String dir = args[0];
    MgfSorter sorter  = new MgfSorter(dir);
    
    try {
      sorter.performSorting();
    }
    catch (Exception e) {
      System.out.println(e.getMessage());
    }
   }
  
  private String dir_;
  
  private MgfSorter(String dir){
    this.dir_ = dir;
  }
  
  private void performSorting() throws Exception{
    File directory = new File(dir_);
    if (!directory.isDirectory()) throw new Exception("The provided path is not a directory");
    for (File file : directory.listFiles()){
      if (!file.isFile() || !file.getAbsolutePath().endsWith(".mgf")) continue;
      sortMgfFile(file);
    }
  }

  @SuppressWarnings("unchecked")
  private void sortMgfFile(File inFile) throws Exception{
    LineNumberReader reader = new LineNumberReader(new FileReader(inFile));
    String outFile = inFile.getAbsolutePath()+"-new";
    BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outFile));
    try{
      String line;
      List<MzInt> values = new ArrayList<MzInt>();
      while ((line=reader.readLine()) != null) {
        if (isDataLine(line)){
          StringTokenizer tokenizer = new StringTokenizer(line);
          if (tokenizer.countTokens()!=2) throw new Exception("There is a problem at line: "+reader.getLineNumber());
          String mz = tokenizer.nextToken();
          String intensity = tokenizer.nextToken();
          try{ Float.parseFloat(mz);} catch (NumberFormatException nfx){throw new Exception("The m/z value at line "+reader.getLineNumber()+" is not float format!");};
          try{ Float.parseFloat(mz);} catch (NumberFormatException nfx){throw new Exception("The intensity value at line "+reader.getLineNumber()+" is not float format!");};
          values.add(new MzInt(mz,intensity));
        }else{
          if (values.size()>0){
            Collections.sort(values, new GeneralComparator("at.tugraz.genome.MgfSorter$MzInt", "getMz", "java.lang.Float"));
            String toWrite = "";
            for (MzInt value : values){
              toWrite += value.getMzString()+" "+value.getIntensity()+"\r\n";
            }
            byte[] data = (toWrite).getBytes();
            out.write(data);
          }       
          values = new ArrayList<MzInt>();
          byte[] toWrite = (line+"\r\n").getBytes();
          out.write(toWrite);
        }
      }
    }finally{
      reader.close();
      out.close();
    }
    replaceOldFile(outFile,inFile.getAbsolutePath());
  }
  
  public class MzInt{
    
    private float mz_;
    private String mzString_;
    private String intensity_;
    
    public MzInt (String mz, String intensity){
      this.mz_ = Float.parseFloat(mz);
      this.mzString_ = mz;
      this.intensity_ = intensity;
    }

    public Float getMz()
    {
      return mz_;
    }

    public String getMzString()
    {
      return mzString_;
    }

    public String getIntensity()
    {
      return intensity_;
    } 
  }
  
  private boolean isDataLine(String line){
    boolean isDataLine = false;
    if (line!=null && line.length()>0){
      char[] chars = line.toCharArray();
      boolean foundEmptySpace = false;
      boolean foundOtherChar = false;
      for (char one : chars){
        if (Character.isDigit(one) || one=='.') continue;
        else if (one == ' '){
          foundEmptySpace = true;
          continue;
        }else{
          foundOtherChar = true;
          break;
        }
      }
      if (foundEmptySpace && !foundOtherChar) isDataLine=true;
    }
    return isDataLine;
  }
  
  private void replaceOldFile(String newFile , String oldFile){
    File old = new File(oldFile);
    old.delete();
    File newF = new File(newFile);
    newF.renameTo(old);
  }
}
