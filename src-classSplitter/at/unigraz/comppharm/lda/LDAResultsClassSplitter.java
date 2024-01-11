/*
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2021 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.unigraz.comppharm.lda;

import java.io.File;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.export.QuantificationResultExporter;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LDAResultsClassSplitter
{
  /* the directory that contains the LDA files to be split*/
  File dir_;
  
  /**
   * constructor for the LDA files splitter
   * @param dir directory containing the LDA files to be split
   * @throws Exception general Exception
   */
  public LDAResultsClassSplitter(File dir) throws Exception {
    this.dir_ = dir;
    splitTheFiles();
  }
  
  /**
   * starts the splitting process
   * @throws Exception general Exception
   */
  private void splitTheFiles() throws Exception {
    QuantificationResult quantRes;
    QuantificationResult splitRes;
    File classDir;
    String splitResultFile;
    Hashtable<String,Vector<LipidParameterSet>> identificationClass;
    for (File file : dir_.listFiles()) {
      if (file.isFile() && file.getName().endsWith(".xlsx")) {
        quantRes = LDAResultReader.readResultFile(file.getAbsolutePath(),  new Hashtable<String,Boolean>());
        for (String lClass : quantRes.getIdentifications().keySet()) {
          classDir = new File(dir_.getAbsolutePath()+File.separator+lClass);
          if (!classDir.exists())
            classDir.mkdir();
          splitResultFile = new String(classDir.getAbsolutePath()+File.separator+file.getName().substring(0,file.getName().length()-".xlsx".length())+"_"+lClass+".xlsx");
          identificationClass = new Hashtable<String,Vector<LipidParameterSet>>();
          identificationClass.put(lClass, quantRes.getIdentifications().get(lClass));
          splitRes = new QuantificationResult(identificationClass, quantRes.getConstants(), quantRes.getMsLevels(), quantRes.getFaHydroxyEncoding(), quantRes.getLcbHydroxyEncoding());
          QuantificationResultExporter.writeResultsToExcel(splitResultFile, splitRes);          
        }
      }
    }
  }
  
  

  public static void main(String[] args)
  {
    File dir = new File(System.getProperty("user.dir"));
    //TODO: the next line is for testing - remove it!
    //dir = new File("C:\\LDA-Collaborations\\Martin\\20210706");
    try {
      new LDAResultsClassSplitter(dir);
    }
    catch (Exception e) {
      System.out.println("ERROR:");
      e.printStackTrace();
    }
  }

}
