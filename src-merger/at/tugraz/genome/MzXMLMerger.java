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

import java.io.File;

import at.tugraz.genome.lda.parser.MzXMLMergerForWaters;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MzXMLMerger
{
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
    File mzXMLFile = null;
    int msLevel = -1;
    mzXMLFile = new File(args[0]);
    if (!mzXMLFile.exists()){
      System.out.println("The file \""+mzXMLFile.getAbsolutePath()+"\" does not exist!");
      System.out.println();
      return;
    }
    if (mzXMLFile.isDirectory()){
      System.out.println("The file \""+mzXMLFile.getAbsolutePath()+"\" is a directory! Only files are allowed!");
      System.out.println();
      return;
    }
    try{msLevel=Integer.parseInt(args[1].replace(",", "."));}catch (NumberFormatException nfx){
      System.out.println("ms-levels is not integer format");
      System.out.println();
      return;
    }
    if (msLevel<1){
      System.out.println("Only ms-levels>=1 are allowed!");
      System.out.println();
      return;      
    }
    MzXMLMergerForWaters merger = new MzXMLMergerForWaters(mzXMLFile.getAbsolutePath(),msLevel);
    try {
      merger.merge();
    }
    catch (Exception e) {
      System.out.println("ERROR:");
      e.printStackTrace();
    }

  }
  
  private static void printUsage()
  {
    System.out.println("\nMzXML-merger for several MS-levels:");
    System.out.println();
    System.out.println("mzXMLMerger $mzXMLFile $ms-levels");
  }

  
}
