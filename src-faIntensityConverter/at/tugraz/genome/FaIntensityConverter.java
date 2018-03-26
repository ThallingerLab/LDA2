/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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
import java.io.FileNotFoundException;

import at.tugraz.genome.lda.exception.ExcelInputFileException;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class FaIntensityConverter
{

  public FaIntensityConverter(String arg){
    convertLDAFilesInCurrentDirectory(arg);
  }
  
  
  public static void main(String[] args)
  {
    String argument = null;
    if (args !=null && args.length>0)
      argument = args[0];
    new FaIntensityConverter(argument);
    
  }

  private void convertLDAFilesInCurrentDirectory(String arg){
    File dir = new File(System.getProperty("user.dir"));
    if (arg!=null && arg.equalsIgnoreCase("-summary")){
      LDAToFASummaryConverter converter = new LDAToFASummaryConverter(dir.getAbsolutePath());
      try {
        converter.convert();
      }
      catch (Exception e) {
        e.printStackTrace();
      }
    }else{
      File[] files = dir.listFiles();
      for (File file : files){
        if (!file.getName().endsWith(".xlsx")) continue;
        LDAToFaConverter converter = new LDAToFaConverter(file.getAbsolutePath());
        try {
          converter.convert();
        } catch (ExcelInputFileException | FileNotFoundException e) {
          System.out.println("---------------------------------------------------");
          System.out.println("There is something wrong with the file: "+file.getName());
          e.printStackTrace();
        }
      }
    }
  }
}
