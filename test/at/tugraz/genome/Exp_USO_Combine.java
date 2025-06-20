/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2024 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.export.QuantificationResultExporter;
import org.apache.commons.math3.util.Pair;

/**
 * This class will allow for appending the entries of a results file to another one.
 * @author Leonida M. Lamp
 *
 */
public class Exp_USO_Combine
{
//	private static final String FOLDER_GENERAL = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_II_b\\Samples_Unlabeled\\";
	private static final String FOLDER_GENERAL = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_II_c_and_III_a_60min\\Samples_Unlabeled\\";
//	private static final String FOLDER_GENERAL = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\NIST_repeat\\30min\\";
	private static final String EXPERIMENT_PATH = FOLDER_GENERAL+"before_comb_USO\\";
	private static final String USO_PATH = FOLDER_GENERAL+"USO\\";
	private static final String OUT_PATH = FOLDER_GENERAL;
	private static final String EXCEL_SUFFIX = ".xlsx";
	private static final String EXPERIMENT_TL_SUFFIX = "_SILDA_60min";
	private static final String USO_TL_SUFFIX = "_USO_IS";
	
//	private static final String EXPERIMENT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\labeled_n7\\corr_no_18_3\\";
//	private static final String USO_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\labeled_n7\\corr_18_3\\";
//	private static final String OUT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\labeled_n7\\";
//	private static final String EXCEL_SUFFIX = ".xlsx";
//	private static final String EXPERIMENT_TL_SUFFIX = "_negative_D";
//	private static final String USO_TL_SUFFIX = "_negative_D_only18_3";


  public static void main(String[] args)
  {
  	try
  	{
  		new Exp_USO_Combine();
  	}
  	catch (Exception ex)
  	{
  		ex.printStackTrace();
  	}

  }

  private Exp_USO_Combine() throws ExcelInputFileException, ExportException
  {
  	File experimentDir = new File(EXPERIMENT_PATH);
    File usoDir = new File(USO_PATH);
  	ArrayList<File> experimentFiles = new ArrayList<File>(Arrays.asList(experimentDir.listFiles()));
  	ArrayList<File> usoFiles = new ArrayList<File>(Arrays.asList(usoDir.listFiles()));
  	ArrayList<Pair<File,File>> filePairs = new ArrayList<Pair<File,File>>();

  	for (File exp : experimentFiles)
  	{

  		String expName = getExcelFileNameMinusSuffix(exp, EXPERIMENT_TL_SUFFIX);
  		if (expName.length()<1) continue;

  		for (File uso : usoFiles)
	  	{
  			String usoName = getExcelFileNameMinusSuffix(uso, USO_TL_SUFFIX);
    		if (usoName.length()<1) continue;
    		if (expName.equalsIgnoreCase(usoName))
    		{
    			filePairs.add(new Pair<File,File>(exp,uso));
    			continue;
    		}
	  	}
  	}	

  	for (Pair<File,File> pair : filePairs)
  	{
  		System.out.println(pair.getKey().getName());
  		Hashtable<String,Boolean> showMods = new Hashtable<String,Boolean>();
      QuantificationResult quantExp = LDAResultReader.readResultFile(pair.getKey().getAbsolutePath(), showMods);
      QuantificationResult quantUSO = LDAResultReader.readResultFile(pair.getValue().getAbsolutePath(), showMods);
      Hashtable<String,Vector<LipidParameterSet>> resultsExp = quantExp.getIdentifications();
      Hashtable<String,Vector<LipidParameterSet>> resultsUSO = quantUSO.getIdentifications();
      for (String lipidClass : resultsExp.keySet())
      {
      	if (resultsUSO.containsKey(lipidClass))
      	{
      		resultsExp.get(lipidClass).addAll(resultsUSO.get(lipidClass));
      	}
      }
      String outPath = String.format("%s%s", OUT_PATH, pair.getKey().getName());

      QuantificationResultExporter.writeResultsToExcel(outPath, quantExp);
  	}

  }

  private String getExcelFileNameMinusSuffix(File file, String suffix)
  {
  	String fileName = "";
  	if (file.toString().endsWith(EXCEL_SUFFIX))
  	{
  		String namePath = file.toPath().getFileName().toString();
  		fileName = namePath.substring(0,namePath.indexOf(suffix));
  	}
  	return fileName;  	
  }

}
