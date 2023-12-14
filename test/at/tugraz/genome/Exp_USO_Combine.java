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
import javafx.util.Pair;

/**
 * This class will allow for appending the entries of a results file to another one.
 * @author Leonida M. Lamp
 *
 */
public class Exp_USO_Combine
{
	private static final String EXPERIMENT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_II_b\\Samples_Unlabeled\\before_comb_USO\\";
	private static final String USO_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_II_b\\Samples_Unlabeled\\USO_PC\\";
	private static final String OUT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_II_b\\Samples_Unlabeled\\";
	private static final String EXCEL_SUFFIX = ".xlsx";
	private static final String EXPERIMENT_TL_SUFFIX = "_SILDA_30min";
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
