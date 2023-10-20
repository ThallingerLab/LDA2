package at.tugraz.genome;

import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.export.QuantificationResultExporter;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class RTFilter
{
	private static final String RESULT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\labeled_n3\\28_D5-18-3(n-3)_327_30min_negative_C.xlsx";
	private static final String LABEL = "C";
	
	public static void main(String[] args)
  {
    // MainFrame frame = new MainFrame(new TestClass(), 1024, 1024);
    filterRT(100.0, 110.0);
  }
	
	/**
	 * 
	 * @param lowerBound 		lowest allowed rt
	 * @param upperBound		highest allowed rt
	 */
	private static void filterRT(double lowerBound, double upperBound)
	{
		try
		{
			QuantificationResult quantRes = LDAResultReader.readResultFile(RESULT_PATH, new Hashtable<String,Boolean>());
			Hashtable<String,Vector<LipidParameterSet>> results = quantRes.getIdentifications();
			Hashtable<String,Vector<LipidParameterSet>> newResults = new Hashtable<String,Vector<LipidParameterSet>>();
			for (String lipidClass : results.keySet())
			{
				newResults.put(lipidClass, new Vector<LipidParameterSet>());
				Vector<LipidParameterSet> params = results.get(lipidClass);
				for (LipidParameterSet param : params)
				{
					if (!param.getNameStringWithoutRt().contains(LABEL))
					{
						newResults.get(lipidClass).add(param);
					}
					else if (param.getPreciseRT() >= lowerBound && param.getPreciseRT() <= upperBound)
					{
						newResults.get(lipidClass).add(param);
					}
				}
			}
			quantRes.setIdentifications(newResults);
			QuantificationResultExporter.writeResultsToExcel(RESULT_PATH, quantRes);
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}
}
