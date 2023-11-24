package at.tugraz.genome;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.QuantificationThread;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.target.export.TargetListExporter;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.QuantVO;
import javafx.util.Pair;

public class TargetListCombiner
{
	static String outPath_ = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_30min\\TL\\Recalibration\\Recalibration_Comparison_Grad_Adj_Final\\Combined_30min\\SILDA_30min_autoCombined.xlsx";
	static String inPath_SILDA_I_a_ = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_30min\\TL\\Recalibration\\Recalibration_Comparison_Grad_Adj_Final\\Combined_30min\\SILDA_I_a_to_II_b.xlsx";
	static String inPath_SILDA_II_a_ = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_30min\\TL\\Recalibration\\Recalibration_Comparison_Grad_Adj_Final\\Combined_30min\\SILDA_II_a_to_II_b.xlsx";
	static String inPath_SILDA_II_b_ = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_30min\\TL\\Recalibration\\Recalibration_Comparison_Grad_Adj_Final\\Combined_30min\\SILDA_II_b.xlsx";
	static String inPath_SILDA_60min_ = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_30min\\TL\\Recalibration\\Recalibration_Comparison_Grad_Adj_Final\\Combined_30min\\SILDA_60min_to_II_b.xlsx";
	static Double rtCutoff_ = 0.1;
	static Double rtCutoff60min_ = 0.15;
	private final static String FOUND_IN = "Found in:";
	static Pair<String,Integer> description_SILDA_I_a_ = new Pair<String,Integer>("SILDA_I_a", 0);
	static Pair<String,Integer> description_SILDA_II_a_ = new Pair<String,Integer>("SILDA_II_a", 1);
	static Pair<String,Integer> description_SILDA_II_b_ = new Pair<String,Integer>("SILDA_II_b", 2);
	static Pair<String,Integer> description_SILDA_60min_ = new Pair<String,Integer>("SILDA_60min", 3);
	
	public static void main(String[] args)
  {
		
		try (	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outPath_));
					XSSFWorkbook workbook = new XSSFWorkbook();)
		{
			
			Vector<?> quantInfo_I_a = QuantificationThread.getCorrectAnalyteSequence(inPath_SILDA_I_a_,false);
			Vector<?> quantInfo_II_a = QuantificationThread.getCorrectAnalyteSequence(inPath_SILDA_II_a_,false);
			Vector<?> quantInfo_II_b = QuantificationThread.getCorrectAnalyteSequence(inPath_SILDA_II_b_,false);
			Vector<?> quantInfo_60min = QuantificationThread.getCorrectAnalyteSequence(inPath_SILDA_60min_,false);
			
			Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects_I_a = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo_I_a.get(4);
			Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects_II_a = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo_II_a.get(4);
			Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects_II_b = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo_II_b.get(4);
			Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects_60min = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo_60min.get(4);
			
			LinkedHashMap<String,Integer> classSequence = (LinkedHashMap<String,Integer>)quantInfo_II_b.get(0);
			LinkedHashMap<String,Vector<String>> analyteSequence = (LinkedHashMap<String,Vector<String>>)quantInfo_II_b.get(1);
			
			for (String cName : classSequence.keySet()) 
	    {
				Hashtable<String,Hashtable<String,QuantVO>> toExport = new Hashtable<String,Hashtable<String,QuantVO>>();
				for (String analyte : analyteSequence.get(cName)) 
		    {
					toExport.put(analyte, new Hashtable<String,QuantVO>());
					Hashtable<String,QuantVO> quantAnalytes_I_a = quantObjects_I_a.get(cName).get(analyte);
					Hashtable<String,QuantVO> quantAnalytes_II_a = quantObjects_II_a.get(cName).get(analyte);
					Hashtable<String,QuantVO> quantAnalytes_II_b = quantObjects_II_b.get(cName).get(analyte);
					Hashtable<String,QuantVO> quantAnalytes_60min = quantObjects_60min.get(cName).get(analyte);
					
					for (String mod : quantAnalytes_II_b.keySet()) 
		      {
						
						QuantVO quant_I_a = quantAnalytes_I_a.get(mod);
						QuantVO quant_II_a = quantAnalytes_II_a.get(mod);
						QuantVO quant_II_b = quantAnalytes_II_b.get(mod);
						QuantVO quant_60min = quantAnalytes_60min.get(mod);
						QuantVO quantToExport = quant_II_b;
						
		        Vector<DoubleBondPositionVO> quantDB_I_a = quant_I_a.getInfoForOmegaAssignment();
		        Vector<DoubleBondPositionVO> quantDB_II_a = quant_II_a.getInfoForOmegaAssignment();
		        Vector<DoubleBondPositionVO> quantDB_II_b = quant_II_b.getInfoForOmegaAssignment();
		        Vector<DoubleBondPositionVO> quantDB_60min = quant_60min.getInfoForOmegaAssignment();
		        
		        Vector<DoubleBondPositionVO> combinedDB = initCombinedDB(quantDB_II_b, description_SILDA_II_b_);
		        if (!quantDB_I_a.isEmpty() || !quantDB_II_a.isEmpty() || !quantDB_60min.isEmpty())
		        {
		        	
		        	findCloseMatches(combinedDB, quantDB_I_a, description_SILDA_I_a_);
		        	findCloseMatches(combinedDB, quantDB_II_a, description_SILDA_II_a_);
		        	findCloseMatches(combinedDB, quantDB_60min, description_SILDA_60min_);
		        	
		        	Collections.sort(combinedDB);
		        }
		        quantToExport.setInfoForOmegaAssignment(combinedDB);
		        toExport.get(analyte).put(mod, quantToExport);
		      }
		    }
				Sheet sheet = workbook.createSheet(cName);
        writeMassListForSheet(sheet, TargetListExporter.getHeaderStyle(workbook), TargetListExporter.getNumberStyle(workbook), 
        		analyteSequence.get(cName), toExport, cName);
	    }
			workbook.write(out);
      System.out.println("workbook written!");
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
  }
	
	@SuppressWarnings("unchecked")
	private static void writeMassListForSheet(
			Sheet sheet, 
			XSSFCellStyle headerStyle, 
			XSSFCellStyle numberStyle,
			Vector<String> analytes, 
			Hashtable<String,Hashtable<String,QuantVO>> quantObjects, 
			String cName)
	{
		TargetListExporter exporter = new TargetListExporter(null,null,null);
		try 
		{
			Vector<Object> elModCharge = exporter.getAvailableElementsAndModificationsPlusCharge(quantObjects);
	    Vector<String> elements = (Vector<String>)elModCharge.get(0);
	    LinkedHashMap<String,String> mods = (LinkedHashMap<String,String>)elModCharge.get(1);
	    Hashtable<String,Integer> modToCharge = (Hashtable<String,Integer>)elModCharge.get(2);
	    
	    List<String> headerTitles = TargetListExporter.createHeaderTitles(elements, mods, modToCharge);
	    headerTitles.add(FOUND_IN);
	    headerTitles.add(description_SILDA_I_a_.getKey());
	    headerTitles.add(description_SILDA_II_a_.getKey());
	    headerTitles.add(description_SILDA_II_b_.getKey());
	    headerTitles.add(description_SILDA_60min_.getKey());
	    
	    if (cName.equals("SM") || cName.equals("Cer"))
	    {
	    	int rowCount = 0;
    		Row outRow = sheet.createRow(rowCount);
        Cell cell = outRow.createCell(9,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("OH-Number: 2");
        cell = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("OH-Range: 2-3");
        cell = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("adductInsensitiveRtFilter");
	    }
	    TargetListExporter.createHeader(sheet, headerTitles, headerStyle);
	    sheet.setColumnWidth(headerTitles.indexOf(FOUND_IN), 45 * 256);
	    
	    ArrayList<String> usedQuantIDs = new ArrayList<String>();
	    int firstModColumn = headerTitles.indexOf(headerTitles.stream().filter((s)->(s).startsWith("mass(form[")).findFirst().get());
	    int rowCount = TargetListExporter.HEADER_ROW+1;
	    for (String analyte : analytes) 
	    {
	    	Cell cell;
	    	
	    	Hashtable<String,QuantVO> quantAnalytes = quantObjects.get(analyte);
	      Hashtable<String,Integer> formula = new Hashtable<String,Integer>();
	      Vector<DoubleBondPositionVO> doubleBondPositionVOs = new Vector<DoubleBondPositionVO>(); //they will be identical for all mods..
	      for (String mod : mods.keySet()) 
	      {
	        QuantVO quant = quantAnalytes.get(mod);
	        doubleBondPositionVOs = quant.getInfoForOmegaAssignment();
	      }
	      
	      if (!doubleBondPositionVOs.isEmpty())
	      	Collections.sort(doubleBondPositionVOs);
	      
	      for (int i=-1;i<doubleBondPositionVOs.size();i++)
	      {
	      	int modCount = 0;
	      	Row row = sheet.createRow(rowCount);
	      	for (String mod : mods.keySet()) 
		      {
		      	QuantVO quant = quantAnalytes.get(mod);
		      	
		      	//the following lines only affect Cer, SM (classes with hydroxylation sites).
		      	if ((i<0 && usedQuantIDs.contains(exporter.computeQuantID(quant, false)) ||
		      			(i>-1 && !usedQuantIDs.contains(exporter.computeQuantID(quant, true)))))
		        {
		      		rowCount--;
		        	continue;
		        }
		        usedQuantIDs.add(exporter.computeQuantID(quant, true));
		        usedQuantIDs.add(exporter.computeQuantID(quant, false));
		      	
		      	if (modCount==0) 
		        {   	
		          cell = row.createCell(headerTitles.indexOf(TargetListExporter.HEADER_NAME),HSSFCell.CELL_TYPE_STRING);
		          cell.setCellValue(quant.getCarbons());
		          cell = row.createCell(headerTitles.indexOf(TargetListExporter.HEADER_COLON),HSSFCell.CELL_TYPE_STRING);
		          cell.setCellValue(":");
		          cell = row.createCell(headerTitles.indexOf(TargetListExporter.HEADER_DBS),HSSFCell.CELL_TYPE_NUMERIC);
		          cell.setCellValue(quant.getDbs());
		          if (i>=0)
		          {
		          	cell = row.createCell(headerTitles.indexOf(TargetListExporter.HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS),HSSFCell.CELL_TYPE_NUMERIC);
			          cell.setCellValue(doubleBondPositionVOs.get(i).getDoubleBondPositionsHumanReadable(LipidomicsConstants.EXPORT_ANALYTE_TYPE_CHAIN));
		          }
		          formula = StaticUtils.categorizeFormula(quant.getAnalyteFormula());
		          for (String element : formula.keySet()) {
		            cell = row.createCell(headerTitles.indexOf(element),HSSFCell.CELL_TYPE_NUMERIC);
		            cell.setCellValue(formula.get(element));
		          }
		        }
		        cell = row.createCell(firstModColumn+modCount,HSSFCell.CELL_TYPE_NUMERIC);
		        cell.setCellValue(quant.getAnalyteMass());
		        modCount++;
		        if (i>=0)
	          {
	          	cell = row.createCell(headerTitles.indexOf(TargetListExporter.HEADER_RETENTION_TIME),HSSFCell.CELL_TYPE_NUMERIC);
		          cell.setCellValue(doubleBondPositionVOs.get(i).getExpectedRetentionTime());
		          cell.setCellStyle(numberStyle);
		          
		          cell = row.createCell(headerTitles.indexOf(FOUND_IN), HSSFCell.CELL_TYPE_STRING);
		          cell.setCellValue(getDescription(doubleBondPositionVOs.get(i).getExperimentRTLookup()));
		        	
		          if (doubleBondPositionVOs.get(i).getExperimentRTLookup().keySet().size()>1)
		          {
		          	for (Pair<String,Integer> experiment : doubleBondPositionVOs.get(i).getExperimentRTLookup().keySet())
			          {
			          	cell = row.createCell(headerTitles.indexOf(experiment.getKey()), HSSFCell.CELL_TYPE_NUMERIC);
			          	cell.setCellValue(doubleBondPositionVOs.get(i).getExpectedRetentionTime() - 
			          			doubleBondPositionVOs.get(i).getExperimentRTLookup().get(experiment));
			          	cell.setCellStyle(numberStyle);
			          }
		          }
	          }
		      }
	      	rowCount++;
	      }
	    }
		}
		catch (ChemicalFormulaException ex)
		{
			ex.printStackTrace();
		}
	}
	
	
	private static Vector<DoubleBondPositionVO> initCombinedDB(Vector<DoubleBondPositionVO> target, Pair<String,Integer> description)
	{
		Vector<DoubleBondPositionVO> combinedDB = new Vector<DoubleBondPositionVO>();
		for (DoubleBondPositionVO db : target)
		{
			db.addToExperimentRTLookup(description, db.getExpectedRetentionTime());
			combinedDB.add(db);
		}
		return combinedDB;
	}
	
	
	private static void findCloseMatches(Vector<DoubleBondPositionVO> target, Vector<DoubleBondPositionVO> toAdd, Pair<String,Integer> toAddName)
	{
		Vector<DoubleBondPositionVO> added = new Vector<DoubleBondPositionVO>();
		
		for (DoubleBondPositionVO dbOne : toAdd)
		{
			for (DoubleBondPositionVO dbTwo : target)
			{
				if (dbOne.getEncodedDetailed(true, true).equals(dbTwo.getEncodedDetailed(true, true)))
				{
					float rtDiff = Math.abs(dbOne.getExpectedRetentionTime() - dbTwo.getExpectedRetentionTime());
					if (rtDiff < rtCutoff_ || 
							(toAddName.equals(description_SILDA_60min_) && rtDiff < rtCutoff60min_))
					{
						dbTwo.addToExperimentRTLookup(toAddName, dbOne.getExpectedRetentionTime());
						if (!dbTwo.getExperimentRTLookup().containsKey(description_SILDA_II_b_) && 
								!toAddName.equals(description_SILDA_60min_))
						{
							dbTwo.setExpectedRetentionTime((dbOne.getExpectedRetentionTime()+dbTwo.getExpectedRetentionTime())/2);
						}
						added.add(dbOne);
					}
				}
			}
		}
		toAdd.removeAll(added);
		
		for (DoubleBondPositionVO db : toAdd)
		{
			db.addToExperimentRTLookup(toAddName, db.getExpectedRetentionTime());
			target.add(db);
		}
	}
	
	private static String getDescription(Hashtable<Pair<String,Integer>,Float> experimentRTLookup)
	{
		Vector<Pair<String,Integer>> keys = new Vector<Pair<String,Integer>>(experimentRTLookup.keySet());
		Collections.sort(keys, new ExperimentRankComparator());
		StringBuilder description = new StringBuilder();
		for (int i=0;i<keys.size();i++)
		{
			description.append(keys.get(i).getKey());
			if(i+1<keys.size())
			{
				description.append(", ");
			}
		}
		return description.toString();
	}
	
	static class ExperimentRankComparator implements Comparator<Pair<String,Integer>> {

		@Override
		public int compare(Pair<String,Integer> arg0, Pair<String,Integer> arg1)
		{
			return Integer.compare(arg0.getValue(), arg1.getValue());
		}
	}
	
}
