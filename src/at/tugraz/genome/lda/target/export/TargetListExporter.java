/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2023 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.target.export;

import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.swing.JFrame;

import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFCellStyle;
import org.apache.poi.hssf.usermodel.HSSFFont;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFFont;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.dhatim.fastexcel.reader.CellType;
import org.dhatim.fastexcel.reader.ReadableWorkbook;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.QuantificationThread;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.target.IsotopeLabelVO;
import at.tugraz.genome.lda.target.calibration.CalibrationGraphPanel;
import at.tugraz.genome.lda.target.calibration.RecalibrationRegression;
import at.tugraz.genome.lda.target.experiment.IsotopeEffectRegression;
import at.tugraz.genome.lda.utils.RangeDouble;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.lda.vos.ResultFileVO;
import javafx.util.Pair;


/**
 * 
 * Writes RTDBs either newly generated from experimental data or mapped to chromatographic conditions.
 * 
 * @author Leonida M. Lamp
 *
 */
public class TargetListExporter
{	
	public final static String HEADER_NAME = "Name";
	public final static String HEADER_COLON = "";
	public final static String HEADER_DBS = "dbs";
  public final static String HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS = "mol. species";
  public final static String HEADER_RETENTION_TIME = "tR (min)";
  
  public final static int HEADER_ROW = 1;
  
//  private Hashtable<String, RecalibrationRegression> recalibrationRegressionsForClass_;
//  private RecalibrationRegression recalibrationRegression_;
  private IsotopeEffectRegression isotopeEffectRegression_;
  private Vector<ResultFileVO> resultFileVO_;
  private Vector<IsotopeLabelVO> labels_;
  private String templatePath_;
  private boolean calibrateSeparately_ = false;
  private CalibrationGraphPanel calibrationGraphPanel_ = null;
  
  //TODO: just for development! for comparing target lists - remove later (first cName, second: original rt, third recalibr. doublebondPositionVO)
  Hashtable<String,Set<Pair<Double,DoubleBondPositionVO>>> beforeAfter_ = new Hashtable<String,Set<Pair<Double,DoubleBondPositionVO>>>(); 
  Hashtable<String,Vector<Pair<Double,DoubleBondPositionVO>>> comparisonPairsOfClass_ = new Hashtable<String,Vector<Pair<Double,DoubleBondPositionVO>>>();
  
  
  /**
   * Constructor for RTDB mapping
   * The regressions for RT mapping are retrieved from @param calibrationGraphPanel depending on @param calibrateSeparately.
   * @param templatePath
   * @param calibrateSeparately
   * @param calibrationGraphPanel
   */
  public TargetListExporter(String templatePath, boolean calibrateSeparately, CalibrationGraphPanel calibrationGraphPanel)
  {
  	this.calibrateSeparately_ = calibrateSeparately;
  	this.templatePath_ = templatePath;
  	this.calibrationGraphPanel_ = calibrationGraphPanel;
  	this.beforeAfter_ = new Hashtable<String,Set<Pair<Double,DoubleBondPositionVO>>>(); //TODO: just for development!
  }
  
  
  /**
   * Constructor for RTDB creation with stable isotope labeled standards
   * @param isotopeEffectRegression
   * @param resultFileVO
   * @param labels
   */
	public TargetListExporter(
			IsotopeEffectRegression isotopeEffectRegression, Vector<ResultFileVO> resultFileVO, Vector<IsotopeLabelVO> labels)
	{
		this.isotopeEffectRegression_ = isotopeEffectRegression;
		this.resultFileVO_ = resultFileVO;
		this.labels_ = labels;
	}
	
	/**
	 * 
	 * @param templatePath
	 * @param exportPanel
	 * @throws ExportException
	 */
	protected void export(ExportPanel exportPanel) throws ExportException
	{		
		try (	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(exportPanel.getOutPath()));
					XSSFWorkbook workbook = new XSSFWorkbook();)
		{
			if (isRTMapping())
			{
				this.exportRTMappedDB(workbook);
			}
			else
			{
				this.exportNewRTDB(workbook, exportPanel);
			}
			workbook.write(out);
	    System.out.println("workbook written!");
    } 
		catch (Exception e) 
		{
      throw new ExportException(e);
    } 
	}
	
	/**
	 * Exports a RTDB mapped to provided chromatographic conditions
	 * @param workbook
	 * @throws Exception
	 */
	@SuppressWarnings("unchecked")
	private void exportRTMappedDB(XSSFWorkbook workbook) throws Exception
	{
		Vector<?> quantInfo = QuantificationThread.getCorrectAnalyteSequence(this.getTemplatePath(),false);
		LinkedHashMap<String,Integer> classSequence = (LinkedHashMap<String,Integer>)quantInfo.get(0);
		LinkedHashMap<String,Vector<String>> analyteSequence = (LinkedHashMap<String,Vector<String>>)quantInfo.get(1);
		Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo.get(4);
    
    for (String cName : classSequence.keySet()) 
    {
    	Vector<DoubleBondPositionVO> allLabeledVOsToAdd = new Vector<DoubleBondPositionVO>();
      Sheet sheet = workbook.createSheet(cName);
      writeMassListForSheet(sheet, getHeaderStyle(workbook), getNumberStyle(workbook), analyteSequence.get(cName), quantObjects.get(cName), cName, allLabeledVOsToAdd);
    }
	}
	
	/**
	 * Exports a RTDB generated from experiments
	 * For the calculation of the expected RT, aTIE is computed based on the gradient chosen in @param exportPanel and the TIE parameters identified.
	 * For adding the omega position information, the prefix of the fatty acid chain is considered (can only have one prefix)
	 * Based on the user provided clustering threshold in @param exportPanel, identifications are grouped.
	 * If within a group we have species with only one different deuterium label each and the labels are compatible, add omega pos to all in the cluster 
	 * For more than one deuterium label in one species, the TIE is corrected additively for each label, and if a singly labeled species with matching expected RT exists, 
	 * the species is assigned and given the expected RT of the singly labeled species as that is deemed more reliable.
	 * 
	 * @param workbook
	 * @param exportPanel
	 * @throws Exception
	 */
	@SuppressWarnings("unchecked")
	private void exportNewRTDB(XSSFWorkbook workbook, ExportPanel exportPanel) throws Exception
	{
		Vector<?> quantInfo = QuantificationThread.getCorrectAnalyteSequence(exportPanel.getTemplatePath(),false);
		LinkedHashMap<String,Integer> classSequence = (LinkedHashMap<String,Integer>)quantInfo.get(0);
		LinkedHashMap<String,Vector<String>> analyteSequence = (LinkedHashMap<String,Vector<String>>)quantInfo.get(1);
		Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo.get(4);
		
		GradientAdjuster adjuster = GradientParser.parseGradient(exportPanel.getExportOptions().getSelectedGradient());
		if (adjuster == null && exportPanel.getExportOptions().isGradientSelected()) throw new ExportException("The gradient file could not be read!");
    
    for (String cName : classSequence.keySet()) 
    {
    	Vector<DoubleBondPositionVO> allLabeledVOsToAdd = new Vector<DoubleBondPositionVO>();
    	
    	Double clusteringThreshold = exportPanel.getExportOptions().getSelectedClustering();
  		MolecularSpeciesContainer container = new MolecularSpeciesContainer();
    	fillContainerForClass(cName, container, adjuster);
    	for (String molecularSpecies : container.getSingleLabeledMolecularSpecies())
    	{
    		Vector<DoubleBondPositionVO> singleLabeledVOs = container.getSingleLabeledSpecies(molecularSpecies);
    		Vector<DoubleBondPositionVO> multiLabeledVOs = container.getMultiLabeledSpecies(molecularSpecies); //TODO: here we only take those that were also found as single labels, this can be expanded later
    		
    		//add multi labeled info to single label vos before clustering them.
    		if (multiLabeledVOs != null)
    		{
    			Vector<DoubleBondPositionVO> clusteredMultiLabeledVOs = clusterMolecularSpecies(multiLabeledVOs, true, clusteringThreshold);
    			for (DoubleBondPositionVO voMulti : clusteredMultiLabeledVOs)
    			{
    				double multiRT = voMulti.getExpectedRetentionTime();
    				RangeDouble range = new RangeDouble(multiRT-clusteringThreshold/60.0, multiRT+clusteringThreshold/60.0);
    				for (DoubleBondPositionVO voSingle : singleLabeledVOs)
    				{
    					double singleRT = voSingle.getExpectedRetentionTime();
    					if (range.insideRange(singleRT) && computeCombinedPattern(voMulti, voSingle).equals(voMulti.getPositionAssignmentPattern()))
    					{
    						voSingle.setChainCombination(voMulti.getChainCombination());
    					}
    				}
    			}
    		}
    		
    		Vector<DoubleBondPositionVO> clusteredSingleLabeledVOs = clusterMolecularSpecies(singleLabeledVOs, true, clusteringThreshold);
    		
    		allLabeledVOsToAdd.addAll(clusteredSingleLabeledVOs); 	
    	}
    	
      Sheet sheet = workbook.createSheet(cName);
      writeMassListForSheet(sheet, getHeaderStyle(workbook), getNumberStyle(workbook), analyteSequence.get(cName), quantObjects.get(cName), cName, allLabeledVOsToAdd);
    }
	}
	
	@SuppressWarnings("unchecked")
	private void writeMassListForSheet(
			Sheet sheet, 
			XSSFCellStyle headerStyle, 
			XSSFCellStyle numberStyle,
			Vector<String> analytes, 
			Hashtable<String,Hashtable<String,QuantVO>> quantObjects, 
			String cName,
			Vector<DoubleBondPositionVO> allLabeledVOsToAdd)
	{
		try 
		{
			Vector<Object> elModCharge = getAvailableElementsAndModificationsPlusCharge(quantObjects);
	    Vector<String> elements = (Vector<String>)elModCharge.get(0);
	    LinkedHashMap<String,String> mods = (LinkedHashMap<String,String>)elModCharge.get(1);
	    Hashtable<String,Integer> modToCharge = (Hashtable<String,Integer>)elModCharge.get(2);
	    
	    List<String> headerTitles = createHeaderTitles(elements, mods, modToCharge);
	    
	    //TODO: read this from reference file.
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
	    createHeader(sheet, headerTitles, headerStyle);
	    
	    ArrayList<String> usedQuantIDs = new ArrayList<String>();
	    int firstModColumn = headerTitles.indexOf(headerTitles.stream().filter((s)->(s).startsWith("mass(form[")).findFirst().get());
	    int rowCount = HEADER_ROW+1;
	    for (String analyte : analytes) 
	    {
	    	Cell cell;
	    	
	    	Hashtable<String,QuantVO> quantAnalytes = quantObjects.get(analyte);
	      Hashtable<String,Integer> formula = new Hashtable<String,Integer>();
	      Vector<DoubleBondPositionVO> doubleBondPositionVOs = new Vector<DoubleBondPositionVO>(); //they will be identical for all mods..
	      for (String mod : mods.keySet()) 
	      {
	        QuantVO quant = quantAnalytes.get(mod);
	        if (!isRTMapping() && quant.getDbs() > 0) //only species with double bonds are relevant
	        {
	        	for (DoubleBondPositionVO doubleBondPositionVO : allLabeledVOsToAdd)
		        {
		        	if (quant.getCarbons() == doubleBondPositionVO.getNumberOfCarbons() &&
		        			quant.getDbs() == doubleBondPositionVO.getNumberOfDoubleBonds()) //number of oh etc are not considered on purpose, they should be grouped in the target list
		        	{
		        		quant.addInfoForOmegaAssignment(doubleBondPositionVO); 
		        	}
		        }
	        }
	        doubleBondPositionVOs = quant.getInfoForOmegaAssignment();
	      }
	      
	      if (isRTMapping())
	      {
	      	RecalibrationRegression regression = calibrationGraphPanel_.getRegressionByFields(CalibrationGraphPanel.PLOT_ALL, 
      				calibrateSeparately_ ? calibrationGraphPanel_.getRelevantRegressionName(cName) : CalibrationGraphPanel.PLOT_ALL);
      		boolean noRegressionAvailable = regression == null;
	      	for (int i=0;i<doubleBondPositionVOs.size();i++)
		      {
	      		DoubleBondPositionVO vo = doubleBondPositionVOs.get(i);
	      		if (noRegressionAvailable) 
	      		{
	      			vo.setExpectedRetentionTime((float)-1.0);
	      		}
	      		else
	      		{
	      			try
	      			{
	      				computeExpectedRetentionTime(vo, regression, cName);
	      			}
	      			catch (Exception ex)
	      			{
	      				ex.printStackTrace();
	      			}
	      		}
		      }
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
		      	if ((i<0 && usedQuantIDs.contains(computeQuantID(quant, false)) ||
		      			(i>-1 && !usedQuantIDs.contains(computeQuantID(quant, true)))))
		        {
		      		rowCount--;
		        	continue;
		        }
		        usedQuantIDs.add(computeQuantID(quant, true));
		        usedQuantIDs.add(computeQuantID(quant, false));
		      	
		      	if (modCount==0) 
		        {   	
		          cell = row.createCell(headerTitles.indexOf(HEADER_NAME),HSSFCell.CELL_TYPE_STRING);
		          cell.setCellValue(quant.getCarbons());
		          cell = row.createCell(headerTitles.indexOf(HEADER_COLON),HSSFCell.CELL_TYPE_STRING);
		          cell.setCellValue(":");
		          cell = row.createCell(headerTitles.indexOf(HEADER_DBS),HSSFCell.CELL_TYPE_NUMERIC);
		          cell.setCellValue(quant.getDbs());
		          if (i>=0)
		          {
		          	cell = row.createCell(headerTitles.indexOf(HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS),HSSFCell.CELL_TYPE_NUMERIC);
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
	          	cell = row.createCell(headerTitles.indexOf(HEADER_RETENTION_TIME),HSSFCell.CELL_TYPE_NUMERIC);
		          cell.setCellValue(doubleBondPositionVOs.get(i).getExpectedRetentionTime());
		          cell.setCellStyle(numberStyle);
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
	
	public String computeQuantID(QuantVO quant, boolean includeOH)
	{
		return String.format("%s, %s, %s, %s", 
				quant.getCarbons(), quant.getDbs(), quant.getModName(), includeOH ? quant.getOhNumber() : "");
	}
	
	private void computeExpectedRetentionTime(DoubleBondPositionVO vo, RecalibrationRegression regression, String cName)
	{
		double before = vo.getExpectedRetentionTime();
		try
		{
			double after = regression.getTargetRT(before);
			vo.setExpectedRetentionTime((float)after);
//			System.out.println(String.format("original: %s, recalibrated: %s", before, after));
		}
		catch (OutOfRangeException ex)
		{
			double extrapolatedRT = extrapolateRTOutOfRange(before, regression);
			vo.setExpectedRetentionTime((float)extrapolatedRT);
		}
		
		//TODO: just for development, remove after!
		if (!beforeAfter_.containsKey(cName))
		{
			beforeAfter_.put(cName, ConcurrentHashMap.newKeySet());
		}
		beforeAfter_.get(cName).add(new Pair<Double,DoubleBondPositionVO>(before,vo));
	}
	
	/**
	 * Extrapolates calibrated retention times when the values are out of the range of calibration.
	 * @param before
	 * @param regression
	 * @return
	 */
	private double extrapolateRTOutOfRange(double before, RecalibrationRegression regression)
	{
		Double returnValue = -1.0;
		ArrayList<Pair<Double,Double>> clustered = regression.getClustered();
		Pair<Double,Double> firstPair = null;
		Pair<Double,Double> secondPair = null;
		if (before < clustered.get(0).getKey())
		{
			firstPair = clustered.get(0);
			secondPair = clustered.get(1);
		}
		else if (before > clustered.get(clustered.size()-1).getKey())
		{
			firstPair = clustered.get(clustered.size()-1);
			secondPair = clustered.get(clustered.size()-2);
		}
		if (firstPair != null && secondPair != null)
		{
			double m = (secondPair.getValue() - firstPair.getValue()) / (secondPair.getKey() - firstPair.getKey());
			double b = firstPair.getValue()-m*firstPair.getKey();
			returnValue = before - (m*before+b);
		}
		System.out.println(returnValue);
		return returnValue;
	}
	
	/**
	 * first all vos are sorted according to the labeling pattern (they already all have the same molecular species) into objects
	 * then for the vos within these objects the mean and standard deviation is calculated
	 * if standard deviation is below our threshold, then this is an accepted cluster, we save it as a new vo
	 * when we are done, the final step is to take all obtained doublebondpositionvos and check if another dbpvo with a compatible labeling pattern is within the time threshold, if so, merge
	 * TODO: in this class we get only vos of identical molecular species. this means we can compute permutable positions for this molecular species and use this to compute overlapping patterns!
	 * 
	 * @param doubleBondPositionVOs
	 * @param labeled
	 * @param clusteringThreshold
	 * @return
	 */
	private Vector<DoubleBondPositionVO> clusterMolecularSpecies(Vector<DoubleBondPositionVO> doubleBondPositionVOs, boolean labeled, double clusteringThreshold)
	{
		/**
		 * if we have labeled species:
		 * we first get the positionassignmentpattern for each : hashtable key pattern value Vector of vos
		 * each pattern is clustered
		 * then matching patterns are checked for overlaps, if so the respective vos are combined to one new vo with higher assignment
		 */
		if (labeled && !doubleBondPositionVOs.isEmpty())
		{
			Vector<DoubleBondPositionVO> returnElements = new Vector<DoubleBondPositionVO>();
			Vector<Boolean> permutationPattern = doubleBondPositionVOs.get(0).getPermutationPattern();
			
			Hashtable<Vector<Integer>,Vector<DoubleBondPositionVO>> patternLookup = computePatternLookup(doubleBondPositionVOs);
			for (Vector<Integer> pattern : patternLookup.keySet()) //cluster elements for each pattern
			{
				Vector<DoubleBondPositionVO> averageElements = computeGroups(patternLookup.get(pattern),clusteringThreshold);
				patternLookup.put(pattern, averageElements);
			}
			
			/**
			 * iterate over markedToCombine (they are within threshold and have an overlapping pattern) and combine to a common vo 
			 */
			
			Hashtable<Vector<Integer>,Vector<Vector<Integer>>> overlappingPatterns = computeOverlappingPatterns(patternLookup.keySet(),permutationPattern);
			
			Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>> markedToCombine = markToCombine(patternLookup, overlappingPatterns,clusteringThreshold);
			Vector<Pair<RangeDouble,Vector<DoubleBondPositionVO>>> combineRanges = new Vector<Pair<RangeDouble,Vector<DoubleBondPositionVO>>>();
			
			for (DoubleBondPositionVO doubleBondPositionVO : markedToCombine.keySet())
			{
				boolean isRepresented = false;
				for(Pair<RangeDouble,Vector<DoubleBondPositionVO>> range1 : combineRanges)
				{
					if (range1.getKey().insideRange(doubleBondPositionVO.getExpectedRetentionTime()))
					{
						RangeDouble range2 = new RangeDouble(doubleBondPositionVO.getExpectedRetentionTime()-clusteringThreshold/60.0, doubleBondPositionVO.getExpectedRetentionTime()+clusteringThreshold/60.0);
						range1.getKey().extendToOtherRanges(range2);
						range1.getValue().add(doubleBondPositionVO);
						isRepresented = true;
					}
				}
				if (!isRepresented)
				{
					RangeDouble range2 = new RangeDouble(doubleBondPositionVO.getExpectedRetentionTime()-clusteringThreshold/60.0, doubleBondPositionVO.getExpectedRetentionTime()+clusteringThreshold/60.0);
					Vector<DoubleBondPositionVO> newVector = new Vector<DoubleBondPositionVO>();
					newVector.add(doubleBondPositionVO);
					Pair<RangeDouble,Vector<DoubleBondPositionVO>> newRange = new Pair<RangeDouble,Vector<DoubleBondPositionVO>>(range2,newVector);
					combineRanges.add(newRange);
				}
			}
			
			for (Pair<RangeDouble,Vector<DoubleBondPositionVO>> range : combineRanges)
			{
				boolean combineSpecies = true;
				Vector<Vector<Integer>> added = new Vector<Vector<Integer>>();
				Vector<DoubleBondPositionVO> toBeRemovedFromLookup = new Vector<DoubleBondPositionVO>();
				Vector<DoubleBondPositionVO> toBeAddedToLookup = new Vector<DoubleBondPositionVO>();
				for (DoubleBondPositionVO doubleBondPositionVO : range.getValue())
				{
					Vector<DoubleBondPositionVO> vos = new Vector<DoubleBondPositionVO>();
					vos.addAll(markedToCombine.get(doubleBondPositionVO));
					Set<Vector<Integer>> combinedPatterns = computeCombinedPatterns(doubleBondPositionVO, vos);
					if (combinedPatterns.size() != 1)
					{
						combineSpecies = false;
					}
					else
					{
						Vector<Integer> pattern = combinedPatterns.iterator().next();
						Vector<FattyAcidVO> original = doubleBondPositionVO.getChainCombination(); //original must remain unchanged!
						Vector<FattyAcidVO> combinedFattyAcidVO = new Vector<FattyAcidVO>();
			  		for (int i=0;i<pattern.size();i++)
			  		{
			  			FattyAcidVO fattyAcid = new FattyAcidVO(original.get(i));
			  			fattyAcid.setOmegaPosition(pattern.get(i));
			  			combinedFattyAcidVO.add(fattyAcid);
			  		}
						double expectedRetentionTime = doubleBondPositionVO.getExpectedRetentionTime();
						int count = 1;
						for (DoubleBondPositionVO vo : vos)
						{
							toBeRemovedFromLookup.add(vo);
							expectedRetentionTime+= vo.getExpectedRetentionTime();
							count++;
						}
						DoubleBondPositionVO combined = new DoubleBondPositionVO(
								combinedFattyAcidVO, expectedRetentionTime/count, DoubleBondPositionVO.ACCURACY_LOW, doubleBondPositionVO.getMolecularSpecies());
						toBeAddedToLookup.add(combined);
						toBeRemovedFromLookup.add(doubleBondPositionVO);
					}
					
				}
				if (combineSpecies == true)
				{
					for (DoubleBondPositionVO toBeAdded : toBeAddedToLookup)
					{
						if (!added.contains(toBeAdded.getPositionAssignmentPattern()))
						{
							if (!patternLookup.containsKey(toBeAdded.getPositionAssignmentPattern()))
							{
								patternLookup.put(toBeAdded.getPositionAssignmentPattern(), new Vector<DoubleBondPositionVO>());
							}
							patternLookup.get(toBeAdded.getPositionAssignmentPattern()).add(toBeAdded);
							added.add(toBeAdded.getPositionAssignmentPattern());
						}
					}
					for (DoubleBondPositionVO toBeRemoved : toBeRemovedFromLookup)
					{
						patternLookup.get(toBeRemoved.getPositionAssignmentPattern()).remove(toBeRemoved);
					}
				}
			}
			for (Vector<Integer> pattern : patternLookup.keySet()) 
			{
				returnElements.addAll(patternLookup.get(pattern));
			}
			return returnElements;
		}
		else
		{
			/**
			 * we calculate the clusters, then these are transformed into a Vector of vos and returned
			 */	
			return computeGroups(doubleBondPositionVOs,clusteringThreshold);
		}
	}
	
	
	private Set<Vector<Integer>> computeCombinedPatterns(DoubleBondPositionVO vo1, Vector<DoubleBondPositionVO> doubleBondPositionVOs)
	{
		Set<Vector<Integer>> combinedPatterns = new HashSet<Vector<Integer>>();
		for (DoubleBondPositionVO vo2 : doubleBondPositionVOs)
		{
			Vector<Integer> combinedPattern = computeCombinedPattern(vo1, vo2);
			combinedPatterns.add(combinedPattern);
		}
		return combinedPatterns;
	}
	
	private Vector<Integer> computeCombinedPattern(DoubleBondPositionVO vo1, DoubleBondPositionVO vo2)
	{
		Vector<Integer> combinedPattern = new Vector<Integer>();
		Vector<Integer> pattern1 = vo1.getPositionAssignmentPattern();
		Vector<Integer> pattern2 = vo2.getPositionAssignmentPattern();
		if (vo2.getPermutationPattern().contains(Boolean.TRUE) && pattern2.size() == 2) //TODO: adapt for more than 2 chains
		{
			int temp = pattern2.get(0);
			pattern2.setElementAt(pattern2.get(1), 0);
			pattern2.setElementAt(temp, 1);
		}
		for (int i=0;i<pattern1.size();i++)
		{
			if (pattern1.get(i) == pattern2.get(i)) {combinedPattern.add(pattern1.get(i));}
			else if ((pattern1.get(i) > 0) && (pattern2.get(i) == -1)) {combinedPattern.add(pattern1.get(i));}
			else if ((pattern1.get(i) == -1) && (pattern2.get(i) > 0)) {combinedPattern.add(pattern2.get(i));}
			else {combinedPattern.add(-1);}
		}
		return combinedPattern;
	}
	
	
	private Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>> markToCombine(
			Hashtable<Vector<Integer>,Vector<DoubleBondPositionVO>> patternLookup,
			Hashtable<Vector<Integer>,Vector<Vector<Integer>>> overlappingPatterns,
			Double clusteringThreshold)
	{
		Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>> markedToCombine = new Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>>();
		for (Vector<Integer> pattern : overlappingPatterns.keySet()) //combine overlapping patterns with same retention times
		{
			Vector<DoubleBondPositionVO> patternVOs = patternLookup.get(pattern); //all species matching the pattern exactly
			for (Vector<Integer> overlappingPattern : overlappingPatterns.get(pattern)) //pattern overlap with the queried pattern
			{
				Vector<DoubleBondPositionVO> overlappingPatternVOs = patternLookup.get(overlappingPattern); //species with pattern overlap with the queried pattern
				for (DoubleBondPositionVO vo1 : patternVOs)
				{			
//					System.out.println("VO1 RT: "+vo1.getExpectedRetentionTime()+" pattern: "+vo1.getPositionAssignmentPattern());
					Float retentionTime1 = (float)vo1.getExpectedRetentionTime();
					for (DoubleBondPositionVO vo2 : overlappingPatternVOs)
					{
//						System.out.println("   > VO2 RT: "+vo2.getExpectedRetentionTime()+" pattern: "+vo2.getPositionAssignmentPattern());
						Float retentionTime2 = (float)vo2.getExpectedRetentionTime();
//							System.out.println(String.format("Species1: %s, RT1: %s, Species2: %s, RT2: %s ", 
//									vo1.getDoubleBondPositionsHumanReadable(), retentionTime1, vo2.getDoubleBondPositionsHumanReadable(), retentionTime2));
						if (Math.abs(retentionTime1-retentionTime2) < (clusteringThreshold/60)) //if is within threshold
						{
							
							if (!markedToCombine.containsKey(vo1))
							{
								markedToCombine.put(vo1, new Vector<DoubleBondPositionVO>());
							}
							if (!markedToCombine.get(vo1).contains(vo2)) 
							{
								markedToCombine.get(vo1).add(vo2);
							}
							
							if (!markedToCombine.containsKey(vo2))
							{
								markedToCombine.put(vo2, new Vector<DoubleBondPositionVO>());
							}
							if (!markedToCombine.get(vo2).contains(vo1)) 
							{
								markedToCombine.get(vo2).add(vo1);
							}
						}
					}
				}
			}				
		}
		return removeDuplicates(markedToCombine);
	}
	
	/** 
	 * Removes a collection in markedToCombine, if it is a subcollection of another one.
	 * A partial overlap of these collections is in certain cases needed, but full overlaps need to be merged.
	 * @param markToCombine
	 * @return
	 */
	private Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>> removeDuplicates(Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>> markedToCombine)
	{	
		Set<DoubleBondPositionVO> keySet = markedToCombine.keySet();
		List<DoubleBondPositionVO> toRemoveList = new ArrayList<DoubleBondPositionVO>();
		for (DoubleBondPositionVO key : keySet)
		{
			if (!toRemoveList.contains(key))
			{
				Vector<DoubleBondPositionVO> collection1 = new Vector<DoubleBondPositionVO>();
				collection1.addAll(markedToCombine.get(key));
				collection1.add(key);
				Vector<DoubleBondPositionVO> keyVector = markedToCombine.get(key);
				for (int i=0; i<keyVector.size(); i++)
				{
					if (!keyVector.get(i).equals(key) 
							&& !toRemoveList.contains(keyVector.get(i)) 
							&& markedToCombine.containsKey(keyVector.get(i)))
					{
						Vector<DoubleBondPositionVO> collection2 = new Vector<DoubleBondPositionVO>();
						collection2.addAll(markedToCombine.get(keyVector.get(i)));
						collection2.add(keyVector.get(i));
						if (collection2.containsAll(collection1))
						{
							toRemoveList.add(key);
						}
						else if (collection1.containsAll(collection2))
						{
							toRemoveList.add(keyVector.get(i));
						}
					}
				}
			}
		}
		for (DoubleBondPositionVO toRemove : toRemoveList)
		{
			markedToCombine.remove(toRemove);
		}
		return markedToCombine;
	}
	
	
	/**
	 * Computes strictly just overlapping patterns, whether species can overlap due to their RT is not considered
	 * @param patterns
	 * @param permutationPattern required to find overlaps with permutable patterns (this is the case, when the molecule contains several identical fatty acids)
	 * @return
	 */
	private Hashtable<Vector<Integer>,Vector<Vector<Integer>>> computeOverlappingPatterns(Set<Vector<Integer>> patterns,Vector<Boolean> permutationPattern)
	{
		Hashtable<Vector<Integer>,Vector<Vector<Integer>>> overlappingPatterns = new Hashtable<Vector<Integer>,Vector<Vector<Integer>>>();
		for (Vector<Integer> pattern1 : patterns) //computing overlapping patterns
		{
			overlappingPatterns.put(pattern1, new Vector<Vector<Integer>>());
			for (Vector<Integer> pattern2 : patterns)
			{ 
				boolean isOverlapping = true;
				boolean isIdentical = pattern1.equals(pattern2);
				if (!isIdentical)
				{
					Vector<Vector<Integer>> patternPermutations = computePatternPermutations(pattern2, permutationPattern);
					for (Vector<Integer> pattern3 : patternPermutations)
					{
						isOverlapping = true;
						for (int i=0; i<pattern1.size(); i++)
						{
							if (!(pattern1.get(i).equals(pattern3.get(i)) || pattern1.get(i) == -1 || pattern3.get(i) == -1)) //if neither of the patterns has the differing position unassigned, the pattern cannot overlap
							{
								isOverlapping = false;
							} 
						}
					}
					if (isOverlapping) 
					{
						overlappingPatterns.get(pattern1).add(pattern2);
					}
				}
			}
		}
		return overlappingPatterns;
	}
	
	/**
	 * Computes all possible permutations for a given double bond position pattern and a permutation pattern
	 * TODO: if more than three fatty acid chains are involved, this method needs an extension
	 * @param pattern
	 * @param permutationPattern
	 * @return
	 */
	private Vector<Vector<Integer>> computePatternPermutations(Vector<Integer> pattern, Vector<Boolean> permutationPattern)
	{
		Vector<Vector<Integer>> patternPermutations = new Vector<Vector<Integer>>();
		
		int size = pattern.size();
		patternPermutations.add(new Vector<Integer>(pattern));
		if (!permutationPattern.contains(Boolean.TRUE))
		{
			return patternPermutations;
		}
		else
		{
			for (int i=0;i<size;i++)
			{
				if (permutationPattern.get(i) == Boolean.TRUE)
				{
					if (i<size-1 && permutationPattern.get(i+1) == Boolean.TRUE)
					{
						Vector<Integer> tempPattern = new Vector<Integer>(pattern);
						int temp = tempPattern.get(i);
						tempPattern.setElementAt(tempPattern.get(i+1), i);
						tempPattern.setElementAt(temp, i+1);
						patternPermutations.add(tempPattern);
					}
					else if (size>2 && i==size-1 && permutationPattern.firstElement() == Boolean.TRUE)
					{
						Vector<Integer> tempPattern = new Vector<Integer>(pattern);
						int temp = tempPattern.get(i);
						tempPattern.setElementAt(tempPattern.firstElement(), i);
						tempPattern.setElementAt(temp, 0);
						patternPermutations.add(tempPattern);
					}
				}
			}
		}
		return patternPermutations;
	}
	
	
	public DoubleBondPositionVO combineDoubleBondPositionVO(DoubleBondPositionVO vo1, DoubleBondPositionVO vo2)
  {
  	if (vo1.getMolecularSpecies().equals(vo2.getMolecularSpecies()))
  	{
  		double combinedRetentionTime = (vo1.getExpectedRetentionTime()+vo2.getExpectedRetentionTime())/2;
  		Vector<Integer> pattern1 = vo1.getPositionAssignmentPattern();
  		Vector<Integer> pattern2 = vo2.getPositionAssignmentPattern();
  		if (vo2.getPermutationPattern().contains(Boolean.TRUE) && pattern2.size() == 2) //TODO: adapt for more than 2 chains
			{
				int temp = pattern2.get(0);
				pattern2.setElementAt(pattern2.get(1), 0);
				pattern2.setElementAt(temp, 1);
			}
  		Vector<Integer> combinedPattern = new Vector<Integer>();
  		Vector<FattyAcidVO> combinedFattyAcidVO = new Vector<FattyAcidVO>(vo1.getChainCombination());
  		for (int i=0;i<pattern1.size();i++)
  		{
  			if (pattern1.get(i) == pattern2.get(i)) {combinedPattern.add(pattern1.get(i));}
  			else if ((pattern1.get(i) > 0) && (pattern2.get(i) == -1)) {combinedPattern.add(pattern1.get(i));}
  			else if ((pattern1.get(i) == -1) && (pattern2.get(i) > 0)) {combinedPattern.add(pattern2.get(i));}
  			else {combinedPattern.add(-1);}
  		}
  		for (int i=0;i<combinedPattern.size();i++)
  		{
  			combinedFattyAcidVO.get(i).setOmegaPosition(combinedPattern.get(i));
  		}
  		return new DoubleBondPositionVO(combinedFattyAcidVO, combinedRetentionTime, DoubleBondPositionVO.ACCURACY_LOW, vo1.getMolecularSpecies());
  	}
  	return null;
  }
	
	
	/**
	 * Creates a lookup table. 
	 * Keys: double bond position labeling pattern.
	 * Values: VOs with this labeling pattern.
	 * @param doubleBondPositionVOs
	 * @return
	 */
	private Hashtable<Vector<Integer>,Vector<DoubleBondPositionVO>> computePatternLookup(Vector<DoubleBondPositionVO> doubleBondPositionVOs)
	{
		Hashtable<Vector<Integer>,Vector<DoubleBondPositionVO>> patternLookup = new Hashtable<Vector<Integer>,Vector<DoubleBondPositionVO>>();
		for (DoubleBondPositionVO doubleBondPositionVO : doubleBondPositionVOs) //dividing vos according to their db assignment pattern
		{
			Vector<Integer> pattern = doubleBondPositionVO.getPositionAssignmentPattern();
			if (!patternLookup.containsKey(pattern))
			{
				patternLookup.put(pattern, new Vector<DoubleBondPositionVO>());
			}
			patternLookup.get(pattern).add(doubleBondPositionVO);
		}
		return patternLookup;
	}
	
	/**
	 * Groups identifications based on the provided @param clusteringThreshold
	 * @param doubleBondPositionVOs
	 * @param clusteringThreshold
	 * @return
	 */
	private Vector<DoubleBondPositionVO> computeGroups(Vector<DoubleBondPositionVO> doubleBondPositionVOs, double clusteringThreshold)
	{
		Vector<DoubleBondPositionVO> averageElements = new Vector<DoubleBondPositionVO>();
		//first index of vo, second number of cluster
		Hashtable<Integer,Integer> clusterLookup = new Hashtable<Integer,Integer>();
		//first number of cluster, second cluster
		Hashtable<Integer,Cluster> clusters = new Hashtable<Integer,Cluster>();
		int numClusters = 0;
		for (int i=0; i<doubleBondPositionVOs.size();i++)
		{
			DoubleBondPositionVO vo1 = doubleBondPositionVOs.get(i);
			Float retentionTime1 = (float)vo1.getExpectedRetentionTime();
			if (!clusterLookup.containsKey(i)) //if this vo is in no cluster yet: new cluster
			{
				Cluster cluster = new Cluster();
				cluster.addDoubleBondPositionVO(vo1);
				clusters.put(numClusters, cluster);
				clusterLookup.put(i, numClusters);
				numClusters++;
			}
			for (int j=0; j<doubleBondPositionVOs.size();j++)
			{
				DoubleBondPositionVO vo2 = doubleBondPositionVOs.get(j);
				if (Math.abs(retentionTime1-(float)vo2.getExpectedRetentionTime()) < (clusteringThreshold/60)) //if is within threshold
				{
					if (!clusterLookup.containsKey(j)) //if this vo is in no cluster yet: add to cluster of vo1
					{
						int numCluster = clusterLookup.get(i);
						clusters.get(numCluster).addDoubleBondPositionVO(vo2);
						clusterLookup.put(j, numCluster);
					}
					else 
					{
						int numCluster1 = clusterLookup.get(i);
						int numCluster2 = clusterLookup.get(j);
						if (numCluster1 != numCluster2) //if they are not in the same cluster: merge clusters
						{
							Cluster cluster = new Cluster(clusters.get(numCluster1), clusters.get(numCluster2));
							clusters.put(numCluster1, cluster);
							clusterLookup.forEach((k, v) -> {if (v == numCluster2) {clusterLookup.replace(k, numCluster1);}}); //replace numCluster2
							clusters.remove(numCluster2); //numCluster2 is obsolete
						}
					}
				}
			}
		}	
		clusters.forEach((k, v) -> {averageElements.add(v.getAverageElement());});
		return averageElements;
	}
	
	private void fillContainerForClass(String cName, MolecularSpeciesContainer container, GradientAdjuster adjuster)
	{
		try
		{
			for (ResultFileVO resultFileVO : resultFileVO_)
			{
				QuantificationResult quantRes = resultFileVO.getQuantificationResult();
				Vector<LipidParameterSet> analytes = quantRes.getIdentifications().get(cName);
				
				for (LipidParameterSet analyte : analytes)
				{
					//TODO: if only one chain for a class, MSn is not necessary. However, some methods would have to be rewritten and identifications are usually more accurate, if there is MSn info available
					if (analyte instanceof LipidomicsMSnSet
							&& analyte.getDoubleBonds()>0 )
					{
						LipidomicsMSnSet analyteMSn = (LipidomicsMSnSet) analyte;
						if (analyteMSn.getStatus() < LipidomicsMSnSet.FRAGMENTS_DETECTED) 
						{
							continue; //it does not make sense to write anything less than a molecular species with double bond information to the masslist
						}
						Vector<Pair<String,String>> labeledUnlabeledPairs = analyteMSn.getLabeledUnlabeledPairs();
						
						double expectedRetentionTime = analyte.getPreciseRT();
						boolean viableForCalibration = true;
						boolean singleLabel = true;
						
						//only deuterated species need recalculation of the expected retention time
						if (isotopeEffectRegression_ != null && analyte.getChemicalFormula().contains("D"))
						{
							int numberDeuterium = StaticUtils.categorizeFormula(analyte.getChemicalFormula()).get("D");
							if (numberDeuterium > isotopeEffectRegression_.getMaxNumDeuteriumAllowed())
							{
								singleLabel = false;
								String validCombi = analyteMSn.getValidChainCombinations().get(0);
								for (FattyAcidVO fa : StaticUtils.decodeLipidNamesFromChainCombi(validCombi))
								{
									String prefix = fa.getPrefix();
									if (prefix != null && !prefix.equals(""))
									{
										for (IsotopeLabelVO label : labels_)
										{
											if (label.getLabelId().equals(prefix))
											{
												int num = label.getLabelElements().containsKey("D") ? label.getLabelElements().get("D") : 0;
												if (num > isotopeEffectRegression_.getMaxNumDeuteriumAllowed())
												{
													viableForCalibration = false; //retention times affected by labels outside the calibration curve cannot be recalculated!
												}
												else
												{
													expectedRetentionTime = adjuster != null ? 
															adjuster.getGradientAdjustedValue(isotopeEffectRegression_.getIsotopeEffect(num), expectedRetentionTime) :
															isotopeEffectRegression_.getIsotopeEffect(num);
												}
												continue;
											}
										}
									}
					      }
							}
							else
							{
								expectedRetentionTime = adjuster != null ? 
										adjuster.getGradientAdjustedValue(isotopeEffectRegression_.getIsotopeEffect(numberDeuterium), analyte.getPreciseRT()) :
										isotopeEffectRegression_.getIsotopeEffect(numberDeuterium);
							}
						}					
						
						if (!viableForCalibration) continue;
						
						for (Pair<String,String> pair : labeledUnlabeledPairs) 
						{
							String labeledSpecies = pair.getValue();
							Vector<FattyAcidVO> chains = new Vector<FattyAcidVO>();
			    		String[] splitName = StaticUtils.splitChainCombinationsAtChainSeparators(labeledSpecies);
			    		for (int i=0; i<splitName.length;i++)
			    		{
			    			FattyAcidVO fa = StaticUtils.decodeHumanReadableChain(splitName[i],Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(),false,null);
			    			String prefix = fa.getPrefix();
			    			
			    			if (!prefix.equals("")) //unlabeled species do not need this loop
			    			{
			    				for (IsotopeLabelVO label : labels_)
				    			{
				    				if (label.getLabelId().equals(prefix))
				    				{
				    					fa.setOmegaPosition(label.getOmegaPosition());
				    					fa.setPrefix("");
				    				}
				    			}
			    			}
			    			
			    			chains.add(fa);
			    		}
			    		
			    		if (analyte.getChemicalFormula().contains("D") || analyte.getChemicalFormula().contains("Cc")) //TODO: consider accounting for additional isotopes?
							{
			    			if (singleLabel)
			    			{
			    				container.addSingleLabeledSpecies(
					    				pair.getKey(), 
					    				new DoubleBondPositionVO(chains,expectedRetentionTime,DoubleBondPositionVO.ACCURACY_LOW,pair.getKey()));
			    			}
			    			else
			    			{
			    				container.addMultiLabeledSpecies(
					    				pair.getKey(), 
					    				new DoubleBondPositionVO(chains,expectedRetentionTime,DoubleBondPositionVO.ACCURACY_LOW,pair.getKey()));
			    			}
							}
						}
					}
				}
			}
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}	
	
	
	
	/**
   * Creates a formatted header row with the given header titles in the row number given by HEADER_ROW
   * @param ws Excel worksheet to write the header to
   * @param headerTitles List of header titles
   */
  public static void createHeader(Sheet sheet, List<String> headerTitles, XSSFCellStyle headerStyle) 
  {
  	Row row = sheet.createRow(HEADER_ROW);
  	Cell cell;
  	for (int i=0; i<headerTitles.size(); i++) 
  	{
  		cell = row.createCell(i, HSSFCell.CELL_TYPE_STRING);
  		cell.setCellValue(headerTitles.get(i));
  		cell.setCellStyle(headerStyle);
  	}
    sheet.setColumnWidth(headerTitles.indexOf(HEADER_NAME), 10 * 256);
    sheet.setColumnWidth(headerTitles.indexOf(HEADER_COLON), 2 * 256);
    sheet.setColumnWidth(headerTitles.indexOf(HEADER_DBS), 10 * 256);
  	sheet.setColumnWidth(headerTitles.indexOf(HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS), 30 * 256);
  	sheet.setColumnWidth(headerTitles.indexOf(HEADER_RETENTION_TIME), 10 * 256);
  }
	
	/**
	 * Creates a list of header titles
	 * @param elements
	 * @param mods
	 * @return
	 */
  public static List<String> createHeaderTitles(
  		Vector<String> elements, LinkedHashMap<String,String> mods, Hashtable<String,Integer> modToCharge) 
  {
    List<String> headerTitles = new ArrayList<String>();
    headerTitles.add(HEADER_NAME);
    headerTitles.add(HEADER_COLON);
    headerTitles.add(HEADER_DBS);
    headerTitles.add(HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS);
    
    for (String element : elements) 
    {
    	headerTitles.add(element);
    }
    
    for (String mod : mods.keySet()) 
    {
    	String modHeader = "mass(form["+mods.get(mod)+"] name["+mod+"]";
    	if (modToCharge.get(mod)>1)
        modHeader += " charge="+modToCharge.get(mod);
      modHeader += ")";
    	headerTitles.add(modHeader);
    }
    
    headerTitles.add(HEADER_RETENTION_TIME);
    
    return headerTitles;
  }
	
  public static XSSFCellStyle getNumberStyle(XSSFWorkbook wb)
  {
  	XSSFCellStyle numberStyle = wb.createCellStyle();
  	numberStyle.setDataFormat(2);
  	return numberStyle;
  }
  
	public static XSSFCellStyle getHeaderStyle(XSSFWorkbook wb)
	{
    XSSFCellStyle arial12style = wb.createCellStyle();
    XSSFFont arial12font = wb.createFont();
    arial12font.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);
    arial12font.setFontName("Arial");
    arial12font.setFontHeightInPoints((short)12);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(HSSFCellStyle.ALIGN_CENTER);
    return arial12style;
  }
	
	/**
   * extracts from the original mass list the chemical elements used, the modifications applied and the charge of those modifications
   * @param quantObjects the information form the original Excel mass list; first key: analyte name; second key modification name; value: information about the analyte entry
   * @return vector containing three objects: first object: a vector containing the elements in the order of the Hill notation; second object: a LinkedHashMap containing the modifications as key and the modification formula as value; third object: a hash table; key: modification name; value: charge of the modification 
   * @throws ChemicalFormulaException thrown when there is something wrong with the chemical formula
   */
  public Vector<Object> getAvailableElementsAndModificationsPlusCharge(Hashtable<String,Hashtable<String,QuantVO>> quantObjects) throws ChemicalFormulaException
  {
    Vector<Object> result = new Vector<Object>();
    Set<String> elements = new HashSet<String>();
    LinkedHashMap<String,String> modifications = new LinkedHashMap<String,String>();
    Hashtable<String,Integer> modToCharge = new Hashtable<String,Integer>();
    for (Hashtable<String,QuantVO> quantAnal : quantObjects.values()) {
      for (QuantVO quant : quantAnal.values()) {
        for (String element : StaticUtils.categorizeFormula(quant.getAnalyteFormula()).keySet()) {
          if (!elements.contains(element))
            elements.add(element);
          if (!modifications.containsKey(quant.getModName())) {
            modifications.put(quant.getModName(),StaticUtils.getFormulaInHillNotation(StaticUtils.categorizeFormula(quant.getModFormula()),false));
            modToCharge.put(quant.getModName(), quant.getCharge());
          }
        }
      }
    }
    Vector<String> sorted = new Vector<String>();
    if (elements.contains("C"))
      sorted.add("C");
    if (elements.contains("H"))
      sorted.add("H");
    List<String> otherThanCH = new ArrayList<String>();
    for (String element : elements) {
      if (!element.equalsIgnoreCase("C") && !element.equalsIgnoreCase("H"))
        otherThanCH.add(element);
    }
    Collections.sort(otherThanCH);
    for (String element : otherThanCH)
      sorted.add(element);
    result.add(sorted);
    result.add(modifications);
    result.add(modToCharge);
    return result;
  }
  
  private String getTemplatePath()
  {
  	return this.templatePath_;
  }
  
  private boolean isRTMapping()
	{
		return this.calibrationGraphPanel_ != null;
	}
  
  //TODO: just for development, remove later!
  private Hashtable<String,Set<Pair<Double,DoubleBondPositionVO>>> getBeforeAfter()
	{
		return this.beforeAfter_;
	}
  
  //TODO: just for development, remove later!
  private Hashtable<String,Vector<Pair<Double,DoubleBondPositionVO>>> parseComparisonSheet(org.dhatim.fastexcel.reader.Sheet sheet) throws Exception
  {
  	Hashtable<String,Vector<Pair<Double,DoubleBondPositionVO>>> allowedPairsOfClass = new Hashtable<String,Vector<Pair<Double,DoubleBondPositionVO>>>();
  	List<org.dhatim.fastexcel.reader.Row> rows = null;
  	rows = sheet.read();
  	Integer rowNr = 0;
    org.dhatim.fastexcel.reader.Row headerRow = rows.get(rowNr++);
    List<String> headerTitles = null;
		try (Stream<org.dhatim.fastexcel.reader.Cell> cells = headerRow.stream();) {
			headerTitles = cells.map((c) -> (!(c==null || c.getType().equals(CellType.ERROR)) ? c.getText() : "null")).collect(Collectors.toList());
    }
		if (headerTitles == null) throw new IOException("No headertitles...");
		
    List<org.dhatim.fastexcel.reader.Row> contentRows = rows.subList(rowNr, rows.size());
    
    String lClass = null;
    String molecularSpecies = null;
    float targetRT = 0f;
    float originalRT = 0f;
    int index;
    String rawValue;
    
    for (org.dhatim.fastexcel.reader.Row row : contentRows) {
      List<org.dhatim.fastexcel.reader.Cell> cells = row.stream().filter((c) -> !(c==null || c.getType().equals(CellType.ERROR))).collect(Collectors.toList());
      for (org.dhatim.fastexcel.reader.Cell cell : cells) {
        index = cell.getColumnIndex();
        rawValue = cell.getRawValue();
        
        if (index == headerTitles.indexOf("Lipid Class")) {
        	lClass = rawValue;
        } else if (index == headerTitles.indexOf("Lipid Molecular Species")) {
          molecularSpecies = rawValue;
        } else if (index == headerTitles.indexOf("RT Target DB /min")) {
        	targetRT = Float.parseFloat(rawValue);
        } else if (index == headerTitles.indexOf("RT Original DB /min")) {
        	originalRT = Float.parseFloat(rawValue);
        }
        
      }
      Vector<FattyAcidVO> chainCombination = StaticUtils.decodeFAsFromHumanReadableName(
      		molecularSpecies, Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), false, null);
      DoubleBondPositionVO vo = new DoubleBondPositionVO(chainCombination, targetRT, 0, molecularSpecies);
      if (!allowedPairsOfClass.containsKey(lClass))
      {
      	allowedPairsOfClass.put(lClass, new Vector<Pair<Double,DoubleBondPositionVO>>());
      }
      allowedPairsOfClass.get(lClass).add(new Pair<Double,DoubleBondPositionVO>(new Double(originalRT), vo));
    }
  	return allowedPairsOfClass;
  }
  
  /**
   * 
   * @param targetPath		file created with the new conditions.
   * @param outPath
   * @param outPath				if there is a comparison template to adhere to
   * @throws ExportException
   */
  //TODO: just for development, remove later!
  @SuppressWarnings("unchecked")
	public void exportBeforeAfter(String targetPath, String outPath, String comparisonPath) throws ExportException
  {
  	comparisonPairsOfClass_ = new Hashtable<String,Vector<Pair<Double,DoubleBondPositionVO>>>();
  	String sheetName = "RTDB_A_to_B1";
  	if (comparisonPath != null)
  	{
  		try (InputStream is = new FileInputStream(comparisonPath);
  				ReadableWorkbook wb = new ReadableWorkbook(is);
          Stream<org.dhatim.fastexcel.reader.Sheet> sheets = wb.getSheets();) {
        sheets.forEach((s) -> {
        				if (s.getName().equalsIgnoreCase(sheetName)) 
        				{
        					try {
                  	comparisonPairsOfClass_ = parseComparisonSheet(s);
                  } catch (Exception ex) {
                    new WarningMessage(new JFrame(), "ERROR", ex.getMessage());
                  }
        				}
              });
        
      } catch (IOException ex){
        ex.printStackTrace();
      }
  	}
  	
  	
  	try (	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outPath));
				XSSFWorkbook workbook = new XSSFWorkbook();)
  	{
  		Sheet sheet = workbook.createSheet("beforeAfter");
  		Cell cell;
  		createBeforeAfterTitle(sheet);
  		int rowCount = 1;
  		@SuppressWarnings("rawtypes")
			Vector quantInfo = QuantificationThread.getCorrectAnalyteSequence(targetPath,false);
			LinkedHashMap<String,Integer> classSequence = (LinkedHashMap<String,Integer>)quantInfo.get(0);
			LinkedHashMap<String,Vector<String>> analyteSequence = (LinkedHashMap<String,Vector<String>>)quantInfo.get(1);
			Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo.get(4);
	  	for (String cName : classSequence.keySet())
			{
	  		//Key is the pair from the original target list (orig. rt, recalibrated), value from the target (rtError, target vo).
	  		Hashtable<Pair<Double,DoubleBondPositionVO>, Pair<Double,DoubleBondPositionVO>> matchedRecalibratedToTarget = 
	  				new Hashtable<Pair<Double,DoubleBondPositionVO>, Pair<Double,DoubleBondPositionVO>>();
	  		//To ensure vos are not matched twice.
	  		Hashtable<DoubleBondPositionVO,Pair<Double,Pair<Double,DoubleBondPositionVO>>> matchedTargetToRecalibrated =
	  				new Hashtable<DoubleBondPositionVO,Pair<Double,Pair<Double,DoubleBondPositionVO>>>();
	  		
	  		Vector<Pair<Double,DoubleBondPositionVO>> pairsOfClassPriority = new Vector<Pair<Double,DoubleBondPositionVO>>();
	  		if (!comparisonPairsOfClass_.isEmpty())
	  		{
	  			pairsOfClassPriority = comparisonPairsOfClass_.get(cName);
	  			if (pairsOfClassPriority == null) pairsOfClassPriority = new Vector<Pair<Double,DoubleBondPositionVO>>();
	  		}
	  		
	  		Set<Pair<Double,DoubleBondPositionVO>> beforeAfter = getBeforeAfter().get(cName);
	  		Vector<Pair<Double,DoubleBondPositionVO>> pairsOfClass = new Vector<Pair<Double,DoubleBondPositionVO>>();
	  		if (beforeAfter != null)
	  		{
	  			pairsOfClass = new Vector<Pair<Double,DoubleBondPositionVO>>(beforeAfter);
	  		}
	  		if (!comparisonPairsOfClass_.isEmpty())
	  		{
	  			pairsOfClass = pairsOfClassPriority;
	  		}
	  		if (pairsOfClass.isEmpty()) continue;
	  			
	  		Vector<Object> elModCharge = getAvailableElementsAndModificationsPlusCharge(quantObjects.get(cName));
		    LinkedHashMap<String,String> mods = (LinkedHashMap<String,String>)elModCharge.get(1);
		    
		    for (Pair<Double,DoubleBondPositionVO> pair : pairsOfClass) //iterating over the original TG entries first to limit false matches.
		    {
		    	boolean pairAdded = false;
		    	for (Pair<Double,DoubleBondPositionVO> original : beforeAfter)
		    	{
		    		if (original.getKey().equals(pair.getKey()) && original.getValue().getChainCombination().equals(pair.getValue().getChainCombination()))
		    		{
		    			pair = original;
		    		}
		    	}
		    	for (String analyte : analyteSequence.get(cName)) 
			    {
			    	Hashtable<String,QuantVO> quantAnalytes = quantObjects.get(cName).get(analyte);
			    	for (String mod : mods.keySet()) 
			      {
			    		QuantVO quant = quantAnalytes.get(mod);
			        Vector<DoubleBondPositionVO> doubleBondPositionVOs = quant.getInfoForOmegaAssignment();
			        for (DoubleBondPositionVO vo : doubleBondPositionVOs)
			        {
			        	if (pair.getValue().getDoubleBondPositionsHumanReadable().equals(vo.getDoubleBondPositionsHumanReadable()))
		        		{
			        		double rtError = vo.getExpectedRetentionTime()-pair.getValue().getExpectedRetentionTime();
			        		double existingRtError = Double.MAX_VALUE;
			        		if (matchedRecalibratedToTarget.containsKey(pair))
			        		{
			        			existingRtError = matchedRecalibratedToTarget.get(pair).getKey();
			        		}
			        		else if (matchedTargetToRecalibrated.containsKey(vo))
			        		{
			        			existingRtError = matchedTargetToRecalibrated.get(vo).getKey();
			        		}
			        		else
			        		{
			        			if (Math.abs(rtError) < 1.0)
			        			{
			        				matchedRecalibratedToTarget.put(pair, new Pair<Double,DoubleBondPositionVO>(rtError,vo));
				        			matchedTargetToRecalibrated.put(vo, new Pair<Double,Pair<Double,DoubleBondPositionVO>>(rtError,pair));
				        			pairAdded = true;
			        			}
			        		}
			        		
			        		if (existingRtError < Double.MAX_VALUE)
			        		{
			        			if (Math.abs(rtError) <= Math.abs(existingRtError))
			        			{
			        				if (!matchedRecalibratedToTarget.containsKey(pair)) //remove existing stuff, so we don't get the false matches
					        		{
			        					Pair<Double,DoubleBondPositionVO> toRemove = null;
					        			for (Pair<Double,DoubleBondPositionVO> previous : matchedRecalibratedToTarget.keySet())
					        			{
					        				if (matchedRecalibratedToTarget.get(previous).getValue().equals(vo) && comparisonPairsOfClass_.isEmpty())
					        				{
					        					toRemove = previous;
					        				}
					        			}
					        			if (toRemove != null) matchedRecalibratedToTarget.remove(toRemove);
					        		}
			        				matchedRecalibratedToTarget.put(pair, new Pair<Double,DoubleBondPositionVO>(rtError,vo));
			        				matchedTargetToRecalibrated.put(vo, new Pair<Double,Pair<Double,DoubleBondPositionVO>>(rtError,pair));
			        				pairAdded = true;
			        			}
			        		}
		        		}
			        }
			      }
			    }
		    	if (!pairAdded && !comparisonPairsOfClass_.isEmpty())
			    {
			    	matchedRecalibratedToTarget.put(pair, pair);
			    }
		    }
		    if (!matchedRecalibratedToTarget.isEmpty())
		    {
		    	ArrayList<DoubleBondPositionVO> usedVOs = new ArrayList<DoubleBondPositionVO>();
		    	for (Pair<Double,DoubleBondPositionVO> originalPair : matchedRecalibratedToTarget.keySet())
		    	{
		    		//first double: rtError, second double: originalRT, vo: recalibrated
		    		Pair<Double,DoubleBondPositionVO> targetPair = matchedRecalibratedToTarget.get(originalPair);
		    		if (usedVOs.contains(targetPair.getValue()) && comparisonPairsOfClass_.isEmpty()) continue;
		    		usedVOs.add(targetPair.getValue());
		    		
		    		String molName = targetPair.getValue().getDoubleBondPositionsHumanReadable();
		    		double targetRT = targetPair.getValue().getExpectedRetentionTime();
        		double originalRT = originalPair.getKey();
        		double recalibratedRT = originalPair.getValue().getExpectedRetentionTime();
        		double rtError = targetPair.getKey()>1.0 ? 0.0 : targetPair.getKey();
        		
        		Row row = sheet.createRow(rowCount);
      			cell = row.createCell(0,HSSFCell.CELL_TYPE_STRING);
	          cell.setCellValue(cName);
      			cell = row.createCell(1,HSSFCell.CELL_TYPE_STRING);
	          cell.setCellValue(molName);
	          cell = row.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
	          cell.setCellValue(targetRT);
	          cell = row.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
	          cell.setCellValue(originalRT);
	          cell = row.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
	          cell.setCellValue(recalibratedRT);
	          cell = row.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
	          cell.setCellValue(originalRT-targetRT);
	          cell = row.createCell(6,HSSFCell.CELL_TYPE_NUMERIC);
	          cell.setCellValue(rtError);
	          cell = row.createCell(7,HSSFCell.CELL_TYPE_NUMERIC);
	          cell.setCellValue(Math.abs(rtError));
	          rowCount++;
		    	}
		    }
			}	
	  	workbook.write(out);
	  	System.out.println("recalibration_comparison.xlsx written!");
  	}
  	catch (Exception e) 
		{
      throw new ExportException(e);
    } 
  }
  
  private void createBeforeAfterTitle(Sheet sheet)
  {
  	Row row = sheet.createRow(0);
  	Cell cell = row.createCell(0,HSSFCell.CELL_TYPE_STRING);
    cell.setCellValue("Lipid Class");
    sheet.setColumnWidth(0, 10 * 256);
    cell = row.createCell(1,HSSFCell.CELL_TYPE_STRING);
    cell.setCellValue("Lipid Molecular Species");
    sheet.setColumnWidth(1, 25 * 256);
    cell = row.createCell(2,HSSFCell.CELL_TYPE_STRING);
    cell.setCellValue("Target DB");
    sheet.setColumnWidth(2, 15 * 256);
    cell = row.createCell(3,HSSFCell.CELL_TYPE_STRING);
    cell.setCellValue("Original DB");
    sheet.setColumnWidth(3, 15 * 256);
    cell = row.createCell(4,HSSFCell.CELL_TYPE_STRING);
    cell.setCellValue("Recalibrated DB");
    sheet.setColumnWidth(4, 15 * 256);
    cell = row.createCell(5,HSSFCell.CELL_TYPE_STRING);
    cell.setCellValue("RT differences (original - target)");
    sheet.setColumnWidth(5, 15 * 256);
    cell = row.createCell(6,HSSFCell.CELL_TYPE_STRING);
    cell.setCellValue("RT error");
    sheet.setColumnWidth(6, 15 * 256);
    cell = row.createCell(7,HSSFCell.CELL_TYPE_STRING);
    cell.setCellValue("Abs RT error");
    sheet.setColumnWidth(7, 15 * 256);
  }
  
  
  
  
  private class Cluster
	{
		private Vector<DoubleBondPositionVO> doubleBondPositionVOs_ = new Vector<DoubleBondPositionVO>();
		
		private Cluster()
		{
			this.doubleBondPositionVOs_ = new Vector<DoubleBondPositionVO>(); 
		}
		
		private Cluster(Cluster cluster1, Cluster cluster2)
		{
			this();
			this.doubleBondPositionVOs_.addAll(cluster1.getDoubleBondPositionVOs());
			this.doubleBondPositionVOs_.addAll(cluster2.getDoubleBondPositionVOs());
		}
		
		private void addDoubleBondPositionVO(DoubleBondPositionVO vo)
		{
			this.doubleBondPositionVOs_.add(vo);
		}
		
		private Vector<DoubleBondPositionVO> getDoubleBondPositionVOs()
		{
			return doubleBondPositionVOs_;
		}
		
		private DoubleBondPositionVO getAverageElement() //combines elements of cluster to one vo with average rt
		{
			if (!doubleBondPositionVOs_.isEmpty())
			{
				DoubleBondPositionVO combined = new DoubleBondPositionVO(doubleBondPositionVOs_.get(0));
				float sum = 0;
				int count = 0;
				for (DoubleBondPositionVO vo : doubleBondPositionVOs_) 
				{
					sum += vo.getExpectedRetentionTime();
					count++;
				}
				float average = sum/count;
				combined.setExpectedRetentionTime(average);
				return combined;
			}
			else
			{
				return null;
			}
		}
	}
  
  
  /**
   * Saves species of a class which were detected with a stable isotope label.
   * All retention times must be corrected for the isotope effect on retention time.
   * 
   * @author Leonida M. Lamp
   *
   */
  private class MolecularSpeciesContainer
  {
  	private Hashtable<String, Vector<DoubleBondPositionVO>> singleLabeledSpecies_;
  	private Hashtable<String, Vector<DoubleBondPositionVO>> multiLabeledSpecies_;
  	
  	private MolecularSpeciesContainer()
  	{
  		singleLabeledSpecies_ = new Hashtable<String, Vector<DoubleBondPositionVO>>();
  		multiLabeledSpecies_ = new Hashtable<String, Vector<DoubleBondPositionVO>>();
  	}
  	
  	private void addSingleLabeledSpecies(String molecularSpecies, DoubleBondPositionVO doubleBondPositionVO)
  	{
  		if (!this.singleLabeledSpecies_.containsKey(molecularSpecies))
  		{
  			this.singleLabeledSpecies_.put(molecularSpecies, new Vector<DoubleBondPositionVO>());
  		}
  		this.singleLabeledSpecies_.get(molecularSpecies).add(doubleBondPositionVO);
  	}
  	
  	private void addMultiLabeledSpecies(String molecularSpecies, DoubleBondPositionVO doubleBondPositionVO)
  	{
  		if (!this.multiLabeledSpecies_.containsKey(molecularSpecies))
  		{
  			this.multiLabeledSpecies_.put(molecularSpecies, new Vector<DoubleBondPositionVO>());
  		}
  		this.multiLabeledSpecies_.get(molecularSpecies).add(doubleBondPositionVO);
  	}
  	
  	private Set<String> getSingleLabeledMolecularSpecies()
  	{
  		return singleLabeledSpecies_.keySet();
  	}
  	
  	private Vector<DoubleBondPositionVO> getSingleLabeledSpecies(String molecularSpecies)
  	{
  		if (singleLabeledSpecies_.containsKey(molecularSpecies))
  		{
  			return singleLabeledSpecies_.get(molecularSpecies);
  		}
  		else
  		{
  			return null;
  		}
  	}
  	
  	private Vector<DoubleBondPositionVO> getMultiLabeledSpecies(String molecularSpecies)
  	{
  		if (multiLabeledSpecies_.containsKey(molecularSpecies))
  		{
  			return multiLabeledSpecies_.get(molecularSpecies);
  		}
  		else
  		{
  			return null;
  		}
  	}
  }
  
  
  
  
}





