package at.tugraz.genome.lda.target.export;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.Vector;

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

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.QuantificationThread;
import at.tugraz.genome.lda.Settings;
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
import at.tugraz.genome.lda.utils.Pair;
import at.tugraz.genome.lda.utils.RangeDouble;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.lda.vos.ResultFileVO;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;

/**
 * We do the following: 
 * possibly run the analysis module on our read in files for clustering 
 * then for each species containing a defined label (here we just use the prefix I guess):
 * we write the molecular species according to the defined labels and calculate the RT according to the formula
 * then it's written out; we likely do not need a quant file for ordering because we can just remove the prefix via regex I believe
 *
 * !!!!!!!!!!!! we need a template masslist !!!!!!
 * 
 * with this things changed.. we do not need comparative analysis: we just take our template masslist.. 
 * ...just locate it somewhere in the filesystem and give it a default path for now.
 * 
 * then we go through the class and correct analyte sequence and just write the masslist out. 
 * group hits together if they are below a certain threshold apart
 * otherwise this should be relatively straightforward
 * 
 * we first iterate through the classes: the analytes in correctAnalyteSequence and get the quantObjects for each;
 * we iterate over all resultfiles and labeled species are added as new QuantVO to the analyte sequence;
 * if there is already an identical analyte with just a different RT: maybe just add the RTs of all to the quantVO, or an average if they are close enough?
 * add them all as doublebondpositionVOs, 
 * when we iterate over this for the export we can do a median over identical ones.
 * as a final step we just export the QuantVOs.
 * 
 * 
 * 
 * 
 */
public class TargetListExporter
{	
  private final static String HEADER_NAME = "Name";
  private final static String HEADER_COLON = "";
  private final static String HEADER_DBS = "dbs";
  public final static String HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS = "mol. species";
  private final static String HEADER_RETENTION_TIME = "tR (min)";
  
  private final static int HEADER_ROW = 1;
  private final static double THRESHOLD_FOR_CLUSTERING = 5.0; //TODO: this should be a user decision, as it depends on the length of the chromatography
  
  private Hashtable<String, RecalibrationRegression> recalibrationRegressionsForClass_;
  private RecalibrationRegression recalibrationRegression_;
  private IsotopeEffectRegression isotopeEffectRegression_;
  private Vector<ResultFileVO> resultFileVO_;
  private Vector<IsotopeLabelVO> labels_;
  private String templatePath_;
  Vector<Pair<Double,Pair<String,DoubleBondPositionVO>>> beforeAfter_ = new Vector<Pair<Double,Pair<String,DoubleBondPositionVO>>>(); //TODO: just for development! remove later
  
  /**
   * Constructor for omega MassList Recalibration
   * @param recalibrationRegression
   */
  public TargetListExporter(RecalibrationRegression recalibrationRegression, String templatePath)
	{
		this.recalibrationRegression_ = recalibrationRegression;
		this.isotopeEffectRegression_ = null;
		this.resultFileVO_ = new Vector<ResultFileVO>();
		this.labels_ = new Vector<IsotopeLabelVO>();
		this.templatePath_ = templatePath;
	}
  
  
  public TargetListExporter(Hashtable<String, RecalibrationRegression> recalibrationRegressionsForClass, String templatePath)
	{
  	this.recalibrationRegressionsForClass_ = recalibrationRegressionsForClass;
  	this.isotopeEffectRegression_ = null;
		this.isotopeEffectRegression_ = null;
		this.resultFileVO_ = new Vector<ResultFileVO>();
		this.labels_ = new Vector<IsotopeLabelVO>();
		this.templatePath_ = templatePath;
	}
  
  
	//TODO: create MassList dynamically
  /**
   * Constructor for omega MassList creation with stable isotope labeled standards
   * @param isotopeEffectRegression
   * @param resultIndices
   * @param resultFileVO
   * @param labels
   */
	public TargetListExporter(
			IsotopeEffectRegression isotopeEffectRegression, Vector<ResultFileVO> resultFileVO, Vector<IsotopeLabelVO> labels)
	{
		this.recalibrationRegression_ = null;
		this.isotopeEffectRegression_ = isotopeEffectRegression;
		this.resultFileVO_ = resultFileVO;
		this.labels_ = labels;
	}
	
	
	
	
	/**
	 * How do we do this....
	 * we go through all result files, 
	 * > filter for only analytes with MSn info
	 * 
	 * > Hashtable className, new java class: 
	 * 				saves labeled and unlabeled species of class separately
	 * 				method getlabeled ( position ind. mol species)
	 * 				method getunlabeled ( position ind. mol species) ...probably not needed anymore
	 * 				method getlabeledmolecularspecies
	 * 				2x field hashtable pos ind mol species, vector lipidomicsmsnset
	 * 				
	 * iterate over labeledmolecularspecies.
	 * get labeled -> expected RT is calculated, then species are clustered
	 * get unlabeled -> species are clustered
	 * 
	 * for the calculation of the expected RT, only the number of deuteriums in the chemical formula is considered
	 * for adding of the omega position, the prefix of the fatty acid chain is considered (can only have one prefix)
	 * 
	 * for each labeled cluster 
	 * 		if within our cluster we have deuterated species with more labels but matching expected RT, add omega pos to all in the cluster 
	 * 		(and expected sn position if all species in the cluster have identical ones?), 
	 * 		
	 * 		we then see if within threshold there's a matching partner with same mol. species, if so adjust RT
	 * 		then add the doublebondpositionVOs to the quants
	 * 		then export
	 * 	
	 * @throws ExportException
	 */
	@SuppressWarnings("unchecked")
	protected void export(String templatePath, String outPath) throws ExportException
	{		
		try (	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outPath));
					XSSFWorkbook workbook = new XSSFWorkbook();)
		{
			Vector<?> quantInfo = QuantificationThread.getCorrectAnalyteSequence(templatePath,false);
			LinkedHashMap<String,Integer> classSequence = (LinkedHashMap<String,Integer>)quantInfo.get(0);
			LinkedHashMap<String,Vector<String>> analyteSequence = (LinkedHashMap<String,Vector<String>>)quantInfo.get(1);
			Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo.get(4);
		  
      new ElementConfigParser("elementconfig.xml").parse();
      
      XSSFCellStyle headerStyle = getHeaderStyle(workbook);  
      for (String cName : classSequence.keySet()) 
      {
      	System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+ cName + "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
      	Vector<DoubleBondPositionVO> allLabeledVOsToAdd = new Vector<DoubleBondPositionVO>();
      	if (this.recalibrationRegression_ == null && this.recalibrationRegressionsForClass_ == null)
      	{
      		MolecularSpeciesContainer container = new MolecularSpeciesContainer();
        	fillContainerForClass(cName, container);
        	Set<String> labeledMolecularSpecies = container.getLabeledMolecularSpecies();
        	for (String molecularSpecies : labeledMolecularSpecies)
        	{
        		Vector<DoubleBondPositionVO> labeledVOs = container.getLabeledSpecies(molecularSpecies);
        		Vector<DoubleBondPositionVO> clusteredLabeledVOs = clusterMolecularSpecies(labeledVOs, true);
        		allLabeledVOsToAdd.addAll(clusteredLabeledVOs); 	
        	}
      	}
        Sheet sheet = workbook.createSheet(cName);
        writeMassListForSheet(sheet, headerStyle, analyteSequence.get(cName), quantObjects.get(cName), cName, allLabeledVOsToAdd);
      }
      workbook.write(out);
      System.out.println("workbook written!");
    } 
		catch (Exception e) 
		{
      throw new ExportException(e);
    } 
	}
	
	@SuppressWarnings("unchecked")
	private void writeMassListForSheet(
			Sheet sheet, 
			XSSFCellStyle headerStyle, 
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
	        if (this.recalibrationRegression_ == null && this.recalibrationRegressionsForClass_ == null && quant.getDbs() > 0) //only species with double bonds are relevant
	        {
	        	for (DoubleBondPositionVO doubleBondPositionVO : allLabeledVOsToAdd)
		        {
		        	if (quant.getCarbons() == doubleBondPositionVO.getNumberOfCarbons() &&
		        			quant.getDbs() == doubleBondPositionVO.getNumberOfDoubleBonds() &&
		        			quant.getOhNumber() == doubleBondPositionVO.getNumberOfOH())
		        	{
		        		quant.addInfoForOmegaAssignment(doubleBondPositionVO);
		        	}
		        }
	        }
	        doubleBondPositionVOs = quant.getInfoForOmegaAssignment();
	      }
	      
	      if (this.recalibrationRegression_ != null || this.recalibrationRegressionsForClass_ != null)
	      {
	      	RecalibrationRegression regression = recalibrationRegression_;
      		if (recalibrationRegressionsForClass_ != null)
      		{
      			String key = CalibrationGraphPanel.PLOT_ALL; //lipid classes without a specific regression will be calibrated with the regression for all.
      			if (recalibrationRegressionsForClass_.contains(cName))
      			{
      				key = cName;
      			}
      			regression = recalibrationRegressionsForClass_.get(key);
      		}
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
	      
//	      if (!doubleBondPositionVOs.isEmpty())
//	      {
//	      	System.out.println("before: "+doubleBondPositionVOs.size());
//	      	Collections.sort(doubleBondPositionVOs);
//	      	System.out.println("after: "+doubleBondPositionVOs.size());
//	      }
	      if (!doubleBondPositionVOs.isEmpty())
	      	Collections.sort(doubleBondPositionVOs);
	      
	      
	      
	      for (int i=-1;i<doubleBondPositionVOs.size();i++)
	      {
	      	int modCount = 0;
	      	Row row = sheet.createRow(rowCount);
	      	for (String mod : mods.keySet()) 
		      {
		      	QuantVO quant = quantAnalytes.get(mod);
		      	if (modCount==0) 
		        {   	
		          cell = row.createCell(headerTitles.indexOf(HEADER_NAME),HSSFCell.CELL_TYPE_STRING);
		          cell.setCellValue(quant.getAnalyteName());
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
	
	private void computeExpectedRetentionTime(DoubleBondPositionVO vo, RecalibrationRegression regression, String cName)
	{
		double before = vo.getExpectedRetentionTime();
		try
		{
			double after = regression.getTargetRT(vo.getExpectedRetentionTime());
			vo.setExpectedRetentionTime((float)after);
			System.out.println(String.format("original: %s, recalibrated: %s", before, after));
		}
		catch (OutOfRangeException ex)
		{
			vo.setExpectedRetentionTime((float)-1.0);
		}
		
		beforeAfter_.add(new Pair<Double,Pair<String,DoubleBondPositionVO>>(before,new Pair<String,DoubleBondPositionVO>(cName,vo))); //TODO: just for development, remove after!
	}
	
	
	//first all vos are sorted according to the labeling pattern (they already all have the same molecular species) into objects
	//then for the vos within these objects the mean and standard deviation is calculated
	//if standard deviation is below our threshold, then this is an accepted cluster, we save it as a new vo
	//when we are done, the final step is to take all obtained doublebondpositionvos and check if another dbpvo with a compatible labeling pattern is within the time threshold, if so, merge
	//TODO: in this class we get only vos of identical molecular species. this means we can compute permutable positions for this molecular species and use this to compute overlapping patterns!
	private Vector<DoubleBondPositionVO> clusterMolecularSpecies(Vector<DoubleBondPositionVO> doubleBondPositionVOs, boolean labeled)
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
				Vector<DoubleBondPositionVO> averageElements = computeClusters(patternLookup.get(pattern));
				patternLookup.put(pattern, averageElements);
			}
			
			/**
			 * TODO: for now we do nothing with patterns which could be combined, because I am not sure if it's scientifically accurate to do
			 * plan is:  iterate over markedToCombine (they are within threshold and have an overlapping pattern) and combine to a common vo 
			 */
			
			Hashtable<Vector<Integer>,Vector<Vector<Integer>>> overlappingPatterns = computeOverlappingPatterns(patternLookup.keySet(),permutationPattern);
			
			Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>> markedToCombine = markToCombine(patternLookup, overlappingPatterns);
			Vector<Pair<RangeDouble,Vector<DoubleBondPositionVO>>> combineRanges = new Vector<Pair<RangeDouble,Vector<DoubleBondPositionVO>>>();
			
			
			/**
			 * TODO: I have an interesting bug here: the current settings will print the weirdness out. I do not understand it. Sorry future me.
			 */
			for (DoubleBondPositionVO doubleBondPositionVO : markedToCombine.keySet())
			{
				boolean isRepresented = false;
				for(Pair<RangeDouble,Vector<DoubleBondPositionVO>> range1 : combineRanges)
				{
					if (range1.getKey().insideRange(doubleBondPositionVO.getExpectedRetentionTime()))
					{
						RangeDouble range2 = new RangeDouble(doubleBondPositionVO.getExpectedRetentionTime()-THRESHOLD_FOR_CLUSTERING/60.0, doubleBondPositionVO.getExpectedRetentionTime()+THRESHOLD_FOR_CLUSTERING/60.0);
						range1.getKey().extendToOtherRanges(range2);
						range1.getValue().add(doubleBondPositionVO);
						isRepresented = true;
					}
				}
				if (!isRepresented)
				{
					RangeDouble range2 = new RangeDouble(doubleBondPositionVO.getExpectedRetentionTime()-THRESHOLD_FOR_CLUSTERING/60.0, doubleBondPositionVO.getExpectedRetentionTime()+THRESHOLD_FOR_CLUSTERING/60.0);
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
						float expectedRetentionTime = doubleBondPositionVO.getExpectedRetentionTime();
						int count = 1;
						for (DoubleBondPositionVO vo : vos)
						{
							toBeRemovedFromLookup.add(vo);
//							System.out.println("To be removed vo1: "+vo.getMolecularSpecies()+"        pattern: "+vo.getPositionAssignmentPattern()+"     RT: "+vo.getExpectedRetentionTime());
							expectedRetentionTime+= vo.getExpectedRetentionTime();
							count++;
						}
						DoubleBondPositionVO combined = new DoubleBondPositionVO(combinedFattyAcidVO, expectedRetentionTime/count, 0, doubleBondPositionVO.getMolecularSpecies());
						toBeAddedToLookup.add(combined);
//						System.out.println("To be removed parent: "+doubleBondPositionVO.getMolecularSpecies()+"        pattern: "+doubleBondPositionVO.getPositionAssignmentPattern()+"     RT: "+doubleBondPositionVO.getExpectedRetentionTime());
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
//							System.out.println("To be added: "+toBeAdded.getMolecularSpecies()+"        pattern: "+toBeAdded.getPositionAssignmentPattern()+"     RT: "+toBeAdded.getExpectedRetentionTime());
							added.add(toBeAdded.getPositionAssignmentPattern());
						}
					}
					for (DoubleBondPositionVO toBeRemoved : toBeRemovedFromLookup)
					{
						patternLookup.get(toBeRemoved.getPositionAssignmentPattern()).remove(toBeRemoved);
//						System.out.println("To be removed: "+toBeRemoved.getMolecularSpecies()+"        pattern: "+toBeRemoved.getPositionAssignmentPattern()+"     RT: "+toBeRemoved.getExpectedRetentionTime());
					}
				}
			}
			
			
//			List<DoubleBondPositionVO> blackList = new ArrayList<DoubleBondPositionVO>();
//			Vector<DoubleBondPositionVO> whiteList = new Vector<DoubleBondPositionVO>();
			/**
			 * markedToCombine must be filtered to some degree:
			 * if only 1 match found: combine without issue (includes chains with permuation etc.)
			 * if multiple matches found: check if one of the vos has double bond positions assigned for all chains:
			 * if yes, check if all positions different from -1 are identical: 
			 * 			if yes, combine all
			 * 			if no, all of them should be removed from the data: add them to a list of elements to be removed and check each time if an equivalent vo is contained there.
			 * if no, check total number of assigned positions for individual fatty acid species:
			 * 			if number matches for all, combine without issue (includes chains with permuation etc.)
			 * 			if no, remove species with the offending fatty acid assigned and combine the rest
			 * TODO: this needs to be expanded, not general enough for triglycerides and some other cases
			 */
			
			/**
			 * new approach to deal with markedToCombine:
			 * we do NOT remove any analytes. we consider all possible combinations for a peak to be correct!
			 * 
			 * 
			 */
//			for (DoubleBondPositionVO vo1 : markedToCombine.keySet())
//			{
//				System.out.println("Parent VO pattern: "+vo1.getPositionAssignmentPattern()+" RT: "+vo1.getExpectedRetentionTime());
////				patternLookup.get(vo1.getPositionAssignmentPattern()).remove(vo1); //removing element
//				Vector<DoubleBondPositionVO> toCombine = markedToCombine.get(vo1);
//				if (vo1.areDoubleBondPositionsAssignedForAllChains())
//				{
//					DoubleBondPositionVO combined = new DoubleBondPositionVO(vo1);
//					float expectedRetentionTime = combined.getExpectedRetentionTime();
//					int count = 1;
//					for (DoubleBondPositionVO vo : toCombine)
//					{
//						patternLookup.get(vo.getPositionAssignmentPattern()).remove(vo); //removing element
//						expectedRetentionTime+= vo.getExpectedRetentionTime();
//						count++;
//					}
//					combined.setExpectedRetentionTime(expectedRetentionTime/count);
//					if (!patternLookup.containsKey(combined.getPositionAssignmentPattern()))
//						patternLookup.put(combined.getPositionAssignmentPattern(), new Vector<DoubleBondPositionVO>());
//					patternLookup.get(combined.getPositionAssignmentPattern()).add(combined);
//				}
//				else if (toCombine.size() == 1)
//				{
//					
//					patternLookup.get(toCombine.get(0).getPositionAssignmentPattern()).remove(toCombine.get(0)); //removing element
//					DoubleBondPositionVO combined = combineDoubleBondPositionVO(vo1,toCombine.get(0));
//					if (combined != null)
//					{
//						if (!patternLookup.containsKey(combined.getPositionAssignmentPattern()))
//							patternLookup.put(combined.getPositionAssignmentPattern(), new Vector<DoubleBondPositionVO>());
//						patternLookup.get(combined.getPositionAssignmentPattern()).add(combined);
//					}
//				}
//				
//				toCombine.forEach((n) -> System.out.println((n).getPositionAssignmentPattern()));
//				if(vo1.getMolecularSpecies().equals("16:2_18:1"))
//				{
//					System.out.println("stop?");
//				}
//			}	
			
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
			return computeClusters(doubleBondPositionVOs);
//			System.out.println("next");
//			calculateDistanceMatrix(doubleBondPositionVOs);
		}
	}
	
	/**
	 * returns true if only one possible double bond position combination is possible for the queried data
	 * @param vo1
	 * @param doubleBondPositionVOs
	 * @return
	 */
	private boolean isCombinationStraightforward(Set<Vector<Integer>> combinedPatterns)
	{
		if (combinedPatterns.size() != 1)
		{
			return false;
		}
		return true;
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
	
	
//	private Vector<Integer> combinePattern(Vector<Vector<Integer>> patterns)
//	{
//		Vector<Integer> combinedPattern = new Vector<Integer>();
//		Vector<Integer> temp = new Vector<Integer>();
//		for (Vector<Integer> pattern : patterns)
//		{
//			if (temp.isEmpty())
//			{
//				temp = pattern;
//				continue;
//			}
//			for (int i=0;i<pattern.size();i++)
//			{
//				if ((pattern.get(i) == temp.get(i)) || patter )
//			}
//		}
//		
//		
//	}
	
	
	
	private Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>> markToCombine(
			Hashtable<Vector<Integer>,Vector<DoubleBondPositionVO>> patternLookup,
			Hashtable<Vector<Integer>,Vector<Vector<Integer>>> overlappingPatterns)
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
						if (Math.abs(retentionTime1-retentionTime2) < (THRESHOLD_FOR_CLUSTERING/60)) //if is within threshold
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
	private Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>> markToCombine(
			Hashtable<Vector<Integer>,Vector<DoubleBondPositionVO>> patternLookup,
			Hashtable<Vector<Integer>,Vector<Vector<Integer>>> overlappingPatterns)
	{
		Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>> markedToCombine = new Hashtable<DoubleBondPositionVO,Vector<DoubleBondPositionVO>>();
		try
		{
			for (Vector<Integer> pattern : overlappingPatterns.keySet()) //combine overlapping patterns with same retention times
			{
				Vector<DoubleBondPositionVO> patternVOs = patternLookup.get(pattern); //the species matching the pattern exactly
				for (Vector<Integer> overlappingPattern : overlappingPatterns.get(pattern)) //pattern overlap with the queried pattern
				{
					Vector<DoubleBondPositionVO> overlappingPatternVOs = patternLookup.get(overlappingPattern); //species with pattern overlap with the queried pattern
					if (overlappingPatternVOs!=null)
					{
						for (DoubleBondPositionVO vo1 : patternVOs)
						{
							Set<Vector<Integer>> patternLookupKeySet = patternLookup.keySet();
							for (Vector<Integer> pttrn : patternLookupKeySet)
							{
								if (pttrn.equals(overlappingPattern))
								{
									System.out.println("HI");
									Vector<DoubleBondPositionVO> test = patternLookup.get(pttrn);
									System.out.println("HI!");
								}
							}
							for (Entry<Vector<Integer>,Vector<DoubleBondPositionVO>> dbvo : patternLookup.entrySet())
							{
								Vector<DoubleBondPositionVO> testest = dbvo.getValue();
								Vector<Integer> henlo = dbvo.getKey();
								Vector<DoubleBondPositionVO> testestest = patternLookup.get(henlo);
								System.out.println("HI!");
							}
							
							Float retentionTime1 = (float)vo1.getExpectedRetentionTime();
							for (DoubleBondPositionVO vo2 : overlappingPatternVOs)
							{
								Float retentionTime2 = (float)vo2.getExpectedRetentionTime();
//								System.out.println(String.format("Species1: %s, RT1: %s, Species2: %s, RT2: %s ", 
//										vo1.getDoubleBondPositionsHumanReadable(), retentionTime1, vo2.getDoubleBondPositionsHumanReadable(), retentionTime2));
								if (Math.abs(retentionTime1-retentionTime2) < (THRESHOLD_FOR_CLUSTERING/60)) //if is within threshold
								{
									if (markedToCombine.containsKey(vo2)) 
									{
										if (!markedToCombine.get(vo2).contains(vo1)) 
										{
											markedToCombine.get(vo2).add(vo1);
										}
									} 
									else if (markedToCombine.containsKey(vo1))
									{
										if (!markedToCombine.get(vo1).contains(vo2)) 
										{
											markedToCombine.get(vo1).add(vo2);
										}
									}
									else
									{
										markedToCombine.put(vo1, new Vector<DoubleBondPositionVO>());
										markedToCombine.get(vo1).add(vo2);
									}
								}
							}
						}
					}
					else 
					{
						System.out.println(overlappingPattern);
					}
				}				
			}
		}
		catch (NullPointerException ex)
		{
			ex.printStackTrace();
			System.out.println("HI!");
		}
		
		return markedToCombine;
	}
	*/
	
	
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
  		float combinedRetentionTime = (vo1.getExpectedRetentionTime()+vo2.getExpectedRetentionTime())/2;
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
  		return new DoubleBondPositionVO(combinedFattyAcidVO, combinedRetentionTime, 0, vo1.getMolecularSpecies());
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
	 * just for verifying the clustering is ok
	 * TODO: delete later
	 * @param doubleBondPositionVOs
	 */
	private void calculateDistanceMatrix(Vector<DoubleBondPositionVO> doubleBondPositionVOs)
	{
		Vector<Vector<Float>> distanceMatrix = new Vector<Vector<Float>>();
		for (DoubleBondPositionVO vo1 : doubleBondPositionVOs)
		{
			Float retentionTime1 = (float)vo1.getExpectedRetentionTime();
			Vector<Float> distances = new Vector<Float>();
			for (DoubleBondPositionVO vo2 : doubleBondPositionVOs)
			{
				distances.add(Math.abs(retentionTime1-(float)vo2.getExpectedRetentionTime()));
			}
			distanceMatrix.add(distances);
			System.out.println(distances);
		}
	}
	
	/**
	 * what we do is the following: 
	 * fist make a hashtable of clusters (key index, value cluster)
	 * if a species is not in a cluster yet, add species to new cluster, add it to hashtable with index
	 * else add cluster to hashtable with index
	 * check for each other species that is within threshold if it is in a cluster, if so if it is the same cluster, if not, merge clusters
	 * merging of clusters works as follows: new cluster with elements of both, assign this new cluster to the initial clusters
	 * @param doubleBondPositionVOs
	 */
	private Vector<DoubleBondPositionVO> computeClusters(Vector<DoubleBondPositionVO> doubleBondPositionVOs)
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
				if (Math.abs(retentionTime1-(float)vo2.getExpectedRetentionTime()) < (THRESHOLD_FOR_CLUSTERING/60)) //if is within threshold
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
	 * Check to identify species which may have a more ambiguous retention time, due to coelution
	 */
	private boolean containsMoreThanOneSpecies(Vector<Pair<String,String>> labeledUnlabeledPairs)
	{
		Set<String> uniqueSpecies = new HashSet<String>();
		for (Pair<String,String> pair : labeledUnlabeledPairs)
		{
			uniqueSpecies.add(pair.getKey());
		}
		return uniqueSpecies.size() > 1;
	}
	
	
	
	private void fillContainerForClass(String cName, MolecularSpeciesContainer container)
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
						
						//TODO: only for printing, remove later
//						for (Pair<String,String> pair : labeledUnlabeledPairs)
//						{
//							if (analyte.getChemicalFormula().contains("D"))
//								System.out.println(String.format("class: %s, unlabeled: %s, labeled: %s, RT: %s ", cName, pair.getKey(), pair.getValue(), analyte.getPreciseRT()));
//						}
						
						Hashtable<String,Integer> elements = StaticUtils.categorizeFormula(analyte.getChemicalFormula());
						double expectedRetentionTime = analyte.getPreciseRT();
						
						//only deuterated species need recalculation of the expected retention time
						if (isotopeEffectRegression_ != null && analyte.getChemicalFormula().contains("D"))
						{
							int numberDeuterium = elements.get("D");
							if (numberDeuterium > isotopeEffectRegression_.getMaxNumDeuteriumAllowed())
							{
								continue; //retention times affected by labels outside the calibration curve cannot be recalculated!
							}
							expectedRetentionTime = isotopeEffectRegression_.getRTofUnlabeledSpecies(numberDeuterium, analyte.getPreciseRT());
						}					
						
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
			    		
			    		//all isotope labeled species are stored separate from unlabeled species, as only they get written out in the end
			    		if (analyte.getChemicalFormula().contains("D") || analyte.getChemicalFormula().contains("Cc")) //TODO: consider accounting for additional isotopes?
							{
				    		container.addLabeledSpecies(
				    				pair.getKey(), 
				    				new DoubleBondPositionVO(chains,(float)expectedRetentionTime,0,pair.getKey()));
							}
			    		else 
			    		{
			    			container.addUnlabeledSpecies(
				    				pair.getKey(), 
				    				new DoubleBondPositionVO(chains,(float)expectedRetentionTime,0,pair.getKey()));
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
	
	
//	private Hashtable<String,Vector<LipidParameterSet>> collectLabeledResultsForClass(String cName)
//	{
//		Hashtable<String,Vector<LipidParameterSet>> labeledResultsForClass = new Hashtable<String,Vector<LipidParameterSet>>();
//		for (int index : resultIndices_)
//		{
//			QuantificationResult quantRes = resultFileVO_.get(index).getQuantificationResult();
//			Vector<LipidParameterSet> analytes = quantRes.getIdentifications().get(cName);
//			
//			for (LipidParameterSet analyte : analytes)
//			{
//				if (analyte instanceof LipidomicsMSnSet 
//						&& analyte.getDoubleBonds()>0 
//						&& (analyte.getChemicalFormula().contains("D") || analyte.getChemicalFormula().contains("Cc")))
//				{
//					String name = analyte.getNameStringWithoutRt();
//					if (!labeledResultsForClass.contains(name))
//					{
//						labeledResultsForClass.put(name, new Vector<LipidParameterSet>());
//					}
//					analyte = addDoubleBondPositionVO((LipidomicsMSnSet)analyte);
//					labeledResultsForClass.get(name).add(analyte);
//				}
//			}
//		}
//		return labeledResultsForClass;
//	}
	
	
	//for this task I need the label definitions, then probably replace elements accordingly and create a doublebondpositionvo, which is added to the analyte
	private LipidParameterSet addDoubleBondPositionVO(LipidomicsMSnSet analyte)
	{
		
		return analyte;
	}
	
	
	private Hashtable<String,Hashtable<String,QuantVO>> addDoubleBondPositionsToQuantVOs(
			Hashtable<String,Hashtable<String,QuantVO>> quantObjects, Vector<String> analyteSequence, String cName)
	{
//		this.isotopeEffectRegression_ = isotopeEffectRegression;
		for (ResultFileVO resultFileVO : resultFileVO_)
		{
			QuantificationResult quantRes = resultFileVO.getQuantificationResult();
			Vector<LipidParameterSet> analytes = quantRes.getIdentifications().get(cName);
			
			
		}
		
		
		return quantObjects;
	}
	
	
	
	
	/**
   * Creates a formatted header row with the given header titles in the row number given by HEADER_ROW
   * @param ws Excel worksheet to write the header to
   * @param headerTitles List of header titles
   */
  private static void createHeader(Sheet sheet, List<String> headerTitles, XSSFCellStyle headerStyle) 
  {
  	Row row = sheet.createRow(HEADER_ROW);
  	Cell cell;
  	for (int i=0; i<headerTitles.size(); i++) 
  	{
  		cell = row.createCell(i, HSSFCell.CELL_TYPE_STRING);
  		cell.setCellValue(headerTitles.get(i));
  		cell.setCellStyle(headerStyle);
  	}
  }
	
	/**
	 * Creates a list of header titles
	 * @param elements
	 * @param mods
	 * @return
	 */
  private static List<String> createHeaderTitles(
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
  protected Vector<Object> getAvailableElementsAndModificationsPlusCharge(Hashtable<String,Hashtable<String,QuantVO>> quantObjects) throws ChemicalFormulaException
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
  
  //TODO: just for development, remove later!
  protected Vector<Pair<Double,Pair<String,DoubleBondPositionVO>>> getBeforeAfter()
	{
		return this.beforeAfter_;
	}
  
  public String getTemplatePath()
  {
  	return this.templatePath_;
  }
  
  
  /**
   * Inner helper class
   * 
   * > Hashtable className, new java class: 
	 * 				saves labeled and unlabeled species of class separately
	 * 				method getlabeled ( position ind. mol species)
	 * 				method getunlabeled ( position ind. mol species)
	 * 				method getlabeledmolecularspecies
	 * 				2x field hashtable pos ind mol species, vector lipidomicsmsnset
   */
  private class MolecularSpeciesContainer
  {
  	private Hashtable<String, Vector<DoubleBondPositionVO>> labeledSpecies_;
  	private Hashtable<String, Vector<DoubleBondPositionVO>> unlabeledSpecies_;
  	
  	private MolecularSpeciesContainer()
  	{
  		labeledSpecies_ = new Hashtable<String, Vector<DoubleBondPositionVO>>();
  		unlabeledSpecies_ = new Hashtable<String, Vector<DoubleBondPositionVO>>();
  	}
  	
  	private void addLabeledSpecies(String molecularSpecies, DoubleBondPositionVO doubleBondPositionVO)
  	{
  		if (!this.labeledSpecies_.containsKey(molecularSpecies))
  		{
  			this.labeledSpecies_.put(molecularSpecies, new Vector<DoubleBondPositionVO>());
  		}
  		this.labeledSpecies_.get(molecularSpecies).add(doubleBondPositionVO);
  	}
  	
  	private void addUnlabeledSpecies(String molecularSpecies, DoubleBondPositionVO doubleBondPositionVO)
  	{
  		if (!this.unlabeledSpecies_.containsKey(molecularSpecies))
  		{
  			this.unlabeledSpecies_.put(molecularSpecies, new Vector<DoubleBondPositionVO>());
  		}
  		this.unlabeledSpecies_.get(molecularSpecies).add(doubleBondPositionVO);
  	}
  	
  	private Set<String> getLabeledMolecularSpecies()
  	{
  		return labeledSpecies_.keySet();
  	}
  	
  	private Vector<DoubleBondPositionVO> getLabeledSpecies(String molecularSpecies)
  	{
  		if (labeledSpecies_.containsKey(molecularSpecies))
  		{
  			return labeledSpecies_.get(molecularSpecies);
  		}
  		else
  		{
  			return null;
  		}
  	}
  	
  	private Vector<DoubleBondPositionVO> getUnlabeledSpecies(String molecularSpecies)
  	{
  		if (unlabeledSpecies_.containsKey(molecularSpecies))
  		{
  			return unlabeledSpecies_.get(molecularSpecies);
  		}
  		else
  		{
  			return null;
  		}
  	}
  }
}





