package at.tugraz.genome.lda.target.experiment;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.Vector;

import javax.swing.table.DefaultTableModel;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.target.IsotopeLabelVO;
import at.tugraz.genome.lda.target.JDefaultComponents;
import at.tugraz.genome.lda.target.JOptionPanel;
import at.tugraz.genome.lda.target.LoadingPanel;
import at.tugraz.genome.lda.target.export.ExportOptionsPanel;
import at.tugraz.genome.lda.target.export.ExportPanel;
import at.tugraz.genome.lda.target.export.TargetListExporter;
import at.tugraz.genome.lda.utils.Pair;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.ResultFileVO;


public class ExperimentGraphPanel extends JOptionPanel
{
	private static final long serialVersionUID = 1L;
	
	private static final int EDITABLE_COLUMN = 2;
	private static final String TABLE_FRAME_TITLE = "Computation of the total deuterium isotope effect on retention time.";
	private static final Dimension TABLE_PANEL_DIMENSION = new Dimension(825,200);
	private static final Dimension PLOT_DIMENSION = new Dimension(825,500);
	
	private ExperimentLabelDefinitionPanel labelDefinitionPanel_;
	private LoadingPanel loadingPanel_;
	private JPanel displayPanel_;
	private Vector<MatchedPartnerVO> matchedIsotopologues_;
	private IsotopeEffectRegression isotopeEffectRegression_;
	private TotalIsotopeEffectPlot plot_;
  
  public ExperimentGraphPanel(JDefaultComponents wizardComponents) {
      super(wizardComponents, "Use stable isotope labels specific for \u03C9-C=C positions.");
      this.loadingPanel_ = new LoadingPanel("<html>Computing data, please wait...</html>");
      init(this.loadingPanel_);
  }
  
  protected void init(JPanel panel) 
  {
    this.add(panel, new GridBagConstraints(0, 0, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    this.invalidate();
    this.updateUI();
  }
  
  public void matchIsotopologues(Vector<Pair<ResultFileVO, Integer>> standardsFileVOs, ExperimentLabelDefinitionPanel labelDefinitionPanel) throws LipidCombinameEncodingException
  {
  	Hashtable<String,Vector<LipidomicsMSnSet>> authenticStandards = collectStandards(standardsFileVOs);
  	this.labelDefinitionPanel_ = labelDefinitionPanel;
  	matchedIsotopologues_ = collectMatchedIsotopologues(authenticStandards, labelDefinitionPanel);
  	//fit the curve
  	isotopeEffectRegression_ = new IsotopeEffectRegression(matchedIsotopologues_, labelDefinitionPanel.getAssignedIsotopeLabels());
  }
  
  public void initDataDisplay()
  {
  	cleanPanels();
  	String[] columnNames = { "Authentic standards", "Matched isotopologues", "Choose pairs for the calibration."};
  	Object[][] tableData = generateTableData();
  	getDefaultComponents().updateComponents();
  	init(generateDisplayPanel(columnNames, tableData, EDITABLE_COLUMN));
  }
  
  private void cleanPanels()
  {
  	if (this.loadingPanel_ != null)
  	{
  		this.remove(loadingPanel_);
  		this.loadingPanel_ = null;
  	}
  	if (this.displayPanel_ != null)
  	{
  		this.remove(displayPanel_);
  		this.displayPanel_ = null;
  	}
  }
  
  /**
   * Computes the name for an analyte to be displayed in the table
   * @param match
   * @param isStandard		true if the authentic standard rather then the matching isotopologue should be displayed
   * @return
   */
  private String convertLipidomicsMSnSetToTableDisplayName(MatchedPartnerVO match, boolean isStandard)
  {
  	String displayName = "";
  	LipidomicsMSnSet standard = match.getStandard();
  	if (isStandard)
  	{
  		displayName = String.format("%s %s, RT = %s min", match.getLipidClass(), standard.getOmegaInformation().get(0).getDoubleBondPositionsHumanReadable(), standard.getRt());
  	}
  	else
  	{
  		//for peaks with multiple molecular species the correct one should be displayed
  		LipidomicsMSnSet set = match.getIsotopologue();
  		Vector<String> normed = normChainCombination(standard.getValidChainCombinations().iterator().next(), match.getLabel()); //assuming the standard only consists of 1 molecular species
  		Vector<String> chainCombis = set.getValidChainCombinations();
  		
  		for (String chainCombi : chainCombis)
  		{
  			if (normChainCombination(chainCombi, match.getLabel()).equals(normed))
  			{
  				try
  				{
    				displayName = String.format("%s %s, RT = %s min", match.getLipidClass(), set.getHumanReadableChainCombi(chainCombi, true), set.getRt());
  				}
  				catch (LipidCombinameEncodingException ex)
  				{
  					ex.printStackTrace();
  				}
  			}
  		}
  	}
  	return displayName;
  }
  
  
  protected Object[][] generateTableData()
  {
  	Object[][] tableData = new Object[matchedIsotopologues_.size()][3];
  	int count=0;
  	for (MatchedPartnerVO match : matchedIsotopologues_)
  	{
  		tableData[count][0] = convertLipidomicsMSnSetToTableDisplayName(match, true);
    	tableData[count][1] = convertLipidomicsMSnSetToTableDisplayName(match, false);
    	tableData[count][EDITABLE_COLUMN] = match.isUseForCalibration();
    	count++;
  	}
  	return tableData;
  }
  
  /**
   * 
   * @param columnNames
   * @param tableData
   * @param editableColumn
   * @return
   */
  protected JPanel generateDisplayPanel(String[] columnNames, Object[][] tableData, int editableColumn)
  {
  	displayPanel_ = new JPanel();
  	displayPanel_.setLayout(new GridBagLayout());
    
  	DefaultTableModel model = new DataModel(tableData, columnNames);
  	model.addTableModelListener(new TableModelListener()
  	{
  		@Override
  		public void tableChanged(TableModelEvent e) {
  			MatchedPartnerVO vo = matchedIsotopologues_.get(e.getFirstRow());
  			boolean before = vo.isUseForCalibration();
  			vo.setUseForCalibration(!before);
  			isotopeEffectRegression_ = new IsotopeEffectRegression(matchedIsotopologues_, labelDefinitionPanel_.getAssignedIsotopeLabels());
  			updatePlot();
  	  }
  	});
    JTable standardsTable = new JTable(model);
    standardsTable.setRowHeight(30);

    JScrollPane scrollPane = new JScrollPane(standardsTable);
    scrollPane.setPreferredSize(TABLE_PANEL_DIMENSION);
    
    plot_ = new TotalIsotopeEffectPlot(matchedIsotopologues_, isotopeEffectRegression_, PLOT_DIMENSION);
    
    displayPanel_.add(scrollPane, new GridBagConstraints(0, 0, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    addPlotToDisplayPanel();
    displayPanel_.setBorder(getTitledPanelBorder(TABLE_FRAME_TITLE));
    return displayPanel_;
  }

  private void updatePlot()
  {
  	displayPanel_.remove(plot_);
  	plot_ = new TotalIsotopeEffectPlot(matchedIsotopologues_, isotopeEffectRegression_, PLOT_DIMENSION);
  	addPlotToDisplayPanel();
  	displayPanel_.invalidate();
  	displayPanel_.updateUI();
  }
  
  private void addPlotToDisplayPanel()
  {
  	displayPanel_.add(plot_, new GridBagConstraints(0, 1, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(25, 0, 0, 0), 0, 0));
  }
  
  protected Vector<MatchedPartnerVO> collectMatchedIsotopologues(Hashtable<String,Vector<LipidomicsMSnSet>> authenticStandards, ExperimentLabelDefinitionPanel labelDefinitionPanel)
	{	
  	Vector<MatchedPartnerVO> matched = new Vector<MatchedPartnerVO>();
		for (ResultFileVO resultFileVO : labelDefinitionPanel.getResultFiles())
		{
			QuantificationResult quantRes = resultFileVO.getQuantificationResult();
			//only classes for which there are authentic standards available are parsed.
			for (String lipidClass : authenticStandards.keySet()) 
			{
				Vector<LipidomicsMSnSet> standards = authenticStandards.get(lipidClass);
				Vector<LipidParameterSet> results = quantRes.getIdentifications().get(lipidClass);
				for (LipidParameterSet result : results)
				{
					//only deuterated species detected at MSn level with double bonds are relevant
					if (result instanceof LipidomicsMSnSet && result.getDoubleBonds()>0 && result.getChemicalFormula().contains("D"))
					{
						LipidomicsMSnSet resultMSn = (LipidomicsMSnSet) result;
						for (LipidomicsMSnSet standard : standards) 
						{
							IsotopeLabelVO labelVO = matchLabelToStandard(labelDefinitionPanel.getAssignedIsotopeLabels(), standard);
							if (isEquivalentCompound(standard, resultMSn, labelVO))
							{
								matched.add(new MatchedPartnerVO(lipidClass, standard, resultMSn, isPositionAssignmentMatch(standard, resultMSn, labelVO), labelVO));
							}
						}
					}
				}
			}
		}
		matched.sort(Comparator.comparing(MatchedPartnerVO::isUseForCalibration).reversed());
		return matched;
	}
  
  private boolean isPositionAssignmentMatch(LipidomicsMSnSet compound, LipidomicsMSnSet labeledCompound, IsotopeLabelVO labelVO)
  {
  	boolean matches = false;
  	Set<String> humanReadable1 = compound.getHumanReadableNameSet();
		Set<String> humanReadable2 = labeledCompound.getHumanReadableNameSet();
		for (String name : humanReadable2)
		{
			name = name.replaceAll(labelVO.getLabelId(), "");
			if (humanReadable1.contains(name))
			{
				matches = true;
			}
		}
  	return matches;
  }
  
  //TODO: as of now multiple labeled compounds are not taken into account!!!
	private boolean isEquivalentCompound(LipidomicsMSnSet compound, LipidomicsMSnSet labeledCompound, IsotopeLabelVO labelVO)
	{
		boolean isEquivalent = false;
		//compounds are the same molecular species (prefix needs to be removed for the check)
		if (compound.getNameStringWithoutRt().equals(labeledCompound.getNameStringWithoutRt().replaceAll(labelVO.getLabelId(), "")))
		{
			try
			{
				Hashtable<String,Integer> chemicalFormula1 = StaticUtils.categorizeFormula(compound.getChemicalFormula());
				Hashtable<String,Integer> chemicalFormula2 = replaceIsotopesInFormula(labeledCompound.getChemicalFormula(), labelVO);
				//compounds have the same chemical formula when elements are replaced according to the labelVO => double bond position is correct
				if (StaticUtils.isChemicalFormulaTheSame(chemicalFormula1,chemicalFormula2))
				{
					Vector<String> chainCombis1 = compound.getValidChainCombinations();
					Vector<String> chainCombis2 = labeledCompound.getValidChainCombinations();
					isEquivalent = isEquivalentChainCombinationPresent(chainCombis1, chainCombis2, labelVO);
				}
			}
			catch (ChemicalFormulaException ex)
			{
				ex.printStackTrace();
			}
		}
		return isEquivalent;
	}
	
	private boolean isEquivalentChainCombinationPresent(Vector<String> chainCombis1, Vector<String> labeledChainCombis, IsotopeLabelVO labelVO)
	{
		for (String chainCombi1 : chainCombis1)
		{
			for (String chainCombi2 : labeledChainCombis)
			{
				if (normChainCombination(chainCombi1, labelVO).equals(normChainCombination(chainCombi2, labelVO)))
				{
					return true;
				}
			}
		}
		return false;
	}
	
	/**
	 * Normalizes a chain combination by removing potential present labels, splitting it at the divider and sorting the chains.
	 * @param chainCombination
	 * @param labelVO
	 * @return
	 */
	private Vector<String> normChainCombination(String chainCombination, IsotopeLabelVO labelVO)
	{
		String removedLabels = chainCombination.replaceAll(labelVO.getLabelId(), "");
		Vector<String> split = new Vector<String>(Arrays.asList(removedLabels.split("<->")));
		split.sort(Comparator.nullsFirst(Comparator.comparing(String::length).thenComparing(Comparator.naturalOrder())));
		return split;
	}
	
	
	private Hashtable<String,Integer> replaceIsotopesInFormula(String chemicalFormulaLabeled, IsotopeLabelVO labelVO) throws ChemicalFormulaException
	{
		Hashtable<String,Integer> chemicalFormula = null;
		try
		{
			chemicalFormula = StaticUtils.categorizeFormula(chemicalFormulaLabeled);
			for (String element : labelVO.getLabelElements().keySet()) 
	    {
	      int nr = 0;
	      if (chemicalFormula.containsKey(element))
	        nr = chemicalFormula.get(element);
	      nr -= labelVO.getLabelElements().get(element);
	      if (nr==0)
	      	chemicalFormula.remove(element);
	      else
	      	chemicalFormula.put(element, nr);
	    }
		}
		catch (ChemicalFormulaException ex)
		{
			throw new ChemicalFormulaException(ex);
		}
    return chemicalFormula;
	}
	
  
  private IsotopeLabelVO matchLabelToStandard(Vector<IsotopeLabelVO> labels, LipidParameterSet standard)
	{
		IsotopeLabelVO matchedLabel = null;
		int omegaPosition = -1;
		//if there are more than one, they will either way all have the same assigned double bond position, so get index 0
		DoubleBondPositionVO doubleBondPositionVO = standard.getOmegaInformation().get(0);
		for (FattyAcidVO fattyAcid : doubleBondPositionVO.getChainCombination()) 
		{
			if (fattyAcid.getDoubleBonds()>0)
			{
				omegaPosition = fattyAcid.getOmegaPosition();
				break;
			}
		}
		for (IsotopeLabelVO label : labels)
		{
			if (label.getOmegaPosition() == omegaPosition)
			{
				matchedLabel = label;
				break;
			}
		}
		return matchedLabel;
	}
  
  /**
   * Reads all files in @param resultFiles, of which the index is a key in @param rowsToOmegaPositions 
	 * and assigns the corresponding value in @param rowsToOmegaPositions as omega double bond position 
	 * to all unsaturated fatty acid chains of analytes detected at MSn level in the result file.
	 * 
   * @param standardsFileVOs
   * @return
   * @throws LipidCombinameEncodingException
   */
	private Hashtable<String,Vector<LipidomicsMSnSet>> collectStandards(Vector<Pair<ResultFileVO, Integer>> standardsFileVOs) throws LipidCombinameEncodingException
	{
		Hashtable<String,Vector<LipidomicsMSnSet>> authenticStandards = new Hashtable<String,Vector<LipidomicsMSnSet>>();
		for (Pair<ResultFileVO, Integer> standard : standardsFileVOs)
		{
			QuantificationResult result = standard.getKey().getQuantificationResult();
			Hashtable<String,Vector<LipidParameterSet>> identifications = result.getIdentifications();
      LipidomicsMSnSet mSnSet;
      for (String lipidClass : identifications.keySet())
      {
      	if (!authenticStandards.containsKey(lipidClass))
      	{
      		authenticStandards.put(lipidClass, new Vector<LipidomicsMSnSet>());
      	}
      	for (LipidParameterSet authenticStandard : identifications.get(lipidClass)) 
      	{
      		if (authenticStandard instanceof LipidomicsMSnSet) //only detections at MSn level should be taken into account
      		{
      			mSnSet = (LipidomicsMSnSet)authenticStandard;
      			Set<String> humanReadableNameSet = mSnSet.getHumanReadableNameSet();
      			LinkedHashMap<String,String> lookup = mSnSet.getNameLookupHumReadableToPositionInsensitive();
      			Hashtable<String,FattyAcidVO> involvedFAs = mSnSet.getInvolvedFAs();
      			
      			for (String name : humanReadableNameSet) 
      			{
      				Vector<String> positionInsensitiveChains = StaticUtils.splitChainCombiToEncodedStrings(
      						lookup.get(name),LipidomicsConstants.CHAIN_COMBI_SEPARATOR);
      				Vector<FattyAcidVO> chainCombination = new Vector<FattyAcidVO>();
      				for (String chain : positionInsensitiveChains)
      				{
      					chainCombination.add(involvedFAs.get(chain));
      				}
      				//add omega double bond positions to all chainCombinations
      				for (FattyAcidVO fattyAcid : chainCombination) {
      					if (fattyAcid.getDoubleBonds() > 0)
      					{
      						fattyAcid.setOmegaPosition(standard.getValue());
      					}
      				}
      				DoubleBondPositionVO doubleBondPositionVO = new DoubleBondPositionVO(
      						chainCombination, 
      						Float.parseFloat(authenticStandard.getRt()),
      						DoubleBondPositionVO.ACCURACY_LOW,
      						name);
      				authenticStandard.addOmegaInformation(doubleBondPositionVO);
      				authenticStandards.get(lipidClass).add(mSnSet);
      			}
      		}
      	}
      }
		}
    return authenticStandards;
	}
  
  
  
  @Override
  protected void back() 
  {
  	try {getDefaultComponents().removeOptionPanel(getDefaultComponents().getCurrentPanel());} catch (Exception ex) {}
  	goBack();
  }
  
  @Override
  protected void next()
  {
  	boolean arePairsChosen = false;
  	for (MatchedPartnerVO vo : matchedIsotopologues_)
  	{
  		arePairsChosen = vo.isUseForCalibration() == Boolean.TRUE;
  		if (arePairsChosen) break;
  	}
  	if (arePairsChosen)
  	{
  		goNext();
  		try
  		{
  			ExportPanel panel = (ExportPanel)getDefaultComponents().getCurrentPanel();
  			panel.addExportOptionsPanel(new ExportOptionsPanel(true));
  			panel.setExporter(new TargetListExporter(isotopeEffectRegression_, labelDefinitionPanel_.getResultFiles(), labelDefinitionPanel_.getAssignedIsotopeLabels()));
  		} catch (Exception ex) {}
  	}
  	else
  	{
  		new WarningMessage(new JFrame(), "Warning", "No pair was chosen. This step is necessary to account for the isotope effect on retention time caused by deuterium labels!");
  	}
  }
  
  
  private class DataModel extends DefaultTableModel
  {
		private static final long serialVersionUID = 1L;

		public DataModel(Object[][] data, Object[] columnNames)
  	{
  		super(data, columnNames);
  	}
  	
  	@Override
    public Class<?> getColumnClass(int columnIndex) 
  	{
      if (columnIndex == EDITABLE_COLUMN) 
      {
        return Boolean.class;
      }
      return super.getColumnClass(columnIndex);
    }
  	
  	@Override
    public boolean isCellEditable(int row, int column) 
  	{
        return column == EDITABLE_COLUMN;
    }
  }
}
