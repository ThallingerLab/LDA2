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

package at.tugraz.genome.lda.target.calibration;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JLabel;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.target.JDefaultComponents;
import at.tugraz.genome.lda.target.JOptionPanel;
import at.tugraz.genome.lda.target.LoadingPanel;
import at.tugraz.genome.lda.target.export.ExportPanel;
import at.tugraz.genome.lda.target.export.TargetListExporter;
import javafx.util.Pair;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class CalibrationGraphPanel extends JOptionPanel
{
	private static final long serialVersionUID = 1L;
	
	private LoadingPanel loadingPanel_;
	private JPanel displayPanel_;
	private String[] lipidClasses_ = new String[0];
	private SubgroupDefinitionPanel subGroupDefinitionPanel_;
	private JComboBox<String> classListJComboBox_;
	private JPanel jComboBoxPanel_;
	private JCheckBox classSpecificJCheckBox_;
	private JButton defineSubgroupsJButton_;
	private JSlider granularityJSlider_;
	private Double grouping_;
	private RecalibrationPlot plot_;
	private ArrayList<RecalibrationRegression> regressions_;
	private File originalTargetList_;
	private Double predictionThreshold_;
	public static final String PLOT_ALL = "Combined";
	private static final Dimension PLOT_DIMENSION = new Dimension(825,650);
	private static final String TABLE_FRAME_TITLE = "Recalibrate the \u03C9-C=C target list";
	
	
  
  public CalibrationGraphPanel(JDefaultComponents wizardComponents) {
      super(wizardComponents, "Recalibrate an omega C=C target list to your chromatographic conditions.");
      init(generateLoadingPanel());
  }
  
  protected void init(JPanel panel) 
  {
  	this.grouping_ = 1.0;
    this.add(panel, new GridBagConstraints(0, 0, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    this.invalidate();
    this.updateUI();
  }
  
  public void initDataDisplay()
  {
  	cleanPanels();
  	getDefaultComponents().updateComponents();
  	init(generateDisplayPanel());
  }
  
  protected JPanel generateLoadingPanel()
  {
  	loadingPanel_ = new LoadingPanel("<html>Processing data, please wait...</html>");
  	return loadingPanel_;
  }
  
  protected JPanel generateDisplayPanel()
  {
  	displayPanel_ = new JPanel();
  	displayPanel_.setLayout(new GridBagLayout());
  	
  	initClassSpecificJCheckBox(); //put button here
  	initClassListJComboBox();
  	initGranularityJSlider();
  	initPlot();
    
    displayPanel_.setBorder(getTitledPanelBorder(TABLE_FRAME_TITLE));
    return displayPanel_;
  }
  
  private void initClassSpecificJCheckBox()
  {
  	JPanel panel = new JPanel();
  	panel.setLayout(new GridBagLayout());
  	
  	classSpecificJCheckBox_ = new JCheckBox("Calibrate lipid classes separately.", true); //other code relies on this being initialized with true
  	classSpecificJCheckBox_.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	classSpecificJCheckBox_actionPerformed(e);
		  }
	  });
  	panel.add(classSpecificJCheckBox_, new GridBagConstraints(0, 0, 0, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 5, 0), 0, 0));
  	
  	defineSubgroupsJButton_ = new JButton("Define Subgroups");
  	defineSubgroupsJButton_.setEnabled(classSpecificJCheckBox_.isSelected());
  	defineSubgroupsJButton_.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	defineSubgroupsJButton_actionPerformed(e);
		  }
	  });
  	panel.add(defineSubgroupsJButton_, new GridBagConstraints(0, 1, 0, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
  	
  	addToDisplayPanel(panel, new GridBagConstraints(10, 0, 0, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
  }
  
  private void initClassListJComboBox()
  {
  	jComboBoxPanel_ = new JPanel();
  	jComboBoxPanel_.setLayout(new GridBagLayout());
  	
  	JLabel label = new JLabel("Select displayed lipid class / group: ");
  	jComboBoxPanel_.add(label, new GridBagConstraints(0, 0, 0, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 5, 0), 0, 0));
  	
  	initClassListJComboBox(lipidClasses_);
  	
  	addToDisplayPanel(jComboBoxPanel_, new GridBagConstraints(0, 0, 0, 2, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 15, 5, 5), 0, 0));
  }
  
  private void initClassListJComboBox(String[] toDisplay)
  {
  	classListJComboBox_ = new JComboBox<String>(toDisplay);
  	classListJComboBox_.setSelectedIndex(0);
  	classListJComboBox_.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	classListJComboBox_actionPerformed(e);
		  }
	  });
  	jComboBoxPanel_.add(classListJComboBox_, new GridBagConstraints(0, 1, 0, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
  }
  
  /**
   * Initializes the granularity slider with values from 0 to 100 with the initial value set to 50.
   * When translating these values to a clustering range from 0.0 to 2.0 min, the initial value would here be 1.0.
   */
  private void initGranularityJSlider()
  {
  	JPanel panel = new JPanel();
  	panel.setLayout(new GridBagLayout());
  	JLabel label = new JLabel("Granularity of fit:");
  	panel.add(label, new GridBagConstraints(0, 0, 0, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 5, 0), 0, 0));
  	
  	JSlider slider = new JSlider(0, 100, 50);
  	slider.setMajorTickSpacing(50);
  	slider.setMinorTickSpacing(10);
  	slider.setPaintTicks(true);
  	
  	// Add positions label in the slider
    Hashtable<Integer, JLabel> position = new Hashtable<Integer, JLabel>();
    position.put(0, new JLabel("0%"));
    position.put(50, new JLabel("50%"));
    position.put(100, new JLabel("100%"));
     
    // Set the label to be drawn
    slider.setLabelTable(position);
    slider.setPaintLabels(true);
  	
    slider.addMouseListener(new MouseListener()
    {
			@Override
			public void mouseReleased(MouseEvent e)
			{
				int max = 100;
				Double maxMinutes = 2.0;
				Double minMinutes = 0.2;
				int value = (((JSlider)e.getSource()).getValue());
				int diff = max - value;
				grouping_ = maxMinutes * diff / max;
				if (grouping_ < minMinutes) grouping_ = minMinutes; //to make sure the data points are a monotonic sequence
				System.out.println(grouping_+" value: "+value);
				showViewOfChoice();
			} 
			
			public void mouseClicked(MouseEvent e) {}
			public void mouseEntered(MouseEvent e) {}
			public void mouseExited(MouseEvent e) {}
			public void mousePressed(MouseEvent e) {}
    });
    
  	granularityJSlider_ = slider;
  	panel.add(granularityJSlider_, new GridBagConstraints(0, 1, 0, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
  	
  	addToDisplayPanel(panel, new GridBagConstraints(5, 0, 10, 2, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
  }
  
  /**
   * Initiates the recalibration plot.
   */
  private void initPlot()
  {
  	RecalibrationRegression regression = getRegressionByFields(PLOT_ALL, PLOT_ALL);
  	ArrayList<Pair<Double,Double>> data = regression == null ? new ArrayList<Pair<Double,Double>>() : regression.getDifferences();
  	RecalibrationRegression regressionStandards = getRegressionByFields(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX, PLOT_ALL);
  	ArrayList<Pair<Double,Double>> dataStandards = regressionStandards == null ? new ArrayList<Pair<Double,Double>>() : regressionStandards.getDifferences();
    plot_ = new RecalibrationPlot(data, dataStandards, regression, PLOT_DIMENSION, this, null);
    addPlotToDisplayPanel();
  }
  
  /**
   * Adds the recalibration plot to the display panel. The plot needs to be generated prior to calling this method.
   */
  private void addPlotToDisplayPanel()
  {
  	addToDisplayPanel(plot_, new GridBagConstraints(0, 3, 15, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 5, 5, 5), 0, 0));
  }
  
  /**
   * Adds a component to the display panel with the given constraints.
   * @param component
   * @param constraints
   */
  private void addToDisplayPanel(JComponent component, GridBagConstraints constraints)
  {
  	displayPanel_.add(component, constraints);
  }
  
  /**
   * Removes all panels; needed when going back to the previous JOptionsPanel
   */
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
  	this.subGroupDefinitionPanel_ = null;
  }
  
  /**
   * Shows the plot selected in the class list JComboBox
   */
  public void showViewOfChoice()
  {
  	String selectedItem = (String)classListJComboBox_.getSelectedItem(); 
  	updatePlot(selectedItem);
  }
  
  /**
   * Handles the events called when triggering the class list JComboBox.
   * @param e
   */
  private void classListJComboBox_actionPerformed(ActionEvent e) 
  {
  	showViewOfChoice();
  }
  
  /**
   * Handles the events called when triggering the class specific JCheckBox.
   * @param e
   */
  private void classSpecificJCheckBox_actionPerformed(ActionEvent e)
  {
  	defineSubgroupsJButton_.setEnabled(classSpecificJCheckBox_.isSelected());
  	updateClassListJComboBox();
  }
  
  /**
   * The class list JComboBox contains the lipid classes and names of subgroups that will be used for recalibration; 
   * the combined regression is always selectable.
   * This method adjusts this list depending on user settings.
   */
  public void updateClassListJComboBox()
  {
  	jComboBoxPanel_.remove(classListJComboBox_);
  	if (!classSpecificJCheckBox_.isSelected())
  	{
  		initClassListJComboBox(new String[] {PLOT_ALL});
  	}
  	else if (subGroupDefinitionPanel_ != null)
  	{
  		ArrayList<String> ungroupedLipidClasses = subGroupDefinitionPanel_.getUngroupedLipidClasses();
  		ArrayList<SubGroup> subGroups = subGroupDefinitionPanel_.getDefinedSubgroups();
  		String[] toDisplay = new String[ungroupedLipidClasses.size()+subGroups.size()+1];
  		toDisplay[0] = PLOT_ALL;
  		for (int i=0; i<subGroups.size(); i++)
  		{
  			toDisplay[i+1] = subGroups.get(i).getGroupName();
  		}
  		for (int i=0; i<ungroupedLipidClasses.size();i++)
  		{
  			toDisplay[i+subGroups.size()+1] = ungroupedLipidClasses.get(i);
  		}
  		jComboBoxPanel_.remove(classListJComboBox_);
  		initClassListJComboBox(toDisplay);
  	}
  	else
  	{
  		initClassListJComboBox(lipidClasses_);
  	}
  	displayPanel_.invalidate();
  	displayPanel_.updateUI();
  }
  
  /**
   * Given a lipid class name, this method returns the relevant key to retrieve the correct regression.
   * If no class specific regression is available, the id for the combined regression is returned.
   * @param cName
   * @return
   */
  public String getRelevantRegressionName(String cName)
  {
  	if (subGroupDefinitionPanel_ != null)
  	{
  		for (String name : subGroupDefinitionPanel_.getUngroupedLipidClasses())
  		{
  			if (name.equals(cName)) return cName;
  		}
  		for (SubGroup group : subGroupDefinitionPanel_.getDefinedSubgroups())
  		{
  			for (String name : group.getLipidClasses())
  			{
  				if (name.equals(cName)) return group.getGroupName();
  			}
  		}
  	}
  	else
  	{
  		for (int i=0;i<lipidClasses_.length;i++)
    	{
    		if (cName.equals(lipidClasses_[i])) return cName;
    	}
  	}
  	return PLOT_ALL;
  }
  
  /**
   * Handles the events called when triggering the define subgroups JButton.
   * @param e
   */
  private void defineSubgroupsJButton_actionPerformed(ActionEvent e)
  {
  	subGroupDefinitionPanel_ = new SubgroupDefinitionPanel(
  			subGroupDefinitionPanel_ == null ? new ArrayList<SubGroup>() : subGroupDefinitionPanel_.getDefinedSubgroups(), lipidClasses_, this);
  }
  
  /**
   * Removes all regressions with an id referring to defined sub groups.
   */
  public void removeSubGroupRegressions()
  {
  	ArrayList<RecalibrationRegression> toRemove = new ArrayList<RecalibrationRegression>();
  	for (SubGroup group : subGroupDefinitionPanel_.getDefinedSubgroups())
  	{
  		for (RecalibrationRegression regression : regressions_)
  		{
  			if (regression.getLipidClass().equals(group.getGroupName()))
  			{
  				toRemove.add(regression);
  			}
  		}
  	}
  	regressions_.removeAll(toRemove);
  }
  
  /**
   * Generates regressions for the subgroups.
   * @param subGroups
   */
  public void addSubGroupRegressions(ArrayList<SubGroup> subGroups)
  {
  	for (SubGroup group : subGroups)
  	{
  		ArrayList<Pair<Double,Double>> differencesStandards = new ArrayList<Pair<Double,Double>>();
  		ArrayList<Pair<Double,Double>> differences = new ArrayList<Pair<Double,Double>>();
  		for (String lipidClass : group.getLipidClasses())
  		{
  			RecalibrationRegression regressionStandards = getRegressionByFields(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX, lipidClass);
  			if (regressionStandards != null) differencesStandards.addAll(regressionStandards.getDifferences());
  			RecalibrationRegression regression = getRegressionByFields(PLOT_ALL, lipidClass);
  			if (regression != null) differences.addAll(regression.getDifferences());
  		}
  		RecalibrationRegression regressionStandardsForGroup = new RecalibrationRegression(
  				differencesStandards, grouping_, CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX, group.getGroupName());
  		regressions_.add(regressionStandardsForGroup);
  		RecalibrationRegression regressionForGroup = new RecalibrationRegression(
  				differences, grouping_, PLOT_ALL, group.getGroupName());
  		regressions_.add(regressionForGroup);
  	}
  	fillMissingDataPoints(regressions_);
  }
  
  private void updatePlot(String className)
  {
  	displayPanel_.remove(plot_);
  	for (RecalibrationRegression reg : regressions_)
  	{
  		reg.setGrouping(grouping_);
  	}
  	fillMissingDataPoints(regressions_);
  	RecalibrationRegression regression = getRegressionByFields(PLOT_ALL, className);
  	ArrayList<Pair<Double,Double>> data = regression == null ? new ArrayList<Pair<Double,Double>>() : regression.getDifferences();
  	RecalibrationRegression regressionStandards = getRegressionByFields(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX, className);
  	ArrayList<Pair<Double,Double>> dataStandards = regressionStandards == null ? new ArrayList<Pair<Double,Double>>() : regressionStandards.getDifferences();
  	plot_ = new RecalibrationPlot(data, dataStandards, regression, PLOT_DIMENSION, this, plot_.getXYPlot());
  	addPlotToDisplayPanel();
  	displayPanel_.invalidate();
  	displayPanel_.updateUI();
  }
  
  /**
   * @param originalConditions
   * @param newConditions
   */
  public void parseData(Hashtable<String,ArrayList<File>> originalConditions, Hashtable<String,ArrayList<File>> newConditions) throws ExcelInputFileException
  {
  	this.regressions_ = new ArrayList<RecalibrationRegression>();
  	
  	if (originalConditions.containsKey(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX) &&
  			newConditions.containsKey(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX))
  	{
  		this.regressions_ = parseStandardMix(
  				originalConditions.get(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX), 
  				newConditions.get(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX));
  	}
  	this.regressions_.addAll(parseAllResults(originalConditions, newConditions));
  }
  
  
  /**
   * Parses the standard mix data to obtain all needed information to compute a recalibration regression for each lipid class as well as all calibrants combined.
   * At this point, no check is done regarding whether enough data points are available for each regression. The regression function will be null in this case.
   * @param originalMix
   * @param newMix
   * @return
   * @throws ExcelInputFileException
   */
  private ArrayList<RecalibrationRegression> parseStandardMix(ArrayList<File> originalMix, ArrayList<File> newMix) throws ExcelInputFileException
  {
  	ArrayList<RecalibrationRegression> regressions = new ArrayList<RecalibrationRegression>();
  	ArrayList<IdentificationVO> originalIdentifications = parseResultFiles(originalMix);
		ArrayList<IdentificationVO> newIdentifications = parseResultFiles(newMix);
		ArrayList<MatchedIdentificationVO> matches = computeMatches(originalIdentifications, originalMix.size(), newIdentifications, newMix.size());
		Hashtable<String, ArrayList<Pair<Double,Double>>> differencesForClass = computeDifferencesForClass(matches, 0);
		ArrayList<Pair<Double,Double>> differences = new ArrayList<Pair<Double,Double>>();
		
		lipidClasses_ = new String[differencesForClass.keySet().size()];
		int count = 0;
		for (String lipidClass : differencesForClass.keySet())
		{
			lipidClasses_[count] = lipidClass;
			regressions.add(new RecalibrationRegression(differencesForClass.get(lipidClass), grouping_, CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX, lipidClass));
			differences.addAll(differencesForClass.get(lipidClass));
			count++;
		}
		
		regressions.add(new RecalibrationRegression(differences, grouping_, CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX, PLOT_ALL));
		return regressions;
  }
  
  /**
   * Parses all data and produces combined recalibration regressions for for each lipid class (wherever enough data points are available)
   * as well as fits for all calibrants combined. If not enough combined calibrants are available, an error message is thrown.
   * @param originalConditions
   * @param newConditions
   * @return
   * @throws ExcelInputFileException
   */
  private ArrayList<RecalibrationRegression> parseAllResults(Hashtable<String,ArrayList<File>> originalConditions, Hashtable<String,ArrayList<File>> newConditions) throws ExcelInputFileException
  {
  	ArrayList<RecalibrationRegression> regressions = new ArrayList<RecalibrationRegression>();
  	ArrayList<MatchedIdentificationVO> matches = computeMatchesForDataTypeOther(originalConditions, newConditions);
  	HashSet<String> uniqueLipidClasses = computeCombinedRegressionsPerClass(matches, regressions);
  	fillMissingDataPoints(regressions);
  	ArrayList<String> sortedLipidClasses = new ArrayList<String>(uniqueLipidClasses);
  	Collections.sort(sortedLipidClasses);
  	sortedLipidClasses.add(0, PLOT_ALL); //adding PLOT_ALL at the beginning, ensuring it is always the first element
		lipidClasses_ = sortedLipidClasses.toArray(new String[sortedLipidClasses.size()]);	
  	
  	return regressions;
  }
  
  /**
   * Adds data points of the combined plot to fill in missing data 
   * @param regressions
   */
  private void fillMissingDataPoints(ArrayList<RecalibrationRegression> regressions)
  {
  	RecalibrationRegression regressionCombined = getRegressionByFields(regressions, PLOT_ALL, PLOT_ALL);
  	ArrayList<Pair<Double, Double>> clusteredCombined = regressionCombined.getClustered();
  	Double minKey = regressionCombined.getClustered().get(0).getKey();
  	Double maxKey = clusteredCombined.get(clusteredCombined.size()-1).getKey();
  	for (RecalibrationRegression regression : regressions)
  	{
  		if (regression.getLipidClass().equals(PLOT_ALL) || regression.getDataType().equals(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX)) continue;
  		
  		ArrayList<Pair<Double, Double>> clustered = regression.getClustered();
  		ArrayList<Pair<Double, Double>> clusteredToAdd = new ArrayList<Pair<Double, Double>>();
  		
  		//adding first and last datapoints
  		Double firstKey = clustered.get(0).getKey();
  		Double lastKey = clustered.get(clustered.size()-1).getKey();
  		Double previousKey = firstKey;
  		if (minKey < firstKey-(this.grouping_*0.1))
  		{
  			clustered.add(0, regressionCombined.getClustered().get(0));
  			previousKey = minKey;
  		}
  		if (maxKey > lastKey+(this.grouping_*0.1))
  		{
  			clustered.add(clusteredCombined.get(clusteredCombined.size()-1));
  		}
  		
  		//adding data points in between
  		for (Pair<Double, Double> cluster : clustered)
			{
  			Double currentKey = cluster.getKey();
  			Double diff = currentKey-previousKey; //previous key is always smaller than the current, due to the clusters being ordered
  			if (diff >= (this.grouping_*2))
  			{
  				clusteredToAdd.addAll(getClusterInBetween(clusteredCombined, previousKey+(this.grouping_*0.75), currentKey-(this.grouping_*0.75)));
  			}
  			previousKey = currentKey;
			}
  		
  		clustered.addAll(clusteredToAdd);
  		Collections.sort(clustered, new Comparator<Pair<Double, Double>>() 
  		{
        @Override
        public int compare(Pair<Double, Double> p1, Pair<Double, Double> p2) 
        {
        	return p1.getKey().compareTo(p2.getKey());
        }
  		});
  		regression.initRegression(clustered);
  	}
  }
  
  private ArrayList<Pair<Double, Double>> getClusterInBetween(ArrayList<Pair<Double, Double>> clustered, Double before, Double after)
  {
  	ArrayList<Pair<Double, Double>> clusteredInBetween = new ArrayList<Pair<Double, Double>>();
  	for (Pair<Double, Double> cluster : clustered)
		{
  		if (cluster.getKey() >= before && cluster.getKey() <= after)
  		{
  			clusteredInBetween.add(cluster);
  		}
		}
  	return clusteredInBetween;
  }
  
//  private Double getMean(List<Double> values) 
//  {
//  	Double sum = 0.0;
//    for (Double value : values) 
//    {
//      sum += value;
//    }
//    return (sum / values.size());
//  }
//  
//  private Double getVariance(List<Double> values) 
//  {
//  	Double mean = getMean(values);
//    Double temp = 0.0;
//    for (Double a : values) 
//    {
//      temp += (a - mean) * (a - mean);
//    }
//    return temp / (values.size() - 1.0);
//  }
//  
//  private Double getStdDev(List<Double> values) 
//  {
//    return Math.sqrt(getVariance(values));
//  }
  
  
  
  /**
   * Adds all valid regressions computed from all calibrants of each lipid class to @param regressions.
   * If not enough combined clustered calibrants are available to compute even a combined fit, an error message is shown and the user is redirected.
   * @param matches
   * @param regressions
   * @return A set of all lipid classes for which valid regressions exist.
   */
  private HashSet<String> computeCombinedRegressionsPerClass(ArrayList<MatchedIdentificationVO> matches, ArrayList<RecalibrationRegression> regressions)
  {
    //in case we have valid regressions for additional lipid classes, we need to redefine the lipidClasses_ array
  	HashSet<String> uniqueLipidClasses = new HashSet<String>();
  	Hashtable<String, ArrayList<Pair<Double,Double>>> differencesForClass = computeDifferencesForClass(matches, 1);
  	ArrayList<Pair<Double,Double>> differencesAll = new ArrayList<Pair<Double,Double>>();
		for (String lipidClass : differencesForClass.keySet()) //there might be different lipid classes here than in the standard mix
		{
			addCombinedRegression(differencesForClass.get(lipidClass), lipidClass, uniqueLipidClasses, regressions);
			differencesAll.addAll(differencesForClass.get(lipidClass));
		}
		
		//adding the combined fit, PLOT_ALL is not added to the unique lipid classes
		boolean allSuccessful = addCombinedRegression(differencesAll, PLOT_ALL, new HashSet<String>(), regressions);
		if (!allSuccessful)
		{
			new WarningMessage(new JFrame(), "Error", String.format("Not enough data points are available! Number of data points: %s. Please check your input data and the selected data types!", differencesAll.size()));
  		back(); //nothing can be plotted at this point.
		}
		
		//adding valid regressions that may only exist in the standard mix as PLOT_ALL
		for (int i=0; i<lipidClasses_.length; i++)
		{
			if (!uniqueLipidClasses.contains(lipidClasses_[i]))
			{
				addCombinedRegression(new ArrayList<Pair<Double,Double>>(), lipidClasses_[i], uniqueLipidClasses, regressions);
			}
		}
		return uniqueLipidClasses;
  }
  
  /**
   * For class specific recalibration
   * @param referenceResults
   * @param targetResults
   * @return
   */
  private Hashtable<String, ArrayList<Pair<Double,Double>>> computeDifferencesForClass(ArrayList<MatchedIdentificationVO> matches, int acceptedConfidence)
  {
  	Hashtable<String, ArrayList<Pair<Double,Double>>> differencesForClass = new Hashtable<String, ArrayList<Pair<Double,Double>>>();
  	
  	for (MatchedIdentificationVO match : matches)
		{
  		if (match.getConfidence() < acceptedConfidence) continue;
  		ArrayList<Pair<IdentificationVO,IdentificationVO>> acceptedMatches = new ArrayList<Pair<IdentificationVO,IdentificationVO>>();
  		String lipidClass = match.getHighestConfidencePair().getKey().getLipidClass();
  		if (match.getConfidence() < 2)
  		{
  			acceptedMatches.add(match.getHighestConfidencePair());
  		}
  		else if (match.getConfidence() == 2)
  		{
  			acceptedMatches = match.getAcceptedMatches();
  		}
  		if (!differencesForClass.containsKey(lipidClass))
  		{
  			differencesForClass.put(lipidClass, new ArrayList<Pair<Double,Double>>());
  		}
  		
  		for (Pair<IdentificationVO,IdentificationVO> matchedPair : acceptedMatches)
  		{
  			Double original = matchedPair.getKey().getAverageRT();
  			Double difference = original - matchedPair.getValue().getAverageRT();
  			differencesForClass.get(lipidClass).add(new Pair<Double,Double>(original, difference));
  		}
		}
  	return differencesForClass;
  }
  
  /**
   * Adds a regression consisting of @param differences and the data points from a possibly existing standard mix to @param regressions if the regression is viable (if there are enough data points).
   * If the regression is viable, then the given @param lipidClass is also added to @param uniqueLipidClasses.
   * @param differences
   * @param lipidClass
   * @param uniqueLipidClasses
   * @param regressions
   * @return true if a valid regression was successfully added, otherwise false.
   */
  private boolean addCombinedRegression(ArrayList<Pair<Double,Double>> differences, String lipidClass, HashSet<String> uniqueLipidClasses, ArrayList<RecalibrationRegression> regressions)
  {
		RecalibrationRegression regressionStandards = getRegressionByFields(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX, lipidClass);
		if (regressionStandards != null)
		{
			differences.addAll(regressionStandards.getDifferences());
		}
		RecalibrationRegression regression = new RecalibrationRegression(differences, grouping_, PLOT_ALL, lipidClass);
		if (regression.getFunction() != null)
		{
			regressions.add(new RecalibrationRegression(differences, grouping_, PLOT_ALL, lipidClass));
			uniqueLipidClasses.add(lipidClass);
			return true;
		}
		else
		{
			return false;
		}
  }
  
  
  /**
   * Computes matches for data types other than standard mix
   * @param originalConditions
   * @param newConditions
   * @return
   * @throws ExcelInputFileException
   */
  private ArrayList<MatchedIdentificationVO> computeMatchesForDataTypeOther(Hashtable<String,ArrayList<File>> originalConditions, Hashtable<String,ArrayList<File>> newConditions) throws ExcelInputFileException
  {
  	ArrayList<MatchedIdentificationVO> matches = new ArrayList<MatchedIdentificationVO>();
  	//adding matches for each data type to the same container, the CalibrationFileChooserPanel makes sure all data types are available for both experimental conditions
  	for (String dataType : originalConditions.keySet())
		{
			if (dataType.equals(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX)) continue;
			ArrayList<IdentificationVO> originalIdentifications = parseResultFiles(originalConditions.get(dataType));
			ArrayList<IdentificationVO> newIdentifications = parseResultFiles(newConditions.get(dataType));
			matches.addAll(computeMatches(originalIdentifications, originalConditions.get(dataType).size(), newIdentifications, newConditions.get(dataType).size()));
		}
  	return matches;
  }
  
  
  private ArrayList<MatchedIdentificationVO> computeMatches(ArrayList<IdentificationVO> originalResults, int maxOriginal, ArrayList<IdentificationVO> newResults, int maxNew)
  {
  	ArrayList<MatchedIdentificationVO> matches = new ArrayList<MatchedIdentificationVO>();
  	Hashtable<String, ArrayList<IdentificationVO>> groupedIdentifications = new Hashtable<String, ArrayList<IdentificationVO>>();
  	String originalIdentification = "original";
  	String newIdentification = "new";
  	
  	groupIdentifications(originalResults, groupedIdentifications, originalIdentification);
  	groupIdentifications(newResults, groupedIdentifications, newIdentification);
  	
  	for (String id : groupedIdentifications.keySet())
  	{
  		if (!id.startsWith(originalIdentification)) continue;
  		ArrayList<IdentificationVO> originalIdentifications = groupedIdentifications.get(id);
  		String partnerId = id.replace(originalIdentification, newIdentification);
  		ArrayList<IdentificationVO> newIdentifications = groupedIdentifications.get(partnerId);
  		
  		if (newIdentifications != null)
  		{
  			RecalibrationRegression standardsReg = getRegressionByFields(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX, PLOT_ALL);
  			matches.add(new MatchedIdentificationVO(originalIdentifications, maxOriginal, newIdentifications, maxNew, standardsReg, predictionThreshold_));
  		}
  	}
  	return matches;
  }
  
  
  /**
   * Groups identifications based on lipid class and lipid species.
   * @param identifications
   * @param groupedIdentifications
   * @param origin
   */
  private void groupIdentifications(ArrayList<IdentificationVO> identifications, Hashtable<String, ArrayList<IdentificationVO>> groupedIdentifications, String origin)
  {
  	for (IdentificationVO identification : identifications)
  	{
  		if (!identification.isIdentificationUnambiguous())
  		{
  			//TODO: should this be a warning... or should this not be reported at all?
  			System.out.format("For the lipid class: %s, the lipid species: %s was found more than once in an input file. "
  					+ "The grouping parameter needs to be adjusted to not group these different identifications together. \n", 
  					identification.getLipidClass(), identification.getLipidSpecies());
  		}
  		else
  		{
  			String id = String.format("%s, %s, %s", origin, identification.getLipidClass(), identification.getLipidSpecies());
    		if (!groupedIdentifications.containsKey(id))
    		{
    			groupedIdentifications.put(id, new ArrayList<IdentificationVO>());
    		}
    		groupedIdentifications.get(id).add(identification);
  		}
  	}
  }
  
  private ArrayList<IdentificationVO> parseResultFiles(ArrayList<File> files) throws ExcelInputFileException
  {
  	ArrayList<IdentificationVO> identifications = new ArrayList<IdentificationVO>();
  	for (File file : files)
  	{
  		QuantificationResult quantResOriginal = LDAResultReader.readResultFile(file.getAbsolutePath(), new Hashtable<String,Boolean>());
  		computeIdentificationVOs(file, quantResOriginal, identifications);
  	}
  	return identifications;
  }
	
  /**
   * Groups all analytes of the same lipid class, lipid species and retention time in @param identifications.
   * @param file
   * @param quantRes
   * @param identifications
   */
	private void computeIdentificationVOs(File file, QuantificationResult quantRes, ArrayList<IdentificationVO> identifications)
  {
  	Hashtable<String,Vector<LipidParameterSet>> quantIdentifications = quantRes.getIdentifications();
  	for (String lipidClass : quantIdentifications.keySet())
  	{
  		for (LipidParameterSet param : quantIdentifications.get(lipidClass))
  		{
  			boolean added = false;
  			for (int i = 0; i<identifications.size(); i++)
  			{
  				IdentificationVO identification = identifications.get(i);
  				if (identification.getLipidClass().equals(lipidClass) && identification.getLipidSpecies().equals(param.getNameStringWithoutRt()))
  				{
  					added = identification.addParam(file, param) ? true : added;
  				}
  			}
  			if (added == false)
  			{
  				IdentificationVO identificationVO = new IdentificationVO(file, lipidClass, param);
					identifications.add(identificationVO);
  			}
  		}
  	}
  }
  
	/**
	 * @param dataType						either CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX or PLOT_ALL
	 * @param lipidClass
	 * @return
	 */
  public RecalibrationRegression getRegressionByFields(String dataType, String lipidClass)
  {
  	return getRegressionByFields(regressions_, dataType, lipidClass);
  }
  
  /**
   * @param regressions					the ArrayList to search the desired regression in
   * @param dataType						either CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX or PLOT_ALL
   * @param lipidClass
   * @return
   */
  private RecalibrationRegression getRegressionByFields(ArrayList<RecalibrationRegression> regressions, String dataType, String lipidClass)
  {
  	for (RecalibrationRegression regression : regressions)
  	{
  		if (regression.getDataType().equals(dataType) && regression.getLipidClass().equals(lipidClass))
  		{
  			return regression;
  		}
  	}
  	return null;
  }
  
  public String findLipidClassForDataPoint(Pair<Double,Double> dataPoint)
  {
  	for (int i=0; i<lipidClasses_.length; i++) 
		{
  		if (lipidClasses_[i].equals(PLOT_ALL)) continue;
  		RecalibrationRegression regressionForClass = getRegressionByFields(PLOT_ALL, lipidClasses_[i]);
  		if (regressionForClass.getDifferences().contains(dataPoint))
  		{
  			return lipidClasses_[i];
  		}
		}
  	return null;
  }
  
  private boolean removeDataPointFromRegression(RecalibrationRegression regression, Pair<Double,Double> dataPoint)
  {
  	if (regression.getClusteredWithoutDataPointSize(dataPoint) > 2)
  	{
  		regression.removeDataPoint(dataPoint);
  		return true;
  	}
  	return false;
  }
  
  protected void removeDataPoint(String lipidClass, Pair<Double,Double> dataPoint)
  {
  	RecalibrationRegression regressionAll = getRegressionByFields(PLOT_ALL, PLOT_ALL);
  	RecalibrationRegression regressionClass = getRegressionByFields(PLOT_ALL, lipidClass);
  	RecalibrationRegression regressionSubGroup = getRegressionByFields(PLOT_ALL, getRelevantRegressionName(lipidClass));
  	if (removeDataPointFromRegression(regressionAll, dataPoint)) //assumption that there will not be too few data points in such cases.
  	{
  		if (!regressionClass.equals(regressionSubGroup))
  		{
  			removeDataPointFromRegression(regressionSubGroup, dataPoint); //assumption that there will not be too few data points in such cases.
  		}
  		if (!removeDataPointFromRegression(regressionClass, dataPoint))
  		{
  			new WarningMessage(new JFrame(), "Warning", String.format("Class specific calibration of %s is not possible anymore due to the removed data point!", getRelevantRegressionName(lipidClass)));
    		List<String> lipidClasses = new ArrayList<String>(Arrays.asList(lipidClasses_));
    		lipidClasses.remove(lipidClass);
    		lipidClasses_ = lipidClasses.toArray(new String[lipidClasses.size()]);
    		regressions_.remove(regressionClass);
    		displayPanel_.remove(jComboBoxPanel_);
    		initClassListJComboBox();
  		}
  	}
  	else
  	{
  		new WarningMessage(new JFrame(), "Error", "The fit requires at least 3 data points.");
  	}
  	showViewOfChoice();
  }
  
  @Override
  protected void next()
  {
  	goNext();
		try
		{
			ExportPanel panel = (ExportPanel)getDefaultComponents().getCurrentPanel();
			panel.setExporter(new TargetListExporter(originalTargetList_.getAbsolutePath(), classSpecificJCheckBox_.isSelected(), this));
		} catch (Exception ex) {}
  }
  
  @Override
  protected void back()
  {
  	cleanPanels();
  	init(generateLoadingPanel());
  	goBack();
  }
  
  public void setOriginalTargetList(File originalTargetList)
  {
  	this.originalTargetList_ = originalTargetList;
  }

  
  public void setPredictionThreshold(Double predictionThreshold)
  {
  	this.predictionThreshold_ = predictionThreshold;
  }
}
