/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.table.TableColumnModel;

import at.tugraz.genome.lda.verifier.IntegerMaxVerifier;
import at.tugraz.genome.lda.swing.BatchQuantificationTable;
import at.tugraz.genome.lda.swing.BatchQuantificationTableModel;;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class QuantificationMenu extends JPanel
{
	private static final long serialVersionUID = 1L;
	JTextField selectedMzxmlDirectory_; 
	JTextField selectedQuantDir_; 
	JTextField timeMinusTol_; 
	JTextField timePlusTol_; 
	JTextField cutoff_; 
	JTextField rtShift_; 
	JCheckBox isoValidation_; 
	JTextField amountOfIsotopes_; 
	JTextField amountOfMatchingSearchIsotopes_; 
	JCheckBox searchUnknownTime_; 
	JTextField nrProcessorsChrom_; 
	JTextField nrProcessors_; 
	/** the ion mode for the Alex123 searches*/
	JComboBox<String> ionMode_;
	JButton startQuantification_;
	JPanel quantifyingPanel_;
	
	JLabel quantifyingLabel_;
	JProgressBar progressBar_;
	JLabel spinnerLabel_;
	
	BatchQuantificationTableModel batchQuantTableModel_;
	BatchQuantificationTable batchQuantTable_;
	
	public QuantificationMenu(boolean isBatch, LipidDataAnalyzer parent) 
	{
		this.setLayout(new GridBagLayout());
		
		Integer y=0;
  	JPanel selectionPanel = addJPanelToQuantMenu(++y, this);
  	
    JLabel mzXMLLabel = new JLabel(isBatch ? "Folder with raw file(s): " : "Raw file: ");
    mzXMLLabel.setToolTipText(isBatch ? TooltipTexts.QUANTITATION_BATCH_RAW_FILE : TooltipTexts.QUANTITATION_SINGLE_RAW_FILE);
    selectionPanel.add(mzXMLLabel,new GridBagConstraints(0, y, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    selectedMzxmlDirectory_ = new JTextField(54);
    selectedMzxmlDirectory_.setMinimumSize(selectedMzxmlDirectory_.getPreferredSize());
    selectedMzxmlDirectory_.setToolTipText(isBatch ? TooltipTexts.QUANTITATION_BATCH_RAW_FILE : TooltipTexts.QUANTITATION_SINGLE_RAW_FILE);
    selectionPanel.add(selectedMzxmlDirectory_,new GridBagConstraints(1, y, 6, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JButton jButtonMzXMLOpen = new JButton("Select");
    jButtonMzXMLOpen.addActionListener(parent);
    jButtonMzXMLOpen.setActionCommand(isBatch ? "showMzxmlDirChooser" : "showMzxmlFileChooser");
    jButtonMzXMLOpen.setToolTipText(isBatch ? TooltipTexts.QUANTITATION_BATCH_RAW_FILE : TooltipTexts.QUANTITATION_SINGLE_RAW_FILE);
    selectionPanel.add(jButtonMzXMLOpen,new GridBagConstraints(7, y, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel quantChromLabel = new JLabel(isBatch ? "Folder with mass list or RT-DB file(s): " : "Mass list or RT-DB file: ");
    quantChromLabel.setToolTipText(isBatch ? TooltipTexts.QUANTITATION_BATCH_MASS_LIST : TooltipTexts.QUANTITATION_SINGLE_MASS_LIST);
    
    y++;
    selectionPanel.add(quantChromLabel,new GridBagConstraints(0, y, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JButton jButtonQuantOpen = new JButton("Select");
    jButtonQuantOpen.addActionListener(parent);
    jButtonQuantOpen.setActionCommand(isBatch ? "showQuantDirChooser" : "showQuantFileChooser");
    jButtonQuantOpen.setToolTipText(isBatch ? TooltipTexts.QUANTITATION_BATCH_MASS_LIST : TooltipTexts.QUANTITATION_SINGLE_MASS_LIST);
    selectionPanel.add(jButtonQuantOpen,new GridBagConstraints(7, y, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    selectedQuantDir_ = new JTextField(54);
    selectedQuantDir_.setMinimumSize(selectedQuantDir_.getPreferredSize());
    selectedQuantDir_.setToolTipText(isBatch ? TooltipTexts.QUANTITATION_BATCH_MASS_LIST : TooltipTexts.QUANTITATION_SINGLE_MASS_LIST);
    selectionPanel.add(selectedQuantDir_,new GridBagConstraints(1, y, 6, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    JPanel settingsPanel1 = addJPanelToQuantMenu(++y, this);
    
    JLabel timeMinusTolLabel = new JLabel("Retention time tolerance before reference: ");
    timeMinusTolLabel.setToolTipText(TooltipTexts.QUANTITATION_RET_BEFORE);
    settingsPanel1.add(timeMinusTolLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    timeMinusTol_ = new JTextField(4);
    timeMinusTol_.setMinimumSize(timeMinusTol_.getPreferredSize());
    timeMinusTol_.setText("5");
    timeMinusTol_.setHorizontalAlignment(JTextField.RIGHT);
    timeMinusTol_.setToolTipText(TooltipTexts.QUANTITATION_RET_BEFORE);
    settingsPanel1.add(timeMinusTol_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timeMinusUnit = new JLabel("min, ");
    timeMinusUnit.setToolTipText(TooltipTexts.QUANTITATION_RET_BEFORE);
    settingsPanel1.add(timeMinusUnit,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timePlusTolLabel = new JLabel("after reference: ");
    timePlusTolLabel.setToolTipText(TooltipTexts.QUANTITATION_RET_AFTER);
    settingsPanel1.add(timePlusTolLabel,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    timePlusTol_ = new JTextField(4);
    timePlusTol_.setMinimumSize(timePlusTol_.getPreferredSize());
    timePlusTol_.setText("5");
    timePlusTol_.setHorizontalAlignment(JTextField.RIGHT);
    timePlusTol_.setToolTipText(TooltipTexts.QUANTITATION_RET_AFTER);
    settingsPanel1.add(timePlusTol_,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timePlusUnit = new JLabel("min, ");
    timePlusUnit.setToolTipText(TooltipTexts.QUANTITATION_RET_AFTER);
    settingsPanel1.add(timePlusUnit,new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel rtShiftLabel = new JLabel("RT-shift: ");
    rtShiftLabel.setToolTipText(TooltipTexts.QUANTIFICATION_RET_SHIFT);
    settingsPanel1.add(rtShiftLabel,new GridBagConstraints(6, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    rtShift_ = new JTextField(4);
    rtShift_.setMinimumSize(rtShift_.getPreferredSize());
    rtShift_.setText("0.0");
    rtShift_.setHorizontalAlignment(JTextField.RIGHT);
    rtShift_.setToolTipText(TooltipTexts.QUANTIFICATION_RET_SHIFT);
    settingsPanel1.add(rtShift_,new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel rtShiftUnit = new JLabel("min");
    rtShiftUnit.setToolTipText(TooltipTexts.QUANTIFICATION_RET_SHIFT);
    settingsPanel1.add(rtShiftUnit,new GridBagConstraints(8, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    JPanel settingsPanel2 = addJPanelToQuantMenu(++y, this);
    
    JLabel cutoffLabel = new JLabel("Relative base-peak cutoff: ");
    cutoffLabel.setToolTipText(TooltipTexts.QUANTITATION_CUTOFF);
    settingsPanel2.add(cutoffLabel,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    cutoff_ = new JTextField(4);
    cutoff_.setMinimumSize(cutoff_.getPreferredSize());
    cutoff_.setText(LipidomicsConstants.getBasePeakDefaultCutoff());
    cutoff_.setHorizontalAlignment(JTextField.RIGHT);
    cutoff_.setToolTipText(TooltipTexts.QUANTITATION_CUTOFF);
    settingsPanel2.add(cutoff_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel cutoffUnit = new JLabel("\u2030");
    cutoffUnit.setToolTipText(TooltipTexts.QUANTITATION_CUTOFF);
    settingsPanel2.add(cutoffUnit,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    JPanel settingsPanel3 = addJPanelToQuantMenu(++y, this);
    
    isoValidation_ = new JCheckBox("Isotopic quantitation of ");
    isoValidation_.setSelected(true);
    isoValidation_.addItemListener(parent.new LipidomicsItemListener("endisableAmountOfBatchIsotopes"));
    isoValidation_.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel3.add(isoValidation_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    amountOfIsotopes_ = new JTextField(4);
    amountOfIsotopes_.setMinimumSize(amountOfIsotopes_.getPreferredSize());
    amountOfIsotopes_.setText("2");
    amountOfIsotopes_.setHorizontalAlignment(JTextField.RIGHT);
    amountOfIsotopes_.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel3.add(amountOfIsotopes_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel isotopePeaks = new JLabel("isotopes where ");
    isotopePeaks.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel3.add(isotopePeaks,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    amountOfMatchingSearchIsotopes_ = new JTextField(4);
    amountOfMatchingSearchIsotopes_.setMinimumSize(amountOfMatchingSearchIsotopes_.getPreferredSize());
    amountOfMatchingSearchIsotopes_.setText("1");
    amountOfMatchingSearchIsotopes_.setHorizontalAlignment(JTextField.RIGHT);
    amountOfMatchingSearchIsotopes_.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel3.add(amountOfMatchingSearchIsotopes_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel isotopeMatchPeaks = new JLabel("isotopic peak(s) must match");
    isotopeMatchPeaks.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel3.add(isotopeMatchPeaks,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    JPanel settingsPanel4 = addJPanelToQuantMenu(++y, this);
    
    searchUnknownTime_ = new JCheckBox("Find molecules where retention time is unknown");
    searchUnknownTime_.setSelected(true);
    searchUnknownTime_.setToolTipText(TooltipTexts.QUANTITATION_RET_UNKNOWN);
    settingsPanel4.add(searchUnknownTime_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    JPanel settingsPanel5 = addJPanelToQuantMenu(++y, this);
    
    JLabel maxProcessorsChrom = new JLabel("Processors for file translation:");
    maxProcessorsChrom.setToolTipText(TooltipTexts.CHROM_PROCESSORS);
    settingsPanel5.add(maxProcessorsChrom,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    nrProcessorsChrom_ = new JTextField(3);
    nrProcessorsChrom_.setMinimumSize(nrProcessorsChrom_.getPreferredSize());
    nrProcessorsChrom_.setText(String.valueOf(1));
    nrProcessorsChrom_.setHorizontalAlignment(JTextField.RIGHT);
    nrProcessorsChrom_.setToolTipText(TooltipTexts.CHROM_PROCESSORS);
    nrProcessorsChrom_.setInputVerifier(new IntegerMaxVerifier(true,1,getMaxProcessors()));
    settingsPanel5.add(nrProcessorsChrom_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    JLabel maxProcessors = new JLabel("Processors for quantitation:");
    maxProcessors.setToolTipText(TooltipTexts.QUANTITATION_PROCESSORS);
    settingsPanel5.add(maxProcessors,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    nrProcessors_ = new JTextField(3);
    nrProcessors_.setMinimumSize(nrProcessors_.getPreferredSize());
    nrProcessors_.setText(String.valueOf(this.getAmountOfProcessorsPreferred()));
    nrProcessors_.setHorizontalAlignment(JTextField.RIGHT);
    nrProcessors_.setToolTipText(TooltipTexts.QUANTITATION_PROCESSORS);
    nrProcessors_.setInputVerifier(new IntegerMaxVerifier(true,1,getMaxProcessors()));
    settingsPanel5.add(nrProcessors_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    if (Settings.useAlex()){
      JPanel settingsPanel6 = addJPanelToQuantMenu(++y, this);
      JLabel ionModeLabel = new JLabel("Ion mode:");
      ionModeLabel.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      settingsPanel6.add(ionModeLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      ionMode_ = new JComboBox<String>();
      ionMode_.addItem("+");
      ionMode_.addItem("-");
      ionMode_.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      ionMode_.setFont(LipidDataAnalyzer.SELECT_FIELD_FONT);
      settingsPanel6.add(ionMode_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      JLabel ionModeInfo = new JLabel("(selection required for Alex123 target lists)");
      ionModeInfo.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      settingsPanel6.add(ionModeInfo,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    }
    
    y++;
    startQuantification_ = new JButton("Start Quantitation");
    startQuantification_.addActionListener(parent);
    startQuantification_.setActionCommand(isBatch ? "startBatchQuantification" : "startQuantification");
    startQuantification_.setToolTipText(TooltipTexts.QUANTITATION_BATCH_START);
    this.add(startQuantification_,new GridBagConstraints(0, y, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    
    quantifyingLabel_ = new JLabel("Quantifying");
    quantifyingLabel_.setToolTipText(TooltipTexts.QUANTITATION_STATUS_TEXT);
    progressBar_ = new JProgressBar();
    progressBar_.setMaximum(100);
    progressBar_.setToolTipText(TooltipTexts.QUANTITATION_PROGRESS);
    Icon icon = new ImageIcon(LipidDataAnalyzer.class.getResource("/images/spinner.gif"));
    spinnerLabel_ = new JLabel(icon);
    
    if (isBatch)
    {
    	quantifyingPanel_ = new JPanel();  
      quantifyingPanel_.setLayout(new BorderLayout());
      this.add(quantifyingPanel_,new GridBagConstraints(0, ++y, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    	
    	batchQuantTableModel_ = new BatchQuantificationTableModel();
      batchQuantTable_ = new BatchQuantificationTable(batchQuantTableModel_);

      TableColumnModel tcm = batchQuantTable_.getColumnModel();
      tcm.getColumn(0).setPreferredWidth(15);
      tcm.getColumn(0).setMaxWidth(15);
      tcm.getColumn(1).setPreferredWidth(170);
      tcm.getColumn(2).setPreferredWidth(200);
      tcm.getColumn(3).setPreferredWidth(170);

      JScrollPane scrollPane = new JScrollPane(batchQuantTable_);
      scrollPane.setPreferredSize(new Dimension(650, 200));
      quantifyingPanel_.add(scrollPane,BorderLayout.CENTER);
      JPanel quantProgressPanel = new JPanel();
      quantProgressPanel.setLayout(new GridBagLayout());
      quantifyingPanel_.add(quantProgressPanel,BorderLayout.SOUTH);
      
      quantProgressPanel.add(quantifyingLabel_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));;
      quantProgressPanel.add(progressBar_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      quantProgressPanel.add(spinnerLabel_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    }
    else
    {
    	quantifyingPanel_ = addJPanelToQuantMenu(++y, this);
    	
      quantifyingPanel_.add(quantifyingLabel_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      quantifyingPanel_.add(progressBar_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      quantifyingPanel_.add(spinnerLabel_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    }
    
    
    quantifyingPanel_.setVisible(false);
	}
	
  private JPanel addJPanelToQuantMenu(int y, JPanel quantMenu)
  {
  	JPanel toAdd = new JPanel(); 
  	toAdd.setLayout(new GridBagLayout());
    quantMenu.add(toAdd,new GridBagConstraints(0, y, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    return toAdd;
  }
  
  private int getMaxProcessors(){
    int maxProcessors = Runtime.getRuntime().availableProcessors();
    if (maxProcessors<1) maxProcessors = Integer.MAX_VALUE;
    return maxProcessors;
  }
  
  private int getAmountOfProcessorsPreferred(){
    int maxProcessors = getMaxProcessors();
    if (maxProcessors == Integer.MAX_VALUE) maxProcessors=1;
    int procs = maxProcessors-1;
    if (procs<1) procs=1;
    return procs;
  }

	public JTextField getSelectedMzxmlDirectory()
	{
		return selectedMzxmlDirectory_;
	}

	public void setSelectedMzxmlDirectory(JTextField selectedMzxmlDirectory)
	{
		this.selectedMzxmlDirectory_ = selectedMzxmlDirectory;
	}

	public JTextField getSelectedQuantDir()
	{
		return selectedQuantDir_;
	}

	public void setSelectedQuantDir(JTextField selectedQuantDir)
	{
		this.selectedQuantDir_ = selectedQuantDir;
	}

	public JTextField getTimeMinusTol()
	{
		return timeMinusTol_;
	}

	public void setTimeMinusTol(JTextField timeMinusTol)
	{
		this.timeMinusTol_ = timeMinusTol;
	}

	public JTextField getTimePlusTol()
	{
		return timePlusTol_;
	}

	public void setTimePlusTol(JTextField timePlusTol)
	{
		this.timePlusTol_ = timePlusTol;
	}

	public JTextField getCutoff()
	{
		return cutoff_;
	}

	public void setCutoff(JTextField cutoff)
	{
		this.cutoff_ = cutoff;
	}

	public JTextField getRtShift()
	{
		return rtShift_;
	}

	public void setRtShift(JTextField rtShift)
	{
		this.rtShift_ = rtShift;
	}

	public JCheckBox getIsoValidation()
	{
		return isoValidation_;
	}

	public void setIsoValidation(JCheckBox isoValidation)
	{
		this.isoValidation_ = isoValidation;
	}

	public JTextField getAmountOfIsotopes()
	{
		return amountOfIsotopes_;
	}

	public void setAmountOfIsotopes(JTextField amountOfIsotopes)
	{
		this.amountOfIsotopes_ = amountOfIsotopes;
	}

	public JTextField getAmountOfMatchingSearchIsotopes()
	{
		return amountOfMatchingSearchIsotopes_;
	}

	public void setAmountOfMatchingSearchIsotopes(
			JTextField amountOfMatchingSearchIsotopes)
	{
		this.amountOfMatchingSearchIsotopes_ = amountOfMatchingSearchIsotopes;
	}

	public JCheckBox getSearchUnknownTime()
	{
		return searchUnknownTime_;
	}

	public void setSearchUnknownTime(JCheckBox searchUnknownTime)
	{
		this.searchUnknownTime_ = searchUnknownTime;
	}

	public Integer getNrProcessorsChrom()
	{
		return Integer.parseInt(nrProcessorsChrom_.getText());
	}

	public void setNrProcessorsChrom(Integer nrProcessorsChrom)
	{
		this.nrProcessorsChrom_.setText(String.valueOf(nrProcessorsChrom));
	}

	public Integer getNrProcessors()
	{
		return Integer.parseInt(nrProcessors_.getText());
	}

	public void setNrProcessors(Integer nrProcessors)
	{
		this.nrProcessors_.setText(String.valueOf(nrProcessors));
	}

	public JComboBox<String> getIonMode()
	{
		return ionMode_;
	}

	public void setIonMode(JComboBox<String> ionMode)
	{
		this.ionMode_ = ionMode;
	}

	public JButton getStartQuantification()
	{
		return startQuantification_;
	}

	public void setStartQuantification(JButton startQuantification)
	{
		this.startQuantification_ = startQuantification;
	}

	public JPanel getQuantifyingPanel()
	{
		return quantifyingPanel_;
	}

	public void setQuantifyingPanel(JPanel quantifyingPanel)
	{
		this.quantifyingPanel_ = quantifyingPanel;
	}

	public JLabel getQuantifyingLabel()
	{
		return quantifyingLabel_;
	}

	public void setQuantifyingLabel(JLabel quantifyingLabel)
	{
		this.quantifyingLabel_ = quantifyingLabel;
	}

	public JProgressBar getProgressBar()
	{
		return progressBar_;
	}

	public void setProgressBar(JProgressBar progressBar)
	{
		this.progressBar_ = progressBar;
	}

	public JLabel getSpinnerLabel()
	{
		return spinnerLabel_;
	}

	public void setSpinnerLabel(JLabel spinnerLabel)
	{
		this.spinnerLabel_ = spinnerLabel;
	}

	public BatchQuantificationTableModel getBatchQuantTableModel()
	{
		return batchQuantTableModel_;
	}

	public void setBatchQuantTableModel(
			BatchQuantificationTableModel batchQuantTableModel)
	{
		this.batchQuantTableModel_ = batchQuantTableModel;
	}

	public BatchQuantificationTable getBatchQuantTable()
	{
		return batchQuantTable_;
	}

	public void setBatchQuantTable(BatchQuantificationTable batchQuantTable)
	{
		this.batchQuantTable_ = batchQuantTable;
	}
  
	/**
	 * Sets the batch quantification table (and table model) for the single quant menu from the batch
	 * quant menu. This is required for switching from the single quantification menu to the batch one
	 * (e.g. for SciEX files).
	 * @param theBatchQuantMenu
	 */
	public void setBatchQuantTableForSingleQuant(QuantificationMenu theBatchQuantMenu)
	{
		this.batchQuantTableModel_ = theBatchQuantMenu.getBatchQuantTableModel();
		this.batchQuantTable_ = theBatchQuantMenu.getBatchQuantTable();
	}
  
}
