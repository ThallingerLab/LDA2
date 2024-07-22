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

package at.tugraz.genome.lda.masslist;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.border.TitledBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.AbstractTableModel;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.msn.parser.FALibParser;
import at.tugraz.genome.lda.msn.parser.SPBLibParser;
import at.tugraz.genome.lda.target.JOptionPanel;
import at.tugraz.genome.lda.target.JTargetFileWizard;
import at.tugraz.genome.lda.target.LoadingPanel;
import at.tugraz.genome.lda.vos.AdductVO;
import javafx.util.Pair;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class MassListCreatorPanel extends JPanel 
{
	private static final long serialVersionUID = 1L;
	
	public final static String CHAIN_LIST_FOLDER = "./fattyAcids/";
	public final static String CHAIN_LIST_SUFFIX = ".xlsx";
	private final static String COMMAND_ADDUCT_NAME = "adductName";
	private final static String COMMAND_ADDUCT_FORMULA = "adductFormula";
	private final static String COMMAND_ADDUCT_CHARGE= "adductCharge";
	private final static String COMMAND_ADDUCT_EXPORT= "adductExport";
	private final static String COMMAND_ADDUCT_DELETE= "adductDelete";
	
	private final static String COMMAND_CLASS_NAME = "className";
	private final static String COMMAND_CLASS_FA_CHAIN_NUM= "classFAChainNum";
	private final static String COMMAND_CLASS_LCB_CHAIN_NUM= "classLCBChainNum";
	private final static String COMMAND_CLASS_FORMULA= "classFormula";
	private final static String COMMAND_CLASS_FA_CHAIN_LIST= "classFAChainList";
	private final static String COMMAND_CLASS_LCB_CHAIN_LIST= "classLCBChainList";
	private final static String COMMAND_CLASS_ADDUCT_LIST= "classAdductList";
	private final static String COMMAND_CLASS_CHAIN_C_MIN= "classChainMin";
	private final static String COMMAND_CLASS_DB_MIN= "classDBMin";
	private final static String COMMAND_CLASS_RT_MIN= "classRtMin";
	private final static String COMMAND_CLASS_OH= "classOH";
	private final static String COMMAND_CLASS_OH_MIN= "classOhMin";
	private final static String COMMAND_CLASS_ADDUCT_INSENSITIVE_RT_FILTER= "classRTFilter";
	private final static String COMMAND_CLASS_PICK_BEST= "classPickBest";
	private final static String COMMAND_CLASS_EXPORT= "classExport";
	private final static String COMMAND_CLASS_DELETE= "classDelete";
	
	private final static String OUT_OPEN= "outOpen";
	
	public final static String EXPORT_OPTION_NEG = "negative ion mode";
	public final static String EXPORT_OPTION_POS = "positive ion mode";
	public final static String EXPORT_OPTION_BOTH = "both ion modes";
	private final static String[] EXPORT_OPTIONS_ION_MODE = new String[] {EXPORT_OPTION_NEG, EXPORT_OPTION_POS, EXPORT_OPTION_BOTH};
	public final static String EXPORT_FORMAT_LDA = "LDA Mass List";
	public final static String EXPORT_FORMAT_COMPOUND_DISCOVERER = "Compound Discoverer Format";
	private final static String[] EXPORT_OPTIONS_FORMAT = new String[] {EXPORT_FORMAT_LDA, EXPORT_FORMAT_COMPOUND_DISCOVERER};
	private final static String COMMAND_EXPORT_ION_MODE = "exportIonMode";
	private final static String COMMAND_EXPORT_FORMAT = "exportFormat";
	private final static String EXPORT = "export";
	
	
	private final static String PLACEHOLDER_OUT_PATH = "Enter path and file name for the mass list export.";
	
	private final static int PREFERRED_DISPLAY_COMPONENT_WIDTH = 150;
	private final static int PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH = 50;
	
	
	private ArrayList<AdductVO> allDefinedAdducts_;
	private ArrayList<LipidClassVO> allDefinedLipidClasses_;
	private String[] faChainListNames_;
	private String[] lcbChainListNames_;
	private String[] adductListNames_;
	private AdductVO selectedAdduct_;
	private AdductVO tempAdduct_;
	private LipidClassVO selectedClass_;
	private LipidClassVO tempClass_;
	private JPanel displayPanel_;
	private AdductsTable adductsTable_;
	private JPanel lipidClassPanel_;
	private LipidClassTable lipidClassTable_;
	private JPanel adductPanel_;
	private JTextField outTextField_;
	private Path previousSelection_ = null;
	private JComboBox<String> exportedIonMode_;
	private JComboBox<String> exportedFormat_;
	private JTextField numberLCBChainField_;
	private JComboBox<String> lcbChainList_;
	private JTextField ohField_;
	
	public MassListCreatorPanel()
	{
		displayPanel_ = new JPanel();
		displayPanel_.setLayout(new GridBagLayout());
		try
		{
			allDefinedAdducts_ = (new AdductParser()).parse();
			adductListNames_ = getAdductNames();
			allDefinedLipidClasses_ = (new LipidClassParser(allDefinedAdducts_)).parse();
			faChainListNames_ = getChainListNames(true);
			lcbChainListNames_ = getChainListNames(false);
			selectedClass_ = allDefinedLipidClasses_.get(0);
			tempClass_ = new LipidClassVO(selectedClass_);
			selectedAdduct_ = allDefinedAdducts_.get(0);
			tempAdduct_ = new AdductVO(selectedAdduct_);		
			lipidClassTable_ = new LipidClassTable("Defined lipid (sub)classes", getLipidClassNames(), selectedClass_.getLipidClass());
			displayPanel_.add(lipidClassTable_, getDefaultGridBagConstraints(0,0, GridBagConstraints.EAST, 1, 1));
			adductsTable_ = new AdductsTable("Defined adducts", getAdductNames(), selectedAdduct_.getAdductName());
			displayPanel_.add(adductsTable_, getDefaultGridBagConstraints(0,1, GridBagConstraints.EAST, 1, 1));
			lipidClassPanel_ = getLipidClassPanel(selectedClass_);
			displayPanel_.add(lipidClassPanel_, getDefaultGridBagConstraints(1,0, GridBagConstraints.EAST, 1, 1));
			adductPanel_ = getAdductPanel(selectedAdduct_);
			displayPanel_.add(adductPanel_, getDefaultGridBagConstraints(1,1, GridBagConstraints.EAST, 1, 1));
			outTextField_ = instantiatePlaceholderJTextField(PLACEHOLDER_OUT_PATH, PLACEHOLDER_OUT_PATH, 825);
			displayPanel_.add(getOutPathPanel(outTextField_), getDefaultGridBagConstraints(0,2, GridBagConstraints.CENTER, 2, 1));
			displayPanel_.add(instantiateJButton(EXPORT, "Export", true, TooltipTexts.MASSLIST_GENERAL_EXPORT), 
					getDefaultGridBagConstraints(0,3, GridBagConstraints.CENTER, 2, 1, new Insets(10, 10, 10, 10)));
			
			
	  	this.add(displayPanel_);
		}
		catch (Exception ex)
		{
			this.add(new JLabel("Mass list creation interface unavailable due to invalid and/or missing content in the subfolder ./massListCreation or ./fattyAcids"));
			this.add(new JLabel("Error message: "+ex.getMessage()));
		}
	}
	
	private LipidClassVO getVOfromName(String lipidClassName)
	{
		for (LipidClassVO vo : allDefinedLipidClasses_)
		{ 
			if (vo.getLipidClass().equals(lipidClassName))
			{
				return vo;
			}
		}
		return null;
	}
	
	private void refreshLipidClassScrollPane(String selectedLipidClassName) throws FileNotFoundException, IOException, ChemicalFormulaException
	{
		allDefinedLipidClasses_ = (new LipidClassParser(allDefinedAdducts_)).parse();
		faChainListNames_ = getChainListNames(true);
		lcbChainListNames_ = getChainListNames(false);
		LipidClassVO vo = allDefinedLipidClasses_.get(0);
		if (selectedLipidClassName != null)
		{
			vo = getVOfromName(selectedLipidClassName);
		}
		if (vo != null) 
		{
			selectedClass_ = vo;
			tempClass_ = new LipidClassVO(selectedClass_);
		}
		else
		{
			throw new IOException("Something went wrong with the import.");
		}
		displayPanel_.remove(lipidClassTable_);
		lipidClassTable_ = new LipidClassTable("Defined lipid (sub)classes", getLipidClassNames(), selectedClass_.getLipidClass());
		displayPanel_.add(lipidClassTable_, getDefaultGridBagConstraints(0,0, GridBagConstraints.EAST, 1, 1));
		displayPanel_.invalidate();
		displayPanel_.updateUI();
	}
	
	private void refreshAdductScrollPane(String selectedAddductName) throws FileNotFoundException, IOException, ChemicalFormulaException
	{
		allDefinedAdducts_ = (new AdductParser()).parse();
		adductListNames_ = getAdductNames();
		if (selectedAddductName == null)
		{
			selectedAdduct_ = allDefinedAdducts_.get(0);
			tempAdduct_ = new AdductVO(selectedAdduct_);
		}
		else
		{
			for (AdductVO vo : allDefinedAdducts_)
			{ 
				if (vo.getAdductName().equals(selectedAddductName))
				{
					selectedAdduct_ = vo;
					tempAdduct_ = new AdductVO(selectedAdduct_);
				}
			}
		}
		displayPanel_.remove(adductsTable_);
		adductsTable_ = new AdductsTable("Defined adducts", getAdductNames(), selectedAdduct_.getAdductName());
		displayPanel_.add(adductsTable_, getDefaultGridBagConstraints(0,1, GridBagConstraints.EAST, 1, 1));
		displayPanel_.invalidate();
		displayPanel_.updateUI();
		refreshLipidClassScrollPane(selectedClass_.getLipidClass());
		refreshLipidClassPanel();
	}
	
	private void refreshAdductPanel() throws IOException
	{
		displayPanel_.remove(adductPanel_);
		adductPanel_ = getAdductPanel(selectedAdduct_);
		displayPanel_.add(adductPanel_, getDefaultGridBagConstraints(1,1, GridBagConstraints.EAST, 1, 1));
		displayPanel_.invalidate();
		displayPanel_.updateUI();
	}
	
	private void refreshLipidClassPanel() throws IOException
	{
		displayPanel_.remove(lipidClassPanel_);
		lipidClassPanel_ = getLipidClassPanel(selectedClass_);
		displayPanel_.add(lipidClassPanel_, getDefaultGridBagConstraints(1,0, GridBagConstraints.EAST, 1, 1));
		displayPanel_.invalidate();
		displayPanel_.updateUI();
	}
	
	private JPanel getOutPathPanel(JTextField outPathField)
	{
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		JLabel exportIonModeLabel = new JLabel("Export relevant adducts for: ");
		exportIonModeLabel.setToolTipText(TooltipTexts.MASSLIST_GENERAL_ION_MODE);
		panel.add(exportIonModeLabel, getDefaultGridBagConstraints(0,0, GridBagConstraints.WEST, 1, 1));
		exportedIonMode_ = instantiateJComboBox(COMMAND_EXPORT_ION_MODE, EXPORT_OPTIONS_ION_MODE, 0);
		exportedIonMode_.setToolTipText(TooltipTexts.MASSLIST_GENERAL_ION_MODE);
		panel.add(exportedIonMode_, getDefaultGridBagConstraints(1,0, GridBagConstraints.WEST, 1, 1));
		
		JLabel exportFormatLabel = new JLabel("Export file format: ");
		exportFormatLabel.setToolTipText(TooltipTexts.MASSLIST_GENERAL_ION_MODE);
		panel.add(exportFormatLabel, getDefaultGridBagConstraints(0,1, GridBagConstraints.WEST, 1, 1));
		exportedFormat_ = instantiateJComboBox(COMMAND_EXPORT_FORMAT, EXPORT_OPTIONS_FORMAT, 0);
		exportedFormat_.setToolTipText(TooltipTexts.MASSLIST_GENERAL_FORMAT);
		panel.add(exportedFormat_, getDefaultGridBagConstraints(1,1, GridBagConstraints.WEST, 1, 1));
		
		outPathField.setToolTipText(TooltipTexts.MASSLIST_GENERAL_BROWSE_FIELD);
		panel.add(outPathField, getDefaultGridBagConstraints(0,2, GridBagConstraints.WEST, 2, 1));
  	JButton outPathButton = instantiateJButton(OUT_OPEN, "Browse", true, TooltipTexts.MASSLIST_GENERAL_BROWSE);
		panel.add(outPathButton, getDefaultGridBagConstraints(2,2, GridBagConstraints.EAST, 1, 1, new Insets(10,10,10,10)));
		TitledBorder border = JOptionPanel.getTitledPanelBorder("Mass list export");
		panel.setBorder(border);
		panel.setPreferredSize(new Dimension(975,150));
		return panel;
	}
	
	private JPanel getAdductPanel(AdductVO vo)
	{
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		JTextField nameField = instantiateJTextField(COMMAND_ADDUCT_NAME, vo.getAdductName());
		addLabeledTextField(panel, 0, new JLabel("Adduct name: "), nameField, TooltipTexts.MASSLIST_ADDUCT_NAME);
		JTextField formulaField = instantiateJTextField(COMMAND_ADDUCT_FORMULA, vo.getFormulaString());
		addLabeledTextField(panel, 1, new JLabel("Chemical fomula: "), formulaField, TooltipTexts.MASSLIST_ADDUCT_FORMULA);
		JTextField chargeField = instantiateJTextField(COMMAND_ADDUCT_CHARGE, String.valueOf(vo.getCharge()));
		addLabeledTextField(panel, 2, new JLabel("Charge: "), chargeField, TooltipTexts.MASSLIST_ADDUCT_CHARGE);
		JPanel buttonPanel = instantiateJButtonPanel(COMMAND_ADDUCT_DELETE, COMMAND_ADDUCT_EXPORT, "Delete", "Override", "Save New",
				TooltipTexts.MASSLIST_ADDUCT_DELETE, TooltipTexts.MASSLIST_ADDUCT_OVERRIDE, TooltipTexts.MASSLIST_ADDUCT_SAVE_NEW);
		panel.add(buttonPanel, getDefaultGridBagConstraints(0,13, GridBagConstraints.CENTER, 5, 1));
		TitledBorder border = JOptionPanel.getTitledPanelBorder("Display / edit currently selected adduct definition");
		panel.setBorder(border);
		panel.setPreferredSize(new Dimension(550,175));
		return panel;
	}
	
	private JPanel getLipidClassPanel(LipidClassVO vo) throws IOException
	{
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		addLabeledTextField(panel, 0, new JLabel("Lipid (sub)class name: "),
				instantiateJTextField(COMMAND_CLASS_NAME, vo.getLipidClass()), TooltipTexts.MASSLIST_CLASS_NAME);
		addLabeledTextField(panel, 1, new JLabel("Chemical formula without chains: "), 
				instantiateJTextField(COMMAND_CLASS_FORMULA, vo.getHeadGroupFormulaString()), TooltipTexts.MASSLIST_CLASS_FORMULA);
		JComboBox<String> faChainList = instantiateJComboBox(COMMAND_CLASS_FA_CHAIN_LIST, faChainListNames_, findSelectedChainListIndex(vo.getFAChainList()));
		
		numberLCBChainField_ = instantiateJTextField(COMMAND_CLASS_LCB_CHAIN_NUM, String.valueOf(vo.getNumberOfLCBChains()));
		addNumberOfChainSelection(panel, 2, PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH,
				instantiateJTextField(COMMAND_CLASS_FA_CHAIN_NUM, String.valueOf(vo.getNumberOfFAChains())), numberLCBChainField_, 
				TooltipTexts.MASSLIST_CLASS_NUM_CHAIN);	
		addLabeledComboBox(panel, 3, new JLabel("Selected FA chain list: "), faChainList, TooltipTexts.MASSLIST_CLASS_FA_LIST);
		lcbChainList_ = instantiateJComboBox(COMMAND_CLASS_LCB_CHAIN_LIST, lcbChainListNames_, findSelectedChainListIndex(vo.getLCBChainList()));
		lcbChainList_.setEnabled(vo.getNumberOfLCBChains() > 0);
		addLabeledComboBox(panel, 4, new JLabel("Selected SPB chain list: "), lcbChainList_,  TooltipTexts.MASSLIST_CLASS_SPB_LIST);
		addLabeledJListScrollPane(panel, 5, new JLabel("Selected adducts (hold CNTR for multiple selection): "), 
				instantiateJListScrollPane(COMMAND_CLASS_ADDUCT_LIST, adductListNames_, vo, findSelectedAdductListIndices(vo)),
				TooltipTexts.MASSLIST_CLASS_ADDUCT_SELECTION);
		
		ohField_ = instantiateJTextField(COMMAND_CLASS_OH, String.valueOf(vo.getOhNumber()));
		ohField_.setEnabled(vo.getNumberOfLCBChains() > 0);
		addLabeledTextField(panel, 6, new JLabel("Sphingolipid chain OH number: "), ohField_, TooltipTexts.MASSLIST_CLASS_OH_NUMBER);
		addLabeledRange(panel, 7, new JLabel("Chain OH / oxidation range: "), 
				instantiateJTextFieldRange(COMMAND_CLASS_OH_MIN, String.valueOf(vo.getOhRangeFrom()), String.valueOf(vo.getOhRangeTo()), 
				PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH), TooltipTexts.MASSLIST_CLASS_OH_RANGE);
		
		addLabeledRange(panel, 8, new JLabel("Total number of chain C atoms: "), 
				instantiateJTextFieldRange(COMMAND_CLASS_CHAIN_C_MIN, String.valueOf(vo.getMinChainC()), String.valueOf(vo.getMaxChainC()), 
				PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH), TooltipTexts.MASSLIST_CLASS_C_RANGE);
		addLabeledRange(panel, 9, new JLabel("Total number of chain double bonds: "), 
				instantiateJTextFieldRange(COMMAND_CLASS_DB_MIN, String.valueOf(vo.getMinChainDB()), String.valueOf(vo.getMaxChainDB()), 
				PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH), TooltipTexts.MASSLIST_CLASS_DB_RANGE);
		addLabeledRange(panel, 10, new JLabel("Retention time (RT) range in minutes: "), 
				instantiateJTextFieldRange(COMMAND_CLASS_RT_MIN, String.valueOf(vo.getRtRangeFrom()), String.valueOf(vo.getRtRangeTo()), 
				PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH), TooltipTexts.MASSLIST_CLASS_RT_RANGE);
//		addLabeledRange(panel, 11, new JLabel("Oxidized lipid Ox range (optional): "), 
//				instantiateJTextFieldRange(COMMAND_CLASS_OX_MIN, String.valueOf(vo.getOxRangeFrom()), String.valueOf(vo.getOxRangeTo()), PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH));
		JCheckBox rtFilter = instantiateCheckBox(COMMAND_CLASS_ADDUCT_INSENSITIVE_RT_FILTER, vo.isAdductInsensitiveRtFilter());
		addLabeledCheckBox(panel, 11, new JLabel("Enable adduct insensitive RT filter: "), rtFilter, TooltipTexts.MASSLIST_CLASS_RT_FILTER);
		JCheckBox pickBest = instantiateCheckBox(COMMAND_CLASS_PICK_BEST, vo.isPickBestMatchBySpectrumCoverage());
		addLabeledCheckBox(panel, 12, new JLabel("Pick best match by spectrum coverage: "), pickBest, TooltipTexts.MASSLIST_CLASS_PICK_BEST);
		JPanel buttonPanel = instantiateJButtonPanel(COMMAND_CLASS_DELETE, COMMAND_CLASS_EXPORT, "Delete", "Override", "Save New",
				TooltipTexts.MASSLIST_CLASS_DELETE, TooltipTexts.MASSLIST_CLASS_OVERRIDE, TooltipTexts.MASSLIST_CLASS_SAVE_NEW);
		panel.add(buttonPanel, getDefaultGridBagConstraints(0,13, GridBagConstraints.CENTER, 5, 1));
		
		TitledBorder border = JOptionPanel.getTitledPanelBorder("Display / edit currently selected lipid class definition");
		panel.setBorder(border);
		panel.setPreferredSize(new Dimension(550,500));
		return panel;
	}
	
	private JPanel instantiateJButtonPanel(String actionCommandDelete, String actionCommandExport, 
			String textDelete, String textOverride, String textSaveNew, 
			String tooltipsDelete, String tooltipsOverride, String tooltipsSaveNew)
	{
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		JButton buttonDelete = instantiateJButton(actionCommandDelete, textDelete, false, tooltipsDelete);
		JButton buttonOverride = instantiateJButton(actionCommandExport, textOverride, true, tooltipsOverride);
		JButton buttonSaveNew = instantiateJButton(actionCommandExport, textSaveNew, false, tooltipsSaveNew);
		panel.add(buttonDelete, getDefaultGridBagConstraints(0,0, GridBagConstraints.WEST, 1, 1, new Insets(10, 10, 10, 10)));
		panel.add(buttonOverride, getDefaultGridBagConstraints(1,0, GridBagConstraints.WEST, 1, 1, new Insets(10, 10, 10, 10)));
		panel.add(buttonSaveNew, getDefaultGridBagConstraints(2,0, GridBagConstraints.EAST, 1, 1, new Insets(10, 10, 10, 10)));
		return panel;
	}
	
	private JButton instantiateJButton(String actionCommand, String text, boolean isOverride, String tooltips)
	{
		JButton button = new JButton(text);
		button.setToolTipText(tooltips);
		button.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	jButtonExecuter(actionCommand, isOverride);
		  }
	  });
		return button;
	}
	
	private JScrollPane instantiateJListScrollPane(String actionCommand, String[] entries, LipidClassVO vo, int[] indices) throws IOException
	{
		JList<String> jList = new JList<String>(entries);
		jList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		jList.addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent arg0)
			{
				jListScrollChangeExecuter(actionCommand, jList);
			}
	  });
		jList.setSelectedIndices(indices);
		JScrollPane scrollPane = new JScrollPane(jList);
		scrollPane.setPreferredSize(new Dimension(PREFERRED_DISPLAY_COMPONENT_WIDTH,50));
		return scrollPane;
	}
	
	private JComboBox<String> instantiateJComboBox(String actionCommand, String[] entries, int index)
	{
		JComboBox<String> jComboBox = new JComboBox<String>(entries);
		jComboBox.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	jComboBoxChangeExecuter(actionCommand, jComboBox);
		  }
	  });
		jComboBox.setSelectedIndex(index);
		jComboBox.setPreferredSize(new Dimension(PREFERRED_DISPLAY_COMPONENT_WIDTH,20));
		return jComboBox;
	}
	
	private int findSelectedChainListIndex(String chainListName) throws IOException
	{
		for (int i=0; i<faChainListNames_.length; i++)
		{
			if (faChainListNames_[i].equalsIgnoreCase(chainListName))
			{
				return i;
			}
		}
		for (int i=0; i<lcbChainListNames_.length; i++)
		{
			if (lcbChainListNames_[i].equalsIgnoreCase(chainListName))
			{
				return i;
			}
		}
		return 0;
	}
	
	private int[] findSelectedAdductListIndices(LipidClassVO vo) throws IOException
	{
		ArrayList<Integer> indices = new ArrayList<Integer>();
		for (int i=0; i<adductListNames_.length; i++)
		{
			for (AdductVO adduct : vo.getAdducts())
			{
				if (adductListNames_[i].equalsIgnoreCase(adduct.getAdductName()))
				{
					indices.add(i);
				}
			}
		}
		if (indices.isEmpty())
		{
			throw new IOException(String.format("The file defining the lipid class '%s' does not contain any defined adducts!", vo.getLipidClass()));
		}
		int[] arr = new int[indices.size()];
		for (int i=0; i<indices.size(); i++) arr[i] = indices.get(i);
		return arr;
	}
	
	private JCheckBox instantiateCheckBox(String actionCommand, boolean selected)
	{
		JCheckBox checkBox = new JCheckBox();
		checkBox.addActionListener(new ActionListener() {
			@Override
		  public void actionPerformed(ActionEvent e) 
		  {
		  	jCheckBoxChangeExecuter(actionCommand, checkBox);
		  }
	  });
		checkBox.setSelected(selected);
		return checkBox;
	}
	
	/**
	 * Creates a new JTextField with the provided placeholder and a focuslistener.
	 * @param actionCommand
	 * @param text
	 * @param width
	 * @return
	 */
	private JTextField instantiatePlaceholderJTextField(String actionCommand, String text, Integer width)
	{
		JTextField field = new JTextField();
		field.setPreferredSize(new Dimension(width,20));
		field.setText(text);
		field.setForeground(Color.GRAY);
		field.addFocusListener(new FocusListener() {
		    @Override
		    public void focusGained(FocusEvent e) {
		        if (field.getText().equals(text)) {
		        	  field.setText("");
		        	  field.setForeground(Color.BLACK);
		        }
		    }
		    @Override
		    public void focusLost(FocusEvent e) {
		        if (field.getText().isEmpty()) {
		        	  field.setForeground(Color.GRAY);
		        	  field.setText(text);
		        }
		    }
		    });
		return field;
	}
	
	private JTextField instantiateJTextField(String actionCommand, String text)
	{
		return instantiateJTextField(actionCommand, text, PREFERRED_DISPLAY_COMPONENT_WIDTH);
	}
	
	private Pair<JTextField,JTextField> instantiateJTextFieldRange(String actionCommand, String textFrom, String textTo, Integer width)
	{
		JTextField textFieldFrom = new JTextField(textFrom);
		textFieldFrom.setPreferredSize(new Dimension(width,20));
		setDefaultTextFieldBorder(textFieldFrom);
		JTextField textFieldTo = new JTextField(textTo);
		textFieldTo.setPreferredSize(new Dimension(width,20));
		setDefaultTextFieldBorder(textFieldTo);
		
		DocumentListener listener = new DocumentListener(){
			@Override
			public void insertUpdate(DocumentEvent e){textFieldChangeExecuterRange(actionCommand, textFieldFrom, textFieldTo);}
			@Override
			public void removeUpdate(DocumentEvent e){textFieldChangeExecuterRange(actionCommand, textFieldFrom, textFieldTo);}
			@Override
			public void changedUpdate(DocumentEvent e){textFieldChangeExecuterRange(actionCommand, textFieldFrom, textFieldTo);}
    };
    
		textFieldFrom.getDocument().addDocumentListener(listener);
		textFieldTo.getDocument().addDocumentListener(listener);
		return new Pair<JTextField,JTextField>(textFieldFrom,textFieldTo);
	}
	
	private JTextField instantiateJTextField(String actionCommand, String text, Integer width)
	{
		JTextField textField = new JTextField(text);
		textField.setPreferredSize(new Dimension(width,20));
		setDefaultTextFieldBorder(textField);
		textField.getDocument().addDocumentListener(new DocumentListener(){
			@Override
			public void insertUpdate(DocumentEvent e){textFieldChangeExecuter(actionCommand, textField);}
			@Override
			public void removeUpdate(DocumentEvent e){textFieldChangeExecuter(actionCommand, textField);}
			@Override
			public void changedUpdate(DocumentEvent e){textFieldChangeExecuter(actionCommand, textField);}
    });
		return textField;
	}
	
	private void addLabeledJListScrollPane(JPanel panel, Integer yPos, JLabel label, JScrollPane jScrollPane, String tooltips)
	{
		label.setToolTipText(tooltips);
		jScrollPane.setToolTipText(tooltips);
		panel.add(label, getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(jScrollPane, getDefaultGridBagConstraints(2,yPos, GridBagConstraints.CENTER, 3, 1));
	}
	
	private void addNumberOfChainSelection(JPanel panel, Integer yPos, Integer width,
			JTextField numberFAChainField, JTextField numberSPBChainField, String tooltips)
	{
		JLabel labelNum = new JLabel("Number of chains: ");
		labelNum.setToolTipText(tooltips);
		JLabel labelFA = new JLabel("FA: ");
		labelFA.setToolTipText(tooltips);
		JLabel labelSPB = new JLabel("LCB: ");
		labelSPB.setToolTipText(tooltips);
		numberFAChainField.setPreferredSize(new Dimension(width,20));
		numberFAChainField.setToolTipText(tooltips);
		numberSPBChainField.setPreferredSize(new Dimension(width,20));
		numberSPBChainField.setToolTipText(tooltips);
		
		panel.add(labelNum, getDefaultGridBagConstraints(0,2, GridBagConstraints.WEST, 1, 1));
		panel.add(labelFA, getDefaultGridBagConstraints(1,2, GridBagConstraints.EAST, 1, 1));
		panel.add(numberFAChainField, getDefaultGridBagConstraints(2,2, GridBagConstraints.WEST, 1, 1));
		panel.add(labelSPB, getDefaultGridBagConstraints(3,2, GridBagConstraints.EAST, 1, 1));
		panel.add(numberSPBChainField, getDefaultGridBagConstraints(4,2, GridBagConstraints.EAST, 1, 1));
	}
	
	private void addLabeledComboBox(JPanel panel, Integer yPos, JLabel label, JComboBox<String> comboBox, String tooltips)
	{
		label.setToolTipText(tooltips);
		comboBox.setToolTipText(tooltips);
		panel.add(label, getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(comboBox, getDefaultGridBagConstraints(2,yPos, GridBagConstraints.CENTER, 3, 1));
	}
	
	private void addLabeledCheckBox(JPanel panel, Integer yPos, JLabel label, JCheckBox checkBox, String tooltips)
	{
		label.setToolTipText(tooltips);
		checkBox.setToolTipText(tooltips);
		panel.add(label, getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(checkBox, getDefaultGridBagConstraints(2,yPos, GridBagConstraints.CENTER, 3, 1));
	}
	
	private void addLabeledTextField(JPanel panel, Integer yPos, JLabel label, JTextField textField, String tooltips)
	{
		label.setToolTipText(tooltips);
		textField.setToolTipText(tooltips);
		panel.add(label, getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(textField, getDefaultGridBagConstraints(2,yPos, GridBagConstraints.EAST, 3, 1));
	}
	
	private void addLabeledRange(JPanel panel, Integer yPos, JLabel label, Pair<JTextField,JTextField> range, String tooltips)
	{
		label.setToolTipText(tooltips);
		range.getKey().setToolTipText(tooltips);
		range.getValue().setToolTipText(tooltips);
		JLabel minField = new JLabel("From: ");
		minField.setToolTipText(tooltips);
		JLabel maxField = new JLabel("To: ");
		maxField.setToolTipText(tooltips);
		panel.add(label, getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 1, 1));
		panel.add(minField, getDefaultGridBagConstraints(1,yPos, GridBagConstraints.EAST, 1, 1));
		panel.add(range.getKey(), getDefaultGridBagConstraints(2,yPos, GridBagConstraints.WEST, 1, 1));
		panel.add(maxField, getDefaultGridBagConstraints(3,yPos, GridBagConstraints.EAST, 1, 1));
		panel.add(range.getValue(), getDefaultGridBagConstraints(4,yPos, GridBagConstraints.EAST, 1, 1));
	}
	
	private GridBagConstraints getDefaultGridBagConstraints(int column, int row, int orientation, int width, int height, Insets insets)
	{
		return new GridBagConstraints(
				column, 
				row, 
				width, 
				height, 
				0.0, 
				0.0,
				orientation, 
				GridBagConstraints.NONE, 
				insets, 
				0, 
				5);
	}
	
	private GridBagConstraints getDefaultGridBagConstraints(int column, int row, int orientation, int width, int height)
	{
		return getDefaultGridBagConstraints(column, row, orientation, width, height, new Insets(2, 3, 2, 3));
	}
	
	private void textFieldChangeExecuterRange(String actionCommand, JTextField textfieldFrom, JTextField textfieldTo)
	{
		switch (actionCommand)
		{
			case COMMAND_CLASS_CHAIN_C_MIN:
				tempClass_.setMinChainC(-1);
				tempClass_.setMaxChainC(-1);
				try
				{
					int numFrom = Integer.parseInt(textfieldFrom.getText());
					int numTo = Integer.parseInt(textfieldTo.getText());
					
					if (numFrom > 0 && numTo > 0 && numFrom < numTo)
					{
						setDefaultTextFieldBorder(textfieldFrom);
						setDefaultTextFieldBorder(textfieldTo);
					}	
					else
					{
						setWarningTextFieldBorder(textfieldFrom);
						setWarningTextFieldBorder(textfieldTo);
					}
					tempClass_.setMinChainC(numFrom);
					tempClass_.setMaxChainC(numTo);
				}
				catch (NumberFormatException ex) 
				{
					setWarningTextFieldBorder(textfieldFrom);
					setWarningTextFieldBorder(textfieldTo);
				}
				break;	
			case COMMAND_CLASS_DB_MIN:
				tempClass_.setMinChainDB(-1);
				tempClass_.setMaxChainDB(-1);
				try
				{
					int numFrom = Integer.parseInt(textfieldFrom.getText());
					int numTo = Integer.parseInt(textfieldTo.getText());
					
					if ((numFrom < numTo) || (numFrom < 1 && numTo < 1))
					{
						setDefaultTextFieldBorder(textfieldFrom);
						setDefaultTextFieldBorder(textfieldTo);
					}	
					else
					{
						setWarningTextFieldBorder(textfieldFrom);
						setWarningTextFieldBorder(textfieldTo);
					}
					tempClass_.setMinChainDB(numFrom);
					tempClass_.setMaxChainDB(numTo);
				}
				catch (NumberFormatException ex) 
				{
					setWarningTextFieldBorder(textfieldFrom);
					setWarningTextFieldBorder(textfieldTo);
				}
				break;	
			case COMMAND_CLASS_RT_MIN:
				tempClass_.setRtRangeFrom(-1);
				tempClass_.setRtRangeTo(-1);
				try
				{
					int numFrom = Integer.parseInt(textfieldFrom.getText());
					int numTo = Integer.parseInt(textfieldTo.getText());
					
					if ((numFrom < numTo) || (numFrom <= 0 && numTo <= 0))
					{
						setDefaultTextFieldBorder(textfieldFrom);
						setDefaultTextFieldBorder(textfieldTo);
					}	
					else
					{
						setWarningTextFieldBorder(textfieldFrom);
						setWarningTextFieldBorder(textfieldTo);
					}
					tempClass_.setRtRangeFrom(numFrom);
					tempClass_.setRtRangeTo(numTo);
				}
				catch (NumberFormatException ex) 
				{
					setWarningTextFieldBorder(textfieldFrom);
					setWarningTextFieldBorder(textfieldTo);
				}
				break;	

			case COMMAND_CLASS_OH_MIN:
				tempClass_.setOhRangeFrom(-1);
				tempClass_.setOhRangeTo(-1);
				try
				{
					int numFrom = Integer.parseInt(textfieldFrom.getText());
					int numTo = Integer.parseInt(textfieldTo.getText());
					
					if ((tempClass_.getNumberOfLCBChains() > 0 && numFrom > 0 && (numFrom <= tempClass_.getOhNumber() && tempClass_.getOhNumber() <= numTo))
							|| (tempClass_.getNumberOfLCBChains() == 0 && numFrom >= 0 && (numFrom <= numTo)))
					{
						setDefaultTextFieldBorder(textfieldFrom);
						setDefaultTextFieldBorder(textfieldTo);
					}	
					else
					{
						setWarningTextFieldBorder(textfieldFrom);
						setWarningTextFieldBorder(textfieldTo);
					}
					tempClass_.setOhRangeFrom(numFrom);
					tempClass_.setOhRangeTo(numTo);
				}
				catch (NumberFormatException ex) 
				{
					setWarningTextFieldBorder(textfieldFrom);
					setWarningTextFieldBorder(textfieldTo);
				}
				break;
				
			default:
				break;
		}
	}
	
	private void textFieldChangeExecuter(String actionCommand, JTextField textfield)
	{
		switch (actionCommand)
		{
			case COMMAND_ADDUCT_NAME:
				tempAdduct_.setAdductName(textfield.getText());
				if (isAdductNameAvailable(textfield.getText()))
					setDefaultTextFieldBorder(textfield);
				else
					setWarningTextFieldBorder(textfield);
				break;
			case COMMAND_ADDUCT_FORMULA:
				try
				{
					tempAdduct_.setFormulaString(textfield.getText());
					setDefaultTextFieldBorder(textfield);
				}
				catch (ChemicalFormulaException ex) {setWarningTextFieldBorder(textfield);}
				break;
			case COMMAND_ADDUCT_CHARGE:
				try
				{
					tempAdduct_.setCharge(0);
					tempAdduct_.setCharge(Integer.parseInt(textfield.getText()));
					setDefaultTextFieldBorder(textfield);
				}
				catch (NumberFormatException ex) {setWarningTextFieldBorder(textfield);}
				break;	
			case COMMAND_CLASS_NAME:
				tempClass_.setLipidClass(textfield.getText());
				if (isLipidClassNameAvailable(textfield.getText()))
					setDefaultTextFieldBorder(textfield);
				else
					setWarningTextFieldBorder(textfield);
				break;
//			case COMMAND_CLASS_CHAIN_NUM:
//				tempClass_.setNumberOfChains(0);
//				try
//				{
//					int chainNum = Integer.parseInt(textfield.getText());
//					if (chainNum > 0)
//						setDefaultTextFieldBorder(textfield);
//					else
//						setWarningTextFieldBorder(textfield);
//					tempClass_.setNumberOfChains(chainNum);
//				}
//				catch (NumberFormatException ex) {setWarningTextFieldBorder(textfield);}
//				break;
			case COMMAND_CLASS_FA_CHAIN_NUM:
				tempClass_.setNumberOfFAChains(0);
				try
				{
					int chainNum = Integer.parseInt(textfield.getText());
					if (chainNum > 0)
						setDefaultTextFieldBorder(textfield);
					else
						setWarningTextFieldBorder(textfield);
					tempClass_.setNumberOfFAChains(chainNum);
				}
				catch (NumberFormatException ex) {setWarningTextFieldBorder(textfield);}
				break;
			case COMMAND_CLASS_LCB_CHAIN_NUM:
				tempClass_.setNumberOfLCBChains(0);
				try
				{
					int chainNum = Integer.parseInt(textfield.getText());
					if (chainNum >= 0)
						setDefaultTextFieldBorder(textfield);
					else
						setWarningTextFieldBorder(textfield);
					tempClass_.setNumberOfLCBChains(chainNum);
					lcbChainList_.setEnabled(chainNum > 0);
					tempClass_.setOhNumber(chainNum > 0 ? (tempClass_.getOhNumber()>0 ? tempClass_.getOhNumber() : 2) : 0);
					ohField_.setText(String.valueOf(tempClass_.getOhNumber()));
					ohField_.setEnabled(chainNum > 0);
				}
				catch (NumberFormatException ex) {setWarningTextFieldBorder(textfield);}
				break;
			case COMMAND_CLASS_FORMULA:
				try
				{
					tempClass_.setHeadGroupFormulaString(textfield.getText());
					setDefaultTextFieldBorder(textfield);
				}
				catch (ChemicalFormulaException ex) {setWarningTextFieldBorder(textfield);}
				break;
			case COMMAND_CLASS_OH:
				tempClass_.setOhNumber(-1);
				try
				{
					int num = Integer.parseInt(textfield.getText());
					if (num >= 0)
						setDefaultTextFieldBorder(textfield);
					else
						setWarningTextFieldBorder(textfield);
					tempClass_.setOhNumber(num);
				}
				catch (NumberFormatException ex) {setWarningTextFieldBorder(textfield);}
				break;
				
			default:
				break;
		}
	}
	
	private void jListScrollChangeExecuter(String actionCommand, JList<String> jList)
	{
		switch (actionCommand)
		{
			case COMMAND_CLASS_ADDUCT_LIST:
				ArrayList<AdductVO> selectedAdducts = new ArrayList<AdductVO>();
				for (String adductName : jList.getSelectedValuesList())
				{
					for (AdductVO vo : allDefinedAdducts_)
					{
						if (adductName.equals(vo.getAdductName()))
						{
							selectedAdducts.add(vo);
						}
					}
				}
				tempClass_.setAdducts(selectedAdducts);
				break;
			default:
				break;
		}
	}
	
	private void jComboBoxChangeExecuter(String actionCommand, JComboBox<String> jComboBox)
	{
		switch (actionCommand)
		{
			case COMMAND_CLASS_FA_CHAIN_LIST:
				tempClass_.setFaChainList((String)jComboBox.getSelectedItem());
				break;
			case COMMAND_CLASS_LCB_CHAIN_LIST:
				tempClass_.setLCBChainList((String)jComboBox.getSelectedItem());
				break;
			default:
				break;
		}
	}
	
	private void jCheckBoxChangeExecuter(String actionCommand, JCheckBox checkBox)
	{
		switch (actionCommand)
		{
			case COMMAND_CLASS_ADDUCT_INSENSITIVE_RT_FILTER:
				tempClass_.setAdductInsensitiveRtFilter(checkBox.isSelected());
				break;
			case COMMAND_CLASS_PICK_BEST:
				tempClass_.setPickBestMatchBySpectrumCoverage(checkBox.isSelected());
				break;
			default:
				break;
		}
	}
	
	private void jButtonExecuter(String actionCommand, boolean isOverride)
	{
		switch (actionCommand)
		{
			case COMMAND_ADDUCT_DELETE:
				deleteAdduct();
				break;
			case COMMAND_ADDUCT_EXPORT:
				if (isAdductViable())
				{
					if (!isOverride ||
							( isOverride && JOptionPane.showConfirmDialog(new JFrame(),
								String.format("Would you like to override the adduct definition '%s'?", selectedAdduct_.getAdductName()),
							  "Override definition", JOptionPane.YES_NO_OPTION) == 0 ) )
					{
						exportAdduct(isOverride);
					}
				}
				else
				{
					new WarningMessage(new JFrame(), "Error", "The adduct definition contains erroneous user-entries, please correct textfields highlighted in red before exporting.");
				}
				break;
			case COMMAND_CLASS_DELETE:
				deleteLipidClass();
				break;
			case COMMAND_CLASS_EXPORT:
				if (isLipidClassViable())
				{
					if (!isOverride ||
							( isOverride && JOptionPane.showConfirmDialog(new JFrame(),
								String.format("Would you like to override the lipid (sub)class definition '%s'?", selectedClass_.getLipidClass()),
							  "Override definition", JOptionPane.YES_NO_OPTION) == 0 ) )
					{
						exportLipidClass(isOverride);
					}
				}	
				else
				{
					new WarningMessage(new JFrame(), "Error", "The lipid class definition contains erroneous user-entries, please correct textfields highlighted in red before exporting and make sure at least one adduct is selected.");
				}
				break;
			case OUT_OPEN:
				selectOutPath();
				break;
			case EXPORT:
				String outPath = outTextField_.getText();
				if (outPath.equals(PLACEHOLDER_OUT_PATH) || outPath.length() < 1)
				{
					new WarningMessage(new JFrame(), "Error", "A filepath to write the mass list to must be defined prior to the export!");
				}
				else if (lipidClassTable_.getSelectedLipidClasses().isEmpty())
				{
					new WarningMessage(new JFrame(), "Error", "Please select at least one lipid class prior to the export!");
				}
				else
				{
					if (!outPath.endsWith(".xlsx"))
					{
						outPath += ".xlsx";
						outTextField_.setText(outPath);
					}
					MassListExporter exporter = new MassListExporter(outPath, lipidClassTable_.getSelectedLipidClasses(), 
							(String)exportedIonMode_.getSelectedItem(), (String)exportedFormat_.getSelectedItem());
					
					StringBuilder builder = new StringBuilder();
					builder.append("<html>Writing the mass list to the specified file.<br>");
					builder.append("Please wait... </html>");
					
					LoadingPanel waitPanel = new LoadingPanel(builder.toString());
					
					Thread thread = new Thread(new Runnable()
		  		{
		  			public void run()
		  			{
		  				updateUI(displayPanel_, waitPanel);
		  				exporter.export();
		  				updateUI(waitPanel, displayPanel_);
		  			}
		  		});
		    	thread.start(); 
		    	
				}
				break;
			default:
				break;
		}
	}
	
	private void updateUI(JPanel current, JPanel update)
	{
		this.remove(current);
		this.add(update);
		this.invalidate();
		this.updateUI();
	}
	
	/**
	 * Sets the text of the provided JTextField to the selected file or folder path depending on the selection mode.
	 */
	private void selectOutPath()
	{
		JFileChooser chooser = new JFileChooser();
		chooser.setFileFilter(new FileNameExtensionFilter("Only .xlsx", "xlsx"));
		chooser.setPreferredSize(JTargetFileWizard.DEFAULT_FILE_CHOOSER_DIMENSION);
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
		chooser.setDialogTitle("Select the path for the mass list export");
		if (previousSelection_ != null)
		{
			chooser.setCurrentDirectory(previousSelection_.getParent().toFile());
		}
		int val = chooser.showOpenDialog(new JFrame());
		if (val == JFileChooser.APPROVE_OPTION)
		{
			String text = chooser.getSelectedFile().getAbsolutePath();
			previousSelection_ = Paths.get(text);
			outTextField_.setText(text);
			outTextField_.setForeground(Color.BLACK);
		}
		else
		{
			return;
		}
	}
	
	private boolean isLipidClassOxDefinitionViable()
	{
		//default
		if (tempClass_.getOhNumber() == 0 && tempClass_.getOhRangeFrom() == 0 && tempClass_.getOhRangeTo() == 0)
		{
			return true;
		}
		//sphingolipids
		else if ((tempClass_.getOhNumber() > 0 && (tempClass_.getOhRangeFrom() > 0 && tempClass_.getOhRangeFrom() < tempClass_.getOhRangeTo())
				&& (tempClass_.getOhRangeFrom() <= tempClass_.getOhNumber() && tempClass_.getOhRangeTo() >= tempClass_.getOhNumber())))
		{
			return true;
		}
		//oxidized lipids
		else if ((tempClass_.getOhNumber() == 0 && (tempClass_.getOhRangeFrom() >= 0 && tempClass_.getOhRangeFrom() < tempClass_.getOhRangeTo()))) 
		{
			return true;
		}
		return false;
	}
	
	private boolean isLipidClassViable()
	{
		if (isLipidClassNameAvailable(tempClass_.getLipidClass())
				&& (tempClass_.getNumberOfFAChains() + tempClass_.getNumberOfLCBChains() > 0)
				&& tempClass_.getHeadgroupFormula() != null
				&& tempClass_.getMinChainC() > 0 && tempClass_.getMinChainC() < tempClass_.getMaxChainC()
				&& (tempClass_.getMinChainDB() >= 0 && tempClass_.getMinChainDB() < tempClass_.getMaxChainDB() || (tempClass_.getMinChainDB() == 0 && tempClass_.getMaxChainDB() == 0))
				&& ((tempClass_.getRtRangeFrom() >= 0 && tempClass_.getRtRangeFrom() < tempClass_.getRtRangeTo()) || (tempClass_.getRtRangeFrom() < 0 && tempClass_.getRtRangeTo() < 0))
				&& isLipidClassOxDefinitionViable()
				&& !tempClass_.getAdducts().isEmpty()
				)
		{
			return true;
		}
		return false;
	}
	
	/**
	 * Deletes an adduct definition file
	 */
	private void deleteAdduct()
	{
		Object[] options = {"Delete","Cancel"};
		int n = JOptionPane.showOptionDialog(new JFrame(), String.format("Are you sure you want to delete the adduct definition '%s'?", 
				selectedAdduct_.getAdductName()), "Deleting adduct definition", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null,
				options,options[0]);
		
		if (n == 0)
		{
			try 
			{
				String originalFilePath = AdductExporter.buildAdductPath(selectedAdduct_.getFileName());
				File file = new File(originalFilePath);
				file.delete();
				updateAdductForLipidClasses(selectedAdduct_, tempAdduct_, true);
				refreshAdductScrollPane(null);
			}
			catch (IOException | ChemicalFormulaException ex) {
				new WarningMessage(new JFrame(), "Error", "An error occurred. Error message: "+ex.getMessage());
			}
		}
	}
	
	/**
	 * Exports an adduct definition file
	 * @param isOverride	true if the selected adduct definition file should be overriden / replaced with the new one
	 */
	private void exportAdduct(boolean isOverride)
	{
		try
		{
			if (isOverride)
			{
				String originalFilePath = AdductExporter.buildAdductPath(selectedAdduct_.getFileName());
				File file = new File(originalFilePath);
				file.delete();
			}
			String fileName = buildAdductFileName(tempAdduct_.getAdductName());
			tempAdduct_.setFileName(fileName);
			AdductExporter exporter = new AdductExporter(tempAdduct_);
			exporter.export();
			if (isOverride)
			{
				updateAdductForLipidClasses(selectedAdduct_, tempAdduct_, false);
			}
			selectedAdduct_ = tempAdduct_;
			refreshAdductScrollPane(selectedAdduct_.getAdductName());
		}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	private void updateAdductForLipidClasses(AdductVO oldAdduct, AdductVO newAdduct, boolean isDelete)
	{
		for (LipidClassVO vo : allDefinedLipidClasses_)
		{
			boolean isChanged = false;
			ArrayList<AdductVO> newAdductVOList = new ArrayList<AdductVO>();
			for (AdductVO adduct : vo.getAdducts())
			{
				if (adduct.getAdductName().equals(oldAdduct.getAdductName()))
				{
					if (!isDelete)
					{
						newAdductVOList.add(newAdduct);
					}
					isChanged = true;
				}
				else
				{
					newAdductVOList.add(adduct);
				}
			}
			if (isChanged)
			{
				vo.setAdducts(newAdductVOList);
				exportLipidClass(vo, true);
			}
		}
	}
	
	/**
	 * Deletes a lipid class definition file
	 */
	private void deleteLipidClass()
	{
		Object[] options = {"Delete","Cancel"};
		int n = JOptionPane.showOptionDialog(new JFrame(), String.format("Are you sure you want to delete the lipid (sub)class definition '%s'?", 
				selectedClass_.getLipidClass()), "Deleting lipid (sub)class definition", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null,
				options,options[0]);
		
		if (n == 0)
		{
			try 
			{
				String originalFilePath = LipidClassExporter.buildLipidClassPath(selectedClass_.getLipidClass());
				File file = new File(originalFilePath);
				file.delete();
				refreshLipidClassScrollPane(null);
			}
			catch (IOException | ChemicalFormulaException ex) {
				new WarningMessage(new JFrame(), "Error", "An error occurred. Error message: "+ex.getMessage());
			}
		}
	}
	
	/**
	 * Exports a lipid class definition file
	 * @param isOverride	true if the selected adduct definition file should be overriden / replaced with the new one
	 */
	private void exportLipidClass(LipidClassVO toExport, boolean isOverride)
	{
		if (isOverride)
		{
			String originalFilePath = LipidClassExporter.buildLipidClassPath(toExport.getLipidClass());
			File file = new File(originalFilePath);
			file.delete();
		}
		LipidClassExporter exporter = new LipidClassExporter(toExport);
		exporter.export();
	}
	
	/**
	 * Exports a lipid class definition file
	 * @param isOverride	true if the selected adduct definition file should be overriden / replaced with the new one
	 */
	private void exportLipidClass(boolean isOverride)
	{
		try
		{
			if (isOverride)
			{
				String originalFilePath = LipidClassExporter.buildLipidClassPath(selectedClass_.getLipidClass());
				File file = new File(originalFilePath);
				file.delete();
			}
			selectedClass_ = tempClass_;
			LipidClassExporter exporter = new LipidClassExporter(selectedClass_);
			exporter.export();
			refreshLipidClassScrollPane(selectedClass_.getLipidClass());
		}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	private String buildAdductFileName(String adductName)
	{
		return "adduct_"+adductName+AdductParser.ADDUCT_SUFFIX;
	}
	
	private boolean isAdductViable()
	{
		if (isAdductNameAvailable(tempAdduct_.getAdductName())
				&& tempAdduct_.getFormula() != null
				&& tempAdduct_.getCharge() != 0)
		{
			return true;
		}
		return false;
	}
	
	private boolean isAdductNameAvailable(String name)
	{
		ArrayList<AdductVO> other = new ArrayList<AdductVO>(allDefinedAdducts_);
		other.remove(selectedAdduct_);
		for (AdductVO vo : other)
		{
			if (vo.getAdductName().equalsIgnoreCase(name)) return false;
		}
		return true;
	}
	
	private boolean isLipidClassNameAvailable(String name)
	{
		ArrayList<LipidClassVO> other = new ArrayList<LipidClassVO>(allDefinedLipidClasses_);
		other.remove(selectedClass_);
		for (LipidClassVO vo : other)
		{
			if (vo.getLipidClass().equalsIgnoreCase(name)) return false;
		}
		return true;
	}
	
	private void setDefaultTextFieldBorder(JTextField textfield)
	{
		textfield.setBorder(BorderFactory.createLineBorder(Color.darkGray));
	}
	
	private void setWarningTextFieldBorder(JTextField textfield)
	{
		textfield.setBorder(BorderFactory.createLineBorder(Color.red));
	}
	
	private String[] getAdductNames()
	{
		ArrayList<String> names = new ArrayList<String>();
		for (AdductVO vo : allDefinedAdducts_)
		{
			names.add(vo.getAdductName());
		}
		return toArray(names);
	}
	
	private String[] getLipidClassNames()
	{
		ArrayList<String> names = new ArrayList<String>();
		for (LipidClassVO vo : allDefinedLipidClasses_)
		{
			names.add(vo.getLipidClass());
		}
		return toArray(names);
	}
	
	private String[] getChainListNames(boolean isFA) throws IOException
	{
		ArrayList<String> names = new ArrayList<String>();
		File folder = new File(CHAIN_LIST_FOLDER);
		if (!folder.exists())
		{
			throw new IOException(String.format("The chain list folder '%s' does not exist!", CHAIN_LIST_FOLDER));
		}
		File[] fileCandidates = folder.listFiles();
		for (int i=0; i<fileCandidates.length;i++)
		{
			String fileName = fileCandidates[i].getName(); 
			if (fileName.endsWith(CHAIN_LIST_SUFFIX))
			{
				try
				{
					if (isFA && new FALibParser(fileCandidates[i]).isFAFile())
					{
				    names.add(fileName.substring(0, fileName.indexOf(CHAIN_LIST_SUFFIX)));
					}
					else if (!isFA && new SPBLibParser(fileCandidates[i]).isLCBFile())
					{
				    names.add(fileName.substring(0, fileName.indexOf(CHAIN_LIST_SUFFIX)));
					}
				}
				catch (Exception ex) {}
			}
		};
		return toArray(names);
	}
	
	private String[] toArray(ArrayList<String> list)
	{
		String[] arr = new String[list.size()];
		return list.toArray(arr);
	}
	
	private void handleLipidClassSelection(String lipidClass) throws IOException, ChemicalFormulaException
	{
		for (LipidClassVO vo : allDefinedLipidClasses_)
		{
			if (vo.getLipidClass().equals(lipidClass))
			{
				selectedClass_ = vo;
				tempClass_ = new LipidClassVO(selectedClass_);
				refreshLipidClassPanel();
				return;
			}
		}
	}
	
	private void handleAdductSelection(String adduct) throws IOException, ChemicalFormulaException
	{
		for (AdductVO vo : allDefinedAdducts_)
		{
			if (vo.getAdductName().equals(adduct))
			{
				selectedAdduct_ = vo;
				tempAdduct_ = new AdductVO(selectedAdduct_);
				refreshAdductPanel();
				return;
			}
		}
	}
	
	
	private class AdductsTable extends JPanel
  {
		private static final long serialVersionUID = 1L;
  	private static final int COLUMN_NAME= 0;
  	
    private JPanel selectionTablePanel_;
    private JTable displayTable_;
    private JScrollPane scrollPane_;
    
    private AdductsTable(String title, String[] adductNames, String selected)
    {
    	this.setPreferredSize(new Dimension(400,175));
    	selectionTablePanel_ = new JPanel();
    	generateSelectionTablePanel(initializeTableData(adductNames), selected);
    	this.setLayout(new GridBagLayout());
    	this.add(selectionTablePanel_, getDefaultGridBagConstraints(0, 1, GridBagConstraints.CENTER, 1, 1));
    	this.setBorder(JOptionPanel.getTitledPanelBorder(title));
    }
    
    private void generateSelectionTablePanel(Object[][] tableData, String selected)
    {
    	String[] columnNames = { "adduct name" };
    	BooleanTableModel model = new BooleanTableModel(tableData, columnNames);
    	displayTable_ = new JTable(model);
    	scrollPane_ = new JScrollPane(displayTable_);
    	scrollPane_.setPreferredSize(new Dimension(325,125));
    	displayTable_.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    	ListSelectionModel selectionModel = displayTable_.getSelectionModel();
    	selectionModel.setAnchorSelectionIndex(model.indexOf(selected,COLUMN_NAME));
    	selectionModel.setLeadSelectionIndex(model.indexOf(selected,COLUMN_NAME));
    	selectionModel.addListSelectionListener(new ListSelectionListener() {
        public void valueChanged(ListSelectionEvent e) {
        	try
        	{
        		MassListCreatorPanel.this.handleAdductSelection((String)model.getValueAt(((ListSelectionModel)e.getSource()).getMaxSelectionIndex(), COLUMN_NAME));
        	}
        	catch (IOException | ChemicalFormulaException ex)
        	{
        		new WarningMessage(new JFrame(), "Warning", String.format("The definition file for the lipid class '%s' could not be parsed.", 
        				(String)model.getValueAt(((ListSelectionModel)e.getSource()).getMaxSelectionIndex(), COLUMN_NAME)));
        	}
        }
    	});
    	selectionTablePanel_.add(scrollPane_);
    }
    
    private Object[][] initializeTableData(String[] adductNames)
    {
    	Object[][] tableData = new Object[adductNames.length][1];
    	for (int i=0; i<adductNames.length; i++)
	    {
	    	tableData[i][COLUMN_NAME] = adductNames[i];
	    }
    	return tableData;
    }
  }
	
	private class LipidClassTable extends JPanel
  {
		private static final long serialVersionUID = 1L;
  	private static final int COLUMN_NAME= 0;
  	private static final int COLUMN_INCLUDE = 1;
  	
  	private String[] lipidClassNames_;
    private JPanel selectionTablePanel_;
    private JTable displayTable_;
    private JScrollPane scrollPane_;
    private BooleanTableModel model_;
    private HashMap<String,Boolean> lipidClassIncluded_;
    
    private LipidClassTable(String title, String[] lipidClassNames, String selected)
    {
    	this.lipidClassNames_ = lipidClassNames;
    	this.setPreferredSize(new Dimension(400,500));
    	selectionTablePanel_ = new JPanel();
    	generateSelectionTablePanel(initializeTableData(lipidClassNames), selected);
    	this.setLayout(new GridBagLayout());
    	this.add(selectionTablePanel_, getDefaultGridBagConstraints(0, 1, GridBagConstraints.CENTER, 1, 1));
    	this.setBorder(JOptionPanel.getTitledPanelBorder(title));
    }
    
    private void generateSelectionTablePanel(Object[][] tableData, String selected)
    {
    	String[] columnNames = { "lipid class name", "include in mass list" };
    	model_ = new BooleanTableModel(tableData, columnNames);
    	displayTable_ = new JTable(model_);
    	scrollPane_ = new JScrollPane(displayTable_);
    	scrollPane_.setPreferredSize(new Dimension(325,450));
    	displayTable_.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    	ListSelectionModel selectionModel = displayTable_.getSelectionModel();
    	selectionModel.setAnchorSelectionIndex(model_.indexOf(selected,COLUMN_NAME));
    	selectionModel.setLeadSelectionIndex(model_.indexOf(selected,COLUMN_NAME));
    	selectionModel.addListSelectionListener(new ListSelectionListener() {
        public void valueChanged(ListSelectionEvent e) {
        	try
        	{
        		MassListCreatorPanel.this.handleLipidClassSelection((String)model_.getValueAt(((ListSelectionModel)e.getSource()).getMaxSelectionIndex(), COLUMN_NAME));
        	}
        	catch (IOException | ChemicalFormulaException ex)
        	{
        		new WarningMessage(new JFrame(), "Warning", String.format("The definition file for the lipid class '%s' could not be parsed.", 
        				(String)model_.getValueAt(((ListSelectionModel)e.getSource()).getMaxSelectionIndex(), COLUMN_NAME)));
        	}
        }
    	});
    	selectionTablePanel_.add(scrollPane_);
    }
    
    private Object[][] initializeTableData(String[] lipidClassNames)
    {
    	Object[][] tableData = new Object[lipidClassNames.length][2];
    	lipidClassIncluded_ = new HashMap<String,Boolean>();
    	for (int i=0; i<lipidClassNames.length; i++)
	    {
	    	tableData[i][COLUMN_NAME] = lipidClassNames[i];
	      tableData[i][COLUMN_INCLUDE] = Boolean.FALSE;
	      lipidClassIncluded_.put(lipidClassNames[i], true);
	    }
    	return tableData;
    }
    
    private ArrayList<LipidClassVO> getSelectedLipidClasses()
    {
    	ArrayList<LipidClassVO> selected = new ArrayList<LipidClassVO>();
    	for (int i=0; i<lipidClassNames_.length; i++)
    	{
    		if ((Boolean)model_.getValueAt(i, COLUMN_INCLUDE).equals(Boolean.TRUE))
    		{
    			LipidClassVO vo = getVOfromName((String)model_.getValueAt(i, COLUMN_NAME));
    			if (vo != null) 
    				selected.add(vo);
    		}
    	}
    	return selected;
    }
  }
	
	private class BooleanTableModel extends AbstractTableModel {
		
		private static final long serialVersionUID = 1L;
		
		Object tableData_[][];
		String[] columnNames_;
		
		private BooleanTableModel(Object[][] tableData, String[] columnNames)
		{
			this.tableData_ = tableData;
			this.columnNames_ = columnNames;
		}

	  public int getColumnCount() {
	    return columnNames_.length;
	  }

	  public String getColumnName(int column) {
	    return columnNames_[column];
	  }

	  public int getRowCount() {
	    return tableData_.length;
	  }

	  public Object getValueAt(int row, int column) {
	    return tableData_[row][column];
	  }
	  
	  @SuppressWarnings({ "rawtypes", "unchecked" })
		public Class getColumnClass(int column) {
	    return (getValueAt(0, column).getClass());
	  }
	  
	  public int indexOf(Object firstColumnEntry, int columnNum)
	  {
	  	for (int i=0;i<tableData_.length;i++)
	  	{
	  		String data = (String)getValueAt(i, columnNum);
	  		if (data.equals(firstColumnEntry))
				{
	  			return i;
				}
	  	}
	  	return 0;
	  }

	  public void setValueAt(Object value, int row, int column) {
	    tableData_[row][column] = value;
	  }

	  public boolean isCellEditable(int row, int column) {
	    return (column != 0);
	  }
	}
}
