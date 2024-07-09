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

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
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
	
	private final static String COMMAND_CLASS_NAME = "className";
	private final static String COMMAND_CLASS_CHAIN_NUM= "classChainNum";
	private final static String COMMAND_CLASS_FORMULA= "classFormula";
	private final static String COMMAND_CLASS_FA_CHAIN_LIST= "classFAChainList";
	private final static String COMMAND_CLASS_ADDUCT_LIST= "classAdductList";
	private final static String COMMAND_CLASS_CHAIN_C_MIN= "classChainMin";
	private final static String COMMAND_CLASS_DB_MIN= "classDBMin";
	private final static String COMMAND_CLASS_RT_MIN= "classRtMin";
	private final static String COMMAND_CLASS_OH= "classOH";
	private final static String COMMAND_CLASS_OH_MIN= "classOhMin";
	private final static String COMMAND_CLASS_OX_MIN= "classOxMin";
	private final static String COMMAND_CLASS_ADDUCT_INSENSITIVE_RT_FILTER= "classRTFilter";
	private final static String COMMAND_CLASS_PICK_BEST= "classPickBest";
	private final static String COMMAND_CLASS_EXPORT= "classExport";
	
	private final static String OUT_OPEN= "outOpen";
	
	public final static String EXPORT_OPTION_NEG = "negative ion mode";
	public final static String EXPORT_OPTION_POS = "positive ion mode";
	public final static String EXPORT_OPTION_BOTH = "both ion modes";
	private final static String[] EXPORT_OPTIONS = new String[] {EXPORT_OPTION_NEG, EXPORT_OPTION_POS, EXPORT_OPTION_BOTH};
	private final static String COMMAND_EXPORT_OPTIONS = "exportOptions";
	private final static String EXPORT = "export";
	
	
	private final static String PLACEHOLDER_OUT_PATH = "Enter path and file name for the mass list export.";
	
	private final static int PREFERRED_DISPLAY_COMPONENT_WIDTH = 150;
	private final static int PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH = 50;
	
	
	private ArrayList<AdductVO> allDefinedAdducts_;
	private ArrayList<LipidClassVO> allDefinedLipidClasses_;
	private String[] faChainListNames_;
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
	private JComboBox<String> exportOptions_;
	
	public MassListCreatorPanel()
	{
		displayPanel_ = new JPanel();
		displayPanel_.setLayout(new GridBagLayout());
		try
		{
			allDefinedAdducts_ = (new AdductParser()).parse();
			adductListNames_ = getAdductNames();
			allDefinedLipidClasses_ = (new LipidClassParser(allDefinedAdducts_)).parse();
			faChainListNames_ = getFAChainListNames();
			selectedClass_ = allDefinedLipidClasses_.get(0);
			tempClass_ = new LipidClassVO(selectedClass_);
			selectedAdduct_ = allDefinedAdducts_.get(0);
			tempAdduct_ = new AdductVO(selectedAdduct_);		
			lipidClassTable_ = new LipidClassTable("Defined lipid classes", getLipidClassNames(), selectedClass_.getLipidClass());
			displayPanel_.add(lipidClassTable_, getDefaultGridBagConstraints(0,0, GridBagConstraints.EAST, 1, 1));
			adductsTable_ = new AdductsTable("Defined adducts", getAdductNames(), selectedAdduct_.getAdductName());
			displayPanel_.add(adductsTable_, getDefaultGridBagConstraints(0,1, GridBagConstraints.EAST, 1, 1));
			lipidClassPanel_ = getLipidClassPanel(selectedClass_);
			displayPanel_.add(lipidClassPanel_, getDefaultGridBagConstraints(1,0, GridBagConstraints.EAST, 1, 1));
			adductPanel_ = getAdductPanel(selectedAdduct_);
			displayPanel_.add(adductPanel_, getDefaultGridBagConstraints(1,1, GridBagConstraints.EAST, 1, 1));
			outTextField_ = instantiatePlaceholderJTextField(PLACEHOLDER_OUT_PATH, PLACEHOLDER_OUT_PATH, 825);
			displayPanel_.add(getOutPathPanel(outTextField_), getDefaultGridBagConstraints(0,2, GridBagConstraints.CENTER, 2, 1));
			displayPanel_.add(instantiateJButton(EXPORT, "Export", true), getDefaultGridBagConstraints(0,3, GridBagConstraints.CENTER, 2, 1, new Insets(10, 10, 10, 10)));
			
			
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
		faChainListNames_ = getFAChainListNames();
		LipidClassVO vo = getVOfromName(selectedLipidClassName);
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
		lipidClassTable_ = new LipidClassTable("Defined lipid classes", getLipidClassNames(), selectedClass_.getLipidClass());
		displayPanel_.add(lipidClassTable_, getDefaultGridBagConstraints(0,0, GridBagConstraints.EAST, 1, 1));
		displayPanel_.invalidate();
		displayPanel_.updateUI();
	}
	
	private void refreshAdductScrollPane(String selectedAddductName) throws FileNotFoundException, IOException, ChemicalFormulaException
	{
		allDefinedAdducts_ = (new AdductParser()).parse();
		adductListNames_ = getAdductNames();
		for (AdductVO vo : allDefinedAdducts_)
		{ 
			if (vo.getAdductName().equals(selectedAddductName))
			{
				selectedAdduct_ = vo;
				tempAdduct_ = new AdductVO(selectedAdduct_);
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
		
		panel.add(new JLabel("Export relevant adducts for (based on sign of adduct charge): "), getDefaultGridBagConstraints(0,0, GridBagConstraints.WEST, 1, 1));
		exportOptions_ = instantiateJComboBox(COMMAND_EXPORT_OPTIONS, EXPORT_OPTIONS, 0);
		panel.add(exportOptions_, getDefaultGridBagConstraints(1,0, GridBagConstraints.WEST, 1, 1));
		
		panel.add(outPathField, getDefaultGridBagConstraints(0,1, GridBagConstraints.WEST, 2, 1));
  	JButton outPathButton = instantiateJButton(OUT_OPEN, "Browse", true);
		panel.add(outPathButton, getDefaultGridBagConstraints(2,1, GridBagConstraints.EAST, 1, 1, new Insets(10,10,10,10)));
		TitledBorder border = JOptionPanel.getTitledPanelBorder("Mass list export");
		panel.setBorder(border);
		panel.setPreferredSize(new Dimension(975,125));
		return panel;
	}
	
	private JPanel getAdductPanel(AdductVO vo)
	{
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		JTextField nameField = instantiateJTextField(COMMAND_ADDUCT_NAME, vo.getAdductName());
		addLabeledTextField(panel, 0, new JLabel("Adduct name: "), nameField);
		JTextField formulaField = instantiateJTextField(COMMAND_ADDUCT_FORMULA, vo.getFormulaString());
		addLabeledTextField(panel, 1, new JLabel("Chemical fomula: "), formulaField);
		JTextField chargeField = instantiateJTextField(COMMAND_ADDUCT_CHARGE, String.valueOf(vo.getCharge()));
		addLabeledTextField(panel, 2, new JLabel("Charge: "), chargeField);
		JPanel buttonPanel = instantiateJButtonPanel(COMMAND_ADDUCT_EXPORT, "Override this definition", "Save as new adduct definition");
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
		JTextField nameField = instantiateJTextField(COMMAND_CLASS_NAME, vo.getLipidClass());
		addLabeledTextField(panel, 0, new JLabel("Lipid class name: "), nameField);
		JTextField numberChainField = instantiateJTextField(COMMAND_CLASS_CHAIN_NUM, String.valueOf(vo.getNumberOfChains()));
		addLabeledTextField(panel, 1, new JLabel("Number of FA and/or LCB chains: "), numberChainField);
		JTextField chemicalFormula = instantiateJTextField(COMMAND_CLASS_FORMULA, vo.getHeadGroupFormulaString());
		addLabeledTextField(panel, 2, new JLabel("Chemical formula without chains: "), chemicalFormula);
		JComboBox<String> faChainList = instantiateJComboBox(COMMAND_CLASS_FA_CHAIN_LIST, faChainListNames_, findSelectedFAChainListIndex(vo));
		addLabeledComboBox(panel, 3, new JLabel("Selected FA chain list: "), faChainList);
		JScrollPane adductList = instantiateJListScrollPane(COMMAND_CLASS_ADDUCT_LIST, adductListNames_, vo, findSelectedAdductListIndices(vo));
		addLabeledJListScrollPane(panel, 4, new JLabel("Selected adducts (hold CNTR for multiple selection): "), adductList);
		addLabeledRange(panel, 5, new JLabel("Total number of C atoms in the chains: "), 
				instantiateJTextFieldRange(COMMAND_CLASS_CHAIN_C_MIN, String.valueOf(vo.getMinChainC()), String.valueOf(vo.getMaxChainC()), PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH));
		addLabeledRange(panel, 6, new JLabel("Total number of double bonds (C=C) in the chains: "), 
				instantiateJTextFieldRange(COMMAND_CLASS_DB_MIN, String.valueOf(vo.getMinChainDB()), String.valueOf(vo.getMaxChainDB()), PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH));
		addLabeledRange(panel, 7, new JLabel("Retention time (RT) range in minutes (optional): "), 
				instantiateJTextFieldRange(COMMAND_CLASS_RT_MIN, String.valueOf(vo.getRtRangeFrom()), String.valueOf(vo.getRtRangeTo()), PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH));
		JTextField ohField = instantiateJTextField(COMMAND_CLASS_OH, String.valueOf(vo.getOhNumber()));
		addLabeledTextField(panel, 8, new JLabel("Sphingolipid OH number (optional): "), ohField);
		addLabeledRange(panel, 9, new JLabel("Sphingolipid OH range (optional): "), 
				instantiateJTextFieldRange(COMMAND_CLASS_OH_MIN, String.valueOf(vo.getOhRangeFrom()), String.valueOf(vo.getOhRangeTo()), PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH));
		addLabeledRange(panel, 10, new JLabel("Oxidized lipid Ox range (optional): "), 
				instantiateJTextFieldRange(COMMAND_CLASS_OX_MIN, String.valueOf(vo.getOxRangeFrom()), String.valueOf(vo.getOxRangeTo()), PREFERRED_DISPLAY_COMPONENT_SMALLER_WIDTH));
		JCheckBox rtFilter = instantiateCheckBox(COMMAND_CLASS_ADDUCT_INSENSITIVE_RT_FILTER, vo.isAdductInsensitiveRtFilter());
		addLabeledCheckBox(panel, 11, new JLabel("Enable adduct insensitive RT filter: "), rtFilter);
		JCheckBox pickBest = instantiateCheckBox(COMMAND_CLASS_PICK_BEST, vo.isPickBestMatchBySpectrumCoverage());
		addLabeledCheckBox(panel, 12, new JLabel("Pick best match by spectrum coverage: "), pickBest);
		JPanel buttonPanel = instantiateJButtonPanel(COMMAND_CLASS_EXPORT, "Override this definition", "Save as new lipid class definition");
		panel.add(buttonPanel, getDefaultGridBagConstraints(0,13, GridBagConstraints.CENTER, 5, 1));
		
		TitledBorder border = JOptionPanel.getTitledPanelBorder("Display / edit currently selected lipid class definition");
		panel.setBorder(border);
		panel.setPreferredSize(new Dimension(550,500));
		return panel;
	}
	
	private JPanel instantiateJButtonPanel(String actionCommand, String textOverride, String textSaveNew)
	{
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		JButton buttonOverride = instantiateJButton(actionCommand, textOverride, true);
		JButton buttonSaveNew = instantiateJButton(actionCommand, textSaveNew, false);
		panel.add(buttonOverride, getDefaultGridBagConstraints(0,0, GridBagConstraints.WEST, 1, 1, new Insets(10, 10, 10, 10)));
		panel.add(buttonSaveNew, getDefaultGridBagConstraints(1,0, GridBagConstraints.EAST, 1, 1, new Insets(10, 10, 10, 10)));
		return panel;
	}
	
	private JButton instantiateJButton(String actionCommand, String text, boolean isOverride)
	{
		JButton button = new JButton(text);
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
	
	private int findSelectedFAChainListIndex(LipidClassVO vo) throws IOException
	{
		for (int i=0; i<faChainListNames_.length; i++)
		{
			if (faChainListNames_[i].equalsIgnoreCase(vo.getFaChainList()))
			{
				return i;
			}
		}
		throw new IOException(String.format("The file defining the lipid class '%s' contains a FA chain list name that does not exist!", vo.getLipidClass()));
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
	
	private void addLabeledJListScrollPane(JPanel panel, Integer yPos, JLabel label, JScrollPane jScrollPane)
	{
		panel.add(label, getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(jScrollPane, getDefaultGridBagConstraints(2,yPos, GridBagConstraints.CENTER, 3, 1));
	}
	
	private void addLabeledComboBox(JPanel panel, Integer yPos, JLabel label, JComboBox<String> comboBox)
	{
		panel.add(label, getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(comboBox, getDefaultGridBagConstraints(2,yPos, GridBagConstraints.CENTER, 3, 1));
	}
	
	private void addLabeledCheckBox(JPanel panel, Integer yPos, JLabel label, JCheckBox checkBox)
	{
		panel.add(label, getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(checkBox, getDefaultGridBagConstraints(2,yPos, GridBagConstraints.CENTER, 3, 1));
	}
	
	private void addLabeledTextField(JPanel panel, Integer yPos, JLabel label, JTextField textField)
	{
		panel.add(label, getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(textField, getDefaultGridBagConstraints(2,yPos, GridBagConstraints.EAST, 3, 1));
	}
	
	private void addLabeledRange(JPanel panel, Integer yPos, JLabel label, Pair<JTextField,JTextField> range)
	{
		JLabel minField = new JLabel("From: ");
		JLabel maxField = new JLabel("To: ");
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
					
					if ((numFrom < numTo) || (numFrom < 0 && numTo < 0))
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
					tempClass_.setOhRangeFrom(numFrom);
					tempClass_.setOhRangeTo(numTo);
				}
				catch (NumberFormatException ex) 
				{
					setWarningTextFieldBorder(textfieldFrom);
					setWarningTextFieldBorder(textfieldTo);
				}
				break;
				
			case COMMAND_CLASS_OX_MIN:
				tempClass_.setOxRangeFrom(-1);
				tempClass_.setOxRangeTo(-1);
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
					tempClass_.setOxRangeFrom(numFrom);
					tempClass_.setOxRangeTo(numTo);
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
			case COMMAND_CLASS_CHAIN_NUM:
				tempClass_.setNumberOfChains(0);
				try
				{
					int chainNum = Integer.parseInt(textfield.getText());
					if (chainNum > 0)
						setDefaultTextFieldBorder(textfield);
					else
						setWarningTextFieldBorder(textfield);
					tempClass_.setNumberOfChains(chainNum);
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
			case COMMAND_ADDUCT_EXPORT:
				if (isAdductViable())
				{
					try
					{
						exportAdduct(isOverride);
						refreshAdductScrollPane(selectedAdduct_.getAdductName());
					}
					catch (IOException | ChemicalFormulaException ex)
					{
						new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
					}
				}
				else
				{
					new WarningMessage(new JFrame(), "Error", "The adduct definition contains erroneous user-entries, please correct textfields highlighted in red before exporting.");
				}
				break;
			case COMMAND_CLASS_EXPORT:
				if (isLipidClassViable())
				{
					if (isLipidClassOxDefinitionViable())
					{
						try
						{
							exportLipidClass(isOverride);
							refreshLipidClassScrollPane(selectedClass_.getLipidClass());
						}
						catch (IOException | ChemicalFormulaException ex)
						{
							new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
						}
					}
					else
					{
						new WarningMessage(new JFrame(), "Error", "The lipid class definition may only contain a definition for the Sphingolipid OH number and OH range (with OH number being within the OH range), OR the oxidized lipid ox range. Please correct this before the export");
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
					MassListExporter exporter = new MassListExporter(outPath, lipidClassTable_.getSelectedLipidClasses(), (String)exportOptions_.getSelectedItem());
					
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
		if (tempClass_.getOhNumber() == 0 && tempClass_.getOhRangeFrom() == 0 && tempClass_.getOhRangeTo() == 0 && tempClass_.getOxRangeFrom() == 0 && tempClass_.getOxRangeTo() == 0)
		{
			return true;
		}
		else if ((tempClass_.getOxRangeFrom() == 0 && tempClass_.getOxRangeTo() == 0)
				&& (tempClass_.getOhNumber() > 0 && (tempClass_.getOhRangeFrom() >= 0 && tempClass_.getOhRangeFrom() < tempClass_.getOhRangeTo())
				&& (tempClass_.getOhRangeFrom() <= tempClass_.getOhNumber() && tempClass_.getOhRangeTo() >= tempClass_.getOhNumber())))
		{
			return true;
		}
		else if ( (tempClass_.getOhNumber() == 0 && tempClass_.getOhRangeFrom() == 0 && tempClass_.getOhRangeTo() == 0)
				&& tempClass_.getOxRangeFrom() >= 0 && tempClass_.getOxRangeFrom() < tempClass_.getOxRangeTo())
		{
			return true;
		}
		return false;
	}
	
	private boolean isLipidClassViable()
	{
		if (isLipidClassNameAvailable(tempClass_.getLipidClass())
				&& tempClass_.getNumberOfChains() != 0
				&& tempClass_.getHeadgroupFormula() != null
				&& tempClass_.getMinChainC() > 0 && tempClass_.getMinChainC() < tempClass_.getMaxChainC()
				&& (tempClass_.getMinChainDB() >= 0 && tempClass_.getMinChainDB() < tempClass_.getMaxChainDB() || (tempClass_.getMinChainDB() == 0 && tempClass_.getMaxChainDB() == 0))
				&& ((tempClass_.getRtRangeFrom() >= 0 && tempClass_.getRtRangeFrom() < tempClass_.getRtRangeTo()) || (tempClass_.getRtRangeFrom() < 0 && tempClass_.getRtRangeTo() < 0))
				&& tempClass_.getOhNumber()>=0
				&& (tempClass_.getOhRangeFrom() >= 0 && tempClass_.getOhRangeFrom() < tempClass_.getOhRangeTo() || (tempClass_.getOhRangeFrom() == 0 && tempClass_.getOhRangeTo() == 0))
				&& (tempClass_.getOxRangeFrom() >= 0 && tempClass_.getOxRangeFrom() < tempClass_.getOxRangeTo() || (tempClass_.getOxRangeFrom() == 0 && tempClass_.getOxRangeTo() == 0))
				&& !tempClass_.getAdducts().isEmpty()
				)
		{
			return true;
		}
		return false;
	}
	
	/**
	 * Exports an adduct definition file
	 * @param isOverride	true if the selected adduct definition file should be overriden / replaced with the new one
	 */
	private void exportAdduct(boolean isOverride)
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
			updateAdductForLipidClasses(selectedAdduct_, tempAdduct_);
		}
		selectedAdduct_ = tempAdduct_;
	}
	
	private void updateAdductForLipidClasses(AdductVO oldAdduct, AdductVO newAdduct)
	{
		for (LipidClassVO vo : allDefinedLipidClasses_)
		{
			boolean isChanged = false;
			ArrayList<AdductVO> newAdductVOList = new ArrayList<AdductVO>();
			for (AdductVO adduct : vo.getAdducts())
			{
				if (adduct.getAdductName().equals(oldAdduct.getAdductName()))
				{
					newAdductVOList.add(newAdduct);
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
		if (isOverride)
		{
			String originalFilePath = LipidClassExporter.buildLipidClassPath(selectedClass_.getLipidClass());
			File file = new File(originalFilePath);
			file.delete();
		}
		selectedClass_ = tempClass_;
		LipidClassExporter exporter = new LipidClassExporter(selectedClass_);
		exporter.export();
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
	
	private String[] getFAChainListNames() throws IOException
	{
		ArrayList<String> names = new ArrayList<String>();
		File folder = new File(CHAIN_LIST_FOLDER);
		if (!folder.exists())
		{
			throw new IOException(String.format("The FA chain list folder '%s' does not exist!", CHAIN_LIST_FOLDER));
		}
		File[] fileCandidates = folder.listFiles();
		for (int i=0; i<fileCandidates.length;i++)
		{
			String fileName = fileCandidates[i].getName(); 
			if (fileName.endsWith(CHAIN_LIST_SUFFIX))
			{
				names.add(fileName.substring(0, fileName.indexOf(CHAIN_LIST_SUFFIX)));
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
