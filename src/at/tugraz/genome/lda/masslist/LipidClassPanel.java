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

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.vos.AdductVO;
import org.apache.commons.math3.util.Pair;


/**
 * 
 *@author Leonida M. Lamp 
 * 
 */
public class LipidClassPanel extends JPanel
{
	private static final long serialVersionUID = 1L;
	protected final static String COMMAND_CLASS_NAME = "className";
	protected final static String COMMAND_CLASS_FA_CHAIN_NUM= "classFAChainNum";
	protected final static String COMMAND_CLASS_LCB_CHAIN_NUM= "classLCBChainNum";
	protected final static String COMMAND_CLASS_FORMULA= "classFormula";
	protected final static String COMMAND_CLASS_FA_CHAIN_LIST= "classFAChainList";
	protected final static String COMMAND_CLASS_LCB_CHAIN_LIST= "classLCBChainList";
	protected final static String COMMAND_CLASS_ADDUCT_LIST= "classAdductList";
	protected final static String COMMAND_CLASS_CHAIN_C= "classChainMin";
	protected final static String COMMAND_CLASS_DB= "classDBMin";
	protected final static String COMMAND_CLASS_RT= "classRtMin";
	protected final static String COMMAND_CLASS_OH= "classOhMin";
	protected final static String COMMAND_CLASS_ADDUCT_INSENSITIVE_RT_FILTER= "classRTFilter";
	protected final static String COMMAND_CLASS_PICK_BEST= "classPickBest";
	protected final static String COMMAND_CLASS_EXPORT= "classExport";
	protected final static String COMMAND_CLASS_DELETE= "classDelete";
	protected final static String COMMAND_ALL_OVERRIDE_ALL= "overrideAll";
	protected final static String COMMAND_ALL_CHAIN_NUM = "allChainNum";
	protected final static String COMMAND_ALL_FA_CHAIN_LIST = "allFAChainList";
	protected final static String COMMAND_ALL_LCB_CHAIN_LIST = "allLCBChainList";
	protected final static String COMMAND_ALL_ADDUCT_LIST= "allAdductList";
	protected final static String COMMAND_ALL_CHAIN_C= "allChainMin";
	protected final static String COMMAND_ALL_DB= "allDBMin";
	protected final static String COMMAND_ALL_RT= "allRtMin";
	protected final static String COMMAND_ALL_OH= "allOhMin";
	protected final static String COMMAND_ALL_ADDUCT_INSENSITIVE_RT_FILTER= "allRTFilter";
	protected final static String COMMAND_ALL_PICK_BEST= "allPickBest";
	protected final static String OPTION_EXPORT_ALL= "All listed in table";
	protected final static String OPTION_EXPORT_SELECTED= "Selected only";
	private MassListCreatorPanel parent_;
	private String[] faChainListNames_;
	private String[] lcbChainListNames_;
	private JTextField numberLCBChainField_;
	private JComboBox<String> lcbChainList_;
	private LipidClassVO selectedClass_;
	private LipidClassVO tempClass_;
	private boolean isSingleEdit_;
	private JComboBox<String> exportOption_;
	
	
	public LipidClassPanel(MassListCreatorPanel parent, boolean isSingleEdit, String[] faChainListNames, String[] lcbChainListNames) throws IOException
	{
		this.parent_ = parent;
		this.selectedClass_ = isSingleEdit ? parent.getSelectedClass() : new LipidClassVO(faChainListNames[0],lcbChainListNames[0]);
		this.tempClass_ = isSingleEdit ? parent.getTempClass() : new LipidClassVO(faChainListNames[0],lcbChainListNames[0]);
		this.faChainListNames_ = faChainListNames;
		this.lcbChainListNames_ = lcbChainListNames;
		this.isSingleEdit_ = isSingleEdit;
		init();
	}
	
	private void init() throws IOException
	{
		this.setLayout(new GridBagLayout());
		int y=0;
		if (isSingleEdit_)
		{
			
			addLabeledTextField(this, y++, new JLabel("Lipid (sub)class name: "),
					instantiateJTextField(COMMAND_CLASS_NAME, selectedClass_.getLipidClass()), TooltipTexts.MASSLIST_CLASS_NAME);
			
			addLabeledTextField(this, y++, new JLabel("Chemical formula without chains: "), 
					instantiateJTextField(COMMAND_CLASS_FORMULA, selectedClass_.getHeadGroupFormulaString()), TooltipTexts.MASSLIST_CLASS_FORMULA);
			
		}
		else
		{
			exportOption_ = parent_.instantiateJComboBox(new String[] {OPTION_EXPORT_ALL,OPTION_EXPORT_SELECTED}, 0, isSingleEdit_);
			addLabeledComboBox(this, y++, new JLabel("Lipid (sub)classes to export: "), exportOption_,  TooltipTexts.MASSLIST_CLASS_EXPORT_OPTION);
			JSeparator sep = new JSeparator();
			sep.setPreferredSize(new Dimension(520,10));
			this.add(sep, new GridBagConstraints(0,y++,6,1,1.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.VERTICAL,new Insets(10,0,0,0),0,2));
		}
		if (!isSingleEdit_) addApplyJButton(y,COMMAND_ALL_CHAIN_NUM,TooltipTexts.MASSLIST_CLASS_APPLY);
		numberLCBChainField_ = instantiateJTextField(COMMAND_CLASS_LCB_CHAIN_NUM, String.valueOf(selectedClass_.getNumberOfLCBChains()));
		addNumberOfChainSelection(this, y++, parent_.getPreferredDisplayComponentWidthSmaller(isSingleEdit_ ? 0 : 25),
				instantiateJTextField(COMMAND_CLASS_FA_CHAIN_NUM, String.valueOf(selectedClass_.getNumberOfFAChains())), numberLCBChainField_, 
				TooltipTexts.MASSLIST_CLASS_NUM_CHAIN);	
		
		if (!isSingleEdit_) addApplyJButton(y,COMMAND_ALL_FA_CHAIN_LIST,TooltipTexts.MASSLIST_CLASS_APPLY);
		JComboBox<String> faChainList = instantiateJComboBox(COMMAND_CLASS_FA_CHAIN_LIST, faChainListNames_, 
				findSelectedChainListIndex(selectedClass_.getFAChainList()));
		addLabeledComboBox(this, y++, new JLabel("Selected FA chain list: "), faChainList, TooltipTexts.MASSLIST_CLASS_FA_LIST);
		
		if (!isSingleEdit_) addApplyJButton(y,COMMAND_ALL_LCB_CHAIN_LIST,TooltipTexts.MASSLIST_CLASS_APPLY);
		lcbChainList_ = instantiateJComboBox(COMMAND_CLASS_LCB_CHAIN_LIST, lcbChainListNames_, findSelectedChainListIndex(selectedClass_.getLCBChainList()));
		lcbChainList_.setEnabled(selectedClass_.getNumberOfLCBChains() > 0);
		addLabeledComboBox(this, y++, new JLabel("Selected SPB chain list: "), lcbChainList_,  TooltipTexts.MASSLIST_CLASS_SPB_LIST);
		
		if (!isSingleEdit_) addApplyJButton(y,COMMAND_ALL_ADDUCT_LIST,TooltipTexts.MASSLIST_CLASS_APPLY);
		addLabeledJListScrollPane(this, y++, new JLabel("Selected adducts (hold CNTR for multiselection): "), 
				instantiateJListScrollPane(
						COMMAND_CLASS_ADDUCT_LIST, parent_.getAdductCreatorPanel().getAdductListNames(),
						!isSingleEdit_ && selectedClass_.getAdducts().isEmpty() ? new int[] {} : parent_.findSelectedAdductListIndices(selectedClass_)),
				TooltipTexts.MASSLIST_CLASS_ADDUCT_SELECTION);
		
		if (!isSingleEdit_) addApplyJButton(y,COMMAND_ALL_OH,TooltipTexts.MASSLIST_CLASS_APPLY);
		addLabeledRange(this, y++, new JLabel("Chain OH / oxidation range: "), 
				instantiateJTextFieldRange(COMMAND_CLASS_OH, String.valueOf(selectedClass_.getOhRangeFrom()), String.valueOf(selectedClass_.getOhRangeTo()), 
				parent_.getPreferredDisplayComponentWidthSmaller(isSingleEdit_ ? 0 : 25)), TooltipTexts.MASSLIST_CLASS_OH_RANGE);
		
		if (!isSingleEdit_) addApplyJButton(y,COMMAND_ALL_CHAIN_C,TooltipTexts.MASSLIST_CLASS_APPLY);
		addLabeledRange(this, y++, new JLabel("Total number of chain C atoms: "), 
				instantiateJTextFieldRange(COMMAND_CLASS_CHAIN_C, String.valueOf(selectedClass_.getMinChainC()), String.valueOf(selectedClass_.getMaxChainC()), 
			  parent_.getPreferredDisplayComponentWidthSmaller(isSingleEdit_ ? 0 : 25)), TooltipTexts.MASSLIST_CLASS_C_RANGE);
		
		if (!isSingleEdit_) addApplyJButton(y,COMMAND_ALL_DB,TooltipTexts.MASSLIST_CLASS_APPLY);
		addLabeledRange(this, y++, new JLabel("Total number of chain double bonds: "), 
				instantiateJTextFieldRange(COMMAND_CLASS_DB, String.valueOf(selectedClass_.getMinChainDB()), String.valueOf(selectedClass_.getMaxChainDB()), 
				parent_.getPreferredDisplayComponentWidthSmaller(isSingleEdit_ ? 0 : 25)), TooltipTexts.MASSLIST_CLASS_DB_RANGE);
		
		if (!isSingleEdit_) addApplyJButton(y,COMMAND_ALL_RT,TooltipTexts.MASSLIST_CLASS_APPLY);
		addLabeledRange(this, y++, new JLabel("Retention time (RT) range in minutes: "), 
				instantiateJTextFieldRange(COMMAND_CLASS_RT, String.valueOf(selectedClass_.getRtRangeFrom()), String.valueOf(selectedClass_.getRtRangeTo()), 
				parent_.getPreferredDisplayComponentWidthSmaller(isSingleEdit_ ? 0 : 25)), TooltipTexts.MASSLIST_CLASS_RT_RANGE);
		
		if (!isSingleEdit_) addApplyJButton(y,COMMAND_ALL_ADDUCT_INSENSITIVE_RT_FILTER,TooltipTexts.MASSLIST_CLASS_APPLY);
		JCheckBox rtFilter = instantiateCheckBox(COMMAND_CLASS_ADDUCT_INSENSITIVE_RT_FILTER, selectedClass_.isAdductInsensitiveRtFilter());
		addLabeledCheckBox(this, y++, new JLabel("Enable adduct insensitive RT filter: "), rtFilter, TooltipTexts.MASSLIST_CLASS_RT_FILTER);
		
		if (!isSingleEdit_) addApplyJButton(y,COMMAND_ALL_PICK_BEST,TooltipTexts.MASSLIST_CLASS_APPLY);
		JCheckBox pickBest = instantiateCheckBox(COMMAND_CLASS_PICK_BEST, selectedClass_.isPickBestMatchBySpectrumCoverage());
		addLabeledCheckBox(this, y++, new JLabel("Pick best match by spectrum coverage: "), pickBest, TooltipTexts.MASSLIST_CLASS_PICK_BEST);
		JPanel buttonPanel;
		if (isSingleEdit_)
			buttonPanel = instantiateJButtonPanel(COMMAND_CLASS_DELETE, COMMAND_CLASS_EXPORT, "Delete", "Override", "Save New",
					TooltipTexts.MASSLIST_CLASS_DELETE, TooltipTexts.MASSLIST_CLASS_OVERRIDE, TooltipTexts.MASSLIST_CLASS_SAVE_NEW);
		else
			buttonPanel = instantiateJButtonPanel(COMMAND_ALL_OVERRIDE_ALL,"Apply all",TooltipTexts.MASSLIST_CLASS_OVERRIDE_ALL);
		
		this.add(buttonPanel, parent_.getDefaultGridBagConstraints(0,y, GridBagConstraints.CENTER, 6, 1));
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
		panel.add(buttonDelete, parent_.getDefaultGridBagConstraints(0,0, GridBagConstraints.WEST, 1, 1, new Insets(10, 10, 10, 10)));
		panel.add(buttonOverride, parent_.getDefaultGridBagConstraints(1,0, GridBagConstraints.WEST, 1, 1, new Insets(10, 10, 10, 10)));
		panel.add(buttonSaveNew, parent_.getDefaultGridBagConstraints(2,0, GridBagConstraints.EAST, 1, 1, new Insets(10, 10, 10, 10)));
		return panel;
	}
	
	private JPanel instantiateJButtonPanel(String actionCommand,String text,String tooltips)
	{
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		JButton button = instantiateJButton(actionCommand, text, true, tooltips);
		panel.add(button, parent_.getDefaultGridBagConstraints(0,0, GridBagConstraints.CENTER, 1, 1, new Insets(10, 10, 10, 10)));
		return panel;
	}
	
	private void addApplyJButton(int yPos, String actionCommand, String tooltips)
	{
		JButton button = new JButton("Apply");
		button.setToolTipText(tooltips);
		button.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	jButtonExecuter(actionCommand);
		  }
	  });
		this.add(button, parent_.getDefaultGridBagConstraints(5,yPos, GridBagConstraints.CENTER, 1, 1));
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
		jComboBox.setPreferredSize(new Dimension(parent_.getPreferredDisplayComponentWidth(isSingleEdit_ ? 0 : 50),20));
		return jComboBox;
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
	
	private JScrollPane instantiateJListScrollPane(String actionCommand, String[] entries, int[] indices) throws IOException
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
		scrollPane.setPreferredSize(new Dimension(parent_.getPreferredDisplayComponentWidth(isSingleEdit_ ? 0 : 50),100));
		return scrollPane;
	}
	
	private JTextField instantiateJTextField(String actionCommand, String text)
	{
		return instantiateJTextField(actionCommand, text, parent_.getPreferredDisplayComponentWidth(isSingleEdit_ ? 0 : 50));
	}
	
	private Pair<JTextField,JTextField> instantiateJTextFieldRange(String actionCommand, String textFrom, String textTo, Integer width)
	{
		JTextField textFieldFrom = new JTextField(textFrom);
		textFieldFrom.setPreferredSize(new Dimension(width,20));
		parent_.setDefaultTextFieldBorder(textFieldFrom);
		JTextField textFieldTo = new JTextField(textTo);
		textFieldTo.setPreferredSize(new Dimension(width,20));
		parent_.setDefaultTextFieldBorder(textFieldTo);
		
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
		parent_.setDefaultTextFieldBorder(textField);
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
		panel.add(label, parent_.getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(jScrollPane, parent_.getDefaultGridBagConstraints(2,yPos, GridBagConstraints.CENTER, 3, 1));
	}
	
	private void addNumberOfChainSelection(JPanel panel, Integer yPos, Integer width,
			JTextField numberFAChainField, JTextField numberSPBChainField, String tooltips)
	{
		JLabel labelNum = new JLabel("Number of chains: ");
		labelNum.setToolTipText(tooltips);
		JLabel labelFA = new JLabel("FA: ");
		labelFA.setToolTipText(tooltips);
		JLabel labelSPB = new JLabel("SPB: ");
		labelSPB.setToolTipText(tooltips);
		numberFAChainField.setPreferredSize(new Dimension(width,20));
		numberFAChainField.setToolTipText(tooltips);
		numberSPBChainField.setPreferredSize(new Dimension(width,20));
		numberSPBChainField.setToolTipText(tooltips);
		
		panel.add(labelNum, parent_.getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 1, 1));
		panel.add(labelFA, parent_.getDefaultGridBagConstraints(1,yPos, GridBagConstraints.EAST, 1, 1));
		panel.add(numberFAChainField, parent_.getDefaultGridBagConstraints(2,yPos, GridBagConstraints.WEST, 1, 1));
		panel.add(labelSPB, parent_.getDefaultGridBagConstraints(3,yPos, GridBagConstraints.EAST, 1, 1));
		panel.add(numberSPBChainField, parent_.getDefaultGridBagConstraints(4,yPos, GridBagConstraints.EAST, 1, 1));
	}
	
	private void addLabeledComboBox(JPanel panel, Integer yPos, JLabel label, JComboBox<String> comboBox, String tooltips)
	{
		label.setToolTipText(tooltips);
		comboBox.setToolTipText(tooltips);
		panel.add(label, parent_.getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(comboBox, parent_.getDefaultGridBagConstraints(2,yPos, GridBagConstraints.CENTER, 3, 1));
	}
	
	private void addLabeledCheckBox(JPanel panel, Integer yPos, JLabel label, JCheckBox checkBox, String tooltips)
	{
		label.setToolTipText(tooltips);
		checkBox.setToolTipText(tooltips);
		panel.add(label, parent_.getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(checkBox, parent_.getDefaultGridBagConstraints(2,yPos, GridBagConstraints.CENTER, 3, 1));
	}
	
	private void addLabeledTextField(JPanel panel, Integer yPos, JLabel label, JTextField textField, String tooltips)
	{
		label.setToolTipText(tooltips);
		textField.setToolTipText(tooltips);
		panel.add(label, parent_.getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(textField, parent_.getDefaultGridBagConstraints(2,yPos, GridBagConstraints.EAST, 3, 1));
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
		panel.add(label, parent_.getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 1, 1));
		panel.add(minField, parent_.getDefaultGridBagConstraints(1,yPos, GridBagConstraints.EAST, 1, 1));
		panel.add(range.getKey(), parent_.getDefaultGridBagConstraints(2,yPos, GridBagConstraints.WEST, 1, 1));
		panel.add(maxField, parent_.getDefaultGridBagConstraints(3,yPos, GridBagConstraints.EAST, 1, 1));
		panel.add(range.getValue(), parent_.getDefaultGridBagConstraints(4,yPos, GridBagConstraints.EAST, 1, 1));
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
			case COMMAND_CLASS_DELETE:
				synchronizeCurrentLipidClass();
				parent_.deleteLipidClass();
				break;
			case COMMAND_CLASS_EXPORT:
				parent_.synchronizeCurrentLipidClass(this);
				if (parent_.isLipidClassViable())
				{
					if (!isOverride ||
							( isOverride && JOptionPane.showConfirmDialog(new JFrame(),
								String.format("Would you like to override the lipid (sub)class definition '%s'?", selectedClass_.getLipidClass()),
							  "Override definition", JOptionPane.YES_NO_OPTION) == 0 ) )
					{
						parent_.exportLipidClass(isOverride);
					}
				}	
				else
				{
					new WarningMessage(new JFrame(), "Error", "The lipid class definition contains erroneous user-input, please correct textfields highlighted in red before exporting and make sure at least one adduct is selected.");
				}
				break;
			case COMMAND_ALL_OVERRIDE_ALL:
				if (parent_.isGeneralDefinitionViable(tempClass_))
				{
					parent_.exportAll(tempClass_, parent_.findLipidClassesToExport((String)exportOption_.getSelectedItem()));
				}
				else
				{
					new WarningMessage(new JFrame(), "Error", "The lipid class definition contains erroneous user-input, please correct textfields highlighted in red before exporting and make sure at least one adduct is selected.");
				}
				break;
			default:
				break;
		}
	}
	
	private void jButtonExecuter(String actionCommand)
	{
		ArrayList<LipidClassVO> toExport = parent_.findLipidClassesToExport((String)exportOption_.getSelectedItem());
		switch (actionCommand)
		{
			case COMMAND_ALL_CHAIN_NUM:
				if (parent_.isChainNumViable(tempClass_))
				{
					parent_.exportChainNumbersToSelectedClasses(tempClass_.getNumberOfFAChains(), tempClass_.getNumberOfLCBChains(), toExport);
				}
				else
				{
					new WarningMessage(new JFrame(), "Error", 
							"The chain number definition contains erroneous user-input, please correct textfields highlighted in red before exporting.");
				}
				break;
			case COMMAND_ALL_FA_CHAIN_LIST:
				parent_.exportFAChainListToSelectedClasses(tempClass_.getFAChainList(), toExport);
				break;
			case COMMAND_ALL_LCB_CHAIN_LIST:
				parent_.exportSPBChainListToSelectedClasses(tempClass_.getLCBChainList(), toExport);
				break;
			case COMMAND_ALL_ADDUCT_LIST:
				if (parent_.isAdductListViable(tempClass_))
				{
					parent_.exportAdductListToSelectedClasses(tempClass_.getAdducts(), toExport);
				}
				else
				{
					new WarningMessage(new JFrame(), "Error", 
							"You must select at least one Adduct to export!");
				}
				break;
			case COMMAND_ALL_CHAIN_C:
				if (parent_.isCNumViable(tempClass_))
				{
					parent_.exportCNumToSelectedClasses(tempClass_.getMinChainC(), tempClass_.getMaxChainC(), toExport);
				}
				else
				{
					new WarningMessage(new JFrame(), "Error", 
							"The definition of chain C atom numbers contains erroneous user-input, please correct textfields highlighted in red before exporting.");
				}
				break;
			case COMMAND_ALL_DB:
				if (parent_.isDBNumViable(tempClass_))
				{
					parent_.exportDBNumToSelectedClasses(tempClass_.getMinChainDB(), tempClass_.getMaxChainDB(), toExport);
				}
				else
				{
					new WarningMessage(new JFrame(), "Error", 
							"The definition of chain double bond numbers contains erroneous user-input, please correct textfields highlighted in red before exporting.");
				}
				break;
			case COMMAND_ALL_RT:
				if (parent_.isRTViable(tempClass_))
				{
					parent_.exportRTToSelectedClasses(tempClass_.getRtRangeFrom(), tempClass_.getRtRangeTo(), toExport);
				}
				else
				{
					new WarningMessage(new JFrame(), "Error", 
							"The retention time definition contains erroneous user-input, please correct textfields highlighted in red before exporting.");
				}
				break;
			case COMMAND_ALL_OH:
				if (parent_.isLipidClassOxDefinitionViable(tempClass_))
				{
					parent_.exportOxNumToSelectedClasses(tempClass_.getOhRangeFrom(), tempClass_.getOhRangeTo(), toExport);
				}
				else
				{
					new WarningMessage(new JFrame(), "Error", 
							"The definition of chain double bond numbers contains erroneous user-input, please correct textfields highlighted in red before exporting.");
				}
				break;
			case COMMAND_ALL_ADDUCT_INSENSITIVE_RT_FILTER:
				parent_.exportPickToSelectedClasses(tempClass_.isAdductInsensitiveRtFilter(), toExport);
				break;
			case COMMAND_ALL_PICK_BEST:
				parent_.exportPickToSelectedClasses(tempClass_.isPickBestMatchBySpectrumCoverage(), toExport);
				break;
			default:
				break;
		}
	}
	
	private void textFieldChangeExecuterRange(String actionCommand, JTextField textfieldFrom, JTextField textfieldTo)
	{
		switch (actionCommand)
		{
			case COMMAND_CLASS_CHAIN_C:
				tempClass_.setMinChainC(-1);
				tempClass_.setMaxChainC(-1);
				try
				{
					int numFrom = Integer.parseInt(textfieldFrom.getText());
					int numTo = Integer.parseInt(textfieldTo.getText());
					
					if (numFrom > 0 && numTo > 0 && numFrom < numTo)
					{
						parent_.setDefaultTextFieldBorder(textfieldFrom);
						parent_.setDefaultTextFieldBorder(textfieldTo);
					}	
					else
					{
						parent_.setWarningTextFieldBorder(textfieldFrom);
						parent_.setWarningTextFieldBorder(textfieldTo);
					}
					tempClass_.setMinChainC(numFrom);
					tempClass_.setMaxChainC(numTo);
				}
				catch (NumberFormatException ex) 
				{
					parent_.setWarningTextFieldBorder(textfieldFrom);
					parent_.setWarningTextFieldBorder(textfieldTo);
				}
				break;	
			case COMMAND_CLASS_DB:
				tempClass_.setMinChainDB(-1);
				tempClass_.setMaxChainDB(-1);
				try
				{
					int numFrom = Integer.parseInt(textfieldFrom.getText());
					int numTo = Integer.parseInt(textfieldTo.getText());
					
					if ((numFrom < numTo) || (numFrom < 1 && numTo < 1))
					{
						parent_.setDefaultTextFieldBorder(textfieldFrom);
						parent_.setDefaultTextFieldBorder(textfieldTo);
					}	
					else
					{
						parent_.setWarningTextFieldBorder(textfieldFrom);
						parent_.setWarningTextFieldBorder(textfieldTo);
					}
					tempClass_.setMinChainDB(numFrom);
					tempClass_.setMaxChainDB(numTo);
				}
				catch (NumberFormatException ex) 
				{
					parent_.setWarningTextFieldBorder(textfieldFrom);
					parent_.setWarningTextFieldBorder(textfieldTo);
				}
				break;	
			case COMMAND_CLASS_RT:
				tempClass_.setRtRangeFrom(-1);
				tempClass_.setRtRangeTo(-1);
				try
				{
					int numFrom = Integer.parseInt(textfieldFrom.getText());
					int numTo = Integer.parseInt(textfieldTo.getText());
					
					if ((numFrom < numTo) || (numFrom <= 0 && numTo <= 0))
					{
						parent_.setDefaultTextFieldBorder(textfieldFrom);
						parent_.setDefaultTextFieldBorder(textfieldTo);
					}	
					else
					{
						parent_.setWarningTextFieldBorder(textfieldFrom);
						parent_.setWarningTextFieldBorder(textfieldTo);
					}
					tempClass_.setRtRangeFrom(numFrom);
					tempClass_.setRtRangeTo(numTo);
				}
				catch (NumberFormatException ex) 
				{
					parent_.setWarningTextFieldBorder(textfieldFrom);
					parent_.setWarningTextFieldBorder(textfieldTo);
				}
				break;	

			case COMMAND_CLASS_OH:
				tempClass_.setOhRangeFrom(0);
				tempClass_.setOhRangeTo(0);
				try
				{
					int numFrom = Integer.parseInt(textfieldFrom.getText());
					int numTo = Integer.parseInt(textfieldTo.getText());
					
					if (numFrom >= 0 && numTo > 0 && numFrom <= numTo)
					{
						parent_.setDefaultTextFieldBorder(textfieldFrom);
						parent_.setDefaultTextFieldBorder(textfieldTo);
					}	
					else
					{
						parent_.setWarningTextFieldBorder(textfieldFrom);
						parent_.setWarningTextFieldBorder(textfieldTo);
					}
					tempClass_.setOhRangeFrom(numFrom);
					tempClass_.setOhRangeTo(numTo);
				}
				catch (NumberFormatException ex) 
				{
					parent_.setWarningTextFieldBorder(textfieldFrom);
					parent_.setWarningTextFieldBorder(textfieldTo);
				}
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
					for (AdductVO vo : parent_.getAdductCreatorPanel().getAllDefinedAdducts())
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
	
	private void textFieldChangeExecuter(String actionCommand, JTextField textfield)
	{
		switch (actionCommand)
		{
			case COMMAND_CLASS_NAME:
				tempClass_.setLipidClass(textfield.getText());
				if (parent_.isLipidClassNameAvailable(textfield.getText()))
					parent_.setDefaultTextFieldBorder(textfield);
				else
					parent_.setWarningTextFieldBorder(textfield);
				break;
			case COMMAND_CLASS_FA_CHAIN_NUM:
				tempClass_.setNumberOfFAChains(0);
				try
				{
					int chainNum = Integer.parseInt(textfield.getText());
					if (chainNum > 0)
						parent_.setDefaultTextFieldBorder(textfield);
					else
						parent_.setWarningTextFieldBorder(textfield);
					tempClass_.setNumberOfFAChains(chainNum);
				}
				catch (NumberFormatException ex) {parent_.setWarningTextFieldBorder(textfield);}
				break;
			case COMMAND_CLASS_LCB_CHAIN_NUM:
				tempClass_.setNumberOfLCBChains(0);
				try
				{
					int chainNum = Integer.parseInt(textfield.getText());
					if (chainNum >= 0)
						parent_.setDefaultTextFieldBorder(textfield);
					else
						parent_.setWarningTextFieldBorder(textfield);
					tempClass_.setNumberOfLCBChains(chainNum);
					lcbChainList_.setEnabled(chainNum > 0);
				}
				catch (NumberFormatException ex) {parent_.setWarningTextFieldBorder(textfield);}
				break;
			case COMMAND_CLASS_FORMULA:
				try
				{
					tempClass_.setHeadGroupFormulaString(textfield.getText());
					parent_.setDefaultTextFieldBorder(textfield);
				}
				catch (ChemicalFormulaException ex) {parent_.setWarningTextFieldBorder(textfield);}
				break;
				
			default:
				break;
		}
	}

	public LipidClassVO getSelectedClass()
	{
		return selectedClass_;
	}

	public LipidClassVO getTempClass()
	{
		return tempClass_;
	}
	
	void synchronizeCurrentLipidClass()
	{
		this.selectedClass_ = parent_.getSelectedClass();
		this.tempClass_ = parent_.getTempClass();
	}
	
	
}
