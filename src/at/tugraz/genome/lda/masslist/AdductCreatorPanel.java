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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
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

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.masslist.MassListCreatorPanel.BooleanTableModel;
import at.tugraz.genome.lda.target.JOptionPanel;
import at.tugraz.genome.lda.vos.AdductVO;

/**
 * 
 * @author Leonida M. Lamp
 * 
 */
public class AdductCreatorPanel extends JFrame
{
	private static final long serialVersionUID = 1L;
	private final static String COMMAND_ADDUCT_NAME = "adductName";
	private final static String COMMAND_ADDUCT_FORMULA = "adductFormula";
	private final static String COMMAND_ADDUCT_CHARGE= "adductCharge";
	private final static String COMMAND_ADDUCT_EXPORT= "adductExport";
	private final static String COMMAND_ADDUCT_DELETE= "adductDelete";
	private MassListCreatorPanel parent_;
	private JPanel displayPanel_;
	private AdductsTable adductsTable_;
	private JPanel adductPanel_;
	private ArrayList<AdductVO> allDefinedAdducts_;
	private String[] adductListNames_;
	private AdductVO selectedAdduct_;
	private AdductVO tempAdduct_;
	
	protected AdductCreatorPanel(MassListCreatorPanel parent)
	{
		this.setTitle("Edit list of adducts");
		this.setSize(new Dimension(1050,350));
		this.setVisible(false);
		this.parent_ = parent;
		displayPanel_ = new JPanel();
		displayPanel_.setLayout(new GridBagLayout());
		try
		{
			allDefinedAdducts_ = (new AdductParser()).parse();
			adductListNames_ = getAdductNames();
			selectedAdduct_ = allDefinedAdducts_.get(0);
			tempAdduct_ = new AdductVO(selectedAdduct_);		
			adductsTable_ = new AdductsTable("Defined adducts", getAdductNames(), selectedAdduct_.getAdductName());
			addAdductsTableToDisplayPanel();
			adductPanel_ = getAdductPanel(selectedAdduct_);
			addAdductsPanelToDisplayPanel();
			
	  	this.add(displayPanel_);
		}
		catch (Exception ex)
		{
			this.add(new JLabel("Mass list creation interface unavailable due to invalid and/or missing content in the subfolder ./massListCreation/adducts"));
			this.add(new JLabel("Error message: "+ex.getMessage()));
		}
	}
	
	private void addAdductsTableToDisplayPanel()
	{
		displayPanel_.add(adductsTable_, parent_.getDefaultGridBagConstraints(0,0, GridBagConstraints.EAST, 1, 1));
	}
	
	private void addAdductsPanelToDisplayPanel()
	{
		displayPanel_.add(adductPanel_, parent_.getDefaultGridBagConstraints(1,0, GridBagConstraints.EAST, 1, 1));
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
		addAdductsTableToDisplayPanel();
		displayPanel_.invalidate();
		displayPanel_.updateUI();
		parent_.reloadLipidClassScrollPane();
		parent_.refreshLipidClassPanel();
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
		panel.add(buttonPanel, parent_.getDefaultGridBagConstraints(0,13, GridBagConstraints.CENTER, 5, 1));
		TitledBorder border = JOptionPanel.getTitledPanelBorder("Display / edit currently selected adduct definition");
		panel.setBorder(border);
		panel.setMinimumSize(new Dimension(550,175));
		panel.setPreferredSize(new Dimension(550,250));
		return panel;
	}
	
	private void addLabeledTextField(JPanel panel, Integer yPos, JLabel label, JTextField textField, String tooltips)
	{
		label.setToolTipText(tooltips);
		textField.setToolTipText(tooltips);
		panel.add(label, parent_.getDefaultGridBagConstraints(0,yPos, GridBagConstraints.WEST, 2, 1));
		panel.add(textField, parent_.getDefaultGridBagConstraints(2,yPos, GridBagConstraints.EAST, 3, 1));
	}
	
	private void refreshAdductPanel() throws IOException
	{
		displayPanel_.remove(adductPanel_);
		adductPanel_ = getAdductPanel(selectedAdduct_);
		addAdductsPanelToDisplayPanel();
		displayPanel_.invalidate();
		displayPanel_.updateUI();
	}
	
	private String[] toArray(ArrayList<String> list)
	{
		String[] arr = new String[list.size()];
		return list.toArray(arr);
	}
	
	protected void handleAdductSelection(String adduct) throws IOException, ChemicalFormulaException
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
	
	private JTextField instantiateJTextField(String actionCommand, String text)
	{
		return instantiateJTextField(actionCommand, text, parent_.getPreferredDisplayComponentWidth(0));
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
	
	private void textFieldChangeExecuter(String actionCommand, JTextField textfield)
	{
		switch (actionCommand)
		{
			case COMMAND_ADDUCT_NAME:
				tempAdduct_.setAdductName(textfield.getText());
				if (isAdductNameAvailable(textfield.getText()))
					parent_.setDefaultTextFieldBorder(textfield);
				else
					parent_.setWarningTextFieldBorder(textfield);
				break;
			case COMMAND_ADDUCT_FORMULA:
				try
				{
					tempAdduct_.setFormulaString(textfield.getText());
					parent_.setDefaultTextFieldBorder(textfield);
				}
				catch (ChemicalFormulaException ex) {parent_.setWarningTextFieldBorder(textfield);}
				break;
			case COMMAND_ADDUCT_CHARGE:
				try
				{
					tempAdduct_.setCharge(0);
					tempAdduct_.setCharge(Integer.parseInt(textfield.getText()));
					parent_.setDefaultTextFieldBorder(textfield);
				}
				catch (NumberFormatException ex) {parent_.setWarningTextFieldBorder(textfield);}
				break;
			default:
				break;
		}
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
			default:
				break;
		}
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
				parent_.updateAdductForLipidClasses(selectedAdduct_, tempAdduct_, true);
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
				parent_.updateAdductForLipidClasses(selectedAdduct_, tempAdduct_, false);
			}
			selectedAdduct_ = tempAdduct_;
			refreshAdductScrollPane(selectedAdduct_.getAdductName());
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
	
	protected ArrayList<AdductVO> getAllDefinedAdducts()
	{
		return allDefinedAdducts_;
	}

	protected String[] getAdductListNames()
	{
		return adductListNames_;
	}

	protected String[] getAdductNames()
	{
		ArrayList<String> names = new ArrayList<String>();
		for (AdductVO vo : allDefinedAdducts_)
		{
			names.add(vo.getAdductName());
		}
		return toArray(names);
	}
	
	protected void open()
	{
		this.setVisible(true);
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
    	this.setMinimumSize(new Dimension(400,175));
    	this.setPreferredSize(new Dimension(400,250));
    	selectionTablePanel_ = new JPanel();
    	generateSelectionTablePanel(initializeTableData(adductNames), selected);
    	this.setLayout(new GridBagLayout());
    	this.add(selectionTablePanel_, parent_.getDefaultGridBagConstraints(0, 1, GridBagConstraints.CENTER, 1, 1));
    	this.setBorder(JOptionPanel.getTitledPanelBorder(title));
    }
    
    private void generateSelectionTablePanel(Object[][] tableData, String selected)
    {
    	String[] columnNames = { "adduct name" };
    	BooleanTableModel model = parent_.new BooleanTableModel(tableData, columnNames);
    	displayTable_ = new JTable(model);
    	scrollPane_ = new JScrollPane(displayTable_);
    	scrollPane_.setMinimumSize(new Dimension(325,125));
    	scrollPane_.setPreferredSize(new Dimension(325,175));
    	displayTable_.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    	ListSelectionModel selectionModel = displayTable_.getSelectionModel();
    	selectionModel.setAnchorSelectionIndex(model.indexOf(selected,COLUMN_NAME));
    	selectionModel.setLeadSelectionIndex(model.indexOf(selected,COLUMN_NAME));
    	selectionModel.addListSelectionListener(new ListSelectionListener() {
        public void valueChanged(ListSelectionEvent e) {
        	try
        	{
        		handleAdductSelection((String)model.getValueAt(((ListSelectionModel)e.getSource()).getMaxSelectionIndex(), COLUMN_NAME));
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
}


