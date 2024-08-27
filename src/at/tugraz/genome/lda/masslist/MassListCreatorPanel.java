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
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.border.TitledBorder;
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

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class MassListCreatorPanel extends JPanel 
{
	private static final long serialVersionUID = 1L;
	private static final int LIPID_CLASS_HEIGHT = 675;
	public final static String CHAIN_LIST_FOLDER = "./fattyAcids/";
	public final static String CHAIN_LIST_SUFFIX = ".xlsx";
	
	final static String COMMAND_EDIT_ADDUCT = "editAdduct";
	
	private final static String OUT_OPEN= "outOpen";
	
	public final static String EXPORT_OPTION_NEG = "negative ion mode";
	public final static String EXPORT_OPTION_POS = "positive ion mode";
	public final static String EXPORT_OPTION_BOTH = "both ion modes";
	private final static String[] EXPORT_OPTIONS_ION_MODE = new String[] {EXPORT_OPTION_NEG, EXPORT_OPTION_POS, EXPORT_OPTION_BOTH};
	public final static String EXPORT_FORMAT_LDA = "LDA Mass List";
	public final static String EXPORT_FORMAT_LONG_LIST = "Merged list (long)";
	public final static String EXPORT_FORMAT_SHORT_LIST = "Merged list (short)";
	private final static String[] EXPORT_OPTIONS_FORMAT = new String[] {EXPORT_FORMAT_LDA, EXPORT_FORMAT_LONG_LIST, EXPORT_FORMAT_SHORT_LIST};
	private final static String EXPORT = "export";
	
	
	private final static String PLACEHOLDER_OUT_PATH = "Enter path and file name for the mass list export.";
	
	private LipidClassPanel lipidClassPanelSingle_;
	private LipidClassPanel lipidClassPanelMulti_;
	private AdductCreatorPanel adductCreatorPanel_;
	private ArrayList<AdductVO> allDefinedAdducts_;
	private ArrayList<LipidClassVO> allDefinedLipidClasses_;
	private String[] faChainListNames_;
	private String[] lcbChainListNames_;
	private String[] adductListNames_;
	private LipidClassVO selectedClass_;
	private LipidClassVO tempClass_;
	private JPanel displayPanel_;
	private JTabbedPane lipidClassPane_;
	private LipidClassTable lipidClassTable_;

	private JTextField outTextField_;
	private Path previousSelection_ = null;
	private JComboBox<String> exportedIonMode_;
	private JComboBox<String> exportedFormat_;
//	private JTextField numberLCBChainField_;
//	private JComboBox<String> lcbChainList_;
	
	public MassListCreatorPanel()
	{
		displayPanel_ = new JPanel();
		displayPanel_.setLayout(new GridBagLayout());
		try
		{
			adductCreatorPanel_ = new AdductCreatorPanel(this);
			allDefinedAdducts_ = adductCreatorPanel_.getAllDefinedAdducts();
			adductListNames_ = adductCreatorPanel_.getAdductNames();
			allDefinedLipidClasses_ = (new LipidClassParser(allDefinedAdducts_)).parse();
			faChainListNames_ = getChainListNames(true);
			lcbChainListNames_ = getChainListNames(false);
			selectedClass_ = allDefinedLipidClasses_.get(0);
			tempClass_ = new LipidClassVO(selectedClass_);	
			lipidClassTable_ = new LipidClassTable("Defined lipid (sub)classes", getLipidClassNames(), selectedClass_.getLipidClass());
			displayPanel_.add(lipidClassTable_, getDefaultGridBagConstraints(0,0, GridBagConstraints.EAST, 1, 1));
//			adductsTable_ = new AdductsTable("Defined adducts", getAdductNames(), selectedAdduct_.getAdductName());
//			displayPanel_.add(adductsTable_, getDefaultGridBagConstraints(0,1, GridBagConstraints.EAST, 1, 1));
			JPanel lipidClassPanel = initLipidClassPanel();
			displayPanel_.add(lipidClassPanel, getDefaultGridBagConstraints(1,0, GridBagConstraints.EAST, 1, 1));
//			adductPanel_ = getAdductPanel(selectedAdduct_);
//			displayPanel_.add(adductPanel_, getDefaultGridBagConstraints(1,1, GridBagConstraints.EAST, 1, 1));
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
	
	void reloadLipidClassScrollPane() throws FileNotFoundException, IOException, ChemicalFormulaException
	{
		reloadLipidClassScrollPane(selectedClass_.getLipidClass());
	}
	
	private void reloadLipidClassScrollPane(String selectedLipidClassName) throws FileNotFoundException, IOException, ChemicalFormulaException
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
		lipidClassTable_ = new LipidClassTable("Defined lipid (sub)classes", getLipidClassNames(), selectedClass_.getLipidClass());
		refreshLipidClassScrollPane();
	}
	
	private void refreshLipidClassScrollPane()
	{
		displayPanel_.remove(lipidClassTable_);
		displayPanel_.add(lipidClassTable_, getDefaultGridBagConstraints(0,0, GridBagConstraints.EAST, 1, 1));
		displayPanel_.invalidate();
		displayPanel_.updateUI();
	}
	
	void refreshLipidClassPanel() throws IOException
	{
		int index = lipidClassPane_.indexOfComponent(lipidClassPanelSingle_);
		lipidClassPane_.remove(lipidClassPanelSingle_);
//		displayPanel_.remove(lipidClassPanel_);
		lipidClassPanelSingle_ = new LipidClassPanel(this, true, faChainListNames_, lcbChainListNames_);
//		displayPanel_.add(lipidClassPanel_, getDefaultGridBagConstraints(1,0, GridBagConstraints.EAST, 1, 1));
		lipidClassPane_.add(lipidClassPanelSingle_, index);
		lipidClassPane_.setTitleAt(index, "Edit selected lipid (sub)class");
		lipidClassPane_.setSelectedIndex(index);
		lipidClassPane_.invalidate();
		lipidClassPane_.updateUI();
	}
	
	private JPanel getOutPathPanel(JTextField outPathField)
	{
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		JLabel exportIonModeLabel = new JLabel("Export relevant adducts for: ");
		exportIonModeLabel.setToolTipText(TooltipTexts.MASSLIST_GENERAL_ION_MODE);
		panel.add(exportIonModeLabel, getDefaultGridBagConstraints(0,0, GridBagConstraints.WEST, 1, 1));
		exportedIonMode_ = instantiateJComboBox(EXPORT_OPTIONS_ION_MODE, 0, false);
		exportedIonMode_.setToolTipText(TooltipTexts.MASSLIST_GENERAL_ION_MODE);
		panel.add(exportedIonMode_, getDefaultGridBagConstraints(1,0, GridBagConstraints.WEST, 1, 1));
		
		JLabel exportFormatLabel = new JLabel("Export file format: ");
		exportFormatLabel.setToolTipText(TooltipTexts.MASSLIST_GENERAL_ION_MODE);
		panel.add(exportFormatLabel, getDefaultGridBagConstraints(0,1, GridBagConstraints.WEST, 1, 1));
		exportedFormat_ = instantiateJComboBox(EXPORT_OPTIONS_FORMAT, 0, false);
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
	
	private JPanel initLipidClassPanel() throws IOException
	{
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		lipidClassPane_ = new JTabbedPane();
		panel.add(lipidClassPane_, new GridBagConstraints(0,0,1,1,1.0,1.0,GridBagConstraints.CENTER,GridBagConstraints.BOTH,new Insets(10,10,10,10),0,5));
		lipidClassPanelSingle_ = new LipidClassPanel(this, true, faChainListNames_, lcbChainListNames_);
		lipidClassPane_.add(lipidClassPanelSingle_, "Edit selected lipid (sub)class");
		lipidClassPanelMulti_ = new LipidClassPanel(this, false, faChainListNames_, lcbChainListNames_);
		lipidClassPane_.add(lipidClassPanelMulti_, "Edit multiple lipid (sub)classes");
		
		panel.add(instantiateJButton(COMMAND_EDIT_ADDUCT, "Edit list of adducts", true, TooltipTexts.MASSLIST_EDIT_ADDUCT), 
				getDefaultGridBagConstraints(0,1, GridBagConstraints.CENTER, 1, 1, new Insets(10, 10, 10, 10)));
		
		TitledBorder border = JOptionPanel.getTitledPanelBorder("Display / edit lipid class definitions");
		panel.setBorder(border);
		panel.setPreferredSize(new Dimension(565,LIPID_CLASS_HEIGHT));
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
	
	JComboBox<String> instantiateJComboBox(String[] entries, int index, boolean isSingle)
	{
		JComboBox<String> jComboBox = new JComboBox<String>(entries);
		jComboBox.setSelectedIndex(index);
		jComboBox.setPreferredSize(new Dimension(getPreferredDisplayComponentWidth(isSingle ? 0 : 50),20));
		return jComboBox;
	}
	
	int[] findSelectedAdductListIndices(LipidClassVO vo) throws IOException
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
	
	GridBagConstraints getDefaultGridBagConstraints(int column, int row, int orientation, int width, int height, Insets insets)
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
	
	GridBagConstraints getDefaultGridBagConstraints(int column, int row, int orientation, int width, int height)
	{
		return getDefaultGridBagConstraints(column, row, orientation, width, height, new Insets(2, 5, 2, 5));
	}
	
	private void jButtonExecuter(String actionCommand, boolean isOverride)
	{
		switch (actionCommand)
		{
			case COMMAND_EDIT_ADDUCT:
				adductCreatorPanel_.open();
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
	
	ArrayList<LipidClassVO> findLipidClassesToExport(String option)
	{
		if (option.equalsIgnoreCase(LipidClassPanel.OPTION_EXPORT_ALL))
			return allDefinedLipidClasses_;
		else
			return lipidClassTable_.getSelectedLipidClasses();
	}
	
	boolean isChainNumViable(LipidClassVO vo)
	{
		return vo.getNumberOfFAChains()>=0 && vo.getNumberOfLCBChains()>=0
				&& (vo.getNumberOfFAChains()+vo.getNumberOfLCBChains()>0);
	}
	
	/**
	 * Exports @param toExport with adapted numbers of @param faChains (FA) and @param sbpchains (SPB) chains
	 */
	void exportChainNumbersToSelectedClasses(int faChains, int spbchains, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setNumberOfFAChains(faChains);
			vo.setNumberOfLCBChains(spbchains);
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	/**
	 * Exports @param toExport with the adapted of @param faChainList (FA chain list)
	 */
	void exportFAChainListToSelectedClasses(String faChainList, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setFaChainList(faChainList);
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	/**
	 * Exports @param toExport with the adapted value of @param spbChainList (SPB chain list)
	 */
	void exportSPBChainListToSelectedClasses(String spbChainList, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setLCBChainList(spbChainList);
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	boolean isAdductListViable(LipidClassVO vo)
	{
		return !vo.getAdducts().isEmpty();
	}
	
	/**
	 * Exports @param toExport with the adapted value of @param adductList (the list of adducts)
	 */
	void exportAdductListToSelectedClasses(ArrayList<AdductVO> adductList, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setAdducts(adductList);
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	boolean isCNumViable(LipidClassVO vo)
	{
		return vo.getMinChainC() > 0 && vo.getMinChainC() < vo.getMaxChainC();
	}
	
	/**
	 * Exports @param toExport with the adapted value of @param min (minimum number of C atoms) and @param max (maximum number of C atoms)
	 */
	void exportCNumToSelectedClasses(int min, int max, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setMinChainC(min);
			vo.setMaxChainC(max);
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	boolean isDBNumViable(LipidClassVO vo)
	{
		return vo.getMinChainDB() >= 0 && vo.getMinChainDB() < vo.getMaxChainDB() 
				|| (vo.getMinChainDB() == 0 && vo.getMaxChainDB() == 0);
	}
	
	/**
	 * Exports @param toExport with the adapted value of @param min (minimum number of DB) and @param max (maximum number of DB)
	 */
	void exportDBNumToSelectedClasses(int min, int max, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setMinChainDB(min);
			vo.setMaxChainDB(max);
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	boolean isRTViable(LipidClassVO vo)
	{
		return (vo.getRtRangeFrom() >= 0 && vo.getRtRangeFrom() < vo.getRtRangeTo()) 
				|| (vo.getRtRangeFrom() < 0 && vo.getRtRangeTo() < 0);
	}
	
	/**
	 * Exports @param toExport with the adapted value of @param min (minimum RT) and @param max (maximum RT)
	 */
	void exportRTToSelectedClasses(double min, double max, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setRtRangeFrom(min);
			vo.setRtRangeTo(max);
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	boolean isLipidClassOxDefinitionViable(LipidClassVO vo)
	{
		return (vo.getOhRangeFrom() >= 0 && vo.getOhRangeFrom() <= vo.getOhRangeTo())
				|| (vo.getOhRangeFrom() == 0 && vo.getOhRangeFrom() == 0);
	}
	
	/**
	 * Exports @param toExport with the adapted value of @param min (minimum oxidation) and @param max (maximum oxidation)
	 */
	void exportOxNumToSelectedClasses(int min, int max, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setOhRangeFrom(min);
			vo.setOhRangeTo(max);
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	/**
	 * Exports @param toExport with the adapted value of @param filter (adduct insensitive retention time filter)
	 */
	void exportFilterToSelectedClasses(boolean filter, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setAdductInsensitiveRtFilter(filter);
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	/**
	 * Exports @param toExport with the adapted value of @param pick (pick best match by spectrum coverage)
	 */
	void exportPickToSelectedClasses(boolean pick, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setPickBestMatchBySpectrumCoverage(pick);
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	/**
	 * To be used when overriding a list of lipid class definitions
	 * @param vo
	 * @return
	 */
	boolean isGeneralDefinitionViable(LipidClassVO vo)
	{
		return isChainNumViable(vo)
				&& isCNumViable(vo)
				&& isDBNumViable(vo)
				&& isRTViable(vo)
				&& isLipidClassOxDefinitionViable(vo)
				&& isAdductListViable(vo);
	}
	
	/**
	 * Exports @param toExport with the adapted values taken from @param template
	 */
	void exportAll(LipidClassVO template, ArrayList<LipidClassVO> toExport)
	{
		for (LipidClassVO vo : toExport)
		{
			vo.setNumberOfFAChains(template.getNumberOfFAChains());
			vo.setNumberOfLCBChains(template.getNumberOfLCBChains());
			vo.setFaChainList(template.getFAChainList());
			vo.setLCBChainList(template.getLCBChainList());
			vo.setAdducts(template.getAdducts());
			vo.setMinChainC(template.getMinChainC());
			vo.setMaxChainC(template.getMaxChainC());
			vo.setMinChainDB(template.getMinChainDB());
			vo.setMaxChainDB(template.getMaxChainDB());
			vo.setRtRangeFrom(template.getRtRangeFrom());
			vo.setRtRangeTo(template.getRtRangeTo());
			vo.setOhRangeFrom(template.getOhRangeFrom());
			vo.setOhRangeTo(template.getOhRangeTo());
			vo.setAdductInsensitiveRtFilter(template.isAdductInsensitiveRtFilter());
			vo.setPickBestMatchBySpectrumCoverage(template.isPickBestMatchBySpectrumCoverage());
			exportLipidClass(vo, true);
		}
		try {reloadLipidClassScrollPane(selectedClass_.getLipidClass());}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	
	
	boolean isLipidClassViable()
	{
		return isLipidClassNameAvailable(tempClass_.getLipidClass())
				&& isChainNumViable(tempClass_)
				&& tempClass_.getHeadgroupFormula() != null
				&& isCNumViable(tempClass_)
				&& isDBNumViable(tempClass_)
				&& isRTViable(tempClass_)
				&& isLipidClassOxDefinitionViable(tempClass_)
				&& isAdductListViable(tempClass_);
	}
	
	void updateAdductForLipidClasses(AdductVO oldAdduct, AdductVO newAdduct, boolean isDelete)
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
	
	
	void synchronizeCurrentLipidClass(LipidClassPanel caller)
	{
		this.selectedClass_ = caller.getSelectedClass();
		this.tempClass_ = caller.getTempClass();;
	}
	
	/**
	 * Deletes a lipid class definition file
	 */
	void deleteLipidClass()
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
				reloadLipidClassScrollPane(null);
			}
			catch (IOException | ChemicalFormulaException ex) {
				new WarningMessage(new JFrame(), "Error", "An error occurred. Error message: "+ex.getMessage());
			}
		}
	}
	
	/**
	 * Exports a lipid class definition file. Called only for updating adduct lists.
	 * @param isOverride	true if the selected adduct definition file should be overriden / replaced with the new one
	 */
	void exportLipidClass(LipidClassVO toExport, boolean isOverride)
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
	void exportLipidClass(boolean isOverride)
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
			reloadLipidClassScrollPane(selectedClass_.getLipidClass());
		}
		catch (IOException | ChemicalFormulaException ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export. Error message: "+ex.getMessage());
		}
	}
	
	boolean isLipidClassNameAvailable(String name)
	{
		ArrayList<LipidClassVO> other = new ArrayList<LipidClassVO>(allDefinedLipidClasses_);
		other.remove(selectedClass_);
		for (LipidClassVO vo : other)
		{
			if (vo.getLipidClass().equalsIgnoreCase(name)) return false;
		}
		return true;
	}
	
	void setDefaultTextFieldBorder(JTextField textfield)
	{
		textfield.setBorder(BorderFactory.createLineBorder(Color.darkGray));
	}
	
	void setWarningTextFieldBorder(JTextField textfield)
	{
		textfield.setBorder(BorderFactory.createLineBorder(Color.red));
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
    
    @Override
    public Dimension getPreferredSize()
    {
    	return new Dimension(400,LIPID_CLASS_HEIGHT);
    }
    
    private LipidClassTable(String title, String[] lipidClassNames, String selected)
    {
    	this.lipidClassNames_ = lipidClassNames;
    	selectionTablePanel_ = generateSelectionTablePanel(initializeTableData(lipidClassNames), selected);
    	this.setLayout(new GridBagLayout());
    	this.add(selectionTablePanel_, new GridBagConstraints(0,0,1,1,1.0,1.0,GridBagConstraints.CENTER,GridBagConstraints.BOTH,new Insets(10,10,0,10),0,5));
    	this.setBorder(JOptionPanel.getTitledPanelBorder(title));
    }
    
    private JPanel generateSelectionTablePanel(Object[][] tableData, String selected)
    {
    	JPanel panel = new JPanel();
    	panel.setLayout(new GridBagLayout());
    	String[] columnNames = { "lipid class name", "select" };
    	model_ = new BooleanTableModel(tableData, columnNames);
    	displayTable_ = new JTable(model_);
    	displayTable_.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    	displayTable_.getColumnModel().getColumn(0).setPreferredWidth(250);
    	scrollPane_ = new JScrollPane(displayTable_);
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
    	panel.add(scrollPane_, new GridBagConstraints(0,0,1,1,1.0,1.0,GridBagConstraints.CENTER,GridBagConstraints.BOTH,new Insets(10,10,10,10),0,5));
    	JButton button = new JButton("Invert selection");
    	button.setToolTipText(TooltipTexts.MASSLIST_CLASS_SELECT_ALL);
  		button.addActionListener(new ActionListener() {
  		  public void actionPerformed(ActionEvent e) 
  		  {
  		  	invertSelection();
  		  }
  	  });
  		panel.add(button, getDefaultGridBagConstraints(0,1, GridBagConstraints.CENTER, 1, 1, new Insets(10, 10, 10, 10)));
    	return panel;
    }
    
    private void invertSelection()
    {
    	for (int i=0; i<lipidClassNames_.length; i++)
    	{
    		boolean value = (Boolean)model_.getValueAt(i, COLUMN_INCLUDE);
    		model_.setValueAt(!value, i, COLUMN_INCLUDE);
    	}
    	refreshLipidClassScrollPane();
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
	
	class BooleanTableModel extends AbstractTableModel {
		
		private static final long serialVersionUID = 1L;
		
		Object tableData_[][];
		String[] columnNames_;
		
		BooleanTableModel(Object[][] tableData, String[] columnNames)
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

	public AdductCreatorPanel getAdductCreatorPanel()
	{
		return this.adductCreatorPanel_;
	}

	public LipidClassVO getSelectedClass()
	{
		return selectedClass_;
	}
	
	public LipidClassVO getTempClass()
	{
		return tempClass_;
	}
	
	int getPreferredDisplayComponentWidth(int offset)
	{
		return 180-offset;
	}
	
	int getPreferredDisplayComponentWidthSmaller(int offset)
	{
		return 65-offset;
	}
}
