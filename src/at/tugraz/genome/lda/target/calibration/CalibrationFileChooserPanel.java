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
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;

import javax.swing.JComboBox;
import javax.swing.DefaultCellEditor;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.filechooser.FileNameExtensionFilter;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.target.JDefaultComponents;
import at.tugraz.genome.lda.target.JOptionPanel;
import at.tugraz.genome.lda.target.JTargetFileWizard;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.verifier.DoubleVerifier;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class CalibrationFileChooserPanel extends JOptionPanel implements ActionListener
{
	private static final long serialVersionUID = 1L;
	
	private JTextField targetFileField_;
	private SelectionTable inputPanelOriginal_;
	private SelectionTable inputPanelNew_;
	private JTextField predictionThresholdField_;
	
	private Path previousSelection_ = null;
	
	private static final Double DEFAULT_THRESHOLD = 0.25;
	private static final String PLACEHOLDER_PREFIX = "Enter ";
	private static final String BROWSE = "Browse";
	private static final String COMMAND_OPEN_TARGET_FILE = "Open target file";
	public static final String DATA_TYPE_STANDARD_MIX = "standards";
	private static final Dimension ENTER_FIELD_DIMENSION_MIN = new Dimension(300,15);
	private static final Dimension ENTER_FIELD_DIMENSION = new Dimension(750,30);
	private static final Dimension BUTTON_PANEL_DIMENSION = new Dimension(825,35);
	private static final Dimension TABLE_PANEL_DIMENSION = new Dimension(825,200);
  
  
  public CalibrationFileChooserPanel(JDefaultComponents wizardComponents) {
      super(wizardComponents, "Calibrate a C=C target file to your chromatographic conditions.");
      init();
  }
  
  private void init() 
  {
  	this.setLayout(new GridBagLayout());
  	
  	targetFileField_ = instantiateJTextField(PLACEHOLDER_PREFIX + "path and file name of the original C=C target file.");
  	JButton targetFileButton = instantiateJButton(COMMAND_OPEN_TARGET_FILE, BROWSE);
  	JPanel targetFilePanel = instantiatePanel(targetFileField_, targetFileButton, "Original C=C target file");

  	inputPanelOriginal_ = new SelectionTable("Original chromatographic conditions.");
  	inputPanelNew_ = new SelectionTable("New chromatographic conditions.");
  	
  	this.add(targetFilePanel, 
  			new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        , GridBagConstraints.CENTER, GridBagConstraints.NONE
        , new Insets(5, 5, 5, 5), 0, 0));
  	this.add(inputPanelOriginal_,
        new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        , GridBagConstraints.CENTER, GridBagConstraints.BOTH
        , new Insets(5, 5, 5, 5), 0, 0));
  	this.add(inputPanelNew_,
        new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        , GridBagConstraints.CENTER, GridBagConstraints.BOTH
        , new Insets(5, 5, 5, 5), 0, 0));
  	this.add(instantiateThresholdPanel(), new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
        , GridBagConstraints.EAST, GridBagConstraints.NONE
        , new Insets(5, 5, 5, 5), 0, 0));
  }
  
  private JPanel instantiateThresholdPanel()
  {
  	JPanel panel = new JPanel();
  	panel.setLayout(new GridBagLayout());
  	predictionThresholdField_ = new JTextField(DEFAULT_THRESHOLD.toString(), 5);
  	predictionThresholdField_.setInputVerifier(new DoubleVerifier(true,true));
  	predictionThresholdField_.setHorizontalAlignment(JTextField.RIGHT);
  	panel.add(new JLabel("Maximum allowed deviation from standards-curve: "), new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
  	panel.add(predictionThresholdField_, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
  	panel.add(new JLabel(" min"), new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
  	return panel;
  }
  
  /**
	 * Creates a new JTextField with the provided placeholder and a focuslistener.
	 * @param placeholder
	 * @return
	 */
	private JTextField instantiateJTextField(String placeholder)
	{
		JTextField field = new JTextField();
		field.setMinimumSize(ENTER_FIELD_DIMENSION_MIN);
		field.setPreferredSize(ENTER_FIELD_DIMENSION);
		field.setText(placeholder);
		field.setForeground(Color.GRAY);
		field.addFocusListener(new FocusListener() {
		    @Override
		    public void focusGained(FocusEvent e) {
		        if (field.getText().equals(placeholder)) {
		        	  field.setText("");
		        	  field.setForeground(Color.BLACK);
		        }
		    }
		    @Override
		    public void focusLost(FocusEvent e) {
		        if (field.getText().isEmpty()) {
		        	  field.setForeground(Color.GRAY);
		        	  field.setText(placeholder);
		        }
		    }
		    });
		return field;
	}
	
	private JButton instantiateJButton(String command, String title)
	{
		JButton button = new JButton(title);
		button.addActionListener(this);
		button.setActionCommand(command);
		return button;
	}
	
	@Override
	public void actionPerformed(ActionEvent arg0)
	{
		if (arg0.getActionCommand().equals(COMMAND_OPEN_TARGET_FILE))
		{
			FileNameExtensionFilter filter = new FileNameExtensionFilter("Only .xlsx", "xlsx");
			selectPath(JFileChooser.FILES_ONLY, this.targetFileField_, filter, "Select the C=C target file (.xlsx file)");
		}
	} 
	
	/**
	 * Sets the text of the provided JTextField to the selected file or folder path depending on the selection mode.
	 * @param selectionMode
	 * @param field
	 * @param filter
	 */
	private void selectPath(int selectionMode, JTextField field, FileNameExtensionFilter filter, String title)
	{
		JFileChooser chooser = new JFileChooser();
		chooser.setFileFilter(filter);
		chooser.setPreferredSize(JTargetFileWizard.DEFAULT_FILE_CHOOSER_DIMENSION);
		chooser.setFileSelectionMode(selectionMode);
		chooser.setDialogTitle(title);
		if (previousSelection_ != null)
		{
			chooser.setCurrentDirectory(previousSelection_.getParent().toFile());
		}
		int val = chooser.showOpenDialog(new JFrame());
		if (val == JFileChooser.APPROVE_OPTION)
		{
			String text = chooser.getSelectedFile().getAbsolutePath();
			previousSelection_ = Paths.get(text);
			field.setText(text);
			field.setForeground(Color.BLACK);
		}
		else
		{
			return;
		}
	}
	
	private JPanel instantiatePanel(JTextField textfield, JButton button, String borderTitle)
	{
		JPanel panel = new JPanel();
  	panel.setLayout(new GridBagLayout());
    panel.add(textfield, getDefaultGridBagConstraints(0, 0));
  	panel.add(button, getDefaultGridBagConstraints(7, 0));
  	panel.setBorder(getTitledPanelBorder(borderTitle));
  	return panel;
	}
	
	private GridBagConstraints getDefaultGridBagConstraints(int column, int row)
	{
		return new GridBagConstraints(
				column, 
				row, 
				6, 
				1, 
				0.0, 
				0.0,
				GridBagConstraints.EAST, 
				GridBagConstraints.NONE, 
				new Insets(10, 6, 0, 0), 
				0, 
				10);
	}
  
  @Override
  protected void next() 
  {
  	if (isPlaceholder(targetFileField_) || inputPanelOriginal_.getUniqueFiles().isEmpty() || inputPanelNew_.getUniqueFiles().isEmpty())
  	{
  		new WarningMessage(new JFrame(), "Warning", "Specify paths and file names of the original C=C target file as well as for the reference files of the original and new chromatographic conditions before continuing!");
  	}
  	else
  	{
  		Hashtable<String,ArrayList<File>> originalConditions = inputPanelOriginal_.generateDataTypeToFilesMap();
  		Hashtable<String,ArrayList<File>> newConditions = inputPanelNew_.generateDataTypeToFilesMap();
  		
  		if (!originalConditions.keySet().equals(newConditions.keySet()))
  		{
  			new WarningMessage(new JFrame(), "Warning", "The reference files of the original and new chromatographic conditions must contain files of the same data types to ensure they can be matched!");
  		}
  		else 
  		{
  			if (!originalConditions.keySet().contains(DATA_TYPE_STANDARD_MIX))
    		{
  				//TODO: consider removing this warning or modifying it, currently a standard mix doesn't do all that much for data quality (gives more data points mostly).
    			new WarningMessage(new JFrame(), "Warning", "The reference files do not contain any files of the data type "+DATA_TYPE_STANDARD_MIX+"! The matched chromatographic peaks used for calibration will require more diligent manual curation.");
    		}
  			
  			goNext();
    		getDefaultComponents().disableAllButtons();
    		
    		Thread thread = new Thread(new Runnable()
    		{
    			public void run()
    			{
    				try 
  		  	  {
    					CalibrationGraphPanel panel = (CalibrationGraphPanel)getDefaultComponents().getCurrentPanel();
    					panel.setOriginalTargetList(new File(targetFileField_.getText()));
    					panel.setPredictionThreshold(Double.parseDouble(predictionThresholdField_.getText()));
  		  			panel.parseData(originalConditions, newConditions);
  		  			//during data parsing a redirect back to this panel might occur.
  		  			if (getDefaultComponents().getCurrentPanel() instanceof CalibrationGraphPanel)
  		  			{
  		  				panel.initDataDisplay();
  		  			}
  		  	  }
  		  		catch (Exception ex)
  		  		{
  		  			new WarningMessage(new JFrame(), "Error", "An error occurred: "+ex.getMessage());
  		  			goBack();
  		  		}
    			}
    		});
    		thread.start();
  		}
  	}
  }
  
  /**
   * Checks if the JTextfield has an entry different from the placeholder.
   * @param field
   * @return true if the field contains the placeholder
   */
  private boolean isPlaceholder(JTextField field)
  {
		if (field.getText().startsWith(PLACEHOLDER_PREFIX))
		{
			return true;
		}
		return false;
  }
  
  
  
  
  private class SelectionTable extends JPanel implements ActionListener
  {
		private static final long serialVersionUID = 1L;
		
		private ImageIcon addFilesIcon_ = new ImageIcon(getClass().getResource(
  			"/images/addFiles.gif"));
    private ImageIcon addFolderIcon_ = new ImageIcon(getClass().getResource(
    		"/images/addFromFolder.gif"));
    private ImageIcon removeIcon_ = new ImageIcon(getClass().getResource(
    		"/images/removeFiles.gif"));
  	
  	private static final String COMMAND_ADD_FILES = "Add Files";
  	private static final String COMMAND_ADD_FOLDER = "Add Folder";
  	private static final String COMMAND_REMOVE_FILE = "Remove";
  	private static final String COMMAND_REMOVE_ALL = "Remove all";
  	private static final int COLUMN_FILE_NAME = 0;
  	private static final int COLUMN_DIR_NAME = 1;
  	private static final int COLUMN_DATA_TYPE = 2;
  	
  	private ArrayList<File> uniqueFiles_ = new ArrayList<File>();
    private JPanel selectionTablePanel_;
    private JTable displayTable_;
    private JScrollPane scrollPane_;
    
    private SelectionTable(String title)
    {
    	JPanel buttonPanel = instantiateButtonPanel();
    	selectionTablePanel_ = new JPanel();
    	generateSelectionTablePanel(new String[0][0]);
    	this.setLayout(new GridBagLayout());
      this.add(buttonPanel, getDefaultGridBagConstraints(0, 0));
    	this.add(selectionTablePanel_, getDefaultGridBagConstraints(0, 1));
    	this.setBorder(getTitledPanelBorder(title));
    }
    
    private JPanel instantiateButtonPanel()
    {
    	JPanel buttonPanel = new JPanel();
    	buttonPanel.setLayout(new GridBagLayout()); 
      buttonPanel.setPreferredSize(BUTTON_PANEL_DIMENSION);
      buttonPanel.setMinimumSize(BUTTON_PANEL_DIMENSION);
      
      JButton buttonAddFiles = instantiateJButton(COMMAND_ADD_FILES, COMMAND_ADD_FILES, addFilesIcon_);
      buttonPanel.add(buttonAddFiles,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, -2, 2, 0), 0, 0));
      
      JButton buttonAddFolder = instantiateJButton(COMMAND_ADD_FOLDER, COMMAND_ADD_FOLDER, addFolderIcon_);
      buttonPanel.add(buttonAddFolder,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0));
      
      JButton buttonRemoveFile = instantiateJButton(COMMAND_REMOVE_FILE, COMMAND_REMOVE_FILE, removeIcon_);
      buttonPanel.add(buttonRemoveFile,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0)); 
      
      JButton buttonRemoveAll = instantiateJButton(COMMAND_REMOVE_ALL, COMMAND_REMOVE_ALL, removeIcon_);
      buttonPanel.add(buttonRemoveAll,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0));
      
      return buttonPanel;
    }
    
    private JButton instantiateJButton(String command, String title, ImageIcon imageIcon)
  	{
  		JButton button = new JButton(title, imageIcon);
  		button.addActionListener(this);
  		button.setActionCommand(command);
  		return button;
  	}
    
    private void generateSelectionTablePanel(String[][] tableData)
    {
    	String[] columnNames = { "file name", "directory", "data type" };
    	displayTable_ = new JTable(tableData, columnNames);
    	scrollPane_ = new JScrollPane(displayTable_);
    	displayTable_.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    	scrollPane_.setPreferredSize(TABLE_PANEL_DIMENSION);
    	scrollPane_.setMinimumSize(TABLE_PANEL_DIMENSION);
    	selectionTablePanel_.add(scrollPane_);
    }

		@Override
		public void actionPerformed(ActionEvent arg0)
		{
			FileNameExtensionFilter filter;
			switch (arg0.getActionCommand())
			{
				case COMMAND_ADD_FILES:
					filter = new FileNameExtensionFilter("Only .xlsx", "xlsx");
					selectPath(filter, "Select result files to use as reference for the original chromatographic conditions (.xlsx file)");
					break;
				case COMMAND_ADD_FOLDER:
					filter = new FileNameExtensionFilter("Only .xlsx", "xlsx");
					selectFolder("Select a folder containing result files to use as reference for the original chromatographic conditions (.xlsx file)");
					break;
				case COMMAND_REMOVE_FILE:
					removeFiles();
					break;
				case COMMAND_REMOVE_ALL:
					removeAllFiles();
					break;
				default:
					break;
			}
		}
		
		private void updateSelectionTable(){
	  	selectionTablePanel_.remove(scrollPane_);
	    String[][] tableData = new String[uniqueFiles_.size()][3];
	    int count=0;
	    for (File file : uniqueFiles_){
	      String fileName = StaticUtils.extractFileName(file.getAbsolutePath()); 
	      String dir = StaticUtils.extractDirName(file.getAbsolutePath()); 
	      tableData[count][COLUMN_FILE_NAME] = fileName;
	      tableData[count][COLUMN_DIR_NAME] = dir;
	      tableData[count][COLUMN_DATA_TYPE] = generateDataTypeJComboBox().getItemAt(0);
	      count++;
	    }
	    generateSelectionTablePanel(tableData);
	    displayTable_.getColumnModel().getColumn(COLUMN_DATA_TYPE).setCellEditor(new DefaultCellEditor(generateDataTypeJComboBox()));
	    selectionTablePanel_.invalidate();
	    selectionTablePanel_.updateUI();
	  }
		
		private JComboBox<String> generateDataTypeJComboBox()
		{
			JComboBox<String> comboBox = new JComboBox<String>();
			comboBox.addItem(DATA_TYPE_STANDARD_MIX);
			comboBox.addItem("group 1");
			comboBox.addItem("group 2");
			comboBox.addItem("group 3");
			comboBox.setSelectedIndex(0);
			return comboBox;
		}
		
		private void selectPath(FileNameExtensionFilter filter, String title)
		{
			JFileChooser chooser = new JFileChooser();
			chooser.setFileFilter(filter);
			chooser.setPreferredSize(JTargetFileWizard.DEFAULT_FILE_CHOOSER_DIMENSION);
			chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			chooser.setMultiSelectionEnabled(true);
			chooser.setDialogTitle(title);
			if (previousSelection_ != null)
			{
				chooser.setCurrentDirectory(previousSelection_.getParent().toFile());
			}
			int val = chooser.showOpenDialog(new JFrame());
			if (val == JFileChooser.APPROVE_OPTION)
			{
				Hashtable<String,File> avoidDuplicates = new Hashtable<String,File>();
				for (File file : uniqueFiles_) avoidDuplicates.put(file.getAbsolutePath(), file);
				File[] fileCandidates = chooser.getSelectedFiles();
				for (int i=0; i!=fileCandidates.length;i++){
					previousSelection_ = Paths.get(fileCandidates[i].getAbsolutePath());
	        if (fileCandidates[i].isFile() && !avoidDuplicates.containsKey(fileCandidates[i].getAbsolutePath())){
	          avoidDuplicates.put(fileCandidates[i].getAbsolutePath(),fileCandidates[i]);
	          this.uniqueFiles_.add(fileCandidates[i]);
	        }
	      }
				if (avoidDuplicates.size()>0){
	        updateSelectionTable();
	      }
			}
			else
			{
				return;
			}
		}
		
		private void selectFolder(String title)
		{
			JFileChooser chooser = new JFileChooser();
			chooser.setPreferredSize(JTargetFileWizard.DEFAULT_FILE_CHOOSER_DIMENSION);
			chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			chooser.setDialogTitle(title);
			if (previousSelection_ != null)
			{
				chooser.setCurrentDirectory(previousSelection_.getParent().toFile());
			}
			int val = chooser.showOpenDialog(new JFrame());
			if (val == JFileChooser.APPROVE_OPTION)
			{
				Hashtable<String,File> avoidDuplicates = new Hashtable<String,File>();
				for (File file : uniqueFiles_) avoidDuplicates.put(file.getAbsolutePath(), file);
				
				File folder = new File(chooser.getSelectedFile().getAbsolutePath());
				if (folder.exists() && folder.isDirectory())
				{
					File[] fileCandidates = folder.listFiles();
					for (int i=0; i!=fileCandidates.length;i++)
					{
						if (fileCandidates[i].isFile() && !avoidDuplicates.containsKey(fileCandidates[i].getAbsolutePath()))
						{
							previousSelection_ = Paths.get(fileCandidates[i].getAbsolutePath());
							String fileName = StaticUtils.extractFileName(fileCandidates[i].getAbsolutePath());
							String suffix = fileName.substring(fileName.lastIndexOf(".")+1);
							if (suffix.equalsIgnoreCase("xlsx"))
							{
								avoidDuplicates.put(fileCandidates[i].getAbsolutePath(),fileCandidates[i]);
								this.uniqueFiles_.add(fileCandidates[i]);
							}
						}
					}
				}
				if (avoidDuplicates.size()>0){
	        updateSelectionTable();
	      }
			}
			else
			{
				return;
			}
		}
		
		private void removeFiles()
		{
			if (uniqueFiles_!= null)
			{
				int[] selectedRows = displayTable_.getSelectedRows();
				if (selectedRows.length>0)
				{
					List<Integer> selectedList = new ArrayList<Integer>();
		      for (int i=0;i!=selectedRows.length;i++)
		      {
		        selectedList.add(selectedRows[i]);
		      }
		      Collections.sort(selectedList);
		      ArrayList<File> filesToRemove = new ArrayList<File>();
		      for (int i=(selectedList.size()-1);i!=-1;i--)
		      {
		        filesToRemove.add(this.uniqueFiles_.get((int)selectedList.get(i)));
		        this.uniqueFiles_.remove((int)selectedList.get(i));
		      }
		      updateSelectionTable();
				}
			}
		}
		
		private void removeAllFiles()
	  {
	  	uniqueFiles_ = new ArrayList<File>();
	  	updateSelectionTable();
	  }
		
		private ArrayList<File> getUniqueFiles()
		{
			return this.uniqueFiles_;
		}
		
		private Hashtable<String,ArrayList<File>> generateDataTypeToFilesMap()
		{
			Hashtable<String,ArrayList<File>> dataTypeToFileMap = new Hashtable<String,ArrayList<File>>();
			
			for (int i=0; i<displayTable_.getRowCount(); i++)
			{
				String fileNameRef = (String)displayTable_.getValueAt(i,COLUMN_FILE_NAME);
				String dirRef = (String)displayTable_.getValueAt(i,COLUMN_DIR_NAME);
				for (File file : uniqueFiles_)
				{
					String fileName = StaticUtils.extractFileName(file.getAbsolutePath()); 
		      String dir = StaticUtils.extractDirName(file.getAbsolutePath());
		      if (fileName.equals(fileNameRef) && dir.equals(dirRef))
		      {
		      	String dataType = (String)displayTable_.getValueAt(i,COLUMN_DATA_TYPE);
		      	if (!dataTypeToFileMap.containsKey(dataType))
		      	{
		      		dataTypeToFileMap.put(dataType, new ArrayList<File>());
		      	}
		      	dataTypeToFileMap.get(dataType).add(file);
		      }
				}
			}
			return dataTypeToFileMap;
		}
  }
}
