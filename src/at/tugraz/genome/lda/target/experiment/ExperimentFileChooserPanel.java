package at.tugraz.genome.lda.target.experiment;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.filechooser.FileNameExtensionFilter;
import at.tugraz.genome.lda.target.JDefaultComponents;
import at.tugraz.genome.lda.target.JOptionPanel;
import at.tugraz.genome.lda.target.JTargetFileWizard;
import at.tugraz.genome.lda.utils.StaticUtils;


public class ExperimentFileChooserPanel extends JOptionPanel implements ActionListener
{
	private static final long serialVersionUID = 1L;
	
	private Path previousSelection_ = null;
	
	private static final String COMMAND_ADD_FILES = "Add Files";
	private static final String COMMAND_ADD_FOLDER = "Add Folder";
	private static final String COMMAND_REMOVE_FILE = "Remove";
	private static final String COMMAND_REMOVE_ALL = "Remove all";
	
	private static final Dimension BUTTON_PANEL_DIMENSION = new Dimension(825,35);
	private static final Dimension TABLE_PANEL_DIMENSION = new Dimension(825,300);
	
	private ImageIcon addFilesIcon_ = new ImageIcon(getClass().getResource(
			"/images/addFiles.gif"));
  private ImageIcon addFolderIcon_ = new ImageIcon(getClass().getResource(
  		"/images/addFromFolder.gif"));
  private ImageIcon removeIcon_ = new ImageIcon(getClass().getResource(
  		"/images/removeFiles.gif"));
  
  private Vector<File> uniqueFiles_ = new Vector<File>();
  JPanel selectionTablePanel_;
  JScrollPane scrollPane_;
  JTable displayTable_;
	
	
  
  public ExperimentFileChooserPanel(JDefaultComponents wizardComponents, String title, String frameTitle) {
      super(wizardComponents, title);
      init(frameTitle);
  }
  
  private void init(String frameTitle) 
  {
  	JPanel buttonPanel = instantiateButtonPanel();
  	selectionTablePanel_ = new JPanel();
  	generateSelectionTablePanel(new String[0][0]);
  	
  	JPanel inputPanel = instantiatePanel(buttonPanel, selectionTablePanel_, frameTitle);
  	
    this.add(inputPanel, new GridBagConstraints(0, 0, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
  }
  
  private JPanel instantiatePanel(JPanel buttonPanel, JPanel selectionTablePanel, String borderTitle)
	{
		JPanel panel = new JPanel();
  	panel.setLayout(new GridBagLayout());
    panel.add(buttonPanel, getDefaultGridBagConstraints(0, 0));
  	panel.add(selectionTablePanel, getDefaultGridBagConstraints(0, 1));
  	panel.setBorder(getTitledPanelBorder(borderTitle));
  	return panel;
	}
  
  private JPanel instantiateButtonPanel()
  {
  	JPanel buttonPanel = new JPanel();
  	buttonPanel.setLayout(new GridBagLayout()); 
    buttonPanel.setPreferredSize(BUTTON_PANEL_DIMENSION);
    
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
  
  private void generateSelectionTablePanel(String[][] tableData)
  {
  	String[] columnNames = { "File Name", "Directory" };
  	displayTable_ = new JTable(tableData, columnNames);
  	scrollPane_ = new JScrollPane(displayTable_);
  	displayTable_.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
  	scrollPane_.setPreferredSize(TABLE_PANEL_DIMENSION);
  	selectionTablePanel_.add(scrollPane_);
  }
  
  private void updateSelectionTable(){
  	selectionTablePanel_.remove(scrollPane_);
    String[][] tableData = new String[uniqueFiles_.size()][2];
    int count=0;
    for (File file : uniqueFiles_){
      String fileName = StaticUtils.extractFileName(file.getAbsolutePath()); 
      String dir = StaticUtils.extractDirName(file.getAbsolutePath()); 
      tableData[count][0] = fileName;
      tableData[count][1] = dir;
      count++;
    }
    generateSelectionTablePanel(tableData);
    selectionTablePanel_.invalidate();
    selectionTablePanel_.updateUI();
  }
	
	private JButton instantiateJButton(String command, String title, ImageIcon imageIcon)
	{
		JButton button = new JButton(title, imageIcon);
		button.addActionListener(this);
		button.setActionCommand(command);
		return button;
	}
	
	@Override
	public void actionPerformed(ActionEvent arg0)
	{
		switch (arg0.getActionCommand())
		{
			case COMMAND_ADD_FILES:
				FileNameExtensionFilter filter = new FileNameExtensionFilter("Only .xlsx", "xlsx");
				selectPath(filter, "Select a LDA result file (.xlsx file)");
				break;
			case COMMAND_ADD_FOLDER:
				selectFolder("Select a folder containing LDA result files (.xlsx files)");
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
	      Vector<File> filesToRemove = new Vector<File>();
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
  	uniqueFiles_ = new Vector<File>();
  	updateSelectionTable();
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
  
  protected Vector<File> getFiles()
  {
  	return uniqueFiles_;
  }
}
