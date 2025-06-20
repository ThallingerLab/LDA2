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

package at.tugraz.genome.lda.target.export;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.nio.file.Path;
import java.nio.file.Paths;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.filechooser.FileNameExtensionFilter;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.target.JDefaultComponents;
import at.tugraz.genome.lda.target.JOptionPanel;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class ExportPanel extends JOptionPanel implements ActionListener
{
	private static final long serialVersionUID = 1L;
	
	private TargetListExporter exporter_;
	private ExportOptionsPanel exportOptionsPanel_;
	private JTextField exportField_;
	private JTextField templateField_;
	private JPanel templatePanel_;
	
	private Path previousSelection_ = null;
	private static final String PLACEHOLDER_PREFIX = "Enter ";
	private static final String BROWSE = "Browse";
	private static final String COMMAND_OPEN_EXPORT= "Open export";
	private static final String COMMAND_OPEN_TEMPLATE= "Open template";
	private static final Dimension ENTER_FIELD_DIMENSION_MIN = new Dimension(300,15);
	private static final Dimension ENTER_FIELD_DIMENSION = new Dimension(775,30);
	private static final Dimension DEFAULT_FILE_CHOOSER_DIMENSION = new Dimension(750,750);
  
	
  public ExportPanel(JDefaultComponents wizardComponents) {
    super(wizardComponents, "RT-DB export");
    init();
  }
  
  private void init() 
  {
  	exportField_ = instantiateJTextField(PLACEHOLDER_PREFIX + "path and file name for the new RT-DB.");
  	JButton exportButton = instantiateJButton(COMMAND_OPEN_EXPORT, BROWSE);
  	JPanel exportPanel = instantiatePanel(exportField_, exportButton, "New RT-DB");
  	this.add(exportPanel, 
  			new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        , GridBagConstraints.CENTER, GridBagConstraints.NONE
        , new Insets(100, 0, 0, 0), 0, 0));
  }
  
  public void addNewRTDBCreationOptions(ExportOptionsPanel panel)
  {
  	addTemplatePanel();
  	addExportOptionsPanel(panel);
  }
  
  private void addTemplatePanel()
  {
  	templateField_ = instantiateJTextField(PLACEHOLDER_PREFIX + "path and file name for the template LDA mass list the new RT-DB should be based on.");
  	JButton templateButton = instantiateJButton(COMMAND_OPEN_TEMPLATE, BROWSE);
  	templatePanel_ = instantiatePanel(templateField_, templateButton, "Template LDA mass list");
  	this.add(templatePanel_, 
  			new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        , GridBagConstraints.CENTER, GridBagConstraints.NONE
        , new Insets(0, 0, 0, 0), 0, 0));
  }
  
  private void addExportOptionsPanel(ExportOptionsPanel panel)
  {
  	this.exportOptionsPanel_ = panel;
  	this.add(this.exportOptionsPanel_,
  			new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        , GridBagConstraints.WEST, GridBagConstraints.NONE
        , new Insets(50, 0, 0, 0), 0, 0));
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
		FileNameExtensionFilter filter;
		switch (arg0.getActionCommand())
		{
			case COMMAND_OPEN_TEMPLATE:
				filter = new FileNameExtensionFilter("Only .xlsx", "xlsx");
				selectPath(JFileChooser.FILES_ONLY, this.templateField_, filter, "Select the path and file name for the template LDA mass list (.xlsx file)");
				break;
			case COMMAND_OPEN_EXPORT:
				filter = new FileNameExtensionFilter("Only .xlsx", "xlsx");
				selectPath(JFileChooser.FILES_ONLY, this.exportField_, filter, "Select the path and file name for the new RT-DB (.xlsx file)");
				break;
			default:
				break;
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
		chooser.setPreferredSize(DEFAULT_FILE_CHOOSER_DIMENSION);
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
	
	/**
   * Checks if the template JTextfield has an entry different from the placeholder if it is not null.
   * @return true if the export field contains the placeholder
   */
  public boolean isPlaceholderTemplate()
  {
		if (templateField_ != null && templateField_.getText().startsWith(PLACEHOLDER_PREFIX))
		{
			return true;
		}
		return false;
  }
	
	/**
   * Checks if the export JTextfield has an entry different from the placeholder.
   * @return true if the export field contains the placeholder
   */
  public boolean isPlaceholderExport()
  {
		if (exportField_.getText().startsWith(PLACEHOLDER_PREFIX))
		{
			return true;
		}
		return false;
  }
	
	public void setExporter(TargetListExporter exporter)
	{
		this.exporter_ = exporter;
	}
	
	public String getOutPath()
	{
		return this.exportField_.getText();
	}
	
	public String getTemplatePath()
	{
		return this.templateField_.getText();
	}
	
	public ExportOptionsPanel getExportOptions()
	{
		return this.exportOptionsPanel_;
	}
	
	@Override
  protected void back() 
  {
		if (exportOptionsPanel_ != null)
		{
			this.remove(exportOptionsPanel_);
			exportOptionsPanel_ = null;
		}
		if (templatePanel_ != null)
		{
			this.remove(templatePanel_);
			templatePanel_ = null;
		}
		goBack();
  }
	
	public void export()
	{
		//TODO: write a fully new masslist without requiring the template.
//		String templatePath = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\masslists_all_labels\\negative.xlsx";
		try
		{
//			if (exporter_.getTemplatePath() != null) //TODO: if not writing a fully new target list: add field to allow to enter a template
//			{
//				templatePath = exporter_.getTemplatePath();
//			}
			exporter_.export(this);
			//TODO: this is to create figures to compare target lists; remove after
//			exporter_.exportBeforeAfter(
//					"D:\\Collaborator_Files\\SILDA\\SILDA_final\\Publication\\SI\\SupplementaryData\\Table03_RTDB_B1.xlsx",
//					"D:\\Collaborator_Files\\SILDA\\SILDA_final\\SILDA_30min_final\\Recalibrations\\comparison_RTDB_A_to_B1_USO.xlsx",
//					"D:\\Collaborator_Files\\SILDA\\SILDA_final\\Publication\\SI\\SupplementaryData\\Table06_RT_mapping_plus1837.xlsx"); 
		}
		catch (Exception ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export: "+ex.getMessage());
		}
	}
	
	
}
