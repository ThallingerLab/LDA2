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

package at.tugraz.genome.lda.target.experiment;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.io.File;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.AbstractCellEditor;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.CellEditorListener;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellEditor;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.target.JDefaultComponents;
import at.tugraz.genome.lda.target.JOptionPanel;
import at.tugraz.genome.lda.target.LoadingPanel;
import at.tugraz.genome.lda.verifier.IntegerMaxVerifier;
import at.tugraz.genome.lda.vos.ResultFileVO;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public abstract class ExperimentTableInputPanel extends JOptionPanel
{
	private static final long serialVersionUID = 1L;
	private static final int MINIMUM_ALLOWED_VALUE = 0;
	private static final int MAXIMUM_ALLOWED_VALUE = 49;
	private static final Dimension TABLE_PANEL_DIMENSION = new Dimension(825,300);
	
	private ExperimentFileChooserPanel fileChooserPanel_;
	private String frameTitle_;
	private LoadingPanel loadingPanel_;
	private JPanel displayPanel_;
	private Vector<ResultFileVO> resultFileVO_;
  
	/**
	 * 
	 * @param wizardComponents
	 * @param fileChooserPanel
	 * @param title
	 * @param frameTitles
	 */
  public ExperimentTableInputPanel(JDefaultComponents wizardComponents, ExperimentFileChooserPanel fileChooserPanel, String title, String frameTitle) {
      super(wizardComponents, title);
      this.fileChooserPanel_ = fileChooserPanel;
      this.frameTitle_ = frameTitle;
      initLoadingPanel();
      this.displayPanel_ = null;
      this.resultFileVO_ = null;
      init(this.loadingPanel_);
  }
  
  private void initLoadingPanel()
  {
  	this.loadingPanel_ = new LoadingPanel("<html>Loading data, please wait...</html>");
  }
  
  protected void init(JPanel panel) 
  {
    this.add(panel, new GridBagConstraints(0, 0, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    this.invalidate();
    this.updateUI();
  }
  
  public abstract void initDataDisplay();
  
  protected abstract Object[][] generateTableData();
  
  protected void cleanPanels()
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
  }
  
  protected void reset()
  {
  	this.resultFileVO_ = null;
  	cleanPanels();
  	initLoadingPanel();
  	init(this.loadingPanel_);
  }
  
  /**
   * 
   * @param columnNames
   * @param tableData
   * @param editableColumn
   * @return
   */
  protected JPanel generateDisplayPanel(String[] columnNames, Object[][] tableData, int editableColumn)
  {
  	displayPanel_ = new JPanel();
  	displayPanel_.setLayout(new GridBagLayout());
    
    JTable standardsTable = new JTable(tableData, columnNames)
    { 
			private static final long serialVersionUID = 1L;
			
      @Override
      public boolean isCellEditable (int row, int column)
      {
        return column == editableColumn ? true : false;
      }
    };
    
    //center the editable column
    DefaultTableCellRenderer centerRenderer = new DefaultTableCellRenderer();
    centerRenderer.setHorizontalAlignment( JLabel.CENTER );
    standardsTable.getColumnModel().getColumn(editableColumn).setCellRenderer( centerRenderer );
    
    SpinnerEditor editor = new SpinnerEditor();
    
    //for value changes using the text field
    editor.addCellEditorListener(new CellEditorListener() {
			
			@Override
			public void editingStopped(ChangeEvent e)
			{
				SpinnerEditor source = (SpinnerEditor)e.getSource();
				int value = source.parseTextFieldContent();
				source.setTextFieldContent();
				standardsTable.setValueAt(value, source.getRow(), editableColumn);
				updateValue(source.getRow(), value);
			}

			@Override
			public void editingCanceled(ChangeEvent e)
			{
				System.out.println("Editing canceled!");
			}
			
		});
    standardsTable.getColumnModel().getColumn(editableColumn).setCellEditor(editor);
    standardsTable.setRowHeight(30);

    JScrollPane scrollPane = new JScrollPane(standardsTable);
    scrollPane.setPreferredSize(TABLE_PANEL_DIMENSION);
    
    
    displayPanel_.add(scrollPane, new GridBagConstraints(0, 0, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    displayPanel_.setBorder(getTitledPanelBorder(frameTitle_));
    return displayPanel_;
  }
  
  protected void loadData()
  {
  	Vector<File> files = fileChooserPanel_.getFiles();
		resultFileVO_ = new Vector<ResultFileVO>();
		try
		{
			for (File file : files)
			{
				QuantificationResult quantRes = LDAResultReader.readResultFile(file.getAbsolutePath(), new Hashtable<String,Boolean>());
				resultFileVO_.add(new ResultFileVO(file.getAbsolutePath(), file, quantRes));
			}
		}
		catch (ExcelInputFileException ex)
		{
			new WarningMessage(new JFrame(), "Warning", "An error occured while parsing the provided Excel files: "+ex.getMessage());
			goBack();
		}
  }
  
  protected ExperimentFileChooserPanel getFileChooserPanel()
  {
  	return this.fileChooserPanel_;
  }
  
  protected Vector<ResultFileVO> getResultFiles()
  {
  	return this.resultFileVO_;
  }
  
  protected abstract void updateValue(int row, int value);
  
  /**
	 * Inner helper class for the editable table column using JSpinner
	 */
	private class SpinnerEditor extends AbstractCellEditor implements TableCellEditor
	{
		private static final long serialVersionUID = 1L;
		final JSpinner spinner_ = new JSpinner();
		int row_ = -1;
		
		protected SpinnerEditor()
		{
			spinner_.setModel(new SpinnerNumberModel(0,MINIMUM_ALLOWED_VALUE,MAXIMUM_ALLOWED_VALUE,1));
			JFormattedTextField textField = (JFormattedTextField) spinner_.getEditor().getComponent(0);
			textField.setInputVerifier(new IntegerMaxVerifier(MINIMUM_ALLOWED_VALUE, MAXIMUM_ALLOWED_VALUE));
			
			//for value changes with the spinner arrows
			spinner_.addChangeListener(new ChangeListener() 
			{ 
	      public void stateChanged(ChangeEvent e) 
	      { 
	      	updateValue(row_, (Integer)spinner_.getValue());
	      	update();  
	      }
      }); 
		}
		
		@Override
		public Object getCellEditorValue()
		{
			return spinner_.getValue();
		}

		@Override
		public Component getTableCellEditorComponent(
				JTable table, Object value, boolean isSelected, int row, int column)
		{
			this.row_ = row;
			spinner_.setValue(value);
	    return spinner_;
		}
		
		protected int getRow() {
			return row_;
		}
		
		protected int parseTextFieldContent() 
		{
			JFormattedTextField textField = (JFormattedTextField) spinner_.getEditor().getComponent(0);
			if (isTextFieldContentViable(textField))
			{
				return Integer.parseInt(textField.getText());
			}
			else
			{
				return 0;
			}
			
		}
		
		protected void setTextFieldContent() 
		{
			JFormattedTextField textField = (JFormattedTextField) spinner_.getEditor().getComponent(0);
			int textFieldContent = parseTextFieldContent();
			textField.setText(Integer.toString(textFieldContent));
			spinner_.setValue(textFieldContent);
		}
		
		/**
		 * Checks if the text field content is viable; 
		 * Required because the InputVerifier can't catch the additional handling of the JTextField
		 * @return
		 */
		protected boolean isTextFieldContentViable(JFormattedTextField textField)
		{
			try
			{
        int number = Integer.parseInt(textField.getText());
        if (number>MAXIMUM_ALLOWED_VALUE || number<MINIMUM_ALLOWED_VALUE)
        {
        	return false;
        }
			} 
			catch (NumberFormatException ex) 
			{
				return false;
			}
			return true;
		}
	}
}
