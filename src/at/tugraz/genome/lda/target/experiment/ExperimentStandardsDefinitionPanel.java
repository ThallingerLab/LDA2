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

import java.util.Vector;

import javax.swing.JFrame;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.target.JDefaultComponents;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.utils.UpdatablePair;
import at.tugraz.genome.lda.vos.ResultFileVO;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class ExperimentStandardsDefinitionPanel extends ExperimentTableInputPanel
{
	private static final long serialVersionUID = 1L;
	
	private static final int EDITABLE_COLUMN = 2;
	
	private Vector<UpdatablePair<ResultFileVO, Integer>> positions_;
  
  public ExperimentStandardsDefinitionPanel(JDefaultComponents wizardComponents, ExperimentFileChooserPanel standardsFileChooserPanel) {
      super(wizardComponents, standardsFileChooserPanel,
      		"Use stable isotope labels specific for \u03C9-positions.", "Enter the \u03C9-position associated with each authentic standard file.");
  }
  
  @Override
  public void initDataDisplay()
  {
  	cleanPanels();
  	positions_ = new Vector<UpdatablePair<ResultFileVO, Integer>>();
  	String[] columnNames = { "File Name", "Directory", "Enter \u03C9-position here!"};
  	Object[][] tableData = generateTableData();
  	getDefaultComponents().updateComponents();
  	init(generateDisplayPanel(columnNames, tableData, EDITABLE_COLUMN));
  }
  
  @Override
  protected Object[][] generateTableData()
  {
  	Object[][] tableData = new Object[getResultFiles().size()][3];
  	int count=0;
    for (ResultFileVO resultFileVO : getResultFiles())
    {
    	tableData[count][0] = StaticUtils.extractFileName(resultFileVO.getFileName());
    	tableData[count][1] = StaticUtils.extractDirName(resultFileVO.getFileName());
    	tableData[count][EDITABLE_COLUMN] = 0;
    	count++;
    	
    	positions_.add(new UpdatablePair<>(resultFileVO, 0));
    }
    return tableData;
  }
  
  @Override
  protected void updateValue(int row, int value)
	{
  	positions_.get(row).setValue(value);
	}
  
  private Vector<UpdatablePair<ResultFileVO, Integer>> getAssignedPositions(Vector<UpdatablePair<ResultFileVO, Integer>> positions)
  {
  	Vector<UpdatablePair<ResultFileVO, Integer>> assignedPositions = new Vector<UpdatablePair<ResultFileVO, Integer>>();
  	for (UpdatablePair<ResultFileVO, Integer> position : positions)
  	{
  		if (position.getValue() > 0)
  		{
  			assignedPositions.add(position);
  		}
  	}
  	return assignedPositions;
  }
  
  @Override
  protected void back() 
  {
  	try {getDefaultComponents().removeOptionPanel(getDefaultComponents().getCurrentPanel());} catch (Exception ex) {}
  	goBack();
  }
  
  @Override
  protected void next() 
  {
  	Vector<UpdatablePair<ResultFileVO, Integer>> assignedPositions = getAssignedPositions(positions_);
  	if (assignedPositions.size() > 0)
  	{
  		try
  		{
  			if (!(getDefaultComponents().getNextPanel() instanceof ExperimentGraphPanel))
  			{
  				getDefaultComponents().addOptionPanelAfterCurrent(new ExperimentGraphPanel(getDefaultComponents()));
  			}
  			
  			goNext();
    		getDefaultComponents().disableAllButtons();
    		
    		Thread thread = new Thread(new Runnable()
    		{
    			public void run()
    			{
    				try 
  		  	  {
    					ExperimentGraphPanel panel = (ExperimentGraphPanel) getDefaultComponents().getCurrentPanel();
    	  			ExperimentStandardsFileChooserPanel standardsFileChooserPanel = (ExperimentStandardsFileChooserPanel) getFileChooserPanel();
    	  			ExperimentLabelDefinitionPanel labelDefinitionPanel = (ExperimentLabelDefinitionPanel) standardsFileChooserPanel.getExperimentDefinitionPanel();
    	  			panel.matchIsotopologues(assignedPositions, labelDefinitionPanel);
    	  			panel.initDataDisplay();
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
    	catch (Exception ex) 
    	{
    		new WarningMessage(new JFrame(), "Error", "An error occurred: "+ex.getMessage());
    	}
  	}
  	else
  	{
  		new WarningMessage(new JFrame(), "Warning", "You must specify the \u03C9-position for at least one authentic standard file before continuing!");
  	}
  }
}
