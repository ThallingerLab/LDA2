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

import javax.swing.JFrame;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.target.JDefaultComponents;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class ExperimentStandardsFileChooserPanel extends ExperimentFileChooserPanel
{
	private static final long serialVersionUID = 1L;
	
	private ExperimentTableInputPanel experimentDefinitionPanel_;
	
  
  public ExperimentStandardsFileChooserPanel(JDefaultComponents wizardComponents, ExperimentTableInputPanel experimentDefinitionPanel) {
      super(wizardComponents, "Use stable isotope labels specific for \u03C9-positions.", "Enter data of authentic standards here.");
      this.experimentDefinitionPanel_ = experimentDefinitionPanel;
  }
  
  @Override
  protected void next() 
  {
  	if (getFiles().isEmpty())
  	{
  		new WarningMessage(new JFrame(), "Warning", "Enter data data of authentic standards before continuing!");
  	}
  	else
  	{
  		try
  		{
  			if (!(getDefaultComponents().getNextPanel() instanceof ExperimentStandardsDefinitionPanel))
  			{
  				getDefaultComponents().addOptionPanelAfterCurrent(new ExperimentStandardsDefinitionPanel(getDefaultComponents(), this));
  			}
  		} catch (Exception ex) {}
  		
  		
  		goNext();
  		getDefaultComponents().disableAllButtons();
  		
  		Thread thread = new Thread(new Runnable()
  		{
  			public void run()
  			{
  				try 
		  	  {
  					ExperimentStandardsDefinitionPanel panel = (ExperimentStandardsDefinitionPanel)getDefaultComponents().getCurrentPanel();
		  			panel.loadData();
		  			panel.initDataDisplay();
		  	  }
		  		catch (Exception ex)
		  		{
		  			new WarningMessage(new JFrame(), "Error", "An error occurred: "+ex.getMessage());
		  		}
  			}
  		});
  		thread.start();
  	}
  }
  
  public ExperimentTableInputPanel getExperimentDefinitionPanel()
  {
  	return this.experimentDefinitionPanel_;
  }
}
