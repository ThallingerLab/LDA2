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
public class ExperimentResultFileChooserPanel extends ExperimentFileChooserPanel
{
	private static final long serialVersionUID = 1L;
	
	
  
  public ExperimentResultFileChooserPanel(JDefaultComponents wizardComponents) {
      super(wizardComponents, "Create RT-DB using experimental data of SIL specific for \u03C9-positions.", "Enter LDA result files of SIL experiments here.");
  }
  
  @Override
  protected void next() 
  {
  	if (getFiles().isEmpty())
  	{
  		new WarningMessage(new JFrame(), "Warning", "Enter data of stable isotope labeled experiments before continuing!");
  	}
  	else
  	{
  		goNext();
  		getDefaultComponents().disableAllButtons();
  		
  		Thread thread = new Thread(new Runnable()
  		{
  			public void run()
  			{
  				try 
		  	  {
  					ExperimentLabelDefinitionPanel panel = (ExperimentLabelDefinitionPanel)getDefaultComponents().getCurrentPanel();
		  			panel.loadData();
		  			panel.parseDataForLabels();
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
  }
}
