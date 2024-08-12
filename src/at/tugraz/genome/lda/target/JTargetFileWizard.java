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

package at.tugraz.genome.lda.target;

import java.awt.Dimension;

import at.tugraz.genome.lda.target.calibration.CalibrationFileChooserPanel;
import at.tugraz.genome.lda.target.calibration.CalibrationGraphPanel;
import at.tugraz.genome.lda.target.experiment.ExperimentLabelDefinitionPanel;
import at.tugraz.genome.lda.target.experiment.ExperimentResultFileChooserPanel;
import at.tugraz.genome.lda.target.export.ExportPanel;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class JTargetFileWizard extends JWizardFramework
{
	private static final long serialVersionUID = 1L;
	
	protected static final int PANEL_CHOOSER = 0;
	protected static final int PANEL_CALIBRATION_FILE_CHOOSER = 1;
	protected static final int PANEL_CALIBRATION_GRAPH = 2;
	protected static final int PANEL_EXPERIMENT_FILE_CHOOSER = 1;
	protected static final int PANEL_EXPERIMENT_DEFINITION = 2;
	protected static final int PANEL_EXPORT = 3;
	
	public static final Dimension DEFAULT_FILE_CHOOSER_DIMENSION = new Dimension(750,750);
  

  public JTargetFileWizard()
  {
  	super();
    init();
  }
  
	private void init() 
	{
		JOptionPanel chooserPanel = new ChooserPanel(getDefaultComponents(), this);
		getDefaultComponents().addOptionPanel(PANEL_CHOOSER, chooserPanel);
    updateComponents();
  }
	
	public void initCalibration()
	{
		JOptionPanel calibrationFileChooserPanel = new CalibrationFileChooserPanel(getDefaultComponents());
		JOptionPanel calibrationGraphPanel = new CalibrationGraphPanel(getDefaultComponents());
		JOptionPanel exportPanel = new ExportPanel(getDefaultComponents());
		getDefaultComponents().addOptionPanel(PANEL_CALIBRATION_FILE_CHOOSER, calibrationFileChooserPanel);
		getDefaultComponents().addOptionPanel(PANEL_CALIBRATION_GRAPH, calibrationGraphPanel);
		getDefaultComponents().addOptionPanel(PANEL_EXPORT, exportPanel);
		getDefaultComponents().nextButton_actionPerformed(null);
    updateComponents();
	}
	
	public void initExperiment()
	{
		ExperimentResultFileChooserPanel experimentFileChooserPanel = new ExperimentResultFileChooserPanel(getDefaultComponents());
		JOptionPanel experimentDefinitionPanel = new ExperimentLabelDefinitionPanel(getDefaultComponents(), experimentFileChooserPanel);
		JOptionPanel exportPanel = new ExportPanel(getDefaultComponents());
		getDefaultComponents().addOptionPanel(PANEL_EXPERIMENT_FILE_CHOOSER, experimentFileChooserPanel);
		getDefaultComponents().addOptionPanel(PANEL_EXPERIMENT_DEFINITION, experimentDefinitionPanel);
		getDefaultComponents().addOptionPanel(PANEL_EXPORT, exportPanel);
		getDefaultComponents().nextButton_actionPerformed(null);
    updateComponents();
	}
	
}
