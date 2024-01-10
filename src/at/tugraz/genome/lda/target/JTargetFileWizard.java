package at.tugraz.genome.lda.target;

import java.awt.Dimension;

import at.tugraz.genome.lda.target.calibration.CalibrationFileChooserPanel;
import at.tugraz.genome.lda.target.calibration.CalibrationGraphPanel;
import at.tugraz.genome.lda.target.experiment.ExperimentLabelDefinitionPanel;
import at.tugraz.genome.lda.target.experiment.ExperimentResultFileChooserPanel;
import at.tugraz.genome.lda.target.export.ExportPanel;


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
