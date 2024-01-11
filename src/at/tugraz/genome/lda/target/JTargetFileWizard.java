package at.tugraz.genome.lda.target;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import javax.swing.ButtonGroup;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.target.calibration.CalibrationFileChooserPanel;
import at.tugraz.genome.lda.target.calibration.CalibrationGraphPanel;
import at.tugraz.genome.lda.target.experiment.ExperimentLabelDefinitionPanel;
import at.tugraz.genome.lda.target.experiment.ExperimentResultFileChooserPanel;
import at.tugraz.genome.lda.target.export.ExportPanel;
import at.tugraz.genome.lda.target.experiment.ExperimentGraphPanel;
import at.tugraz.genome.lda.utils.Pair;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.IsotopeLabelVO;
import at.tugraz.genome.lda.vos.ResultFileVO;
import at.tugraz.genome.lda.exception.ExportException;

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
		JOptionPanel exportPanel = new ExportPanel(getDefaultComponents(), this);
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
		JOptionPanel exportPanel = new ExportPanel(getDefaultComponents(), this);
		getDefaultComponents().addOptionPanel(PANEL_EXPERIMENT_FILE_CHOOSER, experimentFileChooserPanel);
		getDefaultComponents().addOptionPanel(PANEL_EXPERIMENT_DEFINITION, experimentDefinitionPanel);
		getDefaultComponents().addOptionPanel(PANEL_EXPORT, exportPanel);
		getDefaultComponents().nextButton_actionPerformed(null);
    updateComponents();
	}
	
}
