package at.tugraz.genome.lda.target.export;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.io.File;
import java.util.LinkedHashMap;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.verifier.DoubleVerifier;


public class ExportOptionsPanel extends JPanel
{
	private static final long serialVersionUID = 1L;
	
	private JComboBox<String> selectedGradient_;
	private JTextField thresholdClustering_;
	private final static String DEFAULT_GRADIENTS_DIR = "gradients";
	private final static String EXCEL_SUFFIX = ".xlsx";
	private final static Double DEFAULT_THRESHOLD = 5.0;
	
	/**
	 * Constructs a panel responsible for export options.
	 * @param isExperiment		true if options relevant for a C=C target list derived from experimental data should be displayed
	 */
	public ExportOptionsPanel(boolean isExperiment)
  {
  	this.setLayout(new GridBagLayout());
  	if (isExperiment)
  	{
  		selectedGradient_ = instantiateSelectedGradient();
  		thresholdClustering_ = instantiateThresholdClustering();
  		
  		this.add(new JLabel("Gradient definition: "), new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(25, 0, 0, 0), 0, 0));
  		this.add(selectedGradient_, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(25, 0, 0, 0), 0, 0));
  		this.add(new JLabel("Time tolerance: "), new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(25, 0, 0, 0), 0, 0));
  		this.add(thresholdClustering_, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(25, 0, 0, 0), 0, 0));
  		this.add(new JLabel(" sec"), new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(25, 0, 0, 0), 0, 0));
  	}
  }
	
	private JTextField instantiateThresholdClustering()
  {
  	JTextField field = new JTextField(DEFAULT_THRESHOLD.toString(), 5);
  	field.setInputVerifier(new DoubleVerifier(true,true));
  	field.setHorizontalAlignment(JTextField.RIGHT);
  	return field;
  }
  
  private JComboBox<String> instantiateSelectedGradient()
  {
  	JComboBox<String> jComboBox = new JComboBox<String>();
  	for (String gradient : getGradientsFiles().keySet())
  	{
  		jComboBox.addItem(gradient);
  	}
  	return jComboBox;
  }
  
  private LinkedHashMap<String,File> getGradientsFiles()
	{
    File gradientsDir = new File(DEFAULT_GRADIENTS_DIR);
    if (!gradientsDir.exists())
      return null;
    File[] files = gradientsDir.listFiles();
    LinkedHashMap<String,File> gradientFiles = new LinkedHashMap<String,File>();
    for (int i=0; i!=files.length;i++){
      String fileName = StaticUtils.extractFileName(files[i].getAbsolutePath());
      if (fileName.endsWith(EXCEL_SUFFIX))
      {
      	String gradientName = fileName.substring(0,fileName.indexOf(EXCEL_SUFFIX));
      	gradientFiles.put(gradientName, files[i]);
      }
    }
    return gradientFiles;
  }
  
  public File getSelectedGradient()
  {
  	LinkedHashMap<String,File> gradientsFiles = getGradientsFiles();
  	return gradientsFiles.get(selectedGradient_.getSelectedItem());
  }
  
  public Double getSelectedClustering()
  {
  	return Double.parseDouble(thresholdClustering_.getText());
  }
	
}
