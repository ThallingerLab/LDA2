package at.tugraz.genome.lda.target.export;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.nio.file.Path;
import java.nio.file.Paths;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;
import javax.swing.filechooser.FileNameExtensionFilter;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.target.JDefaultComponents;
import at.tugraz.genome.lda.target.JOptionPanel;
import at.tugraz.genome.lda.target.JTargetFileWizard;


public class ExportPanel extends JOptionPanel implements ActionListener
{
	private static final long serialVersionUID = 1L;
	private JTargetFileWizard jTargetFileWizard_;
	
	private TargetListExporter exporter_;
	
	private JTextField exportField_;
	
	private Path previousSelection_ = null;
	private static final String PLACEHOLDER_PREFIX = "Enter ";
	private static final String BROWSE = "Browse";
	private static final String COMMAND_OPEN_EXPORT= "Open export";
	private static final Dimension ENTER_FIELD_DIMENSION_MIN = new Dimension(300,15);
	private static final Dimension ENTER_FIELD_DIMENSION = new Dimension(775,30);
	private static final Dimension DEFAULT_FILE_CHOOSER_DIMENSION = new Dimension(750,750);
  
  public ExportPanel(JDefaultComponents wizardComponents, JTargetFileWizard jTargetFileWizard) {
      super(wizardComponents, "Export your new C=C target list.");
      this.jTargetFileWizard_ = jTargetFileWizard;
      init();
  }
  
  private void init() 
  {
  	exportField_ = instantiateJTextField(PLACEHOLDER_PREFIX + "path and file name for the new C=C target file.");
  	JButton exportButton = instantiateJButton(COMMAND_OPEN_EXPORT, BROWSE);
  	
  	JPanel exportPanel = instantiatePanel(exportField_, exportButton, "New C=C target file");
  	
  	this.add(exportPanel, 
  			new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        , GridBagConstraints.CENTER, GridBagConstraints.NONE
        , new Insets(100, 0, 0, 0), 0, 0));
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
			case COMMAND_OPEN_EXPORT:
				filter = new FileNameExtensionFilter("Only .xlsx", "xlsx");
				selectPath(JFileChooser.FILES_ONLY, this.exportField_, filter, "Select the path and file name for the new C=C target file (.xlsx file)");
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
	
	public void export()
	{
		//TODO: write a fully new masslist without requiring the template.
		String templatePath = "D:\\Collaborator_Files\\SILDA\\massLists\\masslists_double_labels\\negative.xlsx";
		try
		{
			if (exporter_.getTemplatePath() != null)
			{
				templatePath = exporter_.getTemplatePath();
			}
			exporter_.export(templatePath, exportField_.getText());
		}
		catch (Exception ex)
		{
			new WarningMessage(new JFrame(), "Error", "An error occurred during the export: "+ex.getMessage());
		}
	}
	
	
}
