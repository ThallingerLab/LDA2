package at.tugraz.genome.lda.target.experiment;
import javax.swing.JFrame;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.target.JDefaultComponents;


public class ExperimentStandardsFileChooserPanel extends ExperimentFileChooserPanel
{
	private static final long serialVersionUID = 1L;
	
	private ExperimentTableInputPanel experimentDefinitionPanel_;
	
  
  public ExperimentStandardsFileChooserPanel(JDefaultComponents wizardComponents, ExperimentTableInputPanel experimentDefinitionPanel) {
      super(wizardComponents, "Use stable isotope labels specific for \u03C9-C=C positions.", "Enter data of authentic standards here.");
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
