package at.tugraz.genome.lda.target.experiment;
import javax.swing.JFrame;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.target.JDefaultComponents;


public class ExperimentResultFileChooserPanel extends ExperimentFileChooserPanel
{
	private static final long serialVersionUID = 1L;
	
	
  
  public ExperimentResultFileChooserPanel(JDefaultComponents wizardComponents) {
      super(wizardComponents, "Use stable isotope labels specific for \u03C9-C=C positions.", "Enter data of stable isotope labeled experiments here.");
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
