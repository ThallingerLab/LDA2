package at.tugraz.genome.lda.target.experiment;
import java.util.Vector;

import javax.swing.JFrame;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.target.JDefaultComponents;
import at.tugraz.genome.lda.utils.Pair;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.ResultFileVO;


public class ExperimentStandardsDefinitionPanel extends ExperimentTableInputPanel
{
	private static final long serialVersionUID = 1L;
	
	private static final int EDITABLE_COLUMN = 2;
	
	private Vector<Pair<ResultFileVO, Integer>> positions_;
  
  public ExperimentStandardsDefinitionPanel(JDefaultComponents wizardComponents, ExperimentFileChooserPanel standardsFileChooserPanel) {
      super(wizardComponents, standardsFileChooserPanel,
      		"Use stable isotope labels specific for \u03C9-C=C positions.", "Enter the \u03C9-C=C position associated with each authentic standard file.");
  }
  
  @Override
  public void initDataDisplay()
  {
  	cleanPanels();
  	positions_ = new Vector<Pair<ResultFileVO, Integer>>();
  	String[] columnNames = { "File Name", "Directory", "Enter \u03C9-C=C position here!"};
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
    	
    	positions_.add(new Pair<>(resultFileVO, 0));
    }
    return tableData;
  }
  
  @Override
  protected void updateValue(int row, int value)
	{
  	positions_.get(row).setValue(value);
	}
  
  private Vector<Pair<ResultFileVO, Integer>> getAssignedPositions(Vector<Pair<ResultFileVO, Integer>> positions)
  {
  	Vector<Pair<ResultFileVO, Integer>> assignedPositions = new Vector<Pair<ResultFileVO, Integer>>();
  	for (Pair<ResultFileVO, Integer> position : positions)
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
  	Vector<Pair<ResultFileVO, Integer>> assignedPositions = getAssignedPositions(positions_);
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
  		new WarningMessage(new JFrame(), "Warning", "You must specify the \u03C9-C=C position for at least one authentic standard file before continuing!");
  	}
  }
}
