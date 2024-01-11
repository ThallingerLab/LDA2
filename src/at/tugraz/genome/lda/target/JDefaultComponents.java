package at.tugraz.genome.lda.target;

import java.awt.CardLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.target.export.ExportPanel;

public class JDefaultComponents
{
	JButton cancelButton_;
	JButton backButton_;
	JButton nextButton_;
	JButton exportButton_;

  List<JOptionPanel> panelList_;
  int currentIndex_;
  JPanel componentsContainer_;
  PropertyChangeSupport propertyChangeListeners_;
  
  protected final static String CURRENT_PANEL_PROPERTY = "currentPanel";
  
  
  protected JDefaultComponents() 
  {
    try {
      init();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  
  private void init() throws Exception 
  {
    this.propertyChangeListeners_ = new PropertyChangeSupport(this);
    
	  cancelButton_ = new JButton();
	  backButton_ = new JButton();
	  nextButton_ = new JButton();
	  exportButton_ = new JButton();
	
	  panelList_ = new ArrayList<JOptionPanel>();
	  currentIndex_ = 0;
	  componentsContainer_ = new JPanel();
	  
	  cancelButton_.setText("Cancel");
	  cancelButton_.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	cancelButton_actionPerformed(e);
		  }
	  });
	  
	  backButton_.setText("< Back");
	  backButton_.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	backButton_actionPerformed(e);
		  }
	  });
	  
	  nextButton_.setText("Next >");
	  nextButton_.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	nextButton_actionPerformed(e);
		  }
	  });
	  
	  exportButton_.setText("Export");
	  exportButton_.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	exportButton_actionPerformed(e);
		  }
	  });
	
	  componentsContainer_.setLayout(new CardLayout());
  }
  
  
  void cancelButton_actionPerformed(ActionEvent e) 
  {
  	try 
  	{
  		if (getOptionPanelList().size()>1)
  		{
  			for (int i = getOptionPanelList().size()-1; i > 0; i--)
  			{
  				removeOptionPanel(i);
  			}
  		}
  		setCurrentIndex(0);
  		getCurrentPanel();
  		updateComponents();
	  } 
  	catch (Exception ex) 
  	{
	  	ex.printStackTrace();
	  }
  }
  
  void backButton_actionPerformed(ActionEvent e) 
  {
  	try 
  	{
  		getCurrentPanel().back();
	  } 
  	catch (Exception ex) 
  	{
	  	ex.printStackTrace();
	  }
  }
  
  void nextButton_actionPerformed(ActionEvent e)
  {
  	try 
	  {
		  getCurrentPanel().next();
	  } 
	  catch (Exception ex) 
	  {
	  	ex.printStackTrace();
	  }
  }
  
  void exportButton_actionPerformed(ActionEvent e)
  {
  	try
  	{
  		JOptionPanel panel = getCurrentPanel();
  		ExportPanel exporter = (ExportPanel) panel;
  		if (exporter.isPlaceholderExport())
    	{
    		new WarningMessage(new JFrame(), "Warning", "The path and file name to write the C=C target list to must be specified before the export!");
    	}
    	else
    	{
    		exporter.export();
    	}
  	} 
  	catch (Exception ex)
  	{
  		ex.printStackTrace();
  	}
  	
  }
  

  protected void addOptionPanel(JOptionPanel panel) 
  {
	  getOptionPanelList().add(panel);
	  componentsContainer_.add(panel,
	  getOptionPanelList().size() - 1 + "");
  }

  protected void addOptionPanel(int index, JOptionPanel panel) 
  {
	  getOptionPanelList().add(index, panel);
	  componentsContainer_.add(panel, index + "", index);
	  if (index < getOptionPanelList().size() - 1) 
	  {
		  for (int i = index + 1; i < getOptionPanelList().size(); i++) 
		  {
			  componentsContainer_.add(
			  (JOptionPanel)getOptionPanelList().get(i), i + "");
			}
	  }
  }

  protected void addOptionPanelAfter(JOptionPanel panelToBePlacedAfter, JOptionPanel panel) 
  {
	  addOptionPanel(
	  getOptionPanelList().indexOf(panelToBePlacedAfter) + 1, panel);
  }

  protected void addOptionPanelBefore(JOptionPanel panelToBePlacedBefore, JOptionPanel panel) 
  {
	  addOptionPanel(
	  getOptionPanelList().indexOf(panelToBePlacedBefore) - 1, panel);
  }

  public void addOptionPanelAfterCurrent(JOptionPanel panel) 
  {
  	addOptionPanel(getCurrentIndex()+1, panel);
  }

  //TODO: make sure this makes sense
  public JOptionPanel removeOptionPanel(JOptionPanel panel) 
  {
	  int index = getOptionPanelList().indexOf(panel);
	  getOptionPanelList().remove(panel);
	  componentsContainer_.remove(panel);
	  for (int i = index; i < getOptionPanelList().size(); i++) 
	  {
		  componentsContainer_.add(
		  (JOptionPanel) getOptionPanelList().get(i), i + "");
	  }
	  return panel;
  }

  // TODO: make sure this makes sense
  protected JOptionPanel removeOptionPanel(int index) 
  {
	  componentsContainer_.remove(index);
	  JOptionPanel panel = (JOptionPanel) getOptionPanelList().remove(index);
	  for (int i = index; i < getOptionPanelList().size(); i++) 
	  {
		  componentsContainer_.add(
		  (JOptionPanel) getOptionPanelList().get(i), i + "");
	  }
	  return panel;
  }

  protected JOptionPanel removeOptionPanelAfter(JOptionPanel panel) 
  {
  	return removeOptionPanel(getOptionPanelList().indexOf(panel) + 1);
  }

  protected JOptionPanel removeOptionPanelBefore(JOptionPanel panel) 
  {
  	return removeOptionPanel(getOptionPanelList().indexOf(panel) - 1);
  }

  protected JOptionPanel getOptionPanel(int index) 
  {
  	return (JOptionPanel) getOptionPanelList().get(index);
  }

  protected int getIndexOfPanel(JOptionPanel panel) 
  {
  	return getOptionPanelList().indexOf(panel);
  }

  protected boolean onLastPanel() 
  {
  	return (getCurrentIndex() == getOptionPanelList().size() - 1);
  }

  public JOptionPanel getCurrentPanel() throws Exception 
  {
  	if (getOptionPanelList().get(currentIndex_) != null) 
	  {
	  	return (JOptionPanel) getOptionPanelList().get(currentIndex_);
	  } 
	  else 
	  {
	  	throw new Exception("Requested panel not in panelList");
	  }
  }
  
  public JOptionPanel getNextPanel() throws Exception
  {
  	if (getOptionPanelList().get(currentIndex_+1) != null) 
	  {
	  	return (JOptionPanel) getOptionPanelList().get(currentIndex_+1);
	  } 
	  else 
	  {
	  	throw new Exception("Requested panel not in panelList");
	  }
  }

  public void updateComponents() 
  {
	  try 
	  {
		  CardLayout cl = (CardLayout) (componentsContainer_.getLayout());
		  cl.show(componentsContainer_, currentIndex_ + "");
		  
		  cancelButton_.setEnabled(true);
	  	backButton_.setEnabled(true);
	  	nextButton_.setEnabled(true);
	  	exportButton_.setEnabled(false);
		
		  if (currentIndex_ == 0)
		  {
		  	cancelButton_.setEnabled(false);
		  	backButton_.setEnabled(false);
		  	nextButton_.setEnabled(false);
		  }
		  else if (currentIndex_ == 1) //going back to the option panel should only work via the cancelButton
		  {
		  	backButton_.setEnabled(false);
		  }
		  else if (onLastPanel()) //the last panel will always have a higher index than 1
		  {
		  	nextButton_.setEnabled(false);
		  	exportButton_.setEnabled(true);
		  }
		  
		  getCurrentPanel().update();
		
		  // inform PropertyChangeListeners
		  PropertyChangeEvent event = new PropertyChangeEvent(this, CURRENT_PANEL_PROPERTY, null,  getCurrentPanel());
		  propertyChangeListeners_.firePropertyChange(event);
	  } 
	  catch (Exception e) 
	  {
	  	e.printStackTrace();
	  }
  }
  
  public void disableAllButtons()
  {
  	cancelButton_.setEnabled(false); //TODO: maybe cancelling should still be possible.
  	backButton_.setEnabled(false);
  	nextButton_.setEnabled(false);
  	exportButton_.setEnabled(false);
  }
  

  protected List<JOptionPanel> getOptionPanelList() 
  {
    return this.panelList_;
  }

  protected void setOptionPanelList(List<JOptionPanel> panelList) 
  {
    this.panelList_ = panelList;
  }

  protected int getCurrentIndex() 
  {
    return currentIndex_;
  }

  protected void setCurrentIndex(int aCurrentIndex) 
  {
    currentIndex_ = aCurrentIndex;
  }

  protected JPanel getComponentsContainer() 
  {
    return componentsContainer_;
  }

  protected void setComponentsContainer(JPanel aComponentsContainer) 
  {
    componentsContainer_ = aComponentsContainer;
  }
  
  protected JButton getCancelButton() 
  {
    return cancelButton_;
  }

  protected void setCancelButton(JButton aCancelButton) 
  {
    cancelButton_ = aCancelButton;
  }

  protected JButton getBackButton() 
  {
    return backButton_;
  }

  protected void setBackButton(JButton aBackButton) 
  {
    backButton_ = aBackButton;
  }

  protected JButton getNextButton() 
  {
    return nextButton_;
  }

  protected void setNextButton(JButton aNextButton) 
  {
    nextButton_ = aNextButton;
  }
  
  protected JButton getExportButton() 
  {
    return exportButton_;
  }

  protected void setExportButton(JButton aExportButton) 
  {
    exportButton_ = aExportButton;
  }

  protected void addPropertyChangeListener(PropertyChangeListener listener) 
  {
    propertyChangeListeners_.addPropertyChangeListener(listener);
  }
           
  protected void removePropertyChangeListener(PropertyChangeListener listener) 
  {
    propertyChangeListeners_.removePropertyChangeListener(listener);
  }
  
}
