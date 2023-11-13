package at.tugraz.genome.lda.target;

import javax.swing.BorderFactory;
import javax.swing.JPanel;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;


public class JOptionPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  
  protected JDefaultComponents defaultComponents_;
  private String panelTitle_;
  
  protected JOptionPanel(JDefaultComponents defaultComponents) 
  {
  	this.defaultComponents_ = defaultComponents;
    this.panelTitle_ = null;
  }
  
  protected JOptionPanel(JDefaultComponents defaultComponents, String title)
  {
  	this.defaultComponents_ = defaultComponents;
    this.panelTitle_ = title;
  }
  
  protected void update() 
  {
  }
  
  protected void next() 
  {
    goNext();
  }
  
  protected void back() 
  {
    goBack();
  }
  
  protected boolean goNext() 
  {
    if (defaultComponents_.getOptionPanelList().size() > defaultComponents_.getCurrentIndex()+1 ) 
    {
      defaultComponents_.setCurrentIndex(defaultComponents_.getCurrentIndex()+1);
      defaultComponents_.updateComponents();
      return true;
    } 
    else 
    {
      return false;
    }
  }
  
  protected boolean goBack() 
  {
    if (defaultComponents_.getCurrentIndex()-1 >= 0) 
    {
      defaultComponents_.setCurrentIndex(defaultComponents_.getCurrentIndex()-1);
      defaultComponents_.updateComponents();
      return true;
    } 
    else 
    {
      return false;
    }
  }
  
  protected JDefaultComponents getDefaultComponents()
  {
    return defaultComponents_;
  }
  
  protected void setDefaultComponents(JDefaultComponents aDefaultComponents)
  {
    defaultComponents_ = aDefaultComponents;
  }
  
  protected String getPanelTitle() 
  {
    return panelTitle_;
  }
  
  protected void setPanelTitle(String title) 
  {
    panelTitle_ = title;
  }
  
  protected void switchPanel(int panelIndex) 
  {
  	defaultComponents_.setCurrentIndex(panelIndex);
  	defaultComponents_.updateComponents();
  }
  
  public static TitledBorder getTitledPanelBorder(String title)
  {
  	Border raisedbevel = BorderFactory.createRaisedBevelBorder();
  	Border loweredbevel = BorderFactory.createLoweredBevelBorder();
  	Border compound = BorderFactory.createCompoundBorder(raisedbevel, loweredbevel);
  	return BorderFactory.createTitledBorder(compound, title);
  }
  
  
}
