package at.tugraz.genome.lda.target;

import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.SwingConstants;

public class JWizardFramework extends JPanel
{
  private static final long serialVersionUID = 1L;
  
  private JLabel panelTitleLabel_;
  private JDefaultComponents defaultComponents_;
  private JPanel buttonPanel_;
  
  public JWizardFramework() 
  {
    init();
  }
  
  private void init() 
  {
    defaultComponents_ = new JDefaultComponents();
    defaultComponents_.addPropertyChangeListener(new PropertyChangeListener() 
    {
      public void propertyChange(PropertyChangeEvent event) {
        setPanelTitle(((JOptionPanel)event.getNewValue()).getPanelTitle());
      }
    });
    
    this.setLayout(new GridBagLayout());
    this.add(createTitlePanel()
    , new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
    , GridBagConstraints.CENTER, GridBagConstraints.BOTH
    , new Insets(5, 5, 5, 5), 0, 0));
    
    this.add(new JSeparator()
    , new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
    , GridBagConstraints.WEST, GridBagConstraints.BOTH,
    new Insets(1, 1, 1, 1), 0, 0));
    
    this.add(defaultComponents_.getComponentsContainer()
    , new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
    , GridBagConstraints.CENTER, GridBagConstraints.BOTH
    , new Insets(0, 0, 0, 0), 0, 0));
    
    this.add(new JSeparator()
    , new GridBagConstraints(0, 3, 1, 1, 1.0, 0.0
    , GridBagConstraints.WEST, GridBagConstraints.BOTH
    , new Insets(1, 1, 1, 1), 0, 0));
    
    this.add(createButtonPanel(),
    new GridBagConstraints(0, 4, 1, 1, 1.0, 0.0
    ,GridBagConstraints.CENTER, GridBagConstraints.NONE,
    new Insets(0,0,0,0), 0, 25));
  }
  
  protected JDefaultComponents getDefaultComponents()
  {
    return defaultComponents_;
  }
  
  protected void setWizardComponents(JDefaultComponents aDefaultComponents)
  {
    defaultComponents_ = aDefaultComponents;
  }
  
  protected void updateComponents() 
  {
  	defaultComponents_.updateComponents();
  }
  
  
  protected void setPanelTitle(String title) 
  {
    panelTitleLabel_.setText(title);
  }
  
  protected JPanel createTitlePanel() 
  {
    JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT));
    panelTitleLabel_ = new JLabel();
    panelTitleLabel_.setFont(new Font("Arial",Font.BOLD, 24));
    panelTitleLabel_.setHorizontalAlignment(SwingConstants.LEADING);
    panel.add(panelTitleLabel_);
    return panel;
  }
  
  protected JPanel createButtonPanel() 
  {
  	buttonPanel_ = new JPanel(new GridBagLayout());
  	buttonPanel_.add(defaultComponents_.getCancelButton(), getButtonPanelConstraints(0, 0, 300));
  	buttonPanel_.add(defaultComponents_.getBackButton(), getButtonPanelConstraints(1, 5, 5));
  	buttonPanel_.add(defaultComponents_.getNextButton(), getButtonPanelConstraints(2, 5, 5));
  	buttonPanel_.add(defaultComponents_.getExportButton(), getButtonPanelConstraints(3, 300, 0));
    return buttonPanel_;
  }
  
  private GridBagConstraints getButtonPanelConstraints(int gridx, int left, int right)
  {
  	return new GridBagConstraints(gridx, 0, 1, 1, 0.0, 0.0,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, left, 0, right), 0, 0);
  }
}
