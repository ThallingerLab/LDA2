package at.tugraz.genome.lda.target;

import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;


public class LoadingPanel extends JPanel
{
	private static final long serialVersionUID = 1L;
	
	private ImageIcon spinnerIcon_ = new ImageIcon(getClass().getResource(
			"/images/BeanEater.gif"));
  
  public LoadingPanel(String message) {
  	JLabel loading = new JLabel(message, spinnerIcon_, JLabel.CENTER);
    loading.setFont(new Font("Arial",Font.BOLD, 24));
    this.setLayout(new GridBagLayout());
    this.add(loading, new GridBagConstraints(0, 0, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
  }
}
