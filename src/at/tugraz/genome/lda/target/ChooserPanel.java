/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2023 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */

package at.tugraz.genome.lda.target;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.Border;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class ChooserPanel extends JOptionPanel
{
	private static final long serialVersionUID = 1L;
	private JTargetFileWizard jTargetFileWizard_;
  
  public ChooserPanel(JDefaultComponents wizardComponents, JTargetFileWizard jTargetFileWizard) {
      super(wizardComponents, ""
      		+ "Create a new RT-DB with information for \u03C9-position identification.");
      this.jTargetFileWizard_ = jTargetFileWizard;
      init();
  }
  
  private void init() 
  {
  	Dimension preferredSizeTextOneLine = new Dimension(825,15);
  	
  	this.setLayout(new GridBagLayout());
  	
  	JPanel calibrationPanel = initCalibrationPanel(preferredSizeTextOneLine);
  	JPanel experimentsPanel = initExperimentsPanel(preferredSizeTextOneLine);
  	
  	calibrationPanel.setBorder(getPanelBorder());
  	experimentsPanel.setBorder(getPanelBorder());
  	
  	this.add(calibrationPanel
        , new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
        , GridBagConstraints.NORTH, GridBagConstraints.BOTH
        , new Insets(25, 25, 25, 25), 0, 0));
  	this.add(experimentsPanel,
        new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
        , GridBagConstraints.SOUTH, GridBagConstraints.BOTH
        , new Insets(25, 25, 25, 25), 0, 0));
  }
  
  private JPanel initCalibrationPanel(Dimension sizeTextOneLine)
  {
  	JPanel calibrationPanel = new JPanel();
  	calibrationPanel.setLayout(new GridBagLayout());
    JLabel title = new JLabel("<html>RT-DB predictor - map an RT-DB to your chromatographic conditions.</html>");
    title.setFont(new Font("Arial",Font.BOLD, 24));
    title.setPreferredSize(new Dimension(sizeTextOneLine.width,sizeTextOneLine.height*4));
    title.setMinimumSize(new Dimension(title.getPreferredSize().width/2, title.getPreferredSize().height));
    JLabel requirementsHeader = new JLabel("<html> Requirements:</html>");
    requirementsHeader.setFont(new Font("Arial",Font.BOLD, 18));
    requirementsHeader.setPreferredSize(new Dimension(sizeTextOneLine.width,sizeTextOneLine.height*2));
    requirementsHeader.setMinimumSize(new Dimension(requirementsHeader.getPreferredSize().width/2, requirementsHeader.getPreferredSize().height));
    JLabel requirements= new JLabel("<html> <ul style=\"list-style-type:circle;\"> "
    		+ "<li> An RT-DB in LDA-format. </li>"
    		+ "<li> Experimental data (analyzed with LDA) of an (internal) standard mix and/or biological samples measured with the exact same chromatographic conditions the RT-DB was created with. </li>"
    		+ "<li> Experimental data (analyzed with LDA) of the same standard mix and/or biological samples measured with your chromatography. </li>"
    		+ "</html>");
    requirements.setFont(new Font("Arial",Font.PLAIN, 12));
    requirements.setPreferredSize(new Dimension(sizeTextOneLine.width,sizeTextOneLine.height*5));
    requirements.setMinimumSize(new Dimension(requirements.getPreferredSize().width/2, requirements.getPreferredSize().height*2));
    
    JButton calibrateButton = new JButton();
  	calibrateButton.setText("Continue");
  	calibrateButton.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	calibrateButton_actionPerformed(e);
		  }
	  });
  	
  	calibrationPanel.add(title,getDefaultConstraints(0,0,GridBagConstraints.WEST,5,0));
    calibrationPanel.add(requirementsHeader,getDefaultConstraints(0,1,GridBagConstraints.WEST,5,0));
    calibrationPanel.add(requirements,getDefaultConstraints(0,2,GridBagConstraints.WEST,5,0));
    calibrationPanel.add(calibrateButton, getDefaultConstraints(0,3,GridBagConstraints.CENTER,0,10));
    
    return calibrationPanel;
  }
  
  private JPanel initExperimentsPanel(Dimension sizeTextOneLine)
  {
  	JPanel experimentsPanel = new JPanel();
  	experimentsPanel.setLayout(new GridBagLayout());
    JLabel title = new JLabel("<html>RT-DB generator - create a new RT-DB.</html>");
    title.setFont(new Font("Arial",Font.BOLD, 24));
    title.setPreferredSize(new Dimension(sizeTextOneLine.width,sizeTextOneLine.height*4));
    title.setMinimumSize(new Dimension(title.getPreferredSize().width/2, title.getPreferredSize().height));
    JLabel requirementsHeader = new JLabel("<html> Requirements:</html>");
    requirementsHeader.setFont(new Font("Arial",Font.BOLD, 18));
    requirementsHeader.setPreferredSize(new Dimension(sizeTextOneLine.width,sizeTextOneLine.height*2));
    requirementsHeader.setMinimumSize(new Dimension(requirementsHeader.getPreferredSize().width/2, requirementsHeader.getPreferredSize().height));
    JLabel requirements= new JLabel("<html> <ul style=\"list-style-type:circle;\"> "
    		+ "<li> Experimental data (analyzed with LDA) of analytes with SIL specific for \u03C9-positions. </li>"
    		+ "<li> If deuterium labels are used: experimental data (analyzed with LDA) of unlabeled authentic standards with known \u03C9-positions. <br>"
    		+ "Ideally at least one standard for each deuterium label should be provided, with the respective labeled isotopologue identified in the SIL data. <br>"
    		+ "(For <sup>13</sup>C labels experimental data of unlabeled authentic standards is not needed.) </li>"
    		+ "</ul> </html>");
    requirements.setFont(new Font("Arial", Font.PLAIN, 12));
    requirements.setPreferredSize(new Dimension(sizeTextOneLine.width,sizeTextOneLine.height*5));
    requirements.setMinimumSize(new Dimension(requirements.getPreferredSize().width/2, requirements.getPreferredSize().height*2));
    
    JButton experimentsButton = new JButton();
  	experimentsButton.setText("Continue");
  	experimentsButton.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	experimentsButton_actionPerformed(e);
		  }
	  });
  	
  	experimentsPanel.add(title,getDefaultConstraints(0,0,GridBagConstraints.WEST,5,0));
    experimentsPanel.add(requirementsHeader,getDefaultConstraints(0,1,GridBagConstraints.WEST,5,0));
    experimentsPanel.add(requirements,getDefaultConstraints(0,2,GridBagConstraints.WEST,5,0));
    experimentsPanel.add(experimentsButton, getDefaultConstraints(0,3,GridBagConstraints.CENTER,0,10));
    
    return experimentsPanel;
  }
  
  private GridBagConstraints getDefaultConstraints(int gridx, int gridy, int anchor, int left, int ipady)
  {
  	return new GridBagConstraints(gridx, gridy, 1, 1, 0.0, 0.0, anchor, 
  			GridBagConstraints.NONE, new Insets(5, left, 5, 0), 0, ipady);
  }
  
  private Border getPanelBorder()
  {
  	Border raisedbevel = BorderFactory.createRaisedBevelBorder();
  	Border loweredbevel = BorderFactory.createLoweredBevelBorder();
  	return BorderFactory.createCompoundBorder(raisedbevel, loweredbevel);
  }
  
  void calibrateButton_actionPerformed(ActionEvent e)
  {
  	jTargetFileWizard_.initCalibration();
  }
  
  void experimentsButton_actionPerformed(ActionEvent e)
  {
  	jTargetFileWizard_.initExperiment();
  }
}
