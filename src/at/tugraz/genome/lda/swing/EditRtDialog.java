/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler 
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

package at.tugraz.genome.lda.swing;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.verifier.DoubleVerifier;

/**
 * displays dialog for editing the retention time
 * 
 * @author Juergen Hartler
 *
 */
public class EditRtDialog extends JDialog implements ActionListener
{

  private static final long serialVersionUID = -7614750501692196104L;
  
  /** the parent action listener*/
  private ActionListener parent_;
  /** the input field for the retention time*/
  private JTextField rTime_;
  
  /**
   * Constructor for a dialog field for changing the retention time of a hit
   * @param rt the retention time
   * @param parent the parent listener that is informed when a button is pressed
   */
  public EditRtDialog(String rt, ActionListener parent){

    this.parent_ = parent;
    this.setLayout(new BorderLayout());
    setLocation(380,240);
    initMainPanel(rt);
    initButtonPanel();
    pack(); 
    setVisible(true);
  }
  
  /**
   * initializes the main panel that contains the input field for the retention time
   * @param rt the old retention time
   */
  private void initMainPanel(String rt){
    JPanel centerPanel = new JPanel();
    centerPanel.setLayout(new GridBagLayout());
    JLabel label = new JLabel("RT: ");
    label.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_RT);
    centerPanel.add(label,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    rTime_ = new JTextField(5);
    rTime_.setHorizontalAlignment(JTextField.RIGHT);
    rTime_.setInputVerifier(new DoubleVerifier(true,true));
    rTime_.setText(rt);
    rTime_.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_RT);
    centerPanel.add(rTime_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    this.add(centerPanel,BorderLayout.CENTER);
  }
  
  /**
   * initializes the panel containing the buttons "OK" and "Cancel"
   */
  private void initButtonPanel(){
    JPanel buttonPanel = new JPanel();
    JButton okButton = new JButton("OK");
    okButton.setActionCommand("AcceptChangedRT");
    okButton.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    buttonPanel.add(okButton);
    okButton.addActionListener(this);
    JButton cancelButton = new JButton("Cancel");
    cancelButton.setActionCommand("DeclineChangedRT");
    cancelButton.setToolTipText(TooltipTexts.CANCEL_GENERAL);
    buttonPanel.add(cancelButton);
    cancelButton.addActionListener(this);
    cancelButton.addActionListener(parent_);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    this.add(buttonPanel,BorderLayout.SOUTH);
    
  }
  
  
  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase("AcceptChangedRT")){
      try{
        Double.parseDouble(rTime_.getText());
        setVisible(false); 
        parent_.actionPerformed(e);
        dispose();      
      } catch (NumberFormatException nfx){}
    }else if (e.getActionCommand().equalsIgnoreCase("DeclineChangedRT")){
      setVisible(false); 
      dispose();      
    }
  }

  /**
   * 
   * @return the entered retention time
   */
  public Double getRt()
  {
    return Double.parseDouble(rTime_.getText());
  }
  
  

}
