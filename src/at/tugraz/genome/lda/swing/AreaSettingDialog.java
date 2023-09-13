/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.verifier.DoubleVerifier;
import at.tugraz.genome.lda.vos.AreaSettingVO;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class AreaSettingDialog extends JDialog implements ActionListener
{

  private static final long serialVersionUID = 6343928984142509434L;

  private JPanel inputPanel_;
  
  private JTextField startTime_;
  private JTextField stopTime_;  
  private JTextField startMz_;
  private JTextField stopMz_;
  private JLabel startMzLabel_;
  private JLabel stopMzLabel_;

  private JCheckBox threeD_;
  private JButton buttonOK_;
  private JButton buttonCancel_;
  
  private ActionListener parent_;
  
  private final static String CHANGE_3D_SELECTION_STATUS = "change3DSelectionStatus";
  
  
  public AreaSettingDialog(JFrame parent, String title, String message, ActionListener parentListener){
    super(parent,title,true);
    this.parent_ = parentListener;
    setLayout(new BorderLayout());
    setLocation(380,240);
    JPanel messagePane = new JPanel();
    messagePane.add(new JLabel(message));
    add(messagePane,BorderLayout.NORTH);
    initInputPanel();
    initButtonPanel();
    this.setVisible(false);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    pack();
  }
  
  private void initInputPanel(){
    inputPanel_ = new JPanel();
    inputPanel_.setLayout(new GridBagLayout());
    
    JLabel label = new JLabel("start-time: ");
    label.setToolTipText("Test");
    label.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_START_TIME);
    inputPanel_.add(label,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    startTime_ = new JTextField(10);
    startTime_.setInputVerifier(new DoubleVerifier());
    startTime_.setHorizontalAlignment(JTextField.RIGHT);
    startTime_.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_START_TIME);
    inputPanel_.add(startTime_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    label = new JLabel("min");
    label.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_START_TIME);
    inputPanel_.add(label,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 1, 0, 1), 0, 0));
    label = new JLabel("stop-time: ");
    label.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_STOP_TIME);  
    inputPanel_.add(label,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    stopTime_ = new JTextField(10);
    stopTime_.setInputVerifier(new DoubleVerifier());
    stopTime_.setHorizontalAlignment(JTextField.RIGHT);
    stopTime_.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_STOP_TIME);  
    inputPanel_.add(stopTime_,new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    label = new JLabel("min");
    label.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_STOP_TIME);
    inputPanel_.add(label,new GridBagConstraints(6, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 1, 0, 1), 0, 0));
    
    startMzLabel_ = new JLabel("start-m/z: ");
    startMzLabel_.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_START_MZ);
    inputPanel_.add(startMzLabel_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    startMz_ = new JTextField(10);
    startMz_.setInputVerifier(new DoubleVerifier());
    startMz_.setHorizontalAlignment(JTextField.RIGHT);
    startMz_.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_START_MZ);
    inputPanel_.add(startMz_,new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    stopMzLabel_ = new JLabel("stop-m/z: ");
    stopMzLabel_.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_STOP_MZ);
    inputPanel_.add(stopMzLabel_,new GridBagConstraints(4, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    stopMz_ = new JTextField(10);
    stopMz_.setInputVerifier(new DoubleVerifier());
    stopMz_.setHorizontalAlignment(JTextField.RIGHT);
    stopMz_.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_STOP_MZ);
    inputPanel_.add(stopMz_,new GridBagConstraints(5, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    this.enableMzSettings(false);
    threeD_ = new JCheckBox("3D"); 
    threeD_.setActionCommand(CHANGE_3D_SELECTION_STATUS);
    threeD_.addActionListener(this);
    threeD_.setSelected(false);
    threeD_.setToolTipText(TooltipTexts.DIALOG_ADD_PEAK_3D);
    inputPanel_.add(threeD_,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    
    add(inputPanel_,BorderLayout.CENTER);
  }
  
  private void initButtonPanel(){
    JPanel buttonPanel = new JPanel();
    buttonPanel.setLayout(new GridBagLayout());
    buttonOK_ = new JButton("OK");
    buttonOK_.addActionListener(this);
    buttonOK_.setActionCommand("acceptAreaSettings");
    buttonOK_.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    buttonPanel.add(buttonOK_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 2, 0, 2), 0, 0));
    buttonCancel_ = new JButton("Cancel");
    buttonCancel_.addActionListener(parent_);
    buttonCancel_.setActionCommand("discardAreaSettings");
    buttonCancel_.setToolTipText(TooltipTexts.CANCEL_GENERAL);
    buttonPanel.add(buttonCancel_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 0, 2), 0, 0));
    add(buttonPanel,BorderLayout.SOUTH);
  }

  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase(CHANGE_3D_SELECTION_STATUS)){
      this.enableMzSettings(threeD_.isSelected());
    } else if (e.getActionCommand().equalsIgnoreCase("acceptAreaSettings")){
      String notFilledOutFields = "";
      if (startTime_.getText()==null || startTime_.getText().length()==0)
        notFilledOutFields+="start-time, ";
      if (stopTime_.getText()==null || stopTime_.getText().length()==0)
        notFilledOutFields+="stop-time, ";
      if (startMz_.getText()==null || startMz_.getText().length()==0)
        notFilledOutFields+="start-m/z, ";
      if (stopMz_.getText()==null || stopMz_.getText().length()==0)
        notFilledOutFields+="stop-m/z, ";
      if (notFilledOutFields.length()==0){
        parent_.actionPerformed(e);
      } else {
        notFilledOutFields = notFilledOutFields.substring(0,notFilledOutFields.length()-2);
        new WarningMessage(new JFrame(), "Error", "You have to fill out the following fields: "+notFilledOutFields+"!");
      }
    }
  }
  
  private void enableMzSettings(boolean enable){
    startMz_.setEnabled(enable);
    stopMz_.setEnabled(enable);
    startMzLabel_.setEnabled(enable);
    stopMzLabel_.setEnabled(enable);
  }
  
  public void setInputFields(boolean use3D, float startTime, float stopTime, float startMz, float stopMz){
    threeD_.setSelected(use3D);
    startTime_.setText(String.valueOf(Calculator.roundFloat(startTime, 3)));
    stopTime_.setText(String.valueOf(Calculator.roundFloat(stopTime, 3)));
    startMz_.setText(String.valueOf(Calculator.roundFloat(startMz, 5)));
    stopMz_.setText(String.valueOf(Calculator.roundFloat(stopMz, 5)));
    enableMzSettings(use3D);
  }
  
  public AreaSettingVO getInputValues(){
    return new AreaSettingVO(threeD_.isSelected(),Double.parseDouble(startTime_.getText()), Double.parseDouble(stopTime_.getText()),
        Double.parseDouble(startMz_.getText()), Double.parseDouble(stopMz_.getText()));
  }
  
}
