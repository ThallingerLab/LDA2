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

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import at.tugraz.genome.lda.TooltipTexts;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class CheckBoxOptionPane extends JDialog implements ActionListener
{

  private static final long serialVersionUID = -5766256473109433281L;
  
  private final int NR_OF_ROWS = 5;
  private int value_ = JOptionPane.NO_OPTION;
  private ArrayList<String> checkBoxStrings_;
  private Hashtable<String,JCheckBox> checkBoxHash_;
  private boolean showOne_;

  /**
   * generic input pane that contains a list of check boxes
   * @param parent the parent JFrame
   * @param title the title of the frame
   * @param messages the messages displayed to the user (each row is a separate line)
   * @param chechboxStrings a list of strings displayed next to the check boxes (each string is an individual check box)
   * @param showOne is false when a single check box should not be displayed
   */
  private CheckBoxOptionPane(JFrame parent, String title, Vector<String> messages, ArrayList<String> chechboxStrings, boolean showOne){
    super(parent, title, true);
    checkBoxStrings_ = chechboxStrings;
    showOne_ = showOne;
    checkBoxHash_ = new Hashtable<String,JCheckBox>();
    Vector<JCheckBox> chechboxes = new Vector<JCheckBox>();
    if ((showOne&&checkBoxStrings_.size()>0)||checkBoxStrings_.size()>1){
      for (String name : checkBoxStrings_){
        JCheckBox box = new JCheckBox(name);
        box.setSelected(true);
        checkBoxHash_.put(name, box);
        chechboxes.add(box);
      }
    }  
    this.setLayout(new GridBagLayout());
    setLocation(380,240);
    value_ = JOptionPane.NO_OPTION;
    int columns = chechboxes.size()/NR_OF_ROWS;
    if (chechboxes.size()%NR_OF_ROWS!=0 || columns==0)columns++;
    int rowCount = 1;
    JPanel textPanel = new JPanel();
    textPanel.setLayout(new GridLayout());
    for (String message : messages) //edit here!!
    {
      if (message!=null&&message.length()>0)
      {
        JLabel titleLabel = new JLabel(message);
        textPanel.add(titleLabel);
      }
    }
    JScrollPane scrollPane = new JScrollPane(textPanel);
    scrollPane.setPreferredSize(new Dimension(950, 250));
    this.add(scrollPane, new GridBagConstraints(0, 0, columns, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
  
//  int rows = moleculeNames_.size()/columns_;
//  if (moleculeNames_.size()>0 && moleculeNames_.size()%columns_!=0)
//    rows++;
    for (int i=0; i!=chechboxes.size(); i++){
      JCheckBox select  = chechboxes.get(i);
      select.setSelected(true);
      int row = i%NR_OF_ROWS+rowCount;
      int column = i/NR_OF_ROWS;
      this.add(select, new GridBagConstraints(column, row, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    }
    int requiredRows = rowCount+chechboxes.size();
    if (requiredRows>1+NR_OF_ROWS)
      requiredRows = 1+NR_OF_ROWS;
    JPanel buttonPanel = new JPanel();
    this.add(buttonPanel, new GridBagConstraints(0, requiredRows+1, columns, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    
    JButton okButton = new JButton("OK");
    okButton.setActionCommand("AcceptSelection");
    okButton.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    buttonPanel.add(okButton);
    okButton.addActionListener(this);
    JButton cancelButton = new JButton("Cancel");
    cancelButton.setActionCommand("DeclineSelection");
    cancelButton.setToolTipText(TooltipTexts.CANCEL_GENERAL);
    buttonPanel.add(cancelButton);
    cancelButton.addActionListener(this);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    pack(); 
    setVisible(true);
  }
  
  /**
   * This static method takes as input parameters the name of the checkboxes for the dialog
   * and returns the names of the checked ones if "OK" has been pressed, for the "Cancel"
   * button an empty vector is returned
   * @param parent the parent frame
   * @param title  the title of the dialog field
   * @param messages the messages in the dialog box (each string is a separate row)
   * @param chechboxStrings the name of the checkboxes that should be displayed in the dialog box
   * @return the checked checkboxes, if cancel is pressed the returned Vector is empty
   */
  
  public static Vector<String> showConfirmDialog(JFrame parent, String title, Vector<String> messages, ArrayList<String> chechboxStrings, boolean showOne){  
    CheckBoxOptionPane dialog = new CheckBoxOptionPane(parent, title, messages, chechboxStrings, showOne);
    return dialog.getSelected();
  }
  
  public Vector<String> getSelected(){
    Vector<String> selected = new Vector<String>();
    if (getValue()==JOptionPane.YES_OPTION){
      if (checkBoxStrings_.size()==1 && !showOne_){
        selected = new Vector<String>(checkBoxStrings_);
      }else{
        for (String name : checkBoxStrings_){
          if (checkBoxHash_.get(name).isSelected()) selected.add(name);
        }
      }
    }
    return selected;
  }

  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase("AcceptSelection")){
      value_ = JOptionPane.YES_OPTION;
      setVisible(false); 
      dispose(); 
    }else if (e.getActionCommand().equalsIgnoreCase("DeclineSelection")){
      value_ = JOptionPane.NO_OPTION;
      setVisible(false); 
      dispose(); 
    }   
  }

  public int getValue()
  {
    return value_;
  }
  
  
}
