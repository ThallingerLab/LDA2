/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import at.tugraz.genome.lda.TooltipTexts;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ResultSelectionSettings extends JDialog implements ActionListener
{

  private static final long serialVersionUID = 8102700619746913853L;
  
  private final int columns_ = 8;

  private Vector<String> moleculeNames_;
  private Hashtable<String,JCheckBox> checkboxes_;
  private JScrollPane scrollPane_;
  
//  private Vector<ActionListener> parents_;
  
  JButton button_;
  boolean isOpen_ = false;
  
  public ResultSelectionSettings(String title, Vector<String> moleculeNames, boolean selected){
  	super(new JFrame(), title);
    moleculeNames_ = moleculeNames;
//    parents_ = new Vector<ActionListener>();
    checkboxes_ = new Hashtable<String,JCheckBox>();
    
    this.setLayout(new BorderLayout());
    setLocation(380,240);
    
    int titleColSpan = moleculeNames.size();
    if (titleColSpan>columns_)
      titleColSpan = columns_;
    
    scrollPane_ = buildJCheckBoxes(selected);
    
    this.add(scrollPane_,BorderLayout.CENTER);
    
    JPanel buttonPanel = new JPanel();
    this.add(buttonPanel,BorderLayout.SOUTH);

    JButton allButton = new JButton("All");
    allButton.setActionCommand("SelectAll");
    allButton.setToolTipText(TooltipTexts.ALL_GENERAL);
    buttonPanel.add(allButton);
    allButton.addActionListener(this);
    JButton noneButton = new JButton("None");
    noneButton.setActionCommand("SelectNone");
    noneButton.setToolTipText(TooltipTexts.NONE_GENERAL);
    buttonPanel.add(noneButton);
    noneButton.addActionListener(this);
    JButton invertButton = new JButton("Invert");
    invertButton.setActionCommand("Invert");
    invertButton.setToolTipText(TooltipTexts.INVERT_GENERAL);
    buttonPanel.add(invertButton);
    invertButton.addActionListener(this);   
    JLabel emptyLabel = new JLabel("   ");
    buttonPanel.add(emptyLabel);
    button_ = new JButton("OK");
    button_.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    button_.setActionCommand("AcceptDisplaySettings");
    buttonPanel.add(button_);
    setVisible(false);
    setDefaultCloseOperation(DO_NOTHING_ON_CLOSE );
    
    this.addWindowListener(new WindowAdapter() {
    	public void windowClosing(WindowEvent we) {
    		close();
    	}
    });
    
    pack(); 
  }
  
  private JScrollPane buildJCheckBoxes(boolean selected)
  {
  	JPanel moleculeSelection = new JPanel(new GridBagLayout());
  	JScrollPane analScrollPane = new JScrollPane(moleculeSelection);
    analScrollPane.setPreferredSize(new Dimension(1000, 600));
  	int rows = moleculeNames_.size()/columns_;
    if (moleculeNames_.size()>0 && moleculeNames_.size()%columns_!=0)
      rows++;
    for (int i=0; i!=moleculeNames_.size(); i++){
      String name = moleculeNames_.get(i);
      JCheckBox select  = new JCheckBox(name);
      select.setSelected(selected);
      int row = i%rows;
      int column = i/rows;
      select.setToolTipText(TooltipTexts.GENERAL_ACCEPT_SINGLE_START+name+TooltipTexts.GENERAL_ACCEPT_SINGLE_END);
      moleculeSelection.add(select, new GridBagConstraints(column, row+1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
      checkboxes_.put(name, select);
    }
    return analScrollPane;
  }
  
  public void close()
  {
  	this.setVisible(false);
  	isOpen_ = false;
  }
  
  public void open()
  {
  	this.setVisible(true);
  	isOpen_ = true;
  }
  
  public boolean isOpen()
  {
  	return this.isOpen_;
  }
  
  public void addActionListener(ActionListener parent){
    button_.addActionListener(parent);
  }
  
  public void removeActionListener(ActionListener parent) {
  	button_.removeActionListener(parent);
  }

  public boolean isSelected(String moleculeName){
    boolean isSelected = false;
    if (checkboxes_.containsKey(moleculeName)&&checkboxes_.get(moleculeName).isSelected())
      isSelected = true;
    return isSelected;
  }
  
  public Vector<String> getSelected(){
    Vector<String> selected = new Vector<String>();
    for (String molName : moleculeNames_){
      if (checkboxes_.containsKey(molName)&&checkboxes_.get(molName).isSelected())
        selected.add(molName);
    }
    return selected;
  }
  
  public void refreshNames(Vector<String> names){
    for (int i=0; i!=names.size();i++){
      String oldName = moleculeNames_.get(i);
      String newName = names.get(i);
      if (!oldName.equalsIgnoreCase(newName)){
        JCheckBox check = checkboxes_.get(oldName);
        check.setText(newName);
        check.invalidate();
        check.updateUI();
        checkboxes_.put(newName,check);
      }
    }
    moleculeNames_ = names;
    invalidate();
    doLayout();
    this.setSize(this.getPreferredSize());
    repaint();
  }
  
  public void updateNames(Vector<String> names, boolean selected)
  {
  	this.remove(scrollPane_);
  	this.checkboxes_.clear();
  	this.moleculeNames_ = names;
  	this.scrollPane_ = buildJCheckBoxes(selected);
  	this.add(scrollPane_,BorderLayout.CENTER);
  	invalidate();
    doLayout();
    this.setSize(this.getPreferredSize());
    repaint();
  }
  

  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase("SelectAll")||
        e.getActionCommand().equalsIgnoreCase("SelectNone")||
        e.getActionCommand().equalsIgnoreCase("Invert")){
      for (String name : checkboxes_.keySet()){
        JCheckBox box = checkboxes_.get(name);
        if (e.getActionCommand().equalsIgnoreCase("SelectAll")) box.setSelected(true);
        if (e.getActionCommand().equalsIgnoreCase("SelectNone")) box.setSelected(false);
        if (e.getActionCommand().equalsIgnoreCase("Invert")) box.setSelected(!box.isSelected());
      }
    }
    
  }
}
