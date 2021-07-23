/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2021 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.IsoLabelInputException;
import at.tugraz.genome.lda.vos.IsotopicLabelVO;


/**
 * The graphical dialog for selecting and entering information about the isotopic labels for the export
 * @author Juergen Hartler
 *
 */
public class OmegaExportDialog extends JDialog implements ActionListener
{

  private static final long serialVersionUID = -3339973452383692517L;

  /** the listener for further processing the export*/
  private ActionListener parent_;
  /** the selected return option*/
  private int value_ = JOptionPane.NO_OPTION;
  /** the selected labels for the export*/
  private JPanel labelSelection_;
  /** the check boxes of the labels*/  
  private Vector<IsotopicLabelSelection> labelInformation_;
  /** button for adding a line for the isotopicLabel */
  private JButton addButton_;
  /** the entered label information*/
  private Vector<IsotopicLabelVO> enteredLabels_ = new Vector<IsotopicLabelVO>();


  /**
   * 
   * @param title title of the select dialog
   * @param potentialLabels automatically extracted information on the isotopic labels
   * @param parent the parent listener to process further actions
   */
  public OmegaExportDialog(String title, List<IsotopicLabelVO> potentialLabels, ActionListener parent) {
    super(new JFrame(), title, true);
    this.parent_ = parent;
    
    labelSelection_ = new JPanel();
    initCheckboxes(labelSelection_, potentialLabels);
 
    this.setLayout(new BorderLayout());
    setLocation(380,240);
    this.add(labelSelection_,BorderLayout.CENTER);
    this.initButtonPanel();
    enteredLabels_ = new Vector<IsotopicLabelVO>();
  }
  
  
  
  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase("AcceptOmegaExport")){
      boolean ok = extractLabelInformation();
      if (ok) {
        value_ = JOptionPane.YES_OPTION;
        setVisible(false); 
        dispose();
        parent_.actionPerformed(e);
      } else 
        value_ = JOptionPane.NO_OPTION;
    }else if (e.getActionCommand().equalsIgnoreCase("DeclineSelection")){
      value_ = JOptionPane.NO_OPTION;
      setVisible(false); 
      dispose(); 
    }else if (e.getActionCommand().equalsIgnoreCase("SelectAll")||
        e.getActionCommand().equalsIgnoreCase("SelectNone")||
        e.getActionCommand().equalsIgnoreCase("Invert")){
      for (IsotopicLabelSelection labelInfo : labelInformation_){
        if (e.getActionCommand().equalsIgnoreCase("SelectAll")) labelInfo.setSelected(true);
        if (e.getActionCommand().equalsIgnoreCase("SelectNone")) labelInfo.setSelected(false);
        if (e.getActionCommand().equalsIgnoreCase("Invert")) labelInfo.invertSelection();
      }
    }else if (e.getActionCommand().startsWith(IsotopicLabelSelection.DELETE_ACTION_PREFIX)) {
      int row = Integer.parseInt(e.getActionCommand().substring(IsotopicLabelSelection.DELETE_ACTION_PREFIX.length()));
      this.labelInformation_.remove(row-1);
      this.labelSelection_.removeAll();
      this.labelSelection_.setVisible(false);
      initCheckboxes(this.labelSelection_,labelInformation_);
      //this.add(labelSelection_,BorderLayout.CENTER);
      this.labelSelection_.setVisible(true);
      this.pack();
      this.repaint();
    }else if (e.getActionCommand().startsWith("Add")) {
      this.labelSelection_.setVisible(false);
      this.labelSelection_.remove(addButton_);
      this.labelInformation_.add(IsotopicLabelSelection.addIsotopicLabelToPanel(new IsotopicLabelVO("", new Hashtable<String,Integer>(), new LinkedHashMap<String,Integer>()),labelSelection_,this.labelInformation_.size()+1,this));
      this.labelSelection_.add(addButton_, new GridBagConstraints(4, this.labelInformation_.size()+2, 2, 1, 0.0, 0.0,
          GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
      this.labelSelection_.setVisible(true);
      this.pack();
      this.repaint();
    }
  }
  
  /**
   * initializes the button panel
   */
  private void initButtonPanel(){
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
    JButton okButton = new JButton("OK");
    okButton.setActionCommand("AcceptOmegaExport");
    okButton.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    buttonPanel.add(okButton);
    okButton.addActionListener(this);
    JButton cancelButton = new JButton("Cancel");
    cancelButton.setActionCommand("DeclineSelection");
    cancelButton.setToolTipText(TooltipTexts.CANCEL_GENERAL);
    buttonPanel.add(cancelButton);
    cancelButton.addActionListener(this);
    setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
    pack(); 
    setVisible(false);

  }
  
  
  /**
   * initializes input panels for entering information about the isotopic labels
   * @param panel the graphical panel
   * @param labels automatically extracted information on the isotopic labels
   */
  private void initCheckboxes(JPanel panel, List<? extends IsotopicLabelVO> labels){
    panel.setLayout(new GridBagLayout());
    labelInformation_ = new Vector<IsotopicLabelSelection>();
    JLabel text;
    
    text = new JLabel("Label");
    panel.add(text, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
        GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    text = new JLabel("\u03C9-Pos.");
    panel.add(text, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
        GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));    
    text = new JLabel("Formula");
    panel.add(text, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0,
        GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    text = new JLabel("rRT-shift");
    panel.add(text, new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0,
        GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    
    for (int i=0; i!=labels.size(); i++){
        labelInformation_.add(IsotopicLabelSelection.addIsotopicLabelToPanel(labels.get(i),panel,i+1,this));      
    }
    int row = labels.size()+1;
    addButton_ = new JButton("Add");
    addButton_.setActionCommand("AddLabel");
    addButton_.setToolTipText(TooltipTexts.ISO_ADD);
    panel.add(addButton_, new GridBagConstraints(4, row, 2, 1, 0.0, 0.0,
        GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    addButton_.addActionListener(this);
  }


  /**
   * 
   * @return selected return option
   */
  public int getValue()
  {
    return value_;
  }
  
  /**
   * verifies if the entered values are conform with the syntax();
   * @return true when the entered values are fit for further processing
   */
  private boolean extractLabelInformation() {
    enteredLabels_ = new Vector<IsotopicLabelVO>();
    boolean allOk = true;
    IsotopicLabelVO label = null;
    for (int i=0; i!=labelInformation_.size(); i++) {
      IsotopicLabelSelection inputField = labelInformation_.get(i);
      try {
        if (!inputField.isSelected())
          continue;
        label = inputField.getLabelInformation();
        enteredLabels_.add(label);
      }
      catch (IsoLabelInputException e) {
        new WarningMessage(new JFrame(), "Error", "Label \""+inputField.getEnteredLabelId()+"\" at line "+(i+1)+": "+e.getMessage());
        allOk = false;
        break;
      }
    }
    return allOk;
  }
  
  /**
   * 
   * @return the final results of the entered labels
   */
  public Vector<IsotopicLabelVO> getEnteredLabelInformation() {
    return enteredLabels_;
  }
  

}
