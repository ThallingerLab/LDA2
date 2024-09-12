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
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.verifier.IntegerVerifier;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ChromExportDialog extends JDialog implements ActionListener
{

  private static final long serialVersionUID = -2009258591034858049L;
  private final int columns_ = 8;
  private final int WIDTH_DEFAULT = 1024;
  private final int HEIGHT_DEFAULT = 768;
  
  private int value_ = JOptionPane.NO_OPTION;
  private JTabbedPane tabs_;
  private ArrayList<String> experiments_;
  private ArrayList<String> analytes_;
  private ArrayList<String> modifications_;
  private Hashtable<String,JCheckBox> expBoxes_;
  private Hashtable<String,JCheckBox> analBoxes_;
  private Hashtable<String,JCheckBox> modBoxes_;
  private JPanel expSelection_;
  private JPanel moleculeSelection_;
  private JScrollPane scrollPane_;
  private JPanel modificationsSelection_;
  private JRadioButton pngRadio_;
  private JRadioButton svgRadio_;
  private ActionListener parent_;
  private JTextField width_;
  private JTextField height_;
  private JButton okButton_;

  public ChromExportDialog(String title, ArrayList<String> experiments, ArrayList<String> analytes, ArrayList<String> modifications, ActionListener parent){
    super(new JFrame(), title, true);
    experiments_ = experiments;
    analytes_ = analytes;
    modifications_ = modifications;
    this.setLayout(new BorderLayout());
    setLocation(380,240);
    tabs_ = new JTabbedPane();
    parent_ = parent;
    
    expSelection_ = new JPanel();
    expBoxes_ = this.initCheckboxes(expSelection_, "Samples", experiments_);
    moleculeSelection_ = new JPanel();
    analBoxes_ = this.initCheckboxes(moleculeSelection_, "Analytes" ,analytes_);
    scrollPane_ = new JScrollPane(moleculeSelection_);
    scrollPane_.setPreferredSize(new Dimension(1000, 600));
    modificationsSelection_ = new JPanel();
    modBoxes_ = this.initCheckboxes(modificationsSelection_, "Modifications",modifications_);
    
    tabs_.addTab("Sample", expSelection_);
    tabs_.setToolTipTextAt(0, TooltipTexts.TABS_CHROMEXPORT_EXPS);
    tabs_.addTab("Analyte", scrollPane_);
    tabs_.setToolTipTextAt(1, TooltipTexts.TABS_CHROMEXPORT_ANALS);
    tabs_.addTab("Modifications", modificationsSelection_);
    tabs_.setToolTipTextAt(2, TooltipTexts.TABS_CHROMEXPORT_MODS);
    this.add(tabs_,BorderLayout.CENTER);
    initButtonPanel();
  }
  
  private Hashtable<String,JCheckBox> initCheckboxes(JPanel panel, String title, ArrayList<String> analytes){
    panel.setLayout(new GridBagLayout());
    Hashtable<String,JCheckBox> boxes = new Hashtable<String,JCheckBox>();
    int titleColSpan = analytes.size();
    if (titleColSpan>columns_)
      titleColSpan = columns_;
    if (title!=null&&title.length()>0){
      JLabel titleLabel = new JLabel(title);
      panel.add(titleLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    }
    
    int rows = analytes.size()/columns_;
    if (analytes.size()>0 && analytes.size()%columns_!=0)
      rows++;
    for (int i=0; i!=analytes.size(); i++){
      String name = analytes.get(i);
      JCheckBox select  = new JCheckBox(name);
      select.setSelected(true);
      int row = i%rows;
      int column = i/rows;
      select.setToolTipText(TooltipTexts.GENERAL_ACCEPT_SINGLE_START+name+TooltipTexts.GENERAL_ACCEPT_SINGLE_END);
      panel.add(select, new GridBagConstraints(column, row+1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
      boxes.put(name, select);
    }
    return boxes;
  }
  
  private void initButtonPanel(){
    JPanel buttonAndRadio = new JPanel();
    buttonAndRadio.setLayout(new GridBagLayout());
    JPanel radioPanel = new JPanel();
    buttonAndRadio.add(radioPanel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 0, 1, 0), 0, 0));
    ButtonGroup pictType = new ButtonGroup();
    pngRadio_ = new JRadioButton("PNG");
    pngRadio_.setToolTipText(TooltipTexts.EXPORT_PNG);
    pngRadio_.setSelected(true);
    pictType.add(pngRadio_);
    radioPanel.add(pngRadio_);
    svgRadio_ = new JRadioButton("SVG");
    svgRadio_.setToolTipText(TooltipTexts.EXPORT_SVG);
    pictType.add(svgRadio_);
    radioPanel.add(svgRadio_);
    JLabel widthLabel = new JLabel("Width: ");
    radioPanel.add(widthLabel);
    width_ = new JTextField(5);
    width_.setHorizontalAlignment(JTextField.RIGHT);
    width_.setInputVerifier(new IntegerVerifier(true));
    width_.setText(String.valueOf(WIDTH_DEFAULT));
    radioPanel.add(width_);
    JLabel heightLabel = new JLabel("Height: ");
    radioPanel.add(heightLabel);
    height_ = new JTextField(5);
    height_.setHorizontalAlignment(JTextField.RIGHT);
    height_.setInputVerifier(new IntegerVerifier(true));
    height_.setText(String.valueOf(HEIGHT_DEFAULT));
    radioPanel.add(height_);
    
    JPanel buttonPanel = new JPanel();
    buttonAndRadio.add(buttonPanel,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 0, 1, 0), 0, 0));
    this.add(buttonAndRadio,BorderLayout.SOUTH);
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
    okButton_ = new JButton("OK");
    okButton_.setActionCommand("AcceptChromExport");
    okButton_.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    buttonPanel.add(okButton_);
    okButton_.addActionListener(this);
    okButton_.addActionListener(parent_);
    JButton cancelButton = new JButton("Cancel");
    cancelButton.setActionCommand("DeclineSelection");
    cancelButton.setToolTipText(TooltipTexts.CANCEL_GENERAL);
    buttonPanel.add(cancelButton);
    cancelButton.addActionListener(this);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    pack(); 
    setVisible(false);
  }

  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase("AcceptChromExport")){
      value_ = JOptionPane.YES_OPTION;
      setVisible(false); 
      dispose(); 
    }else if (e.getActionCommand().equalsIgnoreCase("DeclineSelection")){
      value_ = JOptionPane.NO_OPTION;
      setVisible(false); 
      dispose(); 
    }else if (e.getActionCommand().equalsIgnoreCase("SelectAll")||
        e.getActionCommand().equalsIgnoreCase("SelectNone")||
        e.getActionCommand().equalsIgnoreCase("Invert")){
    	ArrayList<String> names = null;
      Hashtable<String,JCheckBox> boxes = null;
      if (tabs_.getSelectedIndex()==0){
        names = experiments_;
        boxes = expBoxes_;
      } else if (tabs_.getSelectedIndex()==1){
        names = analytes_;
        boxes = analBoxes_;
      } else if (tabs_.getSelectedIndex()==2){
        names = modifications_;
        boxes = modBoxes_;
      }
      if (names==null || boxes==null) return;
      for (String name : names){
        JCheckBox box = boxes.get(name);
        if (e.getActionCommand().equalsIgnoreCase("SelectAll")) box.setSelected(true);
        if (e.getActionCommand().equalsIgnoreCase("SelectNone")) box.setSelected(false);
        if (e.getActionCommand().equalsIgnoreCase("Invert")) box.setSelected(!box.isSelected());
      }
    }
  }
  
  public int getValue()
  {
    return value_;
  }
  
  public void refreshNames(ArrayList<String> names){
    for (int i=0; i!=names.size();i++){
      String oldName = experiments_.get(i);
      String newName = names.get(i);
      if (!oldName.equalsIgnoreCase(newName)){
        JCheckBox check = expBoxes_.get(oldName);
        check.setText(newName);
        check.invalidate();
        check.updateUI();
        expBoxes_.put(newName,check);
      }
    }
    experiments_ = names;
    correctDialogSize();
  }
  
  public void updateAnalyteNames(ArrayList<String> names){
  	this.tabs_.remove(1);
  	this.analytes_ = names;
  	this.moleculeSelection_ = new JPanel();
  	this.analBoxes_ = this.initCheckboxes(moleculeSelection_, "Analytes" ,analytes_);
  	this.scrollPane_ = new JScrollPane(moleculeSelection_);
  	scrollPane_.setPreferredSize(new Dimension(1000, 600));
  	tabs_.insertTab("Analyte", null, scrollPane_, TooltipTexts.TABS_CHROMEXPORT_ANALS, 1);
  	correctDialogSize();
  	invalidate();
    doLayout();
    this.setSize(this.getPreferredSize());
    repaint();
  }
  
  public void checkMolecules(Vector<String> anals, boolean check){
    for (String anal : anals){
      JCheckBox box = analBoxes_.get(anal);
      box.setSelected(check);
      box.invalidate();
      box.updateUI();
    }
    correctDialogSize();
  }
  
  private void correctDialogSize(){
//    expSelection_.invalidate();
//    expSelection_.doLayout();
//    moleculeSelection_.invalidate();
//    moleculeSelection_.doLayout();
//    int width1 = expSelection_.getPreferredSize().width;
//    int height1 = expSelection_.getPreferredSize().height;
//    int width2 = moleculeSelection_.getPreferredSize().width;
//    int height2 = moleculeSelection_.getPreferredSize().height;;
//    int width = width1>width2 ? width1 : width2;
//    int height = height1>height2 ? height1 : height2;
//    System.out.println("Width: "+width+";"+width1+";"+width2);
//    System.out.println("Height: "+height+";"+height1+";"+height2);
//    
//    expSelection_.setSize(width,height);
//    moleculeSelection_.setSize(width,height);
//    expSelection_.repaint();
//    moleculeSelection_.repaint();
//    System.out.println("-------------------------------");
    invalidate();
    doLayout();
    this.setSize(this.getPreferredSize());
    repaint();
    pack();
  }
  
  public boolean isExperimentSelected(String expName){
    if (expBoxes_.containsKey(expName) && expBoxes_.get(expName).isSelected()) return true;
    else return false;
  }
  
  public boolean isMoleculeSelected(String molName){
    if (analBoxes_.containsKey(molName) && analBoxes_.get(molName).isSelected()) return true;
    else return false;
  }
  
  public boolean isModificationSelected(String molName){
    if (modBoxes_.containsKey(molName) && modBoxes_.get(molName).isSelected()) return true;
    else return false;
  }
  
  public String getPictureType(){
    if (pngRadio_.isSelected())return ExportPanel.EXPORT_PNG;
    else return ExportPanel.EXPORT_SVG;
  }
  
  public Dimension getPictureDimension(){
    return new Dimension(Integer.parseInt(this.width_.getText()),Integer.parseInt(this.height_.getText()));
  }
  
  /**
   * Removes parent action listener.
   */
  public void cleanup()
  {
  	okButton_.removeActionListener(parent_);
  	parent_ = null;
  }
}
