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
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.exception.AbsoluteSettingsInputException;
import at.tugraz.genome.lda.verifier.DoubleVerifier;
import at.tugraz.genome.lda.vos.VolumeConcVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class StandardsSettingsPanel extends JPanel implements ActionListener
{
  private static final long serialVersionUID = 1L;
  private String expName_;
  private Vector<String> isNames_;
  private Vector<String> esNames_;
  private StandardsSettingsListener listener_;
  
  private JPanel generalInputPanel_;
  private JPanel specificInputPanel_;
  
  private JTextField dilutionFactor_;
  private PhysicalUnitInput generalExtVolume_;
  private PhysicalUnitInput generalExtConc_;
  private PhysicalUnitInput generalIntVolume_;
  private PhysicalUnitInput generalIntConc_;
  
  private Hashtable<String,PhysicalUnitInput> extVolumeSettings_;
  private Hashtable<String,PhysicalUnitInput> extConcSettings_;
  private Hashtable<String,PhysicalUnitInput> intVolumeSettings_;
  private Hashtable<String,PhysicalUnitInput> intConcSettings_;
  
  
  public StandardsSettingsPanel(String expName, Vector<String> isNames, Vector<String> esNames, StandardsSettingsListener listener){
    expName_ = expName;
    this.isNames_ = isNames;
    this.esNames_ = esNames;
    this.listener_ = listener;
    this.initComponents();
  }
  
  private void initComponents(){
//    setLayout(new BorderLayout());
    generalInputPanel_ = new JPanel();
    this.setLayout(new GridBagLayout());

    JPanel topPanel = new JPanel();
    JLabel label = new JLabel("Dilution factor: ");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_STANDARDS_DILLUTION);
    topPanel.add(label);
    this.add(topPanel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
//    topPanel.add(label,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
//      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    dilutionFactor_ = new JTextField(3);
    dilutionFactor_.setHorizontalAlignment(JTextField.RIGHT);
    dilutionFactor_.setInputVerifier(new DoubleVerifier());
    dilutionFactor_.setToolTipText(TooltipTexts.STATISTICS_ABS_STANDARDS_DILLUTION);
    topPanel.add(dilutionFactor_);
    if (expName_!=null && expName_.length()>0){
      ApplyButton button = new ApplyButton("Apply to all","applySettingsAll"+";dilution");
      button.addActionListener(this);
      button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_ALL);
      topPanel.add(button);
//    generalInputPanel_.add(button,new GridBagConstraints(6, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
      button = new ApplyButton("Apply to group","applySettingsGroup"+";dilution");
      button.addActionListener(this);
      button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_GROUP);
      topPanel.add(button);
    }
//    generalInputPanel_.add(button,new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));

//    topPanel.add(dilutionFactor_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

    
//    this.add(generalInputPanel_/*, BorderLayout.NORTH*/);
    this.add(generalInputPanel_,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    generalInputPanel_.setVisible(true);
    generalInputPanel_.setLayout(new GridBagLayout());
    
    label = new JLabel("Ext volume: ");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_EXT_VOLUME);
    generalInputPanel_.add(label,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    generalExtVolume_ = new PhysicalUnitInput("\u03BC",TooltipTexts.STATISTICS_ABS_EXT_VOLUME);
    generalInputPanel_.add(generalExtVolume_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    label = new JLabel("L");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_EXT_VOLUME);
    generalInputPanel_.add(label,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));  
    label = new JLabel("Ext conc.: ");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_EXT_CONC);
    generalInputPanel_.add(label,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    generalExtConc_ = new PhysicalUnitInput("\u03BC",TooltipTexts.STATISTICS_ABS_EXT_CONC);
    generalInputPanel_.add(generalExtConc_,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    label = new JLabel("mol/L");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_EXT_CONC);
    generalInputPanel_.add(label,new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));
    label = new JLabel("Int volume: ");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_INT_VOLUME);
    generalInputPanel_.add(label,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    generalIntVolume_ = new PhysicalUnitInput("\u03BC",TooltipTexts.STATISTICS_ABS_INT_VOLUME);
    generalInputPanel_.add(generalIntVolume_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    label = new JLabel("L");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_INT_VOLUME);
    generalInputPanel_.add(label,new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));  
    label = new JLabel("Int conc.: ");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_INT_CONC);
    generalInputPanel_.add(label,new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    generalIntConc_ = new PhysicalUnitInput("\u03BC",TooltipTexts.STATISTICS_ABS_INT_CONC);
    generalInputPanel_.add(generalIntConc_,new GridBagConstraints(4, 1, 1, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    label = new JLabel("mol/L");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_INT_CONC);
    generalInputPanel_.add(label,new GridBagConstraints(5, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));
    extVolumeSettings_ = new Hashtable<String,PhysicalUnitInput>();
    extConcSettings_ = new Hashtable<String,PhysicalUnitInput>();
    intVolumeSettings_ = new Hashtable<String,PhysicalUnitInput>();
    intConcSettings_ = new Hashtable<String,PhysicalUnitInput>();
    if (expName_!=null && expName_.length()>0){
      ApplyButton button = new ApplyButton("Apply to all","applySettingsAll"+";generalES");
      button.addActionListener(this);
      button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_ALL);      
      generalInputPanel_.add(button,new GridBagConstraints(6, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
      button = new ApplyButton("Apply to group","applySettingsGroup"+";generalES");
      button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_GROUP);
      button.addActionListener(this);
      generalInputPanel_.add(button,new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
      button = new ApplyButton("Apply to all","applySettingsAll"+";generalIS");
      button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_ALL);
      button.addActionListener(this);
      generalInputPanel_.add(button,new GridBagConstraints(6, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
      button = new ApplyButton("Apply to group","applySettingsGroup"+";generalIS");
      button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_GROUP);
      button.addActionListener(this);
      generalInputPanel_.add(button,new GridBagConstraints(7, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));

    }   

    specificInputPanel_ = new JPanel();  
//    this.add(specificInputPanel_);
    this.add(specificInputPanel_,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    JPanel panelToAdd = specificInputPanel_;
    if ((esNames_.size()+isNames_.size())>2){
     panelToAdd = new JPanel();
     panelToAdd.setLayout(new GridBagLayout());
     panelToAdd.setPreferredSize(new Dimension(350,30*(esNames_.size()+isNames_.size())));
    }else{
      specificInputPanel_.setLayout(new GridBagLayout());
    }
    specificInputPanel_.setVisible(false);
    for (int i=0; i!=esNames_.size(); i++){
      String esName = esNames_.get(i);
      label = new JLabel("Ext "+esName+"-vol");
      label.setToolTipText(TooltipTexts.STATISTICS_ABS_EXT_VOLUME);
      panelToAdd.add(label,new GridBagConstraints(0, i, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
      PhysicalUnitInput volumeSetting  = new PhysicalUnitInput("\u03BC",TooltipTexts.STATISTICS_ABS_EXT_VOLUME);
      panelToAdd.add(volumeSetting,new GridBagConstraints(1, i, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      label = new JLabel("L");
      label.setToolTipText(TooltipTexts.STATISTICS_ABS_EXT_VOLUME);
      panelToAdd.add(label,new GridBagConstraints(2, i, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));
      extVolumeSettings_.put(esName, volumeSetting);
      label = new JLabel("Ext "+esName+"-conc");
      label.setToolTipText(TooltipTexts.STATISTICS_ABS_EXT_CONC);
      panelToAdd.add(label,new GridBagConstraints(3, i, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
      PhysicalUnitInput concSetting = new PhysicalUnitInput("\u03BC",TooltipTexts.STATISTICS_ABS_EXT_CONC);
      panelToAdd.add(concSetting,new GridBagConstraints(4, i, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      label = new JLabel("mol/L");
      label.setToolTipText(TooltipTexts.STATISTICS_ABS_EXT_CONC);
      panelToAdd.add(label,new GridBagConstraints(5, i, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));
      extConcSettings_.put(esName, concSetting);
      if (expName_!=null && expName_.length()>0){
        ApplyButton button = new ApplyButton("Apply to all","applySettingsAll"+";"+esName);
        button.addActionListener(this);
        button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_ALL);
        panelToAdd.add(button,new GridBagConstraints(6, i, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
        button = new ApplyButton("Apply to group","applySettingsGroup"+";"+esName);
        button.addActionListener(this);
        button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_GROUP);
        panelToAdd.add(button,new GridBagConstraints(7, i, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
      }  
    }
    int rowAdd = esNames_.size();
    for (int i=0; i!=isNames_.size(); i++){
      String isName = isNames_.get(i);
      label = new JLabel("Int "+isName+"-vol");
      label.setToolTipText(TooltipTexts.STATISTICS_ABS_INT_VOLUME);
      panelToAdd.add(label,new GridBagConstraints(0, i+rowAdd, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
      PhysicalUnitInput volumeSetting  = new PhysicalUnitInput("\u03BC",TooltipTexts.STATISTICS_ABS_INT_VOLUME);
      panelToAdd.add(volumeSetting,new GridBagConstraints(1, i+rowAdd, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      label.setToolTipText(TooltipTexts.STATISTICS_ABS_INT_VOLUME);
      label = new JLabel("L");
      panelToAdd.add(label,new GridBagConstraints(2, i+rowAdd, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));
      intVolumeSettings_.put(isName, volumeSetting);
      label = new JLabel("Int "+isName+"-conc");
      label.setToolTipText(TooltipTexts.STATISTICS_ABS_INT_CONC);
      panelToAdd.add(label,new GridBagConstraints(3, i+rowAdd, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
      PhysicalUnitInput concSetting = new PhysicalUnitInput("\u03BC",TooltipTexts.STATISTICS_ABS_INT_CONC);
      panelToAdd.add(concSetting,new GridBagConstraints(4, i+rowAdd, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      label = new JLabel("mol/L");
      label.setToolTipText(TooltipTexts.STATISTICS_ABS_INT_CONC);
      panelToAdd.add(label,new GridBagConstraints(5, i+rowAdd, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));
      intConcSettings_.put(isName, concSetting);
      if (expName_!=null && expName_.length()>0){
        ApplyButton button = new ApplyButton("Apply to all","applySettingsAll"+";"+isName);
        button.addActionListener(this);
        button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_ALL);
        panelToAdd.add(button,new GridBagConstraints(6, i+rowAdd, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
        button = new ApplyButton("Apply to group","applySettingsGroup"+";"+isName);
        button.addActionListener(this);
        button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_GROUP);
        panelToAdd.add(button,new GridBagConstraints(7, i+rowAdd, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
      }  
    }
    if ((esNames_.size()+isNames_.size())>2){
      JScrollPane scrollPane = new JScrollPane(panelToAdd);
      scrollPane.setPreferredSize(new Dimension(700, 65));
      specificInputPanel_.add(scrollPane);
    }
      
  }
  
  public void showSpecific(boolean specific){   
    generalInputPanel_.setVisible(!specific);
    specificInputPanel_.setVisible(specific);
    
  }
  
  public void takeOverValues(StandardsSettingsPanel otherPanel, boolean general){
    dilutionFactor_.setText(otherPanel.dilutionFactor_.getText());
    if (general){
      generalExtVolume_.setInputValue(otherPanel.generalExtVolume_);
      generalExtConc_.setInputValue(otherPanel.generalExtConc_);
      generalIntVolume_.setInputValue(otherPanel.generalIntVolume_);
      generalIntConc_.setInputValue(otherPanel.generalIntConc_);
    }else{
      for (String esName : esNames_){
        extVolumeSettings_.get(esName).setInputValue(otherPanel.extVolumeSettings_.get(esName));
        extConcSettings_.get(esName).setInputValue(otherPanel.extConcSettings_.get(esName));
      }
      for (String isName : isNames_){
        intVolumeSettings_.get(isName).setInputValue(otherPanel.intVolumeSettings_.get(isName));
        intConcSettings_.get(isName).setInputValue(otherPanel.intConcSettings_.get(isName));
      }      
    }
  }
  
  public void copyGeneralSettingsToStandards(){
    for (String esName : esNames_){
      extVolumeSettings_.get(esName).setInputValue(generalExtVolume_);
      extConcSettings_.get(esName).setInputValue(generalExtConc_);
    }
    for (String isName : isNames_){
      intVolumeSettings_.get(isName).setInputValue(generalIntVolume_);
      intConcSettings_.get(isName).setInputValue(generalIntConc_);
    }
  }
  
  public void setInputValues(String molName, StandardsSettingsPanel copyFrom){
    if (molName.equalsIgnoreCase("dilution")){
      dilutionFactor_.setText(copyFrom.dilutionFactor_.getText());
    }else if (molName.equalsIgnoreCase("generalES")){
      generalExtVolume_.setInputValue(copyFrom.generalExtVolume_);
      generalExtConc_.setInputValue(copyFrom.generalExtConc_);
    } else if (molName.equalsIgnoreCase("generalIS")){
      generalIntVolume_.setInputValue(copyFrom.generalIntVolume_);
      generalIntConc_.setInputValue(copyFrom.generalIntConc_);
    } else {
      if (extVolumeSettings_.containsKey(molName)){
        extVolumeSettings_.get(molName).setInputValue(copyFrom.extVolumeSettings_.get(molName));
        extConcSettings_.get(molName).setInputValue(copyFrom.extConcSettings_.get(molName));
      }
      if (intVolumeSettings_.containsKey(molName)){
        intVolumeSettings_.get(molName).setInputValue(copyFrom.intVolumeSettings_.get(molName));
        intConcSettings_.get(molName).setInputValue(copyFrom.intConcSettings_.get(molName));
      }
    }
  }

  public void actionPerformed(ActionEvent e)
  {
    String actionCommand = e.getActionCommand();
    if (actionCommand.indexOf(";")!=-1){
      String command = actionCommand.substring(0,actionCommand.indexOf(";"));
      System.out.println("command: "+actionCommand);
      String molecule = actionCommand.substring(actionCommand.indexOf(";")+1);
      if (command.equalsIgnoreCase("applySettingsAll")){
        listener_.applySettingsToAll(expName_, molecule);
      }else if (command.equalsIgnoreCase("applySettingsGroup")){
        listener_.applySettingsToGroup(expName_, molecule);
      }
    }  
  }

  public PhysicalUnitInput getGeneralExtVolume()
  {
    return generalExtVolume_;
  }

  public PhysicalUnitInput getGeneralExtConc()
  {
    return generalExtConc_;
  }

  public PhysicalUnitInput getGeneralIntVolume()
  {
    return generalIntVolume_;
  }

  public PhysicalUnitInput getGeneralIntConc()
  {
    return generalIntConc_;
  }

  public Hashtable<String,PhysicalUnitInput> getExtVolumeSettings()
  {
    return extVolumeSettings_;
  }

  public Hashtable<String,PhysicalUnitInput> getExtConcSettings()
  {
    return extConcSettings_;
  }

  public Hashtable<String,PhysicalUnitInput> getIntVolumeSettings()
  {
    return intVolumeSettings_;
  }

  public Hashtable<String,PhysicalUnitInput> getIntConcSettings()
  {
    return intConcSettings_;
  }
  
  public JTextField getDilutionFactor()
  {
    return dilutionFactor_;
  }

  public void cleanSettings(){
    generalInputPanel_.setVisible(true);
    specificInputPanel_.setVisible(false);
    dilutionFactor_.setText("");
    generalExtVolume_.setDefault();
    generalExtConc_.setDefault();
    generalIntVolume_.setDefault();
    generalIntConc_.setDefault();   
    for (PhysicalUnitInput input: extVolumeSettings_.values()) input.setDefault();
    for (PhysicalUnitInput input: extConcSettings_.values()) input.setDefault();
    for (PhysicalUnitInput input: intVolumeSettings_.values()) input.setDefault();
    for (PhysicalUnitInput input: intConcSettings_.values()) input.setDefault();    
  }
  
  public VolumeConcVO getSettings(boolean isIS)throws AbsoluteSettingsInputException {
    return getSettings("general",isIS);
  }
  
  public VolumeConcVO getSettings(String standName, boolean isIS)throws AbsoluteSettingsInputException {
    VolumeConcVO volVO = null;
    try{
      if (standName.equalsIgnoreCase("general")){
        if (isIS)
          volVO = getSettings(generalIntVolume_, generalIntConc_);
        else
          volVO = getSettings(generalExtVolume_, generalExtConc_);
      }else{
        if (isIS)
          volVO = getSettings(intVolumeSettings_.get(standName), intConcSettings_.get(standName));
        else
          volVO = getSettings(extVolumeSettings_.get(standName), extConcSettings_.get(standName));
      }
    } catch (AbsoluteSettingsInputException ex){
      String standString = "external";
      if (isIS)
        standString = "internal";
      if (!standName.equalsIgnoreCase("general"))
        standString = "standName";
      throw new AbsoluteSettingsInputException("You have to enter values for the "+standString+" standard!");
    }
    return volVO;
  }
  
  private VolumeConcVO getSettings(PhysicalUnitInput volume, PhysicalUnitInput conc) throws AbsoluteSettingsInputException {
    Double volVal = volume.getInputValue();
    if (volVal == null || volVal.doubleValue()<=0)
      throw new  AbsoluteSettingsInputException("The volume is not set");
    Double concVal = conc.getInputValue();
    if (concVal == null || concVal.doubleValue()<=0)
      throw new  AbsoluteSettingsInputException("The concentration is not set");
    return new VolumeConcVO(volVal,concVal);
  }  
}
