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

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Hashtable;
import java.util.Set;
import java.util.Vector;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.analysis.ComparativeNameExtractor;
import at.tugraz.genome.lda.exception.AbsoluteSettingsInputException;
import at.tugraz.genome.lda.vos.LipidClassSettingVO;
import at.tugraz.genome.lda.vos.VolumeConcVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LipidClassSettingsPanel extends JPanel implements ActionListener,StandardsSettingsListener
{
  private static final long serialVersionUID = 4681131520631813127L;
  
  private static final String CHANGE_IS_STATUS = "changeAllSettingsStatus";
  private static final String CHANGE_STAND_STATUS = "changeAllStandardsStatus";

  private JCheckBox allSettingsSame_;
  private JTabbedPane experimentTabs_;
  private JCheckBox allStandardsSame_;
  private JComboBox<String> alternativeClasses_;
  
  private String className_;
  private ComparativeNameExtractor extractor_;
  private StandardsSettingsPanel generalSettingsPanel_;
  private Hashtable<String,StandardsSettingsPanel> standardsSettings_;
  private GroupsPanel groupsPanel_;
  
  private boolean oneTimeExpCopied_;
  private boolean oneTimeStandCopied_;
  
  public LipidClassSettingsPanel(String className, ComparativeNameExtractor extractor, GroupsPanel groupsPanel){
    className_ = className;
    extractor_= extractor;
    groupsPanel_ = groupsPanel;
    oneTimeExpCopied_ = false;
    oneTimeStandCopied_ = false;
    this.initComponents();
  }
  
  private void initComponents()
  {
  	Vector<String> isNames = extractor_.getISNames(className_);
    Vector<String> esNames = extractor_.getESNames(className_);
  	
  	if (isNames.isEmpty() && esNames.isEmpty())
    {
    	JLabel label1 = new JLabel("No standards were detected for this lipid class.");
    	this.add(label1,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    	JLabel label2 = new JLabel("You may choose to use standards from another class: ");
    	this.add(label2,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    	alternativeClasses_ = new JComboBox<String>(extractor_.getAllClassNames());
      alternativeClasses_.setSelectedItem(className_);
      this.add(alternativeClasses_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      
      standardsSettings_ = new Hashtable<String,StandardsSettingsPanel>();
      for (String expName : this.extractor_.getExpNamesInSequence()){
        StandardsSettingsPanel expStandSettings = new StandardsSettingsPanel(expName,isNames,esNames,this);
        standardsSettings_.put(expName, expStandSettings);
      }
    }
  	else
  	{
  		this.setLayout(new GridBagLayout());
      allSettingsSame_  = new JCheckBox("use same settings for all experiments");
      allSettingsSame_.setSelected(true);
      allSettingsSame_.setActionCommand(CHANGE_IS_STATUS);
      allSettingsSame_.addActionListener(this);
      allSettingsSame_.setToolTipText(TooltipTexts.STATISTICS_ABS_STANDARDS_EXPS_SAME);
      this.add(allSettingsSame_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      allStandardsSame_  = new JCheckBox("use same settings for all standards");
      allStandardsSame_.setSelected(true);
      allStandardsSame_.setActionCommand(CHANGE_STAND_STATUS);
      allStandardsSame_.addActionListener(this);
      allStandardsSame_.setToolTipText(TooltipTexts.STATISTICS_ABS_STANDARDS_STDS_SAME);
      this.add(allStandardsSame_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

      generalSettingsPanel_ = new StandardsSettingsPanel(null,isNames,esNames,this);
      this.add(generalSettingsPanel_,new GridBagConstraints(0, 1, 4, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      experimentTabs_ = new JTabbedPane();
      this.add(experimentTabs_,new GridBagConstraints(0, 2, 4, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      standardsSettings_ = new Hashtable<String,StandardsSettingsPanel>();
      int count = 0;
      for (String expName : this.extractor_.getExpNamesInSequence()){
        StandardsSettingsPanel expStandSettings = new StandardsSettingsPanel(expName,isNames,esNames,this);
        standardsSettings_.put(expName, expStandSettings);
        experimentTabs_.addTab(extractor_.getExpNames().get(expName),expStandSettings);
        experimentTabs_.setToolTipTextAt(count, TooltipTexts.STATISTICS_ABS_STANDARDS_TAB_EXP+expName+"</html>");
        count++;
      }
      experimentTabs_.setVisible(false);
  	}
    
  }

  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase(CHANGE_IS_STATUS)){
      this.changeVisibilitySettingsSame();
      this.invalidate();
      this.updateUI();
    }
    if (e.getActionCommand().equalsIgnoreCase(CHANGE_STAND_STATUS)){
      this.changeVisibilityStandardsSame();
      this.invalidate();
      this.updateUI();
    }    
  }
  
  public void changeVisibilitySettingsSame(){
    generalSettingsPanel_.setVisible(allSettingsSame_.isSelected());
    experimentTabs_.setVisible(!allSettingsSame_.isSelected());      
    if (!oneTimeExpCopied_ && !allSettingsSame_.isSelected()){
      for (StandardsSettingsPanel standardPanel:standardsSettings_.values()){
        standardPanel.takeOverValues(generalSettingsPanel_, allStandardsSame_.isSelected());
      }
      oneTimeExpCopied_ = true;
    }    
  }
  
  public void changeVisibilityStandardsSame(){
    generalSettingsPanel_.showSpecific(!allStandardsSame_.isSelected());      
    if (!oneTimeStandCopied_ && !allStandardsSame_.isSelected()){
      if (allSettingsSame_.isSelected()){
        generalSettingsPanel_.copyGeneralSettingsToStandards();
      }else{
        for (StandardsSettingsPanel panel : standardsSettings_.values()){
          panel.copyGeneralSettingsToStandards();
        }
      }      
      oneTimeStandCopied_ = true;
    } 

    for (StandardsSettingsPanel panel : standardsSettings_.values()){
      panel.showSpecific(!allStandardsSame_.isSelected());
    }    
  }
  
  public void applySettingsToAll(String experimentName, String molName)
  {
    StandardsSettingsPanel copyFrom = standardsSettings_.get(experimentName);
    for (String expName : standardsSettings_.keySet()){
      if (!expName.equalsIgnoreCase(experimentName)){
        StandardsSettingsPanel panel = standardsSettings_.get(expName);
        panel.setInputValues(molName,copyFrom);
      }
    }
  }

  public void applySettingsToGroup(String experimentName, String molName)
  {
    StandardsSettingsPanel copyFrom = standardsSettings_.get(experimentName);
    Set<String> exps = groupsPanel_.getExpsOfGroupOneExpBelongsTo(experimentName,this.extractor_.getExpNamesInSequence());
    if (exps.size()>0){
      for (String expName : exps){
        if (!expName.equalsIgnoreCase(experimentName)){
          StandardsSettingsPanel panel = standardsSettings_.get(expName);
          panel.setInputValues(molName,copyFrom);
        }  
      }
    }else{
      new WarningMessage(new JFrame(), "Warning", "This experiment does not belong to any group");
    }
  }

  public JCheckBox getAllSettingsSame()
  {
    return allSettingsSame_;
  }

  public JCheckBox getAllStandardsSame()
  {
    return allStandardsSame_;
  }

  public StandardsSettingsPanel getGeneralSettingsPanel()
  {
    return generalSettingsPanel_;
  }

  public Hashtable<String,StandardsSettingsPanel> getStandardsSettings()
  {
    return standardsSettings_;
  }
  
  public void cleanSettings(){
    oneTimeExpCopied_ = false;
    oneTimeStandCopied_ = false;
    if (allSettingsSame_ != null)
    {
    	allSettingsSame_.setSelected(true);
    }
    if (allStandardsSame_ != null)
    {
    	allStandardsSame_.setSelected(true);
    }
    if (generalSettingsPanel_ != null)
    {
    	generalSettingsPanel_.setVisible(true);
    	generalSettingsPanel_.cleanSettings();
    }
    if (experimentTabs_ != null)
    {
    	experimentTabs_.setVisible(false);
    }
    for (StandardsSettingsPanel standSets : standardsSettings_.values())
      standSets.cleanSettings();
  }
  
  public LipidClassSettingVO getSettings() throws AbsoluteSettingsInputException{
    Hashtable<String,Double> dilutionFactors = new Hashtable<String,Double>();
    Hashtable<String,Hashtable<String,VolumeConcVO>> esStandards = new Hashtable<String,Hashtable<String,VolumeConcVO>>();
    Hashtable<String,Hashtable<String,VolumeConcVO>> isStandards = new Hashtable<String,Hashtable<String,VolumeConcVO>>();
    for (String expName : standardsSettings_.keySet()){
      String valueString = null;
      if (allSettingsSame_ != null && allSettingsSame_.isSelected())
        valueString = generalSettingsPanel_.getDilutionFactor().getText();
      else
        valueString = standardsSettings_.get(expName).getDilutionFactor().getText();
      Double value = 1d;
      if (valueString!=null && valueString.length()>0){
        value = new Double(valueString);
      }
      if (value<1)
        throw new AbsoluteSettingsInputException("The dilution factor must not be smaller than 1");
      dilutionFactors.put(expName, value);
    }
    for (String expName : standardsSettings_.keySet()){
    	if (generalSettingsPanel_ != null)
    	{
    		for (String esName : generalSettingsPanel_.getExtVolumeSettings().keySet()){
          Hashtable<String,VolumeConcVO> expsHash = new Hashtable<String,VolumeConcVO>();
          if (esStandards.containsKey(esName))
            expsHash = esStandards.get(esName);
          VolumeConcVO value = getVolConcVO(expName, esName, false);
          expsHash.put(expName, value);
          esStandards.put(esName, expsHash);
        }
        for (String isName : generalSettingsPanel_.getIntVolumeSettings().keySet()){
          Hashtable<String,VolumeConcVO> expsHash = new Hashtable<String,VolumeConcVO>();
          if (isStandards.containsKey(isName))
            expsHash = isStandards.get(isName);
          VolumeConcVO value = this.getVolConcVO(expName, isName, true);
          expsHash.put(expName, value);
          isStandards.put(isName, expsHash);
        }
    	}
    }
    return new LipidClassSettingVO(dilutionFactors,esStandards,isStandards);
  }
  
  private VolumeConcVO getVolConcVO(String expName, String standName, boolean isIS) throws AbsoluteSettingsInputException{
    VolumeConcVO value = null;
    try{
      if (allSettingsSame_.isSelected()){
        if (allStandardsSame_.isSelected())
          value = generalSettingsPanel_.getSettings(isIS);
        else
          value = generalSettingsPanel_.getSettings(standName,isIS);
      }else{
        if (getAllStandardsSame().isSelected())
          value = standardsSettings_.get(expName).getSettings(isIS);
        else
          value = standardsSettings_.get(expName).getSettings(standName,isIS); 
      }
    }catch (AbsoluteSettingsInputException ex){
      String standType = "external";
      if (isIS)
        standType = "internal";
      if (allSettingsSame_.isSelected()){
        if (allStandardsSame_.isSelected())
          throw new AbsoluteSettingsInputException("You have to enter volume and concentration for your "+standType+" standard!");
        else
          throw new AbsoluteSettingsInputException("You have to enter volume and concentration for your "+standName+" standard!");
      } else {
        if (allStandardsSame_.isSelected())
          throw new AbsoluteSettingsInputException("You have to enter volume and concentration for your "+standType+" standard (e.g. "+expName+" is not filled out)!");
        else
          throw new AbsoluteSettingsInputException("You have to enter volume and concentration for your "+standName+" standard (e.g. "+expName+" is not filled out)!");
      }
    }
    return value;
  }
  
  public String getChosenClass()
  {
  	if (alternativeClasses_ == null)
  	{
  		return className_;
  	}
  	return (String)alternativeClasses_.getSelectedItem();
  }
  
  public boolean areStandardsAvailable()
  {
    return 	!extractor_.getISNames(className_).isEmpty() || 
    				!extractor_.getESNames(className_).isEmpty();
  }
   
}
