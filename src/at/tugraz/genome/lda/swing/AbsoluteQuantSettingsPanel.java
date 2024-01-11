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
import java.util.Hashtable;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.analysis.ComparativeNameExtractor;
import at.tugraz.genome.lda.exception.AbsoluteSettingsInputException;
import at.tugraz.genome.lda.vos.AbsoluteSettingsVO;
import at.tugraz.genome.lda.vos.LipidClassSettingVO;
import at.tugraz.genome.lda.vos.ProbeVolConcVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class AbsoluteQuantSettingsPanel extends JPanel implements ExpVolumeListener
{

  private static final long serialVersionUID = -7145044539526127998L;
  
  private JPanel inputPanel_;
  private JTabbedPane lipidClassesTabs_;
  private JTabbedPane experimentTabs_;
  
  private Hashtable<String,ExpVolumeSettingsPanel> volumeSettings_;
  private Hashtable<String,LipidClassSettingsPanel> classSettings_;
  private GroupsPanel groupsPanel_;
  
  private ComparativeNameExtractor extractor_;
  
  public AbsoluteQuantSettingsPanel(){
    super();
    this.initComponents();    
  }
  
  public AbsoluteQuantSettingsPanel(GroupsPanel groupsPanel){
    this();
    this.groupsPanel_ = groupsPanel;
  }

  
  private void initComponents(){   
    inputPanel_ = new JPanel();
    inputPanel_.setLayout(new GridBagLayout());
    this.add(inputPanel_);  
  }
  
  public void showSettingsPanel(ComparativeNameExtractor extractor){
    this.extractor_ = extractor;
    this.initTabbedPanes();
    this.invalidate();
    this.updateUI();
  }
  
  public void hideSettingsPanel(){
    inputPanel_.removeAll();
    this.invalidate();
    this.updateUI();
  }
  
  
  private void initTabbedPanes(){
    inputPanel_.removeAll();
    lipidClassesTabs_ = new JTabbedPane();
    experimentTabs_ = new JTabbedPane();
    classSettings_ = new Hashtable<String,LipidClassSettingsPanel>();
    volumeSettings_ = new Hashtable<String,ExpVolumeSettingsPanel>();
    int count = 0;
    for (String lipidClass : extractor_.getLipidClasses()){    
      LipidClassSettingsPanel classSettings = new LipidClassSettingsPanel(lipidClass,extractor_,groupsPanel_);
      classSettings_.put(lipidClass, classSettings);   
      lipidClassesTabs_.addTab(lipidClass,classSettings);
      lipidClassesTabs_.setToolTipTextAt(count,TooltipTexts.STATISTICS_ABS_TAB_CLASS+lipidClass+"</html>");
      count++;
    }
    inputPanel_.add(lipidClassesTabs_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    count = 0;
    for (String expName : this.extractor_.getExpNamesInSequence()){
      ExpVolumeSettingsPanel expVolSettings = new ExpVolumeSettingsPanel(expName,this);
      volumeSettings_.put(expName, expVolSettings);
      experimentTabs_.addTab(extractor_.getExpNames().get(expName),expVolSettings);
      experimentTabs_.setToolTipTextAt(count, TooltipTexts.STATISTICS_ABS_TAB_EXP+expName+"</html>");
      count++;
    }
    inputPanel_.add(experimentTabs_,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
  }

  public void actionPerformed(ActionEvent e)
  {
    // TODO Auto-generated method stub
    
  }

  public void applyVolumeSettingsToAll(String expirementName)
  {
    this.applySettingsToAll(expirementName, true);
  }
  public void applyProteinLipidSettingsToAll(String expirementName)
  {
    this.applySettingsToAll(expirementName, false);
  }
  
  public void applyVolumeSettingsToGroup(String expirementName)
  {
    this.applySettingsToGroup(expirementName,true);
  }
  
  public void applyProteinLipidSettingsToGroup(String expirementName)
  {
    this.applySettingsToGroup(expirementName,false);
  }
  
  public void applySettingsToAll(String expirementName, boolean volume)
  {
    ExpVolumeSettingsPanel panel = volumeSettings_.get(expirementName);
    for (String expName : this.extractor_.getExpNamesInSequence()){
      if (!expName.equalsIgnoreCase(expirementName)){
        if (volume)
          volumeSettings_.get(expName).setVolumeInputValues(panel);
        else
          volumeSettings_.get(expName).setProteinLipidContentValues(panel);
      }  
    }
    this.invalidate();
    this.updateUI();
  }
  
  private void applySettingsToGroup(String expirementName,boolean volume)
  {
    ExpVolumeSettingsPanel panel = volumeSettings_.get(expirementName);
    Set<String> exps = groupsPanel_.getExpsOfGroupOneExpBelongsTo(expirementName,this.extractor_.getExpNamesInSequence());
    if (exps.size()>0){
      for (String expName : exps){
        if (!expName.equalsIgnoreCase(expirementName)){
          if (volume)
            volumeSettings_.get(expName).setVolumeInputValues(panel);
          else
            volumeSettings_.get(expName).setProteinLipidContentValues(panel);
        }  
      }
      this.invalidate();
      this.updateUI();
    }else{
      new WarningMessage(new JFrame(), "Warning", "This experiment does not belong to any group");
    }
  }
  
  public Hashtable<String,ExpVolumeSettingsPanel> getVolumeSettings(){
    return this.volumeSettings_;
  }

  public Hashtable<String,LipidClassSettingsPanel> getClassSettings()
  {
    return classSettings_;
  }
  
  public void cleanVolumeSettings(){
    for (ExpVolumeSettingsPanel volSets : volumeSettings_.values()){
      volSets.cleanSettings();
    }
  }
  
  public void cleanClassSettings(){
    for (LipidClassSettingsPanel classSet : classSettings_.values())
      classSet.cleanSettings();
  }
  
  public AbsoluteSettingsVO getSettingsVO() throws AbsoluteSettingsInputException{
    Hashtable<String,ProbeVolConcVO> probeVols = new Hashtable<String,ProbeVolConcVO>();
    Hashtable<String,LipidClassSettingVO> standSet = new Hashtable<String,LipidClassSettingVO>();
    boolean foundProteinConcOnce = false;
    boolean foundNeutralLipidConcOnce = false;
    boolean foundSampleWeightOnce = false;
    boolean foundProteinConcAll = true;
    boolean foundNeutralLipidConcAll = true;
    boolean foundSampleWeightAll = true;
    for (String expName : volumeSettings_.keySet()){
      ProbeVolConcVO probeVol = volumeSettings_.get(expName).getSettingsVO();
      probeVols.put(expName,probeVol);
      if (probeVol.getProteinConc()==null)
        foundProteinConcAll = false;
      else
        foundProteinConcOnce = true;
      if (probeVol.getNeutralLipidConc()==null)
        foundNeutralLipidConcAll = false;
      else
        foundNeutralLipidConcOnce = true;
      if (probeVol.getSampleWeight()==null)
        foundSampleWeightAll = false;
      else
        foundSampleWeightOnce = true;
    }
    if (foundProteinConcOnce && !foundProteinConcAll)
      throw new AbsoluteSettingsInputException("The protein concentration has to be set for all experiments or for none");
    if (foundNeutralLipidConcOnce && !foundNeutralLipidConcAll)
      throw new AbsoluteSettingsInputException("The neutral lipid concentration has to be set for all experiments or for none");
    if (foundSampleWeightOnce && !foundSampleWeightAll)
      throw new AbsoluteSettingsInputException("The sample weight has to be set for all experiments or for none");
    
    Hashtable<String,String> chosenClassLookup = new Hashtable<String,String>();
    for (String className : classSettings_.keySet()){
    	String chosenClass = getChosenClassLookup(className);
      LipidClassSettingVO classVO = classSettings_.get(chosenClass).getSettings();
      standSet.put(className, classVO);
      chosenClassLookup.put(className, chosenClass);
    }
    AbsoluteSettingsVO settingsVO = new AbsoluteSettingsVO(probeVols,standSet,chosenClassLookup);
    return settingsVO;
  }
  
  public String getChosenClassLookup(String className)
  {
  	return classSettings_ == null ? className : classSettings_.get(className).getChosenClass();
  }
  
  
}
