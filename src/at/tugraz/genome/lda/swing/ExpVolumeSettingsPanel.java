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

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.exception.AbsoluteSettingsInputException;
import at.tugraz.genome.lda.vos.ProbeVolConcVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ExpVolumeSettingsPanel extends JPanel implements ActionListener
{
  private static final long serialVersionUID = 9166132661566566048L;

  private String expName_;
  private PhysicalUnitInput probeVolume_;
  private PhysicalUnitInput endVolume_;
  private PhysicalUnitInput sampleWeight_;
  private PhysicalUnitInput proteinConc_;
  private PhysicalUnitInput neutralConc_;
  private ExpVolumeListener parent_;
  
  public ExpVolumeSettingsPanel(String expName, ExpVolumeListener parent){
    expName_ = expName;
    parent_ = parent;
    this.initComponents();
  }
  
  private void initComponents(){
    this.setLayout(new GridBagLayout());
    JLabel label = new JLabel("Sample volume: ");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_SAMPLE_VOLUME);
    this.add(label,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    probeVolume_ = new PhysicalUnitInput("m",TooltipTexts.STATISTICS_ABS_SAMPLE_VOLUME);
    this.add(probeVolume_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    label = new JLabel("L");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_SAMPLE_VOLUME);
    this.add(label,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));  
    label = new JLabel("End volume: ");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_END_VOLUME);
    this.add(label,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    endVolume_ = new PhysicalUnitInput("\u03BC",TooltipTexts.STATISTICS_ABS_END_VOLUME);
    this.add(endVolume_,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    label = new JLabel("L");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_END_VOLUME);
    this.add(label,new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));
    JButton button = new JButton("Apply to all");
    button.addActionListener(this);
    button.setActionCommand("applyVolumeSettingsAll");
    button.setFont(button.getFont().deriveFont(10f));
    button.setMargin(new Insets(1,5,1,5));
    button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_ALL);
    this.add(button,new GridBagConstraints(6, 0, 1, 2, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
    button = new JButton("Apply to group");
    button.addActionListener(this);
    button.setActionCommand("applyVolumeSettingsGroup");
    button.setFont(button.getFont().deriveFont(10f));
    button.setMargin(new Insets(1,5,1,5));
    button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_GROUP);
    this.add(button,new GridBagConstraints(7, 0, 1, 2, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
    label = new JLabel("Sample weight: ");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_SAMPLE_WEIGHT);
    this.add(label,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    sampleWeight_ = new PhysicalUnitInput("m",TooltipTexts.STATISTICS_ABS_SAMPLE_WEIGHT);
    this.add(sampleWeight_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    label = new JLabel("g");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_SAMPLE_WEIGHT);
    this.add(label,new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));  
    
    label = new JLabel("Protein conc.: ");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_PROTEIN_CONC);
    this.add(label,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    proteinConc_ = new PhysicalUnitInput("m",TooltipTexts.STATISTICS_ABS_PROTEIN_CONC);
    this.add(proteinConc_,new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0)); 
    label = new JLabel("g/L");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_PROTEIN_CONC);
    this.add(label,new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));
    label = new JLabel("Neutral lipid con.: ");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_LIPID_CONC);
    this.add(label,new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    neutralConc_ = new PhysicalUnitInput("m",TooltipTexts.STATISTICS_ABS_LIPID_CONC);
    this.add(neutralConc_,new GridBagConstraints(4, 2, 1, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    label = new JLabel("g/L");
    label.setToolTipText(TooltipTexts.STATISTICS_ABS_LIPID_CONC);
    this.add(label,new GridBagConstraints(5, 2, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 8), 0, 0));
    button = new JButton("Apply to all");
    button.addActionListener(this);
    button.setActionCommand("applyProteinLipidSettingsAll");
    button.setFont(button.getFont().deriveFont(10f));
    button.setMargin(new Insets(1,5,1,5));
    button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_ALL);
    this.add(button,new GridBagConstraints(6, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
    button = new JButton("Apply to group");
    button.addActionListener(this);
    button.setActionCommand("applyProteinLipidSettingsGroup");
    button.setFont(button.getFont().deriveFont(10f));
    button.setMargin(new Insets(1,5,1,5));
    button.setToolTipText(TooltipTexts.STATISTICS_ABS_APPLY_TO_GROUP);
    this.add(button,new GridBagConstraints(7, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 0, 0));
  }
  
//  public ExpVolumeSettingVO getEnteredValues(){
//    ExpVolumeSettingVO vo = new ExpVolumeSettingVO(expName_,probeVolume_.getInputValue(), 
//        sampleWeight_.getInputValue(), (String)sampleWeight_.getUnitMagnitude().getSelectedItem(), 
//        endVolume_.getInputValue(), proteinConc_.getInputValue(), neutralConc_.getInputValue());
//    return vo;
//  }

  public void actionPerformed(ActionEvent e)
  {
    String command = e.getActionCommand();
    if (command.equalsIgnoreCase("applyVolumeSettingsAll")){
      parent_.applyVolumeSettingsToAll(expName_);
    }
    if (command.equalsIgnoreCase("applyVolumeSettingsGroup")){
      parent_.applyVolumeSettingsToGroup(expName_);
    }
    if (command.equalsIgnoreCase("applyProteinLipidSettingsAll")){
      parent_.applyProteinLipidSettingsToAll(expName_);
    }
    if (command.equalsIgnoreCase("applyProteinLipidSettingsGroup")){
      parent_.applyProteinLipidSettingsToGroup(expName_);
    }
  }
  
//  public void setVolumeInputValues(ExpVolumeSettingVO settingVO){
//    probeVolume_.setInputValue(settingVO.getProbeVolume());
//    endVolume_.setInputValue(settingVO.getEndVolume());
//  }
  
  public void setVolumeInputValues(ExpVolumeSettingsPanel otherPanel){
    probeVolume_.setInputValue(otherPanel.probeVolume_);
    endVolume_.setInputValue(otherPanel.endVolume_);
    sampleWeight_.setInputValue(otherPanel.sampleWeight_);
  }
  
  public void setProteinLipidContentValues(ExpVolumeSettingsPanel otherPanel){
    proteinConc_.setInputValue(otherPanel.proteinConc_);
    neutralConc_.setInputValue(otherPanel.neutralConc_);
  }

  public PhysicalUnitInput getProbeVolume()
  {
    return probeVolume_;
  }

  public PhysicalUnitInput getEndVolume()
  {
    return endVolume_;
  }

  public PhysicalUnitInput getSampleWeight()
  {
    return sampleWeight_;
  }
  
  public PhysicalUnitInput getProteinConc()
  {
    return proteinConc_;
  }

  public PhysicalUnitInput getNeutralConc()
  {
    return neutralConc_;
  }
  
  public void cleanSettings(){
    probeVolume_.setDefault();
    endVolume_.setDefault();
    sampleWeight_.setDefault();
    proteinConc_.setDefault();
    neutralConc_.setDefault();
  }
  
  public ProbeVolConcVO getSettingsVO() throws AbsoluteSettingsInputException{
    Double probeVol = probeVolume_.getInputValue();
    if (probeVol==null)
      throw new AbsoluteSettingsInputException("For the absolute quantitation the probe volume is mandatory");
    Double endVol = endVolume_.getInputValue();
    if (endVol==null)
      throw new AbsoluteSettingsInputException("For the absolute quantitation the end volume is mandatory");
    return new ProbeVolConcVO(probeVol,endVol,sampleWeight_.getInputValue(),
        proteinConc_.getInputValue(),neutralConc_.getInputValue());
  }
}
