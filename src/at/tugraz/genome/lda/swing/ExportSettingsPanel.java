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
import java.awt.event.ItemEvent;

import javax.swing.ButtonGroup;
import javax.swing.InputVerifier;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.vos.ExportOptionsVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ExportSettingsPanel extends JDialog implements ActionListener
{
  
  private static final long serialVersionUID = 1195564134976647401L;
  private ActionListener parent_;
  
  protected JCheckBox exportDeviation_;
  private JRadioButton standardDeviation_;
  private JTextField multStandardDeviation_;
  private JRadioButton standardErrorMean_;
  private JRadioButton columnAnalyte_;
  private JRadioButton columnExperiment_;
  protected JCheckBox exportRT_;
  protected JCheckBox exportRTDev_;
  /** radio button indicating that the analytes shall be exported on the species level*/
  private JRadioButton speciesLevel_;
  /** radio button indicating that the analytes shall be exported on the chain level*/
  private JRadioButton chainLevel_;
  /** radio button indicating that the analytes shall be exported on the position level*/
  private JRadioButton positionLevel_;

  
  private final static String CHANGE_SELECTION_STATUS = "changeSelectionStatus";
  public final static String CHANGE_RT_SELECTION_STATUS = "changeRTSelectionStatus";
  
  public ExportSettingsPanel(boolean isGrouped, ActionListener parent){
    parent_ = parent;
    
    setLocation(380,240);
    setLayout(new GridBagLayout());
    if (isGrouped){
      exportDeviation_ = new JCheckBox("export deviation value"); 
      exportDeviation_.setActionCommand(CHANGE_SELECTION_STATUS);
      exportDeviation_.addActionListener(this);
      exportDeviation_.setSelected(false);
      exportDeviation_.setEnabled(isGrouped);
      exportDeviation_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_DEVIATION);
      this.add(exportDeviation_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
      ButtonGroup deviationGroup = new ButtonGroup();
      standardDeviation_ = new JRadioButton("standard deviation");
      standardDeviation_.setSelected(true);
      standardDeviation_.setEnabled(false);
      standardDeviation_.addItemListener(new SelectionItemListener("ChangeSD"));
      deviationGroup.add(standardDeviation_);
      standardDeviation_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_STANDARD_DEVIATION);
      this.add(standardDeviation_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
      multStandardDeviation_ = new JTextField(3);
      multStandardDeviation_.setText("1.0");
      multStandardDeviation_.setEnabled(false);
      multStandardDeviation_.setInputVerifier(new DoubleVerifier());
      multStandardDeviation_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_STANDARD_DEVIATION_WHICH);
      this.add(multStandardDeviation_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
      standardErrorMean_ = new JRadioButton("standard error mean");
      standardErrorMean_.setEnabled(false);
      standardErrorMean_.addItemListener(new SelectionItemListener("ChangeSD"));
      deviationGroup.add(standardErrorMean_);
      standardErrorMean_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_STANDARD_ERROR);
      this.add(standardErrorMean_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    }
    ButtonGroup columnGroup = new ButtonGroup();
    columnAnalyte_ = new JRadioButton("analytes in column");
    columnAnalyte_.setSelected(true);
    columnGroup.add(columnAnalyte_);
    columnAnalyte_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_COLUMN_ANALYTE);
    this.add(columnAnalyte_,new GridBagConstraints(0, 2, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));    
    columnExperiment_ = new JRadioButton("experiments in column");
    columnExperiment_.setSelected(false);
    columnGroup.add(columnExperiment_);
    columnExperiment_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_COLUMN_EXPERIMENT);
    this.add(columnExperiment_,new GridBagConstraints(0, 3, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    exportRT_ = new JCheckBox("export retention-time"); 
    exportRT_.setActionCommand(CHANGE_RT_SELECTION_STATUS);
    exportRT_.addActionListener(this);
    exportRT_.setSelected(false);
    exportRT_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_RT);  
    this.add(exportRT_,new GridBagConstraints(0, 4, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    if (isGrouped){
      exportRTDev_ = new JCheckBox("export RT-stdev"); 
      exportRTDev_.addActionListener(this);
      exportRTDev_.setSelected(false);
      exportRTDev_.setEnabled(false);
      exportRTDev_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_RT_SD);  
      this.add(exportRTDev_,new GridBagConstraints(0, 5, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));      
    }
    ButtonGroup speciesGroup = new ButtonGroup();
    speciesLevel_ = new JRadioButton("species level");
    speciesLevel_.setSelected(true);
    speciesGroup.add(speciesLevel_);
    speciesLevel_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_SPECIES);
    this.add(speciesLevel_,new GridBagConstraints(0, 6, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));    
    chainLevel_ = new JRadioButton("chain level");
    chainLevel_.setSelected(false);
    speciesGroup.add(chainLevel_);
    chainLevel_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_CHAIN);
    this.add(chainLevel_,new GridBagConstraints(0, 7, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    positionLevel_ = new JRadioButton("position level");
    positionLevel_.setSelected(false);
    speciesGroup.add(positionLevel_);
    positionLevel_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_POSITION);
    this.add(positionLevel_,new GridBagConstraints(0, 8, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));

    JButton button = new JButton("OK");
    button.setActionCommand("AcceptExportSettings");
    button.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    this.add(button,new GridBagConstraints(0, 9, 3, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    button.addActionListener(parent_);
    setVisible(false);
    setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
    pack(); 
  }

  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase(CHANGE_SELECTION_STATUS)){
      if (exportDeviation_.isSelected()){
        standardDeviation_.setEnabled(true);
        multStandardDeviation_.setEnabled(true);
        standardErrorMean_.setEnabled(true);
      }else{
        standardDeviation_.setEnabled(false);
        multStandardDeviation_.setEnabled(false);
        standardErrorMean_.setEnabled(false);       
      }
    } else if (e.getActionCommand().equalsIgnoreCase(CHANGE_RT_SELECTION_STATUS)){
      if (exportRTDev_!=null){
        if (exportRT_.isSelected() ){
          exportRTDev_.setEnabled(true);
        }else{
          exportRTDev_.setSelected(false);
          exportRTDev_.setEnabled(false);
        }
      }
    }
  }
  
  private boolean verifySDFactInput(){
    try{
      Double.parseDouble(multStandardDeviation_.getText());
      return true;
    } catch (NumberFormatException nfx){
      new WarningMessage(new JFrame(), "Error", "SD value input must be in double format (xxx.xxx) and not "+multStandardDeviation_.getText());
      return false;
    }    
  }
  
  private class SelectionItemListener implements java.awt.event.ItemListener
  {
    private String m_ctrl;
    
    public SelectionItemListener(String ctrl){
      m_ctrl = ctrl;
    }

    public void itemStateChanged(ItemEvent e)  {
      if (m_ctrl.equalsIgnoreCase("ChangeSD")){
        if (e.getStateChange()==ItemEvent.SELECTED){
          if (standardDeviation_.isSelected())
            multStandardDeviation_.setEnabled(true);
          else
            multStandardDeviation_.setEnabled(false);
        }
      }      
    }
  }
  
  private class DoubleVerifier extends InputVerifier{

    public boolean verify(JComponent input)
    {
      return verifySDFactInput();
    } 
  }
  
  public ExportOptionsVO getSettings(){
    int exportType = ExportOptionsVO.EXPORT_NO_DEVIATION;
    String variationValue = null;
    if (this.exportDeviation_!=null && this.exportDeviation_.isSelected()){
      if (standardDeviation_.isSelected()){
        exportType = ExportOptionsVO.EXPORT_SD_DEVIATION;
        variationValue = multStandardDeviation_.getText();
      }else
        exportType = ExportOptionsVO.EXPORT_SD_ERROR;
    }
    boolean exportRTDev = false;
    if (exportRTDev_!=null && exportRTDev_.isSelected())
      exportRTDev = true;
    short speciesType = LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES;
    if (speciesLevel_.isSelected())
      speciesType = LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES;
    else if (chainLevel_.isSelected())
      speciesType = LipidomicsConstants.EXPORT_ANALYTE_TYPE_CHAIN;
    else if (positionLevel_.isSelected())
      speciesType = LipidomicsConstants.EXPORT_ANALYTE_TYPE_POSITION;    
    return new ExportOptionsVO(exportType,variationValue,columnAnalyte_.isSelected(),exportRT_.isSelected(),exportRTDev,6,speciesType);
  }
}
