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

package at.tugraz.genome.lda;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.Border;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;
import javax.swing.filechooser.FileNameExtensionFilter;

import at.tugraz.genome.lda.swing.AbsoluteQuantSettingsPanel;
import at.tugraz.genome.lda.swing.CutoffSettingsPanel;
import at.tugraz.genome.lda.swing.GroupsPanel;
import at.tugraz.genome.lda.verifier.DoubleVerifier;

/**
 * Singleton class for the quant file JPanel
 * 
 * @author Leonida M. Lamp
 *
 */
public class QuantFileJPanel
{
	private static QuantFileJPanel single_instance = null;
	private static final int RESULTS_TABLE_WIDTH = 950;
	
	public static synchronized QuantFileJPanel getInstance()
  {
      if (single_instance == null)
          single_instance = new QuantFileJPanel();

      return single_instance;
  }
	
	
	
	private void createResultsMenu(){
    JPanel resultsMenu = new JPanel();
    resultsMenu.setLayout(new GridBagLayout());
    
    
    JPanel buttonPanel = new JPanel();
    buttonPanel.setLayout(new GridBagLayout()); 
    buttonPanel.setPreferredSize(new Dimension(RESULTS_TABLE_WIDTH, 35));
    
    resultsMenu.add(buttonPanel,new GridBagConstraints(0, 0, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));    
    
//    jButtonResultFilesOpen_ = new JButton("Add Files", addFilesIcon_);
//    jButtonResultFilesOpen_.addActionListener(this);
//    jButtonResultFilesOpen_.setActionCommand("showResultsFilesChooser");
//    jButtonResultFilesOpen_.setToolTipText(TooltipTexts.STATISTICS_ADD_RESULT_FILES);
//    buttonPanel.add(jButtonResultFilesOpen_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, -2, 2, 0), 0, 0));
//    
//    jButtonResultsDirOpen_ = new JButton("Add Dir",addFolderIcon_);
//    jButtonResultsDirOpen_.addActionListener(this);
//    jButtonResultsDirOpen_.setActionCommand("showResultsDirChooser");
//    jButtonResultsDirOpen_.setToolTipText(TooltipTexts.STATISTICS_ADD_DIR);
//    buttonPanel.add(jButtonResultsDirOpen_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0));
//    
//    jButtonResultFilesRemove_ = new JButton("Remove",removeFilesIcon_);
//    jButtonResultFilesRemove_.addActionListener(this);
//    jButtonResultFilesRemove_.setActionCommand("removeResultFiles");
//    jButtonResultFilesRemove_.setToolTipText(TooltipTexts.STATISTICS_REMOVE_SELECTION);
//    buttonPanel.add(jButtonResultFilesRemove_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0));    
//    
//    jButtonResultFilesClean_ = new JButton("Remove all",removeFilesIcon_);
//    jButtonResultFilesClean_.addActionListener(this);
//    jButtonResultFilesClean_.setActionCommand("removeAllResultFiles");
//    jButtonResultFilesClean_.setToolTipText(TooltipTexts.STATISTICS_REMOVE_ALL_SELECTION);
//    buttonPanel.add(jButtonResultFilesClean_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0));
//    
//    jButtonResultAddToGroup_ = new JButton("Add to group",addToGroupIcon_);
//    jButtonResultAddToGroup_.addActionListener(this);
//    jButtonResultAddToGroup_.setActionCommand("addToGroup");
//    jButtonResultAddToGroup_.setToolTipText(TooltipTexts.STATISTICS_ADD_TO_GROUP);
//    buttonPanel.add(jButtonResultAddToGroup_,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0));    
//    
//    jButtonOmegaMassList_ = new JButton("MassList Wizard", omegaIcon_);
//    jButtonOmegaMassList_.addActionListener(this);
//    jButtonOmegaMassList_.setActionCommand("createOmegaMassList");
//    jButtonOmegaMassList_.setToolTipText(TooltipTexts.STATISTICS_CREATE_MASSLIST);
//    buttonPanel.add(jButtonOmegaMassList_,new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 250, 0, 0), 0, 0));
//
//    analysisSelectionTablePanel_ = new JPanel();
//    this.generateResultsAnalysisTablePane(new String[0][0]);
////    resultsMenu_.add(analysisSelectionTablePanel_);
//    resultsMenu.add(analysisSelectionTablePanel_,new GridBagConstraints(0, 1, 5, 1, 0.0, 0.0
//        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));  
//    
//    JPanel statisticalAnalysisSettingsPanel = new JPanel();
//    statisticalAnalysisSettingsPanel.setLayout(new GridBagLayout()); 
//    
//    //frameHeight_+5 is the resulting dimension when analysisTablePane is added to analysisSelectionTablePanel_
//    statisticalAnalysisSettingsPanel.setPreferredSize(new Dimension(FRAME_HEIGHT+5, 350));
//    TitledBorder title;
//    Border loweredEtched = BorderFactory.createEtchedBorder(EtchedBorder.LOWERED);
//    title = BorderFactory.createTitledBorder(loweredEtched,"Settings for Statistical Analysis");
//    statisticalAnalysisSettingsPanel.setBorder(title);
//    
//    groupsPanel_= new GroupsPanel();
//    resultsMenu.add(groupsPanel_ ,new GridBagConstraints(0, 2, 5, 1, 0.0, 0.0
//        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
//    
//    cutoffSettingsPanel_ = new CutoffSettingsPanel(TooltipTexts.STATISTICS_ADD_CUTOFF_SETTINGS);
//    resultsMenu.add(cutoffSettingsPanel_ ,new GridBagConstraints(0, 3, 5, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));  
//    
//    quantSettingsPanel_ = new AbsoluteQuantSettingsPanel (groupsPanel_);
//    resultsMenu.add(quantSettingsPanel_ ,new GridBagConstraints(0, 4, 5, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0)); 
//    
//    resultsMenu.add(statisticalAnalysisSettingsPanel,new GridBagConstraints(0, 5, 5, 1, 0.0, 0.0
//        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(30, 0, 0, 0), 0, 0));
//    
//    JPanel groupRtPanel = new JPanel();
//    groupRtPanel.setLayout(new GridBagLayout());  
//    JLabel sepRtLabel = new JLabel("Show hits with different RT separately: ");
//    sepRtLabel.setToolTipText(TooltipTexts.STATISTICS_SEPARATE_RT);
//    groupRtPanel.add(sepRtLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
  }
	
	
	
}
