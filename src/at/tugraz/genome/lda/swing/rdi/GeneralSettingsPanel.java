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

package at.tugraz.genome.lda.swing.rdi;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.filechooser.FileNameExtensionFilter;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.interfaces.rdi.GeneralSettingsPanelListener;
import at.tugraz.genome.lda.msn.FattyAcidsContainer;
import at.tugraz.genome.lda.msn.FragmentCalculator;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;
import at.tugraz.genome.lda.verifier.GeneralSettingsVerifier;
import at.tugraz.genome.lda.vos.rdi.GeneralSettingsVO;


/**
 * Panel of the RDI for entering the general information of the rules section
 * @author Juergen Hartler
 *
 */
public class GeneralSettingsPanel extends JPanel implements ActionListener
{

  private static final long serialVersionUID = 8237230657153326314L;
  
  /** Field for the amount of chains */  
  private JTextField chainsAmountField_;

  /** Field for the amount of alkyl chains */  
  private JTextField alkylChainsAmountField_;

  /** Field for the amount of alkenyl chains */  
  private JTextField alkenylChainsAmountField_;
  
  /** Field for the chains library */  
  private JTextField chainsLibraryField_;

  /** Field for the carbon atoms rule */  
  private JTextField carbonAtomsRuleField_;
  
  /** Field for the double bond rule */  
  private JTextField doubleBondRuleField_;
  
  /** Combo box for selection if single chain identifications are valid */
  private JComboBox<String> singleChainIdentification_;

  /** Field for the base chain cutoff */  
  private JTextField chainCutoffField_;

  /** Combo box for choosing the unit for the chain cutoff */
  private JComboBox<String> chainCutoffUnitCombo_;
  
  /** Field for the base peak cutoff */  
  private JTextField basePeakCutoffField_;
  
  /** Possibility List for unit combo */
  private String unitList_[] = {"  ","%", "\u2030"};
  
  /** General true/false selction for any combo box */
  private String[] trueFalseSeletion_ = {"true", "false"};

  /** Combo box for choosing the unit for base peak cutoff */
  private JComboBox<String> basePeakUnitCombo_;
  
  /** Field for the spectrum coverage */  
  private JTextField spectrumCoverageField_;

  /** Combo box for choosing the unit for spectrum coverage */
  private JComboBox<String> spectrumCoverageUnitCombo_;

  /** Combo box for the post processing retention time */
  private JComboBox<String> possibibilitiesPostProcessing_;
  
  /** Combo box for the parallel mode retention time */
  private JComboBox<String> possibibilitiesParallelMode_;
  
  /** Field for the retention time max deviation */  
  private JTextField rtMaxDevField_;
  
  /** Field for additional chain positions */  
  private JTextField addChainPositions_;
  
  /** MS1 identification before MSn */
  private final static String MS1_FIRST = "MS1 first";
  /** MSn identification first, then MS1 only for found MSn hits*/
  private final static String MSN_ONLY = "MSn only";
  /** MSn identification first, then MS1 identification for predicted retention times*/
  private final static String MSN_FIRST = "MSn first";
  
  /** Possibility List for MS identification orders */
  private final static String msIdCases_[] = {MS1_FIRST,MSN_ONLY, MSN_FIRST};

  /** Combo box for the parallel mode retention time */
  private JComboBox<String> msIdentificationOrder_;

  
  /** button for saving general setup for rules only */
  private JButton saveGeneralRules_;
  
  /** callback listener for updating depending graphical components */
  private GeneralSettingsPanelListener listener_;
  
  /** the default size of the text fields */
  private int DEFAULT_TEXT_INPUT_SIZE = 12;
  
  /** if default values for the chain cutoff shall be used */
  private final static String CHAIN_CUTOFF_DEFAULT = "default";
  
  /** the displayed text in front of the text field for alkyl chains */
  private final static String ALKYLCHAINS_LABEL = "Alkylchains #";
  /** the displayed text in front of the text field for alkenyl chains */
  private final static String ALKENYLCHAINS_LABEL = "Alkenylchains #";
  /** the displayed text in front of the text field for the chain cutoff */
  private final static String CHAINCUTOFF_LABEL = "Chain Cutoff";
  /** the displayed text in front of the text field for entering the retention time max deviation */
  private final static String RTMAXDEV_LABEL = "Retention Time Max Deviation";
  /** the displayed text in front of the text field for adding additional allowed chain positions */
  private final static String ADD_CHAIN_POSITIONS_LABEL = "Additional allowed positions";
  
  /**
   * constructor of this JPanel
   * @param vo a VO containing all required parameters that have to be set in the input components
   * @param listener the callback listener for updating depending graphical components
   */
  public GeneralSettingsPanel(GeneralSettingsVO vo, GeneralSettingsPanelListener listener) {
    this.listener_ = listener;
    init(vo);
  }
  
  /**
   * initiates the graphical input components of this JPanel
   * @param vo a VO containing all required parameters that have to be set in the input components
   */
  private void init(GeneralSettingsVO vo){
    setLayout(new GridBagLayout());
    
    JLabel spaceClass3 = new JLabel( "    " );
    add(spaceClass3, new GridBagConstraints (0, 0, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
    JLabel chainsAmount = new JLabel( "  Chains #  " );        
    add(chainsAmount, new GridBagConstraints (0, 1, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));  
    chainsAmountField_ = new JTextField();
    add(chainsAmountField_, new GridBagConstraints (1, 1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    chainsAmountField_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    String amountOfChains = null;
    if (vo.getAmountOfChains()!=null)
      amountOfChains = String.valueOf(vo.getAmountOfChains());
    chainsAmountField_.setText(amountOfChains);
    chainsAmountField_.setToolTipText("Enter the amount of chains here!");
    chainsAmountField_.setInputVerifier(new GeneralSettingsVerifier(this));
    chainsAmountField_.setHorizontalAlignment(JTextField.RIGHT);
    
    JLabel alkylChains = new JLabel( "  "+ALKYLCHAINS_LABEL+"  " );
    add(alkylChains, new GridBagConstraints (0, 2, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));  
    alkylChainsAmountField_ = new JTextField();
    add(alkylChainsAmountField_, new GridBagConstraints (1, 2, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    alkylChainsAmountField_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    String amountOfAlkylChains = null;
    if (vo.getAmountOfAlkylChains()!=null)
      amountOfAlkylChains = String.valueOf(vo.getAmountOfAlkylChains());
    alkylChainsAmountField_.setText(amountOfAlkylChains);
    alkylChainsAmountField_.setToolTipText("Enter the amount of alkyl chains here!");
    alkylChainsAmountField_.setInputVerifier(new GeneralSettingsVerifier(this));
    alkylChainsAmountField_.setHorizontalAlignment(JTextField.RIGHT);
    
    JLabel alkenylChains = new JLabel( "  "+ALKENYLCHAINS_LABEL+"  " );
    add(alkenylChains, new GridBagConstraints (0, 3, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));  
    alkenylChainsAmountField_ = new JTextField();
    add(alkenylChainsAmountField_, new GridBagConstraints (1, 3, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    alkenylChainsAmountField_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    String amountOfAlkenylChains = null;
    if (vo.getAmountOfAlkenylChains()!=null)
      amountOfAlkenylChains = String.valueOf(vo.getAmountOfAlkenylChains());
    alkenylChainsAmountField_.setText(amountOfAlkenylChains);
    alkenylChainsAmountField_.setToolTipText("Enter the amount of alkenyl chains here!");
    alkenylChainsAmountField_.setInputVerifier(new GeneralSettingsVerifier(this));
    alkenylChainsAmountField_.setHorizontalAlignment(JTextField.RIGHT);
    
    JLabel chainsLibrary = new JLabel( "  Chain Library  " );        
    add(chainsLibrary, new GridBagConstraints (0, 4, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));                    
    chainsLibraryField_ = new JTextField();
    add(chainsLibraryField_, new GridBagConstraints (1, 4, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    chainsLibraryField_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    chainsLibraryField_.setText(vo.getChainLibrary());
    chainsLibraryField_.setToolTipText("Enter the chain library here!");
    chainsLibraryField_.setInputVerifier(new GeneralSettingsVerifier(this));
    
    JButton selectButton = new JButton( "Select" );
    selectButton.setFont(selectButton.getFont().deriveFont(10f));
    selectButton.setMargin(new Insets(1,5,1,5));
    add(selectButton, new GridBagConstraints (2, 4, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    selectButton.setToolTipText("Opens the file system");
    selectButton.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        JFileChooser chooser = new JFileChooser("fattyAcids");        
        FileNameExtensionFilter filter = new FileNameExtensionFilter(
            "Excel-File", "xlsx");
        chooser.setFileFilter(filter);        
        int returnVal = chooser.showOpenDialog(new JFrame());
        if(returnVal == JFileChooser.APPROVE_OPTION) 
        {           
           chainsLibraryField_.setText( chooser.getSelectedFile().getName() );
        }
      }
    });
    
    JLabel carbonAtomsRule = new JLabel( "  Carbon Atoms Rule  " );        
    add(carbonAtomsRule, new GridBagConstraints (0, 5, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));                  
    carbonAtomsRuleField_ = new JTextField();
    add(carbonAtomsRuleField_, new GridBagConstraints (1, 5, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    carbonAtomsRuleField_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    carbonAtomsRuleField_.setText(vo.getCarbonAtomsRule());
    carbonAtomsRuleField_.setToolTipText("Enter the rule for carbon atoms here!");
    carbonAtomsRuleField_.setInputVerifier(new GeneralSettingsVerifier(this));    
    JButton test1Button = new JButton( "Test" );
    test1Button.setFont(selectButton.getFont().deriveFont(10f));
    test1Button.setMargin(new Insets(1,5,1,5));
    add(test1Button, new GridBagConstraints (2, 5, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    test1Button.setToolTipText("Shows the amount of carbon atoms for the entered rule");
    test1Button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {    
        int  cAtoms = checkCAtomsRule(true);
        if (cAtoms>-1)
          JOptionPane.showMessageDialog(new JFrame(), "Number of Carbon Atoms = " + cAtoms); 
      }
    }); 

    JLabel doubleBondRule = new JLabel( "  Double Bond Rule  " );        
    add(doubleBondRule, new GridBagConstraints (0, 6, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));                    
    doubleBondRuleField_ = new JTextField();
    add(doubleBondRuleField_, new GridBagConstraints (1, 6, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    doubleBondRuleField_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    doubleBondRuleField_.setText(vo.getDoubleBondsRule()); 
    doubleBondRuleField_.setToolTipText("Enter the rule for double bonds here!");
    doubleBondRuleField_.setInputVerifier(new GeneralSettingsVerifier(this));
    
    JButton test2Button = new JButton( "Test" );
    test2Button.setFont(selectButton.getFont().deriveFont(10f));
    test2Button.setMargin(new Insets(1,5,1,5));
    add(test2Button, new GridBagConstraints (2, 6, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    test2Button.setToolTipText("Shows the amount of double bonds for the entered rule");
    test2Button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {    
        int dbs = checkDoubleBondsRule(true);
        if (dbs>-1)
          JOptionPane.showMessageDialog(new JFrame(), "Number of Double Bonds = " + dbs); 
      }     
    }); 
    
    JLabel singleChainIdentification = new JLabel( "  Single Chain Identification " );        
    add(singleChainIdentification, new GridBagConstraints (0, 7, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));    
    singleChainIdentification_ = new JComboBox<String>(this.trueFalseSeletion_);  
    singleChainIdentification_.setToolTipText(TooltipTexts.RDI_GENERAL_SINGLECHAIN);
    if(vo.isAllowSingleChain()) singleChainIdentification_.setSelectedIndex(0);
    else singleChainIdentification_.setSelectedIndex(1);
    add(singleChainIdentification_, new GridBagConstraints (1, 7, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    singleChainIdentification_.setInputVerifier(new GeneralSettingsVerifier(this));

    JLabel chainCutoffRule = new JLabel( "  "+CHAINCUTOFF_LABEL+"  " );        
    add(chainCutoffRule, new GridBagConstraints (0, 8, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));                    
    chainCutoffField_ = new JTextField();
    add(chainCutoffField_, new GridBagConstraints (1, 8, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    chainCutoffField_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    String chainCutoffValue = vo.getChainCutoff();
    String ccMagnitude = " ";
    if (chainCutoffValue!=null && (chainCutoffValue.endsWith("%") || chainCutoffValue.endsWith("\u2030"))){
      ccMagnitude = chainCutoffValue.substring(chainCutoffValue.length()-1);
      chainCutoffValue = chainCutoffValue.substring(0,chainCutoffValue.length()-1);
    }
    if (chainCutoffValue!=null){
      try{
        double chainCutoff = Double.parseDouble(chainCutoffValue);
        if (chainCutoff<0) chainCutoffValue = CHAIN_CUTOFF_DEFAULT;
      } catch (NumberFormatException nfx){
        chainCutoffValue = CHAIN_CUTOFF_DEFAULT;
      }
    } else {
      chainCutoffValue = CHAIN_CUTOFF_DEFAULT;
    }
    if (chainCutoffValue==null || chainCutoffValue.equalsIgnoreCase(CHAIN_CUTOFF_DEFAULT))
      ccMagnitude = "%";
    chainCutoffField_.setText(chainCutoffValue);
    chainCutoffField_.setToolTipText(TooltipTexts.RDI_GENERAL_CHAINCUTOFF);
    chainCutoffField_.setInputVerifier(new GeneralSettingsVerifier(this));
    chainCutoffField_.setHorizontalAlignment(JTextField.RIGHT);    
    chainCutoffUnitCombo_ = new JComboBox<String>(unitList_); 
    chainCutoffUnitCombo_.setSelectedItem(ccMagnitude);
    chainCutoffUnitCombo_.setToolTipText("Choose a unit!");
    add(chainCutoffUnitCombo_, new GridBagConstraints (2, 8, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));  
    chainCutoffUnitCombo_.setInputVerifier(new GeneralSettingsVerifier(this));
        
    // second input column
    JLabel basePeakCutoffRule = new JLabel( "  Base Peak Cutoff  " );        
    add(basePeakCutoffRule, new GridBagConstraints (3, 1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));                    
    basePeakCutoffField_ = new JTextField();
    add(basePeakCutoffField_, new GridBagConstraints (4, 1, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    basePeakCutoffField_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    String basePeakValue = vo.getBasePeakCutoff();
    String bpMagnitude = " ";
    if (basePeakValue.endsWith("%") || basePeakValue.endsWith("\u2030")){
      bpMagnitude = basePeakValue.substring(basePeakValue.length()-1);
      basePeakValue = basePeakValue.substring(0,basePeakValue.length()-1);
    }
    basePeakCutoffField_.setText(basePeakValue);
    basePeakCutoffField_.setToolTipText("Enter the base peak cutoff here!");
    basePeakCutoffField_.setInputVerifier(new GeneralSettingsVerifier(this));
    basePeakCutoffField_.setHorizontalAlignment(JTextField.RIGHT); 
    basePeakUnitCombo_ = new JComboBox<String>(unitList_); 
    basePeakUnitCombo_.setSelectedItem(bpMagnitude);
    basePeakUnitCombo_.setToolTipText("Choose a unit!");
    add(basePeakUnitCombo_, new GridBagConstraints (5, 1, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));  
    basePeakUnitCombo_.setInputVerifier(new GeneralSettingsVerifier(this));
    
    JLabel spectrumCoverageRule = new JLabel( "  Spectrum Coverage  " );        
    add(spectrumCoverageRule, new GridBagConstraints (3, 2, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));                    
    spectrumCoverageField_ = new JTextField();
    add(spectrumCoverageField_, new GridBagConstraints (4, 2, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    spectrumCoverageField_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    String spCovValue = vo.getSpectrumCoverage();
    String spCovMagnitude = " ";
    if (spCovValue!=null && (spCovValue.endsWith("%") || spCovValue.endsWith("\u2030"))){
      spCovMagnitude = spCovValue.substring(spCovValue.length()-1);
      spCovValue = spCovValue.substring(0,spCovValue.length()-1);
    }
    if (spCovValue==null)spCovMagnitude = "%";    
    spectrumCoverageField_.setText(spCovValue);
    spectrumCoverageField_.setToolTipText("Enter the spectrum coverage here!");
    spectrumCoverageField_.setInputVerifier(new GeneralSettingsVerifier(this));
    spectrumCoverageField_.setHorizontalAlignment(JTextField.RIGHT); 
    spectrumCoverageUnitCombo_ = new JComboBox<String>(unitList_); 
    spectrumCoverageUnitCombo_.setSelectedItem(spCovMagnitude);       
    spectrumCoverageUnitCombo_.setToolTipText("Choose a unit!");
    add(spectrumCoverageUnitCombo_, new GridBagConstraints (5, 2, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));  
    spectrumCoverageUnitCombo_.setInputVerifier(new GeneralSettingsVerifier(this));
    
    JLabel retentionTimePostProcessing = new JLabel( "  Retention Time Post Processing " );        
    add(retentionTimePostProcessing, new GridBagConstraints (3, 3, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    String possibibilityListPostProcessing[] = {"true", "false"};
    possibibilitiesPostProcessing_ = new JComboBox<String>(possibibilityListPostProcessing);  
    possibibilitiesPostProcessing_.setToolTipText("Choose the retention time post processing!");
    if(vo.isRtPostProcessing())
    {
      possibibilitiesPostProcessing_.setSelectedIndex(0);
    }
    else  
    {
      possibibilitiesPostProcessing_.setSelectedIndex(1);
    }    
    add(possibibilitiesPostProcessing_, new GridBagConstraints (4, 3, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    
    JLabel retentionTimeParallelModel = new JLabel( "  Retention Time Parallel Series " );        
    add(retentionTimeParallelModel, new GridBagConstraints (3, 4, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));    
    String possibibilitiesListParallelMode[] = {"true", "false"};
    possibibilitiesParallelMode_ = new JComboBox<String>(possibibilitiesListParallelMode);   
    possibibilitiesParallelMode_.setToolTipText("Choose the retention time parallel series!");
    if(vo.isRtParallelSeries())
    {
      possibibilitiesParallelMode_.setSelectedIndex(0);
    }
    else  
    {
      possibibilitiesParallelMode_.setSelectedIndex(1);
    }    
    add(possibibilitiesParallelMode_, new GridBagConstraints (4, 4, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));  
    
    JLabel rtMaxDevLabel = new JLabel( "  "+RTMAXDEV_LABEL+"  " );        
    add(rtMaxDevLabel, new GridBagConstraints (3, 5, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));                    
    rtMaxDevField_ = new JTextField();
    add(rtMaxDevField_, new GridBagConstraints (4, 5, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    rtMaxDevField_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    String rtMaxDevValue = null;
    if (vo.getRtMaxDeviation()!=null) rtMaxDevValue = vo.getRtMaxDeviation().toString();
    rtMaxDevField_.setText(rtMaxDevValue);
    rtMaxDevField_.setToolTipText(TooltipTexts.RDI_GENERAL_RTMAXDEV);
    rtMaxDevField_.setHorizontalAlignment(JTextField.RIGHT); 
    JLabel unitLabel = new JLabel("min");
    add(unitLabel, new GridBagConstraints (5, 5, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));  

    JLabel msIdOrderLabel = new JLabel( "  MS Identification Order " );        
    add(msIdOrderLabel, new GridBagConstraints (3, 6, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));    
    msIdentificationOrder_ = new JComboBox<String>(msIdCases_);   
    msIdentificationOrder_.setToolTipText(TooltipTexts.RDI_GENERAL_MSIDORDER);
    if(vo.getMsIdentificationOrder() == RulesContainer.ORDER_MS1_FIRST)
      msIdentificationOrder_.setSelectedItem(MS1_FIRST);
    else if (vo.getMsIdentificationOrder() == RulesContainer.ORDER_MSN_ONLY)
      msIdentificationOrder_.setSelectedItem(MSN_ONLY);
    else if (vo.getMsIdentificationOrder() == RulesContainer.ORDER_MSN_FIRST)
      msIdentificationOrder_.setSelectedItem(MSN_FIRST);    
    add(msIdentificationOrder_, new GridBagConstraints (4, 6, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));  

    JLabel addChainPositions = new JLabel("  "+ADD_CHAIN_POSITIONS_LABEL+"  ");
    add(addChainPositions, new GridBagConstraints (3, 7, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    addChainPositions_  = new JTextField();
    addChainPositions_.setToolTipText(TooltipTexts.RDI_GENERAL_ADDPOSITIONS);
    add(addChainPositions_, new GridBagConstraints (4, 7, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));  

    addChainPositions_.setColumns( DEFAULT_TEXT_INPUT_SIZE );
    String addPoss = null;
    if (vo.getAddChainPositions()!=null && vo.getAddChainPositions()>0)
      addPoss = String.valueOf(vo.getAddChainPositions());
    addChainPositions_.setText(addPoss);
    addChainPositions_.setInputVerifier(new GeneralSettingsVerifier(this));
    addChainPositions_.setHorizontalAlignment(JTextField.RIGHT);

    
    saveGeneralRules_ = new JButton(" Save General Settings and Continue ");
    saveGeneralRules_.addActionListener(this);
    saveGeneralRules_.setActionCommand("saveGeneralRules");
   
    if(!listener_.areGeneralSettingsDependingFieldsVisible()){
      add(saveGeneralRules_, new GridBagConstraints (0, 9, 6, 1, 0, 0, 
          GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 0, 0, 0), 0, 0 )); 
    }

  }
  
  /**
   * checks the general entries for validity
   * @return true if all of them are OK
   */
  private boolean checkGeneralEntries(boolean verbose){
    int chainsAmount = 0;
    if (this.chainsAmountField_.getText() == null || this.chainsAmountField_.getText().length()<1){
      if (verbose)new WarningMessage(new JFrame(), "Error", "You have to enter a value for the \"Chains #\"");
      chainsAmountField_.requestFocus();
      return false;
    } else {
      boolean inputOK = this.checkIntValueInRange(chainsAmountField_, "Chains #", 1, 10, verbose);
      if (inputOK) chainsAmount = Integer.parseInt(chainsAmountField_.getText());
      else{
        chainsAmountField_.requestFocus();
        return false;
      }
    }
    int alkylChains = 0;
    if (alkylChainsAmountField_ != null && alkylChainsAmountField_.getText().trim().length()>0){
      boolean inputOK = this.checkIntValueInRange(alkylChainsAmountField_, ALKYLCHAINS_LABEL, 0, chainsAmount, verbose);
      if (inputOK) alkylChains = Integer.parseInt(alkylChainsAmountField_.getText());
      else{
        alkylChainsAmountField_.requestFocus();
        return false;
      }
    }
    int alkenylChains = 0;
    if (alkenylChainsAmountField_!=null && alkenylChainsAmountField_.getText().trim().length()>0){
      boolean inputOK = this.checkIntValueInRange(alkenylChainsAmountField_, ALKENYLCHAINS_LABEL, 0, chainsAmount, verbose);
      if (inputOK) alkenylChains = Integer.parseInt(alkenylChainsAmountField_.getText());
      else{
        alkenylChainsAmountField_.requestFocus();
        return false;
      }      
    }
    if ((alkylChains+alkenylChains)>chainsAmount){
      if (verbose)new WarningMessage(new JFrame(), "Error", "The sum of  \""+ALKYLCHAINS_LABEL+"\"+\""+ALKENYLCHAINS_LABEL+"\" must not be greater than \"Chains #\"");
      if (alkenylChains>alkylChains) alkenylChainsAmountField_.requestFocus();
      else alkylChainsAmountField_.requestFocus();
      return false;
    }
    if (this.chainsLibraryField_.getText() == null || this.chainsLibraryField_.getText().length()<1){
      if (verbose) new WarningMessage(new JFrame(), "Error", "You have to enter a value for the \"Chain Library\"");
      chainsLibraryField_.requestFocus();
      return false;
    } else {
      try {
        FattyAcidsContainer.getAllFattyAcidChains(chainsLibraryField_.getText());
      }
      catch (RulesException | NoRuleException | IOException e) {
        if (verbose) new WarningMessage(new JFrame(), "Error", "There is something wrong with your selected chain library: "+e.getMessage());
        chainsLibraryField_.requestFocus();
        return false;
      }
    }
    if (carbonAtomsRuleField_.getText() == null || this.carbonAtomsRuleField_.getText().length()<1){
      if (verbose) new WarningMessage(new JFrame(), "Error", "You have to enter a value for the \"Carbon Atoms Rule\"");
      carbonAtomsRuleField_.requestFocus();
      return false;
    }else{
      if (checkCAtomsRule(verbose)<0){
        carbonAtomsRuleField_.requestFocus();
        return false;
      }
    }
    if (doubleBondRuleField_.getText() == null || this.doubleBondRuleField_.getText().length()<1){
      if (verbose) new WarningMessage(new JFrame(), "Error", "You have to enter a value for the \"Double Bond Rule\"");
      doubleBondRuleField_.requestFocus();
      return false;
    }else{
      if (checkDoubleBondsRule(verbose)<0){
        doubleBondRuleField_.requestFocus();
        return false;
      }
    }
    if (chainCutoffField_.getText()!=null && chainCutoffField_.getText().trim().length()>0 && !chainCutoffField_.getText().trim().equalsIgnoreCase(CHAIN_CUTOFF_DEFAULT)){
      if (!isBetween0And100Percent(chainCutoffField_.getText(),(String)chainCutoffUnitCombo_.getSelectedItem(),CHAINCUTOFF_LABEL,verbose)){
        this.chainCutoffField_.requestFocus();
        return false;
      }
    }
    if (spectrumCoverageField_.getText()!=null && spectrumCoverageField_.getText().trim().length()>0){
      if (!isBetween0And100Percent(spectrumCoverageField_.getText(),(String)this.spectrumCoverageUnitCombo_.getSelectedItem(),"spectrum coverage",verbose)){
        this.spectrumCoverageField_.requestFocus();
        return false;
      }
    }
    if (!isBetween0And100Percent(basePeakCutoffField_.getText(),(String)this.basePeakUnitCombo_.getSelectedItem(),"base peak cutoff",verbose)){
      this.basePeakCutoffField_.requestFocus();
      return false;
    }
    if (rtMaxDevField_.getText()!=null && rtMaxDevField_.getText().trim().length()>0){
      try{
        double value = Double.parseDouble(rtMaxDevField_.getText());
        if (value<=0){
          if (verbose) new WarningMessage(new JFrame(), "Error", "The \""+RTMAXDEV_LABEL+"\" value must be a positive value greather than 0");
          rtMaxDevField_.requestFocus();
          return false;
        }
      } catch (NumberFormatException nfx){
        if (verbose) new WarningMessage(new JFrame(), "Error", "The \""+RTMAXDEV_LABEL+"\" value must be double format");
        rtMaxDevField_.requestFocus();
        return false;
      }
    }
    if (this.addChainPositions_.getText() != null && addChainPositions_.getText().trim().length()>0){
      boolean inputOK = this.checkIntValueInRange(addChainPositions_, ADD_CHAIN_POSITIONS_LABEL, 0, 10, verbose);
      if (!inputOK){
        addChainPositions_.requestFocus();
        return false;
      }
    }
    
    return true;
  }
  
  /**
   * reads from a text field, checks if the input is in integer format, and if the value is inside a definable range
   * @param input the text field to be read
   * @param label the label in front of the text field
   * @param min the range start value
   * @param max the range stop value
   * @param verbose should the errors be displayed on screen
   * @return true if the value is OK and inside the range, false otherwise
   */
  private boolean checkIntValueInRange (JTextField input, String label, int min, int max, boolean verbose){
    try{
      int value = Integer.parseInt(input.getText());
      if (value<min){
        if (verbose)new WarningMessage(new JFrame(), "Error", "The \""+label+"\" value cannot be smaller than "+min);
        return false;
      } else if (value>max){
        if (verbose) new WarningMessage(new JFrame(), "Error", "The \""+label+"\" value cannot be higher than "+max);
        return false;
      } else
        return true;
    } catch (NumberFormatException nfx){
      if (verbose) new WarningMessage(new JFrame(), "Error", "The \""+label+"\" value must be integer format");
      return false;
    }
  }
  
  /**
   * returns the String representation that is combined from the value entering field and selection combo with no magnitude, %, or permille
   * @param input the double value
   * @param magnitude the unit magnifier with no magnitude, %, or permille
   * @return the String representation for fraction, percent and permille input
   */
  private String getFullPercentValue(String input, String magnitude){
    String value = new String(input);
    if (!magnitude.equalsIgnoreCase(" "))
      value += magnitude;
    return value;
  }
  
  /**
   * checks if the entered value is a number, and if the value is between the allowed 0-100% (100% itself is not allowed)
   * @param input the double value as String
   * @param magnitude the unit magnifier with no magnitude, %, or permille
   * @param fieldName the parameter name for which the value is required
   * @param verbose should error messages be printed on screen
   * @return true if the entered values are valid double entries an within the 0-100% range (100% itself is not allowed)
   */
  private boolean isBetween0And100Percent(String input, String magnitude, String fieldName, boolean verbose){
    try {
      FragRuleParser.readPercentPermilleValue(getFullPercentValue(input,magnitude), fieldName, -1);
      return true;
    } catch (RulesException e) {
      if (verbose){
        if (e.getMessage().indexOf("must be double format")!=-1)
          new WarningMessage(new JFrame(), "Error", "The "+fieldName+" value is not a number");
        else if (e.getMessage().indexOf("must not be negative")!=-1)
          new WarningMessage(new JFrame(), "Error", "The "+fieldName+" must not be smaller than 0!");
        else if (e.getMessage().indexOf("must not be bigger than 100%")!=-1)
          new WarningMessage(new JFrame(), "Error", "The "+fieldName+" must not be bigger or equal than 100%!");

      }
      return false;
    }

  }
  
  /**
   * checks if the Java Regular Expression to extract the number of carbon atoms from the lipid species name is valid
   * @return true if the Java Regular Expression to extract the number of carbon atoms from the lipid species name is valid
   */
  private int checkCAtomsRule(boolean verbose){
    try{  
      int value = FragmentCalculator.getIntValueFromParsingRule(carbonAtomsRuleField_.getText(), listener_.getAnalyteName(), listener_.getLipidClassName(), "Carbon Atoms Rule");
      return value;
    } catch (RulesException rx){
      if (verbose) new WarningMessage(new JFrame(), "Error", rx.getMessage());
      return -1;
    }
  }

  /**
   * checks if the Java Regular Expression to extract the number of double bonds from the lipid species name is valid
   * @return true if the Java Regular Expression to extract the number of carbon atoms from the lipid species name is valid
   */
  private int checkDoubleBondsRule(boolean verbose){
    try{  
      int value = FragmentCalculator.getIntValueFromParsingRule(doubleBondRuleField_.getText(), listener_.getAnalyteName(), listener_.getLipidClassName(), "Double Bond Rule");
      return value;
    } catch (RulesException rx){
      if (verbose) new WarningMessage(new JFrame(), "Error", rx.getMessage());
      return -1;
    }
  }
  
  /**
   * checks if the general entries are OK and updates them in the cache
   * @return true if the general entries are OK
   */
  public boolean updateGeneralEntries(){
    return updateGeneralEntries(true);
  }
  
  /**
   * checks if the general entries are OK and updates them in the cache
   * @param verbose should error messages be reported in WarningMessage
   * @return true if the general entries are OK
   */
  public boolean updateGeneralEntries(boolean verbose){
    if (!listener_.areGeneralSettingsDependingFieldsVisible()) return true;
    if (!checkGeneralEntries(verbose)) return false;
    listener_.updateGeneralSettings();
    return true;
  }
  
  public void actionPerformed(ActionEvent e)
  {
    String command = e.getActionCommand();
    if (command.equalsIgnoreCase("saveGeneralRules")){
      if (!checkGeneralEntries(true)) return;
      if (listener_.performStorageOfInitialGeneralSettings()){
        saveGeneralRules_.setVisible(false);
        listener_.refreshGeneralSettingsDependantDisplay();
      }
    }
  }

  /**
   * Checks if a string is numeric
   * @param str
   * @return
   */
  private static boolean isNumeric(String input)  
  {  
    try  
    {  
      Double.parseDouble(input); 
      return true; 
    }  
    catch(NumberFormatException nfe)  
    {  
      return false;  
    }       
  }
  
  /**
   * 
   * @return a VO representation of the entered values
   */
  public GeneralSettingsVO getValues(){
    //TODO: it would be programmatically more correct to perform a checkGeneralEntries(false),
    //      however, I guess this won't work with the current error handling of the RDI
    Integer alkylChains = readIntegerFromTextField(alkylChainsAmountField_);
    Integer alkenylChains = readIntegerFromTextField(alkenylChainsAmountField_);
    boolean isSingleChainAllowed = this.readBooleanComboBox(singleChainIdentification_);
    String chainCutoff = null;
    if (chainCutoffField_.getText()!=null && chainCutoffField_.getText().trim().length()>0 && !chainCutoffField_.getText().trim().equalsIgnoreCase(CHAIN_CUTOFF_DEFAULT)){
      double value = Double.parseDouble(chainCutoffField_.getText());
      if (value>0){
        chainCutoff = chainCutoffField_.getText();
        if (!((String)this.chainCutoffUnitCombo_.getSelectedItem()).endsWith(" "))
          chainCutoff += (String)this.chainCutoffUnitCombo_.getSelectedItem();
      }
    }
    String basePeakCutoff = basePeakCutoffField_.getText();
    if (!((String)this.basePeakUnitCombo_.getSelectedItem()).endsWith(" "))
      basePeakCutoff += (String)this.basePeakUnitCombo_.getSelectedItem();
    String spectrumCoverage =  null;
    if (spectrumCoverageField_.getText()!=null && spectrumCoverageField_.getText().trim().length()>0){
      double value = Double.parseDouble(spectrumCoverageField_.getText());
      if (value>=0){
        spectrumCoverage = spectrumCoverageField_.getText();
        if (!((String)this.spectrumCoverageUnitCombo_.getSelectedItem()).endsWith(" "))
          spectrumCoverage += (String)this.spectrumCoverageUnitCombo_.getSelectedItem();
      }
    }
    boolean rtPostProcessing = this.readBooleanComboBox(possibibilitiesPostProcessing_);
    boolean rtParallelSeries = this.readBooleanComboBox(possibibilitiesParallelMode_);
    Double rtMaxDev = this.readDoubleFromTextField(this.rtMaxDevField_);
    int msIdOrder = RulesContainer.ORDER_MS1_FIRST;
    if (msIdentificationOrder_.getItemAt((msIdentificationOrder_.getSelectedIndex())).equalsIgnoreCase(MSN_ONLY))
      msIdOrder = RulesContainer.ORDER_MSN_ONLY;
    if (msIdentificationOrder_.getItemAt((msIdentificationOrder_.getSelectedIndex())).equalsIgnoreCase(MSN_FIRST))
      msIdOrder = RulesContainer.ORDER_MSN_FIRST;
    Integer addPositions = readIntegerFromTextField(this.addChainPositions_);
    
    //TODO: the number of LCB chains, the LCB chain library, the FaHydroxylationRange and the LcbHydroxylationRange the is not read from the GUI - appropriate input masks have to be implemented
    GeneralSettingsVO settings = new GeneralSettingsVO(new Integer(chainsAmountField_.getText()),
        alkylChains, alkenylChains, null, addPositions, chainsLibraryField_.getText(), null, -1, -1, -1,
        -1, carbonAtomsRuleField_.getText(), this.doubleBondRuleField_.getText(), isSingleChainAllowed,
        chainCutoff, basePeakCutoff, spectrumCoverage, rtPostProcessing, rtParallelSeries, rtMaxDev,
        msIdOrder);
    return settings;
  }
  
  /**
   * reads an value from a text field and returns an integer - if nothing is found "null"
   * @param field the input field
   * @return the integer value if the input text is not empty
   */
  private Integer readIntegerFromTextField(JTextField field){
    Integer value = null;
    if (field.getText()!=null && field.getText().trim().length()>0)
      value = new Integer(field.getText());
    return value;
  }
  
  /**
   * reads an value from a text field and returns a double - if nothing is found "null"
   * @param field the input field
   * @return the double value if the input text is not empty
   */
  private Double readDoubleFromTextField(JTextField field){
    Double value = null;
    if (field.getText()!=null && field.getText().trim().length()>0)
      value = new Double(field.getText());
    return value;
  }
  
  /**
   * returns true/false in Boolean format from a JComboBox
   * @param box the JComboBox
   * @return the boolean value
   */
  private boolean readBooleanComboBox(JComboBox<String> box){
    boolean value = false;
    if (box.getItemAt(box.getSelectedIndex()).equalsIgnoreCase("true")) value = true;    
    return value;
  }
    
  /**
   * Checks if spectrumCoverageField and basePeakCutoffField is a number
   * @return boolean a number
   */
  public boolean checkPercentageFieldsForNumbers(){
    if(isNumeric(basePeakCutoffField_.getText())){
      return true;
    } else {
      JOptionPane.showMessageDialog(new JFrame(), "The Base Peak Cutoff Field is not a number!");
      return false;
    }  
  }
  
}
