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

import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import at.tugraz.genome.lda.verifier.DoubleVerifier;

/**
 * A panel allowing for the selection of display options in the display results tab.
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
public class DisplayTolerancePanel extends JPanel
{
	private static final long serialVersionUID = 1L;
	
	private LipidDataAnalyzer lda_;
	
	private final static Font SMALL_FONT = new Font("Arial",Font.PLAIN,9);
	/** action for setting a locked m/z range for the 3D viewer*/
  public final static String CHANGE_LOCK_MZ_RANGE = "lockMzRange";
  /** the displayed label next to check box for activating the lock mass range*/
  public final static String DISPLAY_LOCK_MZ_TEXT = "Lock m/z range ";
  
  public final static String UPDATE_QUANT_TOL_OF_CURRENTLY_SELECTED = "updateQuantTolOfCurrentlySelected";
  
  public final static String SHOW_MSN_ONLY = "showMSnOnly";
  
  public final static String SHOW_CHAIN_ONLY = "showChainOnly";
  
  public final static String SHOW_MSN_NAMES = "showMSnNames";
  
  public final static String SHOW_2D_CHANGED = "show2dChanged";
  
  public final static String SHOW_OMEGA_NAMES = "showOmegaNames";
  
	
	private JTextField displayMinusTolerance_;
	private JTextField displayPlusTolerance_;
	private JTextField displayRtStart_;
	private JTextField displayRtStop_;
	/** this check box decides whether a static (locked) range is used, or it is adapted to the selected analyte*/
	private JCheckBox lockMzRange_;
	/** text field to set a permanent m/z start value for the 3D viewer*/
	private JTextField displayMzStart_;
	/** text field to set a permanent m/z stop value for the 3D viewer*/
	private JTextField displayMzStop_;
	/** show only names with MSn evidence in the display table or not; oxLipids extension */
	private JCheckBox showMSnEvidence_;
	/** show only names with MSn chain evidence in the display table or not; oxLipids extension */
	private JCheckBox showChainEvidence_;	
	/** show the names in MSn style in the display results or not*/
	private JCheckBox showMSnNames_;
	/** show the 2D chromatogram or not*/
	private JCheckBox show2D_;
	/** show the names in MSn style with available omega double bond assignments in the display results or not*/
	private JCheckBox showOmegaNames_;
	
	
	public DisplayTolerancePanel(LipidDataAnalyzer lda, boolean showOmegaEditingTools)
	{
		this.lda_ = lda;
		
    this.setLayout(new GridBagLayout());
    JLabel diplayTolMinus = new JLabel("- m/z: ");
    diplayTolMinus.setFont(SMALL_FONT);
    diplayTolMinus.setToolTipText(TooltipTexts.DISPLAY_MZ_MINUS);
    this.add(diplayTolMinus,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 2, 0, 0), 0, 0));
    JLabel diplayTolPlus = new JLabel("+ m/z: ");
    diplayTolPlus.setFont(SMALL_FONT);
    diplayTolPlus.setToolTipText(TooltipTexts.DISPLAY_MZ_PLUS);
    this.add(diplayTolPlus,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 2, 0, 0), 0, 0));
    
    displayMinusTolerance_ = new JTextField(2);
    displayMinusTolerance_.setFont(SMALL_FONT);
    displayMinusTolerance_.setText("1.5");
    displayMinusTolerance_.setHorizontalAlignment(JTextField.RIGHT);
    displayMinusTolerance_.setToolTipText(TooltipTexts.DISPLAY_MZ_MINUS);
    displayMinusTolerance_.setInputVerifier(new DoubleVerifier());
    this.add(displayMinusTolerance_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 1, 0, 0), 0, 0));   
    
    displayPlusTolerance_ = new JTextField(2);
    displayPlusTolerance_.setFont(SMALL_FONT);
    displayPlusTolerance_.setText("2.5");
    displayPlusTolerance_.setHorizontalAlignment(JTextField.RIGHT);
    displayPlusTolerance_.setToolTipText(TooltipTexts.DISPLAY_MZ_PLUS);
    displayPlusTolerance_.setInputVerifier(new DoubleVerifier());
    this.add(displayPlusTolerance_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 1, 0, 0), 0, 0));
    JLabel diplayTolUnit1 = new JLabel("[Da]");
    diplayTolUnit1.setFont(SMALL_FONT);
    this.add(diplayTolUnit1,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 1, 0, 0), 0, 0));
    JLabel rtTolMinusText = new JLabel("Start: ");
    rtTolMinusText.setFont(SMALL_FONT);
    rtTolMinusText.setToolTipText(TooltipTexts.DISPLAY_RT_START);
    this.add(rtTolMinusText,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 5, 0, 0), 0, 0));
    JLabel rtTolPlusText = new JLabel("Stop: ");
    rtTolPlusText.setFont(SMALL_FONT);
    rtTolPlusText.setToolTipText(TooltipTexts.DISPLAY_RT_STOP);
    this.add(rtTolPlusText,new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 5, 0, 0), 0, 0));
    
    displayRtStart_ = new JTextField(2);
    displayRtStart_.setFont(SMALL_FONT);
    displayRtStart_.setText("");
    displayRtStart_.setHorizontalAlignment(JTextField.RIGHT);
    displayRtStart_.setToolTipText(TooltipTexts.DISPLAY_RT_START);
    displayRtStart_.setInputVerifier(new DoubleVerifier());
    this.add(displayRtStart_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    
    displayRtStop_ = new JTextField(2);
    displayRtStop_.setFont(SMALL_FONT);
    displayRtStop_.setText("");
    displayRtStop_.setHorizontalAlignment(JTextField.RIGHT);
    displayRtStop_.setToolTipText(TooltipTexts.DISPLAY_RT_STOP);
    displayRtStop_.setInputVerifier(new DoubleVerifier());
    this.add(displayRtStop_,new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    JLabel diplayTolUnit2 = new JLabel("[min]");
    diplayTolUnit2.setFont(SMALL_FONT);
//    diplayTolUnit1.setToolTipText(TooltipTexts.DISPLAY_MZ_MINUS);
    this.add(diplayTolUnit2,new GridBagConstraints(4, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    
    JPanel lockMxSelectionPanel = new JPanel();
    
    lockMzRange_ = new JCheckBox();
    lockMzRange_.setFont(SMALL_FONT);
    lockMzRange_.addItemListener(new PanelItemListener(CHANGE_LOCK_MZ_RANGE));
    lockMzRange_.setToolTipText(TooltipTexts.DISPLAY_LOCK_MZ);
    lockMxSelectionPanel.add(lockMzRange_);
    JLabel lockMz = new JLabel(DISPLAY_LOCK_MZ_TEXT);
    lockMz.setToolTipText(TooltipTexts.DISPLAY_LOCK_MZ);
    lockMz.setFont(SMALL_FONT);
    lockMxSelectionPanel.add(lockMz);
    
    displayMzStart_ = new JTextField(4);
    displayMzStart_.setFont(SMALL_FONT);
    displayMzStart_.setText("");
    displayMzStart_.setHorizontalAlignment(JTextField.RIGHT);
    displayMzStart_.setToolTipText(TooltipTexts.DISPLAY_MZ_MINUS);
    displayMzStart_.setInputVerifier(new DoubleVerifier());
    displayMzStart_.setEnabled(false);
    lockMxSelectionPanel.add(displayMzStart_);
    JLabel mzRangeSign = new JLabel("-");
    lockMxSelectionPanel.add(mzRangeSign);
    
    displayMzStop_ = new JTextField(4);
    displayMzStop_.setFont(SMALL_FONT);
    displayMzStop_.setText("");
    displayMzStop_.setHorizontalAlignment(JTextField.RIGHT);
    displayMzStop_.setToolTipText(TooltipTexts.DISPLAY_MZ_MINUS);
    displayMzStop_.setInputVerifier(new DoubleVerifier());
    displayMzStop_.setEnabled(false);
    lockMxSelectionPanel.add(displayMzStop_);
    
    this.add(lockMxSelectionPanel ,new GridBagConstraints(0, 2, 6, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));

    JButton quantTolUpdate = new JButton("Update");
    quantTolUpdate.addActionListener(lda_);
    quantTolUpdate.setFont(quantTolUpdate.getFont().deriveFont(10f));
    quantTolUpdate.setMargin(new Insets(1,5,1,5));
    quantTolUpdate.setActionCommand(UPDATE_QUANT_TOL_OF_CURRENTLY_SELECTED);
    quantTolUpdate.setToolTipText(TooltipTexts.DISPLAY_UPDATE);
    this.add(quantTolUpdate,new GridBagConstraints(5, 0, 1, 2, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 6, 0, 0), 0, 0));
    JPanel showOptionsPanel = new JPanel();
    this.add(showOptionsPanel,new GridBagConstraints(0, 3, 6, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    //start: added with the oxidized lipids extension
    JPanel showOptionsPanel2 = new JPanel();
    this.add(showOptionsPanel2,new GridBagConstraints(0, 4, 6, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    JPanel showOptionsPanel3 = new JPanel();
    this.add(showOptionsPanel3,new GridBagConstraints(0, 5, 6, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    showMSnEvidence_ = new JCheckBox("");
    showMSnEvidence_.addItemListener(new PanelItemListener(SHOW_MSN_ONLY));
    showOptionsPanel2.add(showMSnEvidence_);
    JLabel showMSnEvidence = new JLabel("Show MSn Only");
    showMSnEvidence.setToolTipText("Show only lipids with MSn headgroup evidence");
    showOptionsPanel2.add(showMSnEvidence);
  	
    showChainEvidence_ = new JCheckBox("");
    showChainEvidence_.setEnabled(false);
    showChainEvidence_.addItemListener(new PanelItemListener(SHOW_CHAIN_ONLY));
    showOptionsPanel3.add(showChainEvidence_);
    JLabel showChainEvidence = new JLabel("Chain Evidence Only");
    showChainEvidence.setToolTipText("Show only lipids with MSn headgroup- and chain evidence");
    showOptionsPanel3.add(showChainEvidence);
    //end: added with the oxidized lipids extension
    
    showMSnNames_ = new JCheckBox();
    showMSnNames_.addItemListener(new PanelItemListener(SHOW_MSN_NAMES));
    showMSnNames_.setToolTipText(TooltipTexts.DISPLAY_SHOW_MSN);
    showOptionsPanel.add(showMSnNames_);
    JLabel showMSn = new JLabel("Show MSn");
    showMSn.setToolTipText(TooltipTexts.DISPLAY_SHOW_MSN);
    showOptionsPanel.add(showMSn);
    
    show2D_ = new JCheckBox();
    show2D_.setSelected(true);
    show2D_.addItemListener(new PanelItemListener(SHOW_2D_CHANGED));
    show2D_.setToolTipText(TooltipTexts.DISPLAY_SHOW_2D);
    showOptionsPanel.add(show2D_);
    JLabel show2d = new JLabel("Show 2D-View");
    show2d.setToolTipText(TooltipTexts.DISPLAY_SHOW_2D);
    showOptionsPanel.add(show2d);
    
    showOmegaNames_ = new JCheckBox();
    if (showOmegaEditingTools) {
      showOmegaNames_.addItemListener(new PanelItemListener(SHOW_OMEGA_NAMES));
      showOmegaNames_.setToolTipText(TooltipTexts.DISPLAY_SHOW_OMEGA);
      JPanel showOmegaOptionPanel = new JPanel();
      showOmegaOptionPanel.add(showOmegaNames_);
      JLabel showOmega = new JLabel("Show assigned \u03C9-DB");
      showOmegaOptionPanel.add(showOmega);
      this.add(showOmegaOptionPanel,new GridBagConstraints(0, 6, 6, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    }
	}
	
	
	public JTextField getDisplayMinusTolerance()
	{
		return displayMinusTolerance_;
	}


	public JTextField getDisplayPlusTolerance()
	{
		return displayPlusTolerance_;
	}


	public JTextField getDisplayRtStart()
	{
		return displayRtStart_;
	}


	public JTextField getDisplayRtStop()
	{
		return displayRtStop_;
	}


	public JCheckBox getLockMzRange()
	{
		return lockMzRange_;
	}


	public JTextField getDisplayMzStart()
	{
		return displayMzStart_;
	}


	public JTextField getDisplayMzStop()
	{
		return displayMzStop_;
	}


	public JCheckBox getShowMSnEvidence()
	{
		return showMSnEvidence_;
	}


	public JCheckBox getShowChainEvidence()
	{
		return showChainEvidence_;
	}


	public JCheckBox getShowMSnNames()
	{
		return showMSnNames_;
	}


	public JCheckBox getShow2D()
	{
		return show2D_;
	}


	public JCheckBox getShowOmegaNames()
	{
		return showOmegaNames_;
	}


	private class PanelItemListener implements ItemListener
  {
    private String m_ctrl;
    
    public PanelItemListener(String ctrl){
      m_ctrl = ctrl;
    }

    public void itemStateChanged(ItemEvent e)  {
      lda_.change2DTypeState(e,m_ctrl);
    }
  }
	
}
