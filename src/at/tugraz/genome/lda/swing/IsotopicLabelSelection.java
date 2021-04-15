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

import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.event.ActionListener;
import java.util.Hashtable;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.IsoLabelInputException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.verifier.DoubleVerifier;
import at.tugraz.genome.lda.verifier.IntegerVerifier;
import at.tugraz.genome.lda.vos.IsotopicLabelVO;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * graphical interface for setting information about the meaning of the isotopic label
 * @author Juergen Hartler
 *
 */
public class IsotopicLabelSelection extends IsotopicLabelVO
{
  
  /** a check box whether to select this label*/
  private JCheckBox selected_ = null;
  /** the encoding for the isotopic label*/
  private JTextField label_ = null;
  /** the omega position*/
  private JTextField omega_ = null;
  /** the difference in chemical elements introduced by the label*/
  private JTextField formula_ = null;
  /** the retention time shift caused by the label*/
  private JTextField rtShift_ = null;
  /** panel containing the omega information*/
  private JPanel omegaPanel_ = null;
  /** the button to delete a label*/
  private JButton delete_ = null;

  /** the delete icon (a red X)*/
  private final ImageIcon deleteIcon_ = new ImageIcon(getClass().getResource(
  "/images/Delete.gif"));
  /** the action prefix for deleting an entry; the prefix is followed by the row number*/
  public final static String DELETE_ACTION_PREFIX = "deleteLabel_";
 

  /**
   * constructor for graphical interface for setting information about the meaning of the isotopic label
   * @param vo value object containing information about the isotopic label
   * @param parentListener the parent listener for action commands
   */
  private IsotopicLabelSelection(IsotopicLabelVO vo, ActionListener parentListener) {
    super(vo);
    selected_  = new JCheckBox();
    selected_.setSelected(true);
    selected_.setToolTipText(TooltipTexts.GENERAL_ACCEPT_SINGLE_START+labelId_+TooltipTexts.GENERAL_ACCEPT_SINGLE_END);

    label_ = new JTextField(3);
    label_.setText(labelId_);
    label_.setToolTipText(TooltipTexts.ISO_LABEL_LABELID);
    omega_ = new JTextField(3);
    if (vo.getOmegaPosition()>0)
      omega_.setText(String.valueOf(vo.getOmegaPosition()));
    else
      omega_.setText("");
    omega_.setInputVerifier(new IntegerVerifier());
    omega_.setToolTipText(TooltipTexts.ISO_LABEL_OMEGA);
    omegaPanel_ = new JPanel();
    omegaPanel_.add(new JLabel("n-"));
    omegaPanel_.add(omega_);
    formula_ = new JTextField(8);
    formula_.setText(StaticUtils.getFormulaInHillNotation_PlusFirst(this.labelElements_, true));
    formula_.setToolTipText(TooltipTexts.ISO_LABEL_FORMULA);
    rtShift_ = new JTextField(6);
    double rtShift = 0;
    if (vo.getRtShift()!=null)
      rtShift = vo.getRtShift().doubleValue();
    rtShift_.setText(Calculator.FormatNumberToString(rtShift,3d));
    rtShift_.setToolTipText(TooltipTexts.ISO_LABEL_RTSHIFT);
    rtShift_.setInputVerifier(new DoubleVerifier(false,false));
    delete_ = new JButton("",deleteIcon_);
    delete_.setIconTextGap(2);
    delete_.setMargin(new Insets(2,2,2,2));
    delete_.addActionListener(parentListener);
    delete_.setToolTipText(TooltipTexts.ISO_DELETE);
    
    //TODO: those if/else is only for testing
/****    if (labelId_.equalsIgnoreCase("A"))
      omega_.setText("9");
    else if (labelId_.equalsIgnoreCase("B"))
      omega_.setText("6");
    else if (labelId_.equalsIgnoreCase("C"))
      omega_.setText("3");
    else if (labelId_.equalsIgnoreCase("D"))
      omega_.setText("7");
    else if (labelId_.equalsIgnoreCase("E"))
      omega_.setText("10");
    else if (labelId_.equalsIgnoreCase("F"))
      omega_.setText("7");*/
  }

  /**
   * @return true if this label is selected
   */
  public boolean isSelected()
  {
    return selected_.isSelected();
  }

  /**
   * selects/deselects this specific label
   * @param select true: select; false: deselect
   */
  public void setSelected(boolean select)
  {
    this.selected_.setSelected(select);
  }
  
  /**
   * inverts the selection
   */
  public void invertSelection() {
    this.selected_.setSelected(!selected_.isSelected());
  }
  
  
  /**
   * constructor for graphical interface for setting information about the meaning of the isotopic label
   * @param vo value object containing information about the isotopic label
   * @param panel the graphical panel where this graphical user interface should be added
   * @param row the row in the panel
   * @param parentListener the parent listener for action commands
   * @return the IsotopicLabelSelection that has been added to the panel
   */
  public static IsotopicLabelSelection addIsotopicLabelToPanel(IsotopicLabelVO vo, JPanel panel, int row, ActionListener parentListener) {
    IsotopicLabelSelection selection = null;
    if (vo instanceof IsotopicLabelSelection)
      selection = (IsotopicLabelSelection)vo;
    else
      selection = new IsotopicLabelSelection(vo,parentListener);
    panel.add(selection.selected_, new GridBagConstraints(0, row, 1, 1, 0.0, 0.0,
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    panel.add(selection.label_, new GridBagConstraints(1, row, 1, 1, 0.0, 0.0,
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0)); 
    panel.add(selection.omegaPanel_, new GridBagConstraints(2, row, 1, 1, 0.0, 0.0,
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    panel.add(selection.formula_, new GridBagConstraints(3, row, 1, 1, 0.0, 0.0,
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    panel.add(selection.rtShift_, new GridBagConstraints(4, row, 1, 1, 0.0, 0.0,
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    selection.delete_.setActionCommand(DELETE_ACTION_PREFIX+String.valueOf(row));
    panel.add(selection.delete_, new GridBagConstraints(5, row, 1, 1, 0.0, 0.0,
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    return selection;
  }
  
  
  /**
   * 
   * @return the currently entered label id
   */
  public String getEnteredLabelId() {
    return label_.getText();
  }
  
  
  /**
   * verifies the input values and returns the corresponding VO
   * @return the value object for the label information
   * @throws IsoLabelInputException when some input is wrong or insufficient
   */
  public IsotopicLabelVO getLabelInformation() throws IsoLabelInputException {
    if (label_.getText()==null || label_.getText().length()==0)
      throw new IsoLabelInputException("the label field must not be empty!");
    if (omega_.getText()==null || omega_.getText().length()==0)
      throw new IsoLabelInputException("the \u03C9 input field must not be empty!");
    if (formula_.getText()==null || formula_.getText().length()==0)
      throw new IsoLabelInputException("the formula field must not be empty!");
    Hashtable<String,Integer> labelElements = null;
    try {labelElements = StaticUtils.categorizeFormula(formula_.getText());
    } catch (ChemicalFormulaException e) {
      throw new IsoLabelInputException(e);
    }
    if (rtShift_.getText()==null || rtShift_.getText().length()==0)
      throw new IsoLabelInputException("the RT-shift field must not be empty!");
    IsotopicLabelVO labelVO = new IsotopicLabelVO (label_.getText(), Integer.parseInt(omega_.getText()), labelElements,Float.parseFloat(rtShift_.getText()),getPrefixes());
    
    return labelVO;
  }

}
