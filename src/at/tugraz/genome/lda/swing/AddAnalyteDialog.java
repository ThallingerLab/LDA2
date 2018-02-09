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

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Hashtable;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.analysis.AnalyteAddRemoveListener;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.verifier.DoubleVerifier;
import at.tugraz.genome.lda.verifier.IntegerVerifier;
import at.tugraz.genome.lda.verifier.MzInputVerifierForTolerance;
import at.tugraz.genome.lda.vos.AddAnalyteVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class AddAnalyteDialog extends JDialog implements ActionListener
{

  private static final long serialVersionUID = -8260667248906145073L;
  
  private JTextField analyteName_;
  private JTextField analyteFormula_;
  private JTextField modName_;
  private JTextField modFormula_;
  private JTextField mzTolerance_;
  private JTextField exactMass_;
  private JTextField charge_;
  private JTextField rTime_;
  private int positionToAddNewAnalyte_;
  private AnalyteAddRemoveListener parentListener_;

  public AddAnalyteDialog(JFrame parent, String title, String message, String analyteName, float mz, String analyteFormula,
      String modName, String modFormula, int positionToAddNewAnalyte, boolean showRt, AnalyteAddRemoveListener parentListener) {
    super(parent, title, true);
    parentListener_ = parentListener;
    positionToAddNewAnalyte_ = positionToAddNewAnalyte;
    setLocation(380,240);
    this.setLayout(new BorderLayout());
    JPanel messagePane = new JPanel();
    messagePane.add(new JLabel(message));
    getContentPane().add(messagePane,BorderLayout.NORTH);
    JPanel inputPane = new JPanel();
    inputPane.setLayout(new GridBagLayout());
    JLabel label = new JLabel("name: ");
    label.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_NAME);
    inputPane.add(label,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    String rt = "";
    if (showRt){
      rt = analyteName.substring(analyteName.lastIndexOf("_")+1);
      analyteName = analyteName.substring(0,analyteName.lastIndexOf("_"));
    }
    analyteName_ = new JTextField(5);
    analyteName_.setText(analyteName);
    analyteName_.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_NAME);
    inputPane.add(analyteName_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    label = new JLabel("formula: ");
    label.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_FORMULA);
    inputPane.add(label,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    analyteFormula_ = new JTextField(10);
    analyteFormula_.setText(analyteFormula);
    analyteFormula_.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_FORMULA);
    inputPane.add(analyteFormula_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));  
    label = new JLabel("modification: ");
    label.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_MOD_NAME);
    inputPane.add(label,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    modName_ = new JTextField(5);
    modName_.setText(modName);
    modName_.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_MOD_NAME);
    inputPane.add(modName_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    label = new JLabel("formula: ");
    label.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_MOD_FORMULA);
    inputPane.add(label,new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    modFormula_ = new JTextField(10);
    modFormula_.setText(modFormula);
    modFormula_.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_MOD_FORMULA);
    inputPane.add(modFormula_,new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    label = new JLabel("m/z tolerance: ");
    label.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_MZ);
    inputPane.add(label,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    mzTolerance_ = new JTextField(5);
    mzTolerance_.setHorizontalAlignment(JTextField.RIGHT);
    mzTolerance_.setInputVerifier(new DoubleVerifier());
    mzTolerance_.setText(String.valueOf(Calculator.roundFloat(LipidomicsConstants.getCoarseChromMzTolerance(mz),5)));
    mzTolerance_.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_MZ);
    inputPane.add(mzTolerance_,new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    label = new JLabel("exact mass: ");
    label.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_MASS);
    inputPane.add(label,new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    exactMass_ = new JTextField(10);
    exactMass_.setHorizontalAlignment(JTextField.RIGHT);
    exactMass_.setInputVerifier(new MzInputVerifierForTolerance(mzTolerance_));
    exactMass_.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_MASS);
    inputPane.add(exactMass_,new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    label = new JLabel("charge: ");
    label.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_CHARGE);
    inputPane.add(label,new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    charge_ = new JTextField(5);
    charge_.setHorizontalAlignment(JTextField.RIGHT);
    charge_.setInputVerifier(new IntegerVerifier());
    charge_.setText("1");
    charge_.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_CHARGE);
    inputPane.add(charge_,new GridBagConstraints(1, 3, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    if (showRt){
      label = new JLabel("RT: ");
      label.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_RT);
      inputPane.add(label,new GridBagConstraints(2, 3, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      rTime_ = new JTextField(5);
      rTime_.setHorizontalAlignment(JTextField.RIGHT);
      rTime_.setInputVerifier(new DoubleVerifier(true,true));
      rTime_.setText(rt);
      rTime_.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_RT);
      inputPane.add(rTime_,new GridBagConstraints(3, 3, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    }
    getContentPane().add(inputPane, BorderLayout.CENTER);    
    JPanel buttonPane = new JPanel();
    JButton cancelButton = new JButton("Cancel");
    buttonPane.add(cancelButton); 
    cancelButton.addActionListener(this);
    cancelButton.setActionCommand("Cancel");
    cancelButton.setToolTipText(TooltipTexts.CANCEL_GENERAL);
    JButton okButton = new JButton("OK"); 
    buttonPane.add(okButton); 
    okButton.addActionListener(this);
    okButton.setActionCommand("OK");
    okButton.setToolTipText(TooltipTexts.DISPLAY_ADD_ANALYTE_ACCEPT);
    getContentPane().add(buttonPane, BorderLayout.SOUTH);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    pack(); 
    setVisible(true);
  }  
  
  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase("OK")){
      String notFilledOutFields = "";
      if (analyteName_.getText()==null || analyteName_.getText().length()==0)
        notFilledOutFields+="name, ";
      if (analyteFormula_.getText()==null || analyteFormula_.getText().length()==0)
        notFilledOutFields+="formula, ";
      if (mzTolerance_.getText()==null || mzTolerance_.getText().length()==0)
        notFilledOutFields+="m/z tolerance, ";
      if (exactMass_.getText()==null || exactMass_.getText().length()==0)
        notFilledOutFields+="exact mass, ";
      if (charge_.getText()==null || charge_.getText().length()==0)
        notFilledOutFields+="charge, ";
      if (notFilledOutFields.length()==0){
        try{
          if (analyteFormula_.getText()!=null && analyteFormula_.getText().length()>0 && checkChemicalFormula(analyteFormula_.getText()) && 
            (modFormula_.getText()==null || modFormula_.getText().length()==0 || checkModFormula(analyteFormula_.getText(),modFormula_.getText()))){
            String modification = "";
            if (modName_.getText()!=null && modName_.getText().length()>0)
              modification = modName_.getText();
            String modFormula = "";
            if (modFormula_.getText()!=null&&modFormula_.getText().length()>0&&checkModFormula(analyteFormula_.getText(),modFormula_.getText()))
              modFormula = modFormula_.getText();
            String rt = "";
            if (this.rTime_!=null && rTime_.getText()!=null && rTime_.getText().length()>0) rt = rTime_.getText().trim();
            AddAnalyteVO addVO = new AddAnalyteVO(analyteName_.getText().trim(),analyteFormula_.getText().replaceAll(" ", "").trim(),
                modification, modFormula.replaceAll(" ", ""),mzTolerance_.getText().trim(),exactMass_.getText().trim(),charge_.getText().trim(),rt);
            setVisible(false);
            dispose();
            parentListener_.addAnalyte(positionToAddNewAnalyte_, addVO);
          }
        }catch(ChemicalFormulaException cfx){
          //Comment: The warning message is already checked in the method
        }
      }else{
        notFilledOutFields = notFilledOutFields.substring(0,notFilledOutFields.length()-2);
        new WarningMessage(new JFrame(), "Error", "You have to fill out the following fields: "+notFilledOutFields+"!");
      }
    } else if (e.getActionCommand().equalsIgnoreCase("Cancel")){
      setVisible(false);
      dispose();
    }
  }
  
  private boolean checkModFormula(String chemicalFormula, String modFormula) throws ChemicalFormulaException{
    Hashtable<String,Integer> formulaElemCard = StaticUtils.categorizeFormula(chemicalFormula);
    Hashtable<String,Integer> modElemCard = StaticUtils.categorizeFormula(modFormula);
    String formulaWithMod = "";
    for (String element : formulaElemCard.keySet()){
      if (formulaWithMod.length()>0) formulaWithMod+=" ";
      int amount = formulaElemCard.get(element);
      if (modElemCard.containsKey(element)) amount+=modElemCard.get(element);
      if (amount<0){
        new WarningMessage(new JFrame(), "Error", "The modification formula "+modFormula+" results in a negative amount of elements for "+element+"!");
        return false;
      }
      formulaWithMod+=element+String.valueOf(amount);
    }
    return checkChemicalFormula(formulaWithMod);
  }

  private boolean checkChemicalFormula(String formula){
    if (formula.contains("-")){
      new WarningMessage(new JFrame(), "Error", "The formula "+formula+" must not contain any negative values!");
      return false;
    }
    char[] formulaChars = formula.toCharArray();
    String formulaToCheck = "";
    boolean isPreviousDigit = false;
    for (int i=0;i!=formulaChars.length;i++){
      char currentChar = formulaChars[i];
      if (isPreviousDigit && !Character.isDigit(currentChar)){
        formulaToCheck+=" ";
      }
      formulaToCheck+=String.valueOf(currentChar);
      isPreviousDigit = Character.isDigit(currentChar);
    }
    ElementConfigParser aaParser = Settings.getElementParser();
    try {
      aaParser.calculateTheoreticalMass(formulaToCheck, false);
      return true;
    }
    catch (SpectrummillParserException e) {
      new WarningMessage(new JFrame(), "Error", "The formula "+formula+" is not OK! "+e.getMessage());
      return false;
    } 
  }

 
}
