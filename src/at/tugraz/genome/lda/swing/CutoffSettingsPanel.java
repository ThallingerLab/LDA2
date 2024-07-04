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

package at.tugraz.genome.lda.swing;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;

import at.tugraz.genome.lda.analysis.ClassNamesExtractor;
import at.tugraz.genome.lda.exception.AbsoluteSettingsInputException;
import at.tugraz.genome.lda.exception.SettingsException;
import at.tugraz.genome.lda.verifier.DoubleVerifier;
import at.tugraz.genome.lda.verifier.IntegerMaxVerifier;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class CutoffSettingsPanel extends JPanel
{

  private static final long serialVersionUID = 7828948004968836687L;
  private JPanel inputPanel_;
  private JTextField maxIsotope_;
  private ClassNamesExtractor extractor_;
  private String tooltip_;
  
  private Hashtable<String,JTextField> inputFields_;
  
  private final static int itemsPerLine_ = 6;

  public CutoffSettingsPanel(String tooltip){
    super();
    this.tooltip_ = tooltip;
    maxIsotope_ = new JTextField(1);
    maxIsotope_.setText("2");
    maxIsotope_.setHorizontalAlignment(JTextField.RIGHT);
    maxIsotope_.setToolTipText(tooltip);
    maxIsotope_.setInputVerifier(new IntegerMaxVerifier(true,0,9));

    this.initComponents();
  }
  
  private void initComponents(){   
    inputPanel_ = new JPanel();
//    inputPanel_.setLayout(new GridBagLayout());
    this.add(inputPanel_);  
  }
  
  public void showSettingsPanel(ClassNamesExtractor extractor){
    this.extractor_ = extractor;
    this.initInputFields();
    this.invalidate();
    this.updateUI();
  }

  private void initInputFields(){
    inputPanel_.removeAll();
    JTabbedPane dummyTab = new JTabbedPane();
    inputPanel_.add(dummyTab);
    
    
    inputFields_ = new Hashtable<String,JTextField>();
    JPanel inputFields = new JPanel();
    inputFields.setLayout(new GridBagLayout());
    dummyTab.addTab("Cutoff-thresholds",inputFields);
    
    Vector<String> lipidClasses = extractor_.getLipidClasses();
//    if (lipidClasses.size()>(itemsPerLine_*2)){
    int lines = lipidClasses.size()/itemsPerLine_;
    if (lipidClasses.size()%itemsPerLine_!=0) lines++;
    inputFields.setPreferredSize(new Dimension(850,30*(lines+1)));
//    }
    JPanel isotopePanel = new JPanel();
    inputFields.add(isotopePanel,new GridBagConstraints(0, 0, itemsPerLine_*3, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
    JLabel label = new JLabel("For cut-off filter, analyte area will be calculated as sum from the isotopes 0 -");
    isotopePanel.add(label);
    isotopePanel.add(maxIsotope_);
    
    JLabel classLabel, unit;
    for (int i=0; i!=lipidClasses.size(); i++){
      int row = i/itemsPerLine_+1;
      int column = (i%itemsPerLine_)*3;
      classLabel = new JLabel(lipidClasses.get(i));
      classLabel.setToolTipText(tooltip_);
      inputFields.add(classLabel,new GridBagConstraints(column, row, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
      JTextField value = new JTextField(3);
      value.setHorizontalAlignment(JTextField.RIGHT);
      value.setInputVerifier(new DoubleVerifier(false,true));
      value.setText("0.0");
      value.setToolTipText(tooltip_);
      inputFields_.put(lipidClasses.get(i), value);
      inputFields.add(value,new GridBagConstraints(column+1, row, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
      unit = new JLabel("\u2030");
      unit.setToolTipText(tooltip_);
      inputFields.add(unit,new GridBagConstraints(column+2, row, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));      
    }
  }

  public void hideSettingsPanel(){
    inputPanel_.removeAll();
    this.invalidate();
    this.updateUI();
  }
  
  public Hashtable<String,Double> getCutoffs() throws AbsoluteSettingsInputException{
    Hashtable<String,Double> cutoffs = new Hashtable<String,Double>();
    for (String lipidClass : extractor_.getLipidClasses()){
      String valueString = getInputValueCorrected(inputFields_.get(lipidClass).getText(),lipidClass);
      double value = Double.parseDouble(valueString);
      value /= 1000d;
      cutoffs.put(lipidClass, value);
    }
    return cutoffs;
  }
  
  public Hashtable<String,String> getCutoffsAsString() throws AbsoluteSettingsInputException{
    Hashtable<String,String> cutoffs = new Hashtable<String,String>();
    for (String lipidClass : extractor_.getLipidClasses()){
      cutoffs.put(lipidClass,getInputValueCorrected(inputFields_.get(lipidClass).getText(),lipidClass));
    }
    return cutoffs;
  }
  
  private String getInputValueCorrected(String input, String lipidClass) throws AbsoluteSettingsInputException{
    try{
      String output = input;
      double value = Double.parseDouble(input);
      if (value<0){
        value *= -1d;
        output = input.substring(1);
      }
      if (value>=1000) throw new AbsoluteSettingsInputException("The cut-off of "+lipidClass+" cannot be greater than 1000\u2030");
      return output;
    } catch (NumberFormatException ex){
      throw new AbsoluteSettingsInputException("The cut-off of the class "+lipidClass+" is not correct! "+ex.getMessage());
    }

  }
  
  public int getMaxIsotope(){
    return Integer.parseInt(maxIsotope_.getText());
  }
  
  public void setMaxIsotope(int maxIsotope) throws SettingsException{
    if (maxIsotope<0) throw new SettingsException("The value must not be negative");
    maxIsotope_.setText(String.valueOf(maxIsotope));
  }
  
  public void setCutoff(String lipidClass, String value) throws NumberFormatException, AbsoluteSettingsInputException, SettingsException{
    double number = Double.parseDouble(value);
    if (number<0) throw new SettingsException("The value must not be negative");
    if (!inputFields_.containsKey(lipidClass)) throw new AbsoluteSettingsInputException("The lipid class "+lipidClass+" is not available as input field");
    inputFields_.get(lipidClass).setText(value);
  }
}
