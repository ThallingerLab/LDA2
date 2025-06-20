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

import java.awt.ComponentOrientation;
import java.awt.Dimension;

import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JTextField;

import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.verifier.DoubleVerifier;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class PhysicalUnitInput extends JPanel
{
  private static final long serialVersionUID = 2647998092139418279L;
  private String preSelected_;
  private JTextField value_;
  private JComboBox<String> unitMagnitude_;
  private String tooltipText_;
    
  public PhysicalUnitInput(){
    this.initComponents();
  }
  
  public PhysicalUnitInput(String preSelected, String tooltipText){
    tooltipText_ = tooltipText;
    this.initComponents();
    preSelected_ = preSelected;
    unitMagnitude_.setSelectedItem(preSelected);
  }

  private void initComponents(){
    value_ = new JTextField(3);
    value_.setHorizontalAlignment(JTextField.RIGHT);
    value_.setInputVerifier(new DoubleVerifier());
    if (tooltipText_!=null&&tooltipText_.length()>0)
      value_.setToolTipText(tooltipText_);
    this.add(value_);
//    Font toolBarFont=new Font("Helvetica",Font.PLAIN,10);
    unitMagnitude_ = new JComboBox<String>(StaticUtils.physicalMagnitudes_);
//    unitMagnitude_.setFont(toolBarFont);
    unitMagnitude_.setComponentOrientation(ComponentOrientation.RIGHT_TO_LEFT);
//    unitMagnitude_.setMinimumSize(new Dimension(1,1));
    unitMagnitude_.setPreferredSize(new Dimension(42,18));
    if (tooltipText_!=null&&tooltipText_.length()>0)
      unitMagnitude_.setToolTipText(tooltipText_);
    this.add(unitMagnitude_);
//    JLabel label = new JLabel(baseUnit_);
//    this.add(label);

  }
  
  public Double getInputValue(){
    Double value = null;
    if (value_.getText()!=null&&value_.getText().length()>0){
      try{
        value = new Double(value_.getText());
        String selectedMagnifier = (String)unitMagnitude_.getSelectedItem();
        value = StaticUtils.getValueDividedByUnit(value,selectedMagnifier);
      } catch (NumberFormatException nfx){
        nfx.printStackTrace();
      }
    }
    return value;
  }
  
  public void setInputValue(PhysicalUnitInput otherInput){
    this.setInputValue(otherInput.value_.getText(),(String)otherInput.unitMagnitude_.getSelectedItem());
  }
  
  public void setInputValue(String value, String magnitude){
    this.value_.setText(value);
    this.unitMagnitude_.setSelectedItem(magnitude);
  }

  public JTextField getValue()
  {
    return value_;
  }

  public JComboBox<String> getUnitMagnitude()
  {
    return unitMagnitude_;
  }
  
  public void setDefault(){
    this.value_.setText("");
    this.unitMagnitude_.setSelectedItem(this.preSelected_);
  }
  
  
//  public void setInputValue(Double value){
//    System.out.println("Setting value: "+value);
//    if (value>=1){
//      unitMagnitude_.setSelectedItem("");
//      value_.setText(String.valueOf(value));
//    } else if (value*1000d >= 1){
//      unitMagnitude_.setSelectedItem("m");
//      value_.setText(String.valueOf(value*1000d));
//    } else if (value*1000000d >= 1){
//      unitMagnitude_.setSelectedItem("\u03BC");
//      value_.setText(String.valueOf(value*1000000d));
//    } else if (value*1000000000d >= 1){
//      unitMagnitude_.setSelectedItem("n");
//      value_.setText(String.valueOf(value*1000000000d));
//    } else if (value*1000000000000d >= 1){
//      unitMagnitude_.setSelectedItem("p");
//      value_.setText(String.valueOf(value*1000000000000d));
//    } else if (value*1000000000000000d >= 1){
//      unitMagnitude_.setSelectedItem("f");
//      value_.setText(String.valueOf(value*1000000000000000d));
//    }else{
//      unitMagnitude_.setSelectedItem("a");
//      value_.setText(String.valueOf(value*1000000000000000000d));     
//    }
//    
//    
//  }
  
}
