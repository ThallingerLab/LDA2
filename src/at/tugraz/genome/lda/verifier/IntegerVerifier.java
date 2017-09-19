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

package at.tugraz.genome.lda.verifier;

import javax.swing.InputVerifier;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JTextField;

import at.tugraz.genome.lda.WarningMessage;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class IntegerVerifier extends InputVerifier
{
  
  private boolean mandatory_;
  
  public IntegerVerifier(){
    this(false);
  }
  
  public IntegerVerifier(boolean mandatory){
    mandatory_ = mandatory;
  }

  public boolean verify(JComponent input)
  {
    String text = ((JTextField)input).getText();
    if (text.contains(",")){
      text = text.replaceAll(",", ".");
      ((JTextField)input).setText(text);
    }  
    if (text!=null&&text.length()>0){
      try{
        int number = Integer.parseInt(text);
        return performOtherChecks(number);
      } catch (NumberFormatException nfx){
        new WarningMessage(new JFrame(), "Error", "The entered value must be integer format");
        return false;
      }
    }else{
      if (mandatory_){
        new WarningMessage(new JFrame(), "Error", "A value must be entered");
        return false;        
      }
    }
    return true;
  }
  
  protected boolean performOtherChecks(int number){
    return true;
  }

}
