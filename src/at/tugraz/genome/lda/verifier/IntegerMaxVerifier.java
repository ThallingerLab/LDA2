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

package at.tugraz.genome.lda.verifier;

import javax.swing.JFrame;

import at.tugraz.genome.lda.WarningMessage;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class IntegerMaxVerifier extends IntegerVerifier
{
  int minValue_;
  int maxValue_;
  
  public IntegerMaxVerifier(int minValue, int maxValue){
    super();
    init(minValue,maxValue);
  }
  
  public IntegerMaxVerifier(boolean mandatory, int minValue, int maxValue){
    super(mandatory);
    init(minValue,maxValue);
  }
  
  private void init (int minValue, int maxValue){
    minValue_ = minValue;
    maxValue_ = maxValue;    
  }

  protected boolean performOtherChecks(int number){
    if (number>maxValue_){
      new WarningMessage(new JFrame(), "Error", "The entered value must be lower or equal "+maxValue_);
      return false;
    }else if (number<minValue_){
      new WarningMessage(new JFrame(), "Error", "The entered value must be greater or equal "+minValue_);
      return false;
    }
    return true;
  }

}
