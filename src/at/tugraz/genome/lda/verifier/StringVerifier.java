/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2023 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp 
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
 * @author Leonida M. Lamp
 *
 */
public class StringVerifier extends InputVerifier
{
	private int minimumTextLength_;
  private int maximumTextLength_;
  private String fieldName_;
  
  public StringVerifier(int min, int max, String fieldName){
  	maximumTextLength_ = max;
  	minimumTextLength_ = min;
  	fieldName_ = fieldName;
  }
  
  @Override
  public boolean verify(JComponent input)
  {
  	String text = ((JTextField)input).getText();
  	if (text.length() < minimumTextLength_)
  	{
  		new WarningMessage(new JFrame(), "Error", 
  				String.format("The text entered in field '%s' has to be at least %s characters long.", fieldName_, minimumTextLength_));
      return false;
  	}
  	else if (text.length() > maximumTextLength_)
  	{
  		new WarningMessage(new JFrame(), "Error", 
  				String.format("The text entered in field '%s' has to be at most %s characters long.", fieldName_, maximumTextLength_));
      return false;
  	}
  	return true;
  }

}
