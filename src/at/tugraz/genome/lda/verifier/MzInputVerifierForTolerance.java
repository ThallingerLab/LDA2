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

import javax.swing.JComponent;
import javax.swing.JTextField;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MzInputVerifierForTolerance extends DoubleVerifier
{
  /** the tolerance input field, that will be updated upon a change of this field*/
  JTextField mzTolerance_;
  
  /**
   * constructor requiring the tolerance input field, that will be updated upon a change of this field
   * @param mzTolerance
   */
  public MzInputVerifierForTolerance(JTextField mzTolerance){
    super();
    this.mzTolerance_ = mzTolerance;
  }
  
  public boolean verify(JComponent input)
  {
    boolean ok = super.verify(input);
    if (ok && LipidomicsConstants.getMzUnit().equalsIgnoreCase(LipidomicsConstants.MZUNIT_PPM)){
      double value = Double.parseDouble(((JTextField)input).getText());
      mzTolerance_.setText(Calculator.FormatNumberToString(LipidomicsConstants.getCoarseChromMzTolerance((float)value), 5));
    }
    return ok;
  }
}
